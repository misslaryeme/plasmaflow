"""
Differential gene expression analysis for PlasmaFlow

This module integrates gene proximity results with differential expression data
to create comprehensive volcano plots and enrichment analyses.
"""

import logging
import os
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from scipy.stats import fisher_exact

from ..config import Config
from ..utils import get_logger

logger = get_logger(__name__)
warnings.filterwarnings('ignore')

class GeneExpressionAnalyzer:
    """Analyzer for differential gene expression near chromatin loops"""
    
    def __init__(self, config: Config):
        """
        Initialize expression analyzer
        
        Args:
            config: PlasmaFlow configuration object
        """
        self.config = config
        
        # Get expression analysis parameters
        self.expression_params = config.genomics.get("expression_analysis", {})
        
        # Default thresholds
        self.fdr_threshold = self.expression_params.get("fdr_threshold", 0.05)
        self.logfc_threshold = self.expression_params.get("logfc_threshold", 1.0)
        
        # Volcano plot parameters
        self.volcano_params = self.expression_params.get("volcano_plot_params", {
            "point_size": 1.0,
            "alpha": 0.6,
            "label_top_genes": 20
        })
        
        # Color schemes
        self.viz_params = config.genomics.get("visualization", {})
        self.category_colors = self.viz_params.get("category_colors", {
            'up': '#E69F00',
            'down': '#56B4E9', 
            'common': '#009E73'
        })
        
        self.multicolor_palette = [
            '#FF0000', '#00FF00', '#0000FF', '#FF00FF', '#FFFF00', '#00FFFF',
            '#FFA500', '#800080', '#FF1493', '#32CD32', '#FF4500', '#1E90FF'
        ]
        
        # Output directories
        self.output_dir = Path(config.output_dir) / "genomics" / "expression"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Key transition genes for plasma cell differentiation
        self.key_transition_genes = [
            "EZH2", "PAX5", "IRF4", "BCL6", "PRDM1", "CIITA", 
            "EGR1", "BAMBI", "FOS", "BACH2", "CXCR4", "ID3", 
            "MYC", "SMAD3", "MCM2", "MCM5", "TFRC", "E2F1", 
            "E2F2", "E2F3", "E2F4", "CD38", "CD27", "IGHG1",
            "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "XBP1",
            "ATF4", "CHOP", "DDIT3", "HSPA5", "CALR", "PDIA4"
        ]
    
    def analyze_differential_expression(
        self,
        de_results_file: Union[str, Path],
        loop_genes_directory: Union[str, Path],
        comparison_name: str = "prePB_vs_memB",
        method_name: str = "diffHic"
    ) -> Dict[str, Any]:
        """
        Perform comprehensive differential expression analysis
        
        Args:
            de_results_file: Path to differential expression results CSV
            loop_genes_directory: Directory containing loop-associated gene files
            comparison_name: Name of comparison (e.g., "prePB_vs_memB")
            method_name: Analysis method name
            
        Returns:
            Dictionary with analysis results
        """
        
        logger.info(f"Starting differential expression analysis for {comparison_name}")
        
        # Step 1: Get gene ID mapping
        gene_mapping = self._get_ensembl_to_symbol_mapping()
        
        # Step 2: Load and process DE results
        de_results = self._load_and_process_de_results(de_results_file, gene_mapping)
        
        # Step 3: Load loop-associated genes
        loop_gene_sets = self._load_loop_associated_genes(loop_genes_directory, comparison_name, method_name)
        
        # Step 4: Map genes and merge with DE data
        merged_results, mapping_stats = self._map_genes_and_merge(loop_gene_sets, de_results)
        
        # Step 5: Calculate enrichment statistics
        enrichment_results = self._calculate_enrichment_statistics(merged_results, de_results)
        
        # Step 6: Create visualizations
        visualization_files = self._create_volcano_plots(
            de_results, merged_results, enrichment_results, comparison_name, method_name
        )
        
        # Step 7: Export results
        export_files = self._export_comprehensive_results(
            merged_results, enrichment_results, de_results, mapping_stats, comparison_name, method_name
        )
        
        logger.info(f"Differential expression analysis completed for {comparison_name}")
        
        return {
            'comparison_name': comparison_name,
            'method_name': method_name,
            'de_results': de_results,
            'merged_results': merged_results,
            'enrichment_results': enrichment_results,
            'visualization_files': visualization_files,
            'export_files': export_files,
            'mapping_statistics': mapping_stats
        }
    
    def _get_ensembl_to_symbol_mapping(self) -> pd.DataFrame:
        """Get ENSEMBL to gene symbol mapping from BioMart"""
        
        logger.info("Fetching ENSEMBL to gene symbol mapping from BioMart...")
        
        try:
            biomart_url = "http://grch37.ensembl.org/biomart/martservice"
            
            query = \"\"\"<?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
                    <Attribute name = "ensembl_gene_id" />
                    <Attribute name = "external_gene_name" />
                    <Attribute name = "gene_biotype" />
                    <Attribute name = "chromosome_name" />
                    <Attribute name = "description" />
                </Dataset>
            </Query>\"\"\"
            
            response = requests.post(biomart_url, data={'query': query}, timeout=30)
            
            if response.status_code != 200:
                raise Exception(f"BioMart HTTP request failed with status {response.status_code}")
            
            if len(response.text.strip()) == 0:
                raise Exception("BioMart returned empty response")
            
            logger.info("Processing BioMart data...")
            
            lines = response.text.strip().split('\\n')
            mapping_data = []
            
            for line in lines:
                if line.strip():
                    parts = line.split('\\t')
                    if len(parts) >= 4:
                        ensembl_id = parts[0]
                        gene_symbol = parts[1]
                        biotype = parts[2]
                        chromosome = parts[3]
                        description = parts[4] if len(parts) > 4 else ""
                        
                        if (chromosome in [str(i) for i in range(1, 23)] + ['X', 'Y'] and
                            gene_symbol and ensembl_id):
                            mapping_data.append({
                                'ensembl_gene_id': ensembl_id,
                                'gene_symbol': gene_symbol,
                                'gene_biotype': biotype,
                                'chromosome': chromosome,
                                'description': description
                            })
            
            if len(mapping_data) == 0:
                raise Exception("No valid gene mappings found in BioMart response")
            
            mapping_df = pd.DataFrame(mapping_data)
            mapping_df = mapping_df.drop_duplicates(subset=['ensembl_gene_id'], keep='first')
            
            logger.info(f"Retrieved {len(mapping_df)} gene mappings from BioMart")
            return mapping_df
            
        except Exception as e:
            logger.error(f"BioMart download failed: {e}")
            # Return empty DataFrame as fallback
            return pd.DataFrame(columns=['ensembl_gene_id', 'gene_symbol', 'gene_biotype', 'chromosome', 'description'])
    
    def _load_and_process_de_results(
        self, 
        de_results_file: Union[str, Path], 
        gene_mapping: pd.DataFrame
    ) -> pd.DataFrame:
        """Load and process differential expression results"""
        
        logger.info("Loading differential expression results...")
        
        try:
            de_df = pd.read_csv(de_results_file)
            logger.info(f"Initial load: {len(de_df)} genes")
            
            # Fix column structure
            if 'Unnamed: 0' in de_df.columns:
                de_df['ensembl_gene_id'] = de_df['Unnamed: 0']
                de_df = de_df.drop('Unnamed: 0', axis=1)
            
            # Handle gene symbol mapping
            if 'name' in de_df.columns:
                de_df['gene_symbol'] = de_df['name']
            
            # Handle European decimal format
            numeric_columns = ['baseMean', 'log2FoldChange', 'padj']
            for col in numeric_columns:
                if col in de_df.columns:
                    de_df[col] = de_df[col].astype(str).str.replace(',', '.').astype(float)
            
            # Fill missing symbols using BioMart
            if len(gene_mapping) > 0:
                empty_symbols = de_df['gene_symbol'].isna().sum()
                if empty_symbols > 0:
                    logger.info(f"Filling {empty_symbols} missing gene symbols...")
                    
                    biomart_dict = dict(zip(gene_mapping['ensembl_gene_id'], gene_mapping['gene_symbol']))
                    
                    for idx, row in de_df[de_df['gene_symbol'].isna()].iterrows():
                        ensembl_id = row['ensembl_gene_id']
                        if pd.notna(ensembl_id) and str(ensembl_id) in biomart_dict:
                            de_df.loc[idx, 'gene_symbol'] = biomart_dict[str(ensembl_id)]
            
            # Clean data
            de_df = de_df.dropna(subset=['ensembl_gene_id'])
            de_df = de_df[de_df['ensembl_gene_id'].astype(str) != 'nan']
            de_df = de_df.drop_duplicates(subset=['ensembl_gene_id'], keep='first')
            
            # Add significance annotations
            de_df['significant'] = (
                (de_df['padj'] < self.fdr_threshold) & 
                (np.abs(de_df['log2FoldChange']) > self.logfc_threshold)
            )
            
            de_df['direction'] = 'NS'
            de_df.loc[(de_df['significant']) & (de_df['log2FoldChange'] > 0), 'direction'] = 'UP_in_treatment'
            de_df.loc[(de_df['significant']) & (de_df['log2FoldChange'] < 0), 'direction'] = 'DOWN_in_treatment'
            
            logger.info(f"DE results processed: {len(de_df)} genes, {de_df['significant'].sum()} significant")
            return de_df
            
        except Exception as e:
            logger.error(f"Error loading DE results: {e}")
            raise
    
    def _load_loop_associated_genes(
        self, 
        loop_genes_directory: Union[str, Path],
        comparison_name: str,
        method_name: str
    ) -> Dict[Tuple[str, str], List[str]]:
        """Load loop-associated genes with cluster support"""
        
        logger.info("Loading loop-associated genes...")
        
        import glob
        import re
        
        loop_gene_sets = {}
        pattern = os.path.join(loop_genes_directory, f"{comparison_name}_{method_name}_*_genes_enhanced.csv")
        csv_files = glob.glob(pattern)
        
        logger.info(f"Found {len(csv_files)} gene files")
        
        for csv_file in csv_files:
            filename = os.path.basename(csv_file)
            
            # Parse filename to extract category and cluster
            category = None
            cluster = None
            
            if filename.endswith('_genes_enhanced.csv'):
                base_name = filename.replace('_genes_enhanced.csv', '')
                parts = base_name.split('_')
                
                if len(parts) >= 5:
                    potential_category = parts[4]
                    
                    if potential_category in ['up', 'down', 'common']:
                        category = potential_category
                        
                        # Check for cluster specification
                        if len(parts) >= 6:
                            cluster_part = parts[5]
                            if cluster_part.startswith('cluster'):
                                cluster_match = re.search(r'cluster(\\d+)', cluster_part)
                                if cluster_match:
                                    cluster = cluster_match.group(1)
            
            if category is None:
                logger.warning(f"Could not parse category from {filename}")
                continue
            
            try:
                df = pd.read_csv(csv_file)
                
                if 'gene_symbol' in df.columns:
                    unique_genes = df['gene_symbol'].dropna().unique()
                    key = (category, cluster)
                    loop_gene_sets[key] = list(unique_genes)
                    
                    cluster_label = f"cluster {cluster}" if cluster else "ALL"
                    logger.info(f"Loaded {category} {cluster_label}: {len(unique_genes)} genes")
                    
            except Exception as e:
                logger.error(f"Error loading {csv_file}: {e}")
        
        logger.info(f"Loaded {len(loop_gene_sets)} gene sets")
        return loop_gene_sets
    
    def _map_genes_and_merge(
        self, 
        loop_gene_sets: Dict[Tuple[str, str], List[str]], 
        de_results: pd.DataFrame
    ) -> Tuple[Dict[Tuple[str, str], pd.DataFrame], Dict[Tuple[str, str], Dict[str, Any]]]:
        """Map loop genes to DE results with comprehensive matching"""
        
        logger.info("Mapping genes to DE results...")
        
        merged_results = {}
        mapping_stats = {}
        
        # Create lookup dictionary
        genes_with_symbols = de_results.dropna(subset=['gene_symbol'])
        unique_symbols = genes_with_symbols.drop_duplicates(subset=['gene_symbol'], keep='first')
        symbol_to_de = {row['gene_symbol']: row for _, row in unique_symbols.iterrows()}
        
        for (category, cluster), gene_symbols in loop_gene_sets.items():
            cluster_label = f"cluster {cluster}" if cluster else "ALL"
            
            mapped_genes = []
            mapping_methods = {'direct_symbol': 0, 'case_insensitive': 0, 'fuzzy_match': 0}
            
            for gene_symbol in gene_symbols:
                gene_data = None
                method_used = None
                
                # Strategy 1: Direct symbol match
                if gene_symbol in symbol_to_de:
                    gene_data = symbol_to_de[gene_symbol].copy()
                    method_used = 'direct_symbol'
                    mapping_methods['direct_symbol'] += 1
                
                # Strategy 2: Case-insensitive match
                elif method_used is None:
                    for de_symbol in symbol_to_de.keys():
                        if isinstance(de_symbol, str) and isinstance(gene_symbol, str):
                            if de_symbol.lower() == gene_symbol.lower():
                                gene_data = symbol_to_de[de_symbol].copy()
                                method_used = 'case_insensitive'
                                mapping_methods['case_insensitive'] += 1
                                break
                
                # Strategy 3: Fuzzy matching
                if method_used is None:
                    variations = [
                        gene_symbol.replace('-', ''),
                        gene_symbol.replace('_', ''),
                        gene_symbol.replace('-', '_'),
                        gene_symbol.replace('_', '-')
                    ]
                    
                    for variation in variations:
                        if variation in symbol_to_de:
                            gene_data = symbol_to_de[variation].copy()
                            method_used = 'fuzzy_match'
                            mapping_methods['fuzzy_match'] += 1
                            break
                
                # If gene was found, add to mapped list
                if gene_data is not None:
                    if isinstance(gene_data, pd.Series):
                        gene_data_dict = gene_data.to_dict()
                    else:
                        gene_data_dict = dict(gene_data)
                    
                    gene_data_dict['original_gene_symbol'] = gene_symbol
                    gene_data_dict['loop_category'] = category
                    gene_data_dict['loop_cluster'] = cluster
                    mapped_genes.append(gene_data_dict)
            
            # Create DataFrame
            if mapped_genes:
                merged_df = pd.DataFrame(mapped_genes).reset_index(drop=True)
            else:
                merged_df = pd.DataFrame()
            
            # Store results
            key = (category, cluster)
            merged_results[key] = merged_df
            
            # Calculate statistics
            total_mapped = len(merged_df)
            mapping_stats[key] = {
                'category': category,
                'cluster': cluster,
                'total_genes': len(gene_symbols),
                'total_mapped': total_mapped,
                'mapping_rate': total_mapped/len(gene_symbols) if len(gene_symbols) > 0 else 0,
                **mapping_methods
            }
            
            logger.info(f"{category} {cluster_label}: {total_mapped}/{len(gene_symbols)} mapped "
                       f"({total_mapped/len(gene_symbols)*100:.1f}%)")
        
        return merged_results, mapping_stats
    
    def _calculate_enrichment_statistics(
        self, 
        merged_results: Dict[Tuple[str, str], pd.DataFrame], 
        de_results: pd.DataFrame
    ) -> Dict[Tuple[str, str], Dict[str, Any]]:
        """Calculate enrichment statistics using Fisher's exact test"""
        
        logger.info("Calculating enrichment statistics...")
        
        total_genes = len(de_results)
        total_de_genes = de_results['significant'].sum()
        background_de_rate = total_de_genes / total_genes
        
        enrichment_results = {}
        
        for key, merged_data in merged_results.items():
            category, cluster = key
            
            if len(merged_data) == 0:
                continue
            
            de_in_category = merged_data['significant'].sum()
            total_in_category = len(merged_data)
            category_de_rate = de_in_category / total_in_category if total_in_category > 0 else 0
            
            # Construct contingency table for Fisher's exact test
            de_not_in_category = total_de_genes - de_in_category
            not_de_in_category = total_in_category - de_in_category
            not_de_not_in_category = total_genes - total_in_category - de_not_in_category
            
            # Fisher's exact test
            if all(val >= 0 for val in [de_not_in_category, not_de_in_category, not_de_not_in_category]):
                contingency_table = [
                    [de_in_category, not_de_in_category],
                    [de_not_in_category, not_de_not_in_category]
                ]
                
                try:
                    odds_ratio, fisher_pvalue = fisher_exact(contingency_table, alternative='greater')
                except:
                    odds_ratio, fisher_pvalue = 1.0, 1.0
            else:
                odds_ratio, fisher_pvalue = 1.0, 1.0
            
            # Count directions
            up_genes = (merged_data['direction'] == 'UP_in_treatment').sum()
            down_genes = (merged_data['direction'] == 'DOWN_in_treatment').sum()
            
            enrichment_results[key] = {
                'category': category,
                'cluster': cluster,
                'total_genes': total_in_category,
                'de_genes': de_in_category,
                'de_rate': category_de_rate,
                'background_de_rate': background_de_rate,
                'fold_enrichment': category_de_rate / background_de_rate if background_de_rate > 0 else 0,
                'odds_ratio': odds_ratio,
                'fisher_pvalue': fisher_pvalue,
                'up_genes': up_genes,
                'down_genes': down_genes
            }
            
            cluster_label = f"cluster {cluster}" if cluster else "ALL"
            logger.info(f"{category} {cluster_label}: {category_de_rate:.3f} DE rate, "
                       f"{enrichment_results[key]['fold_enrichment']:.2f}x enrichment "
                       f"(p={fisher_pvalue:.3e})")
        
        return enrichment_results
    
    def _create_volcano_plots(
        self,
        de_results: pd.DataFrame,
        merged_results: Dict[Tuple[str, str], pd.DataFrame],
        enrichment_results: Dict[Tuple[str, str], Dict[str, Any]],
        comparison_name: str,
        method_name: str
    ) -> Dict[str, Path]:
        """Create comprehensive volcano plots"""
        
        logger.info("Creating volcano plots...")
        
        visualization_files = {}
        
        # Create main categories plot
        main_plot = self._create_main_categories_volcano_plot(
            de_results, merged_results, enrichment_results, comparison_name, method_name
        )
        if main_plot:
            visualization_files['main_categories'] = main_plot
        
        # Create category-specific cluster plots
        for category in ['up', 'down', 'common']:
            cluster_plot = self._create_category_clusters_volcano_plot(
                de_results, merged_results, enrichment_results, category, comparison_name, method_name
            )
            if cluster_plot:
                visualization_files[f'{category}_clusters'] = cluster_plot
        
        return visualization_files
    
    def _create_main_categories_volcano_plot(
        self,
        de_results: pd.DataFrame,
        merged_results: Dict[Tuple[str, str], pd.DataFrame],
        enrichment_results: Dict[Tuple[str, str], Dict[str, Any]],
        comparison_name: str,
        method_name: str
    ) -> Optional[Path]:
        """Create volcano plot for main categories"""
        
        fig, ax = plt.subplots(1, 1, figsize=(16, 12))
        
        # Background genes
        x = -np.log10(de_results['padj'].fillna(1))
        y = de_results['log2FoldChange']
        ax.scatter(x, y, c='lightgray', alpha=0.3, s=3, label='Background genes')
        
        # Main categories
        main_categories = ['up', 'down', 'common']
        main_category_keys = [(category, None) for category in main_categories 
                             if (category, None) in merged_results]
        
        main_category_colors = {
            'up': '#DC143C',      # Crimson
            'down': '#8A2BE2',    # Blue Violet
            'common': '#FFD700'   # Gold
        }
        
        total_highlighted = 0
        legend_entries = []
        
        for key in main_category_keys:
            category, cluster = key
            merged_data = merged_results[key]
            
            if len(merged_data) == 0:
                continue
                
            enrich_data = enrichment_results[key]
            
            # Get genes for this category
            loop_ensembl_ids = set(merged_data['ensembl_gene_id'].unique())
            gene_mask = de_results['ensembl_gene_id'].isin(loop_ensembl_ids)
            highlighted_count = gene_mask.sum()
            total_highlighted += highlighted_count
            
            if highlighted_count > 0:
                color = main_category_colors.get(category, 'gray')
                
                # Adjust point size based on group size
                if highlighted_count > 1000:
                    point_size = 12
                    point_alpha = 0.6
                elif highlighted_count > 500:
                    point_size = 20
                    point_alpha = 0.8
                else:
                    point_size = 30
                    point_alpha = 1.0
                
                ax.scatter(x[gene_mask], y[gene_mask], 
                          c=color, alpha=point_alpha, s=point_size,
                          label=f'{category.upper()} ({highlighted_count} genes)',
                          edgecolors='black', linewidth=0.5)
                
                # Add gene labels for smaller groups
                if highlighted_count <= 100:
                    self._add_gene_labels(ax, de_results[gene_mask], top_n=10)
                
                legend_entries.append((category, color, highlighted_count))
        
        # Add key transition gene labels
        self._add_key_transition_gene_labels(ax, de_results)
        
        # Formatting
        ax.axhline(y=self.logfc_threshold, color='black', linestyle='--', alpha=0.5)
        ax.axhline(y=-self.logfc_threshold, color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=-np.log10(self.fdr_threshold), color='black', linestyle='--', alpha=0.5)
        
        ax.set_xlabel('-log10(adjusted p-value)', fontsize=14)
        ax.set_ylabel('log2(Fold Change)', fontsize=14)
        ax.set_title(f'{comparison_name} - {method_name}: Loop Categories\\n'
                    f'({total_highlighted} genes from differential loops)', fontsize=16)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', framealpha=0.95)
        
        plt.tight_layout()
        
        # Save plot
        plot_file = self.output_dir / f'{comparison_name}_{method_name}_main_categories_volcano.png'
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Main categories volcano plot saved: {plot_file}")
        return plot_file
    
    def _create_category_clusters_volcano_plot(
        self,
        de_results: pd.DataFrame,
        merged_results: Dict[Tuple[str, str], pd.DataFrame],
        enrichment_results: Dict[Tuple[str, str], Dict[str, Any]],
        category: str,
        comparison_name: str,
        method_name: str
    ) -> Optional[Path]:
        """Create volcano plot for clusters within a category"""
        
        # Get cluster keys for this category
        category_keys = [key for key in merged_results.keys() 
                        if key[0] == category and key[1] is not None]
        
        if len(category_keys) < 2:
            return None
        
        fig, ax = plt.subplots(1, 1, figsize=(14, 10))
        
        # Background genes
        x = -np.log10(de_results['padj'].fillna(1))
        y = de_results['log2FoldChange']
        ax.scatter(x, y, c='lightgray', alpha=0.3, s=4, label='Background genes')
        
        category_keys = sorted(category_keys, key=lambda x: int(x[1]) if x[1] is not None else 0)
        
        total_genes = 0
        color_index = 0
        
        for key in category_keys:
            _, cluster = key
            merged_data = merged_results[key]
            
            if len(merged_data) == 0:
                continue
            
            enrich_data = enrichment_results[key]
            
            # Get genes
            loop_ensembl_ids = set(merged_data['ensembl_gene_id'].unique())
            gene_mask = de_results['ensembl_gene_id'].isin(loop_ensembl_ids)
            highlighted_count = gene_mask.sum()
            total_genes += highlighted_count
            
            if highlighted_count > 0:
                color = self.multicolor_palette[color_index % len(self.multicolor_palette)]
                
                ax.scatter(x[gene_mask], y[gene_mask], 
                          c=color, alpha=0.8, s=25,
                          label=f'Cluster {cluster} ({highlighted_count} genes)',
                          edgecolors='black', linewidth=0.3)
                
                # Add gene labels for smaller clusters
                if highlighted_count <= 50:
                    self._add_gene_labels(ax, de_results[gene_mask], top_n=5)
                
                color_index += 1
        
        # Formatting
        ax.axhline(y=self.logfc_threshold, color='black', linestyle='--', alpha=0.5)
        ax.axhline(y=-self.logfc_threshold, color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=-np.log10(self.fdr_threshold), color='black', linestyle='--', alpha=0.5)
        
        ax.set_xlabel('-log10(adjusted p-value)', fontsize=12)
        ax.set_ylabel('log2(Fold Change)', fontsize=12)
        ax.set_title(f'{category.upper()} Clusters - {comparison_name}\\n'
                    f'({total_genes} total genes)', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', framealpha=0.95, fontsize=10)
        
        plt.tight_layout()
        
        # Save plot
        plot_file = self.output_dir / f'{comparison_name}_{method_name}_{category}_clusters_volcano.png'
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"{category} clusters volcano plot saved: {plot_file}")
        return plot_file
    
    def _add_gene_labels(self, ax, gene_subset: pd.DataFrame, top_n: int = 10):
        """Add gene labels for top significant genes"""
        
        if len(gene_subset) == 0:
            return
        
        # Filter for valid genes and calculate combined score
        valid_genes = gene_subset[
            gene_subset['padj'].notna() & 
            (gene_subset['padj'] > 0) &
            gene_subset['log2FoldChange'].notna()
        ].copy()
        
        if len(valid_genes) == 0:
            return
        
        valid_genes['significance_score'] = -np.log10(valid_genes['padj'])
        valid_genes['effect_size'] = np.abs(valid_genes['log2FoldChange'])
        valid_genes['combined_score'] = valid_genes['significance_score'] * valid_genes['effect_size']
        
        # Get top genes
        top_genes = valid_genes.nlargest(top_n, 'combined_score')
        
        for _, row in top_genes.iterrows():
            if pd.notna(row['gene_symbol']):
                x_pos = -np.log10(row['padj'])
                y_pos = row['log2FoldChange']
                
                ax.annotate(row['gene_symbol'], 
                           xy=(x_pos, y_pos),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=8, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.2', 
                                    facecolor='white', 
                                    edgecolor='black',
                                    alpha=0.8),
                           arrowprops=dict(arrowstyle='->', 
                                         color='black',
                                         alpha=0.7,
                                         lw=0.5))
    
    def _add_key_transition_gene_labels(self, ax, de_results: pd.DataFrame):
        """Add labels for key transition genes that are significantly DE"""
        
        # Filter for key transition genes that are significant
        key_genes = de_results[
            de_results['gene_symbol'].isin(self.key_transition_genes)
        ]
        
        significant_key_genes = key_genes[
            (key_genes['padj'] < self.fdr_threshold) & 
            (key_genes['log2FoldChange'].abs() > self.logfc_threshold)
        ]
        
        logger.info(f"Adding labels for {len(significant_key_genes)} significant key transition genes")
        
        for _, row in significant_key_genes.iterrows():
            x_pos = -np.log10(row['padj'])
            y_pos = row['log2FoldChange']
            gene_name = row['gene_symbol']
            
            ax.annotate(gene_name, 
                       xy=(x_pos, y_pos),
                       xytext=(10, 10), textcoords='offset points',
                       fontsize=10, fontweight='bold',
                       color='darkred',
                       bbox=dict(boxstyle='round,pad=0.3',     
                                facecolor='gold',               
                                edgecolor='darkorange',
                                alpha=0.95,                     
                                linewidth=1.2),
                       arrowprops=dict(arrowstyle='->', 
                                     color='darkorange',
                                     alpha=0.9,
                                     linewidth=1.5))
    
    def _export_comprehensive_results(
        self,
        merged_results: Dict[Tuple[str, str], pd.DataFrame],
        enrichment_results: Dict[Tuple[str, str], Dict[str, Any]],
        de_results: pd.DataFrame,
        mapping_stats: Dict[Tuple[str, str], Dict[str, Any]],
        comparison_name: str,
        method_name: str
    ) -> Dict[str, Path]:
        """Export comprehensive analysis results"""
        
        logger.info("Exporting comprehensive results...")
        
        export_files = {}
        
        # Export DE results
        de_file = self.output_dir / f'{comparison_name}_{method_name}_DE_results.csv'
        de_results.to_csv(de_file, index=False)
        export_files['de_results'] = de_file
        
        # Export key transition genes analysis
        key_genes_analysis = self._analyze_key_transition_genes(de_results)
        key_genes_file = self.output_dir / f'{comparison_name}_{method_name}_key_transition_genes.csv'
        key_genes_analysis.to_csv(key_genes_file, index=False)
        export_files['key_transition_genes'] = key_genes_file
        
        # Export summary table
        summary_data = []
        for key, enrich in enrichment_results.items():
            category, cluster = key
            mapping_stat = mapping_stats[key]
            
            cluster_label = f"cluster {cluster}" if cluster else "ALL"
            
            summary_data.append({
                'Category': category,
                'Cluster': cluster_label,
                'Total_Genes': enrich['total_genes'],
                'DE_Genes': enrich['de_genes'],
                'DE_Rate': f"{enrich['de_rate']*100:.1f}%",
                'Up_Genes': enrich['up_genes'],
                'Down_Genes': enrich['down_genes'],
                'Fold_Enrichment': f"{enrich['fold_enrichment']:.2f}x",
                'Fisher_pvalue': f"{enrich['fisher_pvalue']:.2e}",
                'Mapping_Rate': f"{mapping_stat['mapping_rate']*100:.1f}%"
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_file = self.output_dir / f'{comparison_name}_{method_name}_expression_summary.csv'
        summary_df.to_csv(summary_file, index=False)
        export_files['summary'] = summary_file
        
        # Export individual category results
        for key, merged_data in merged_results.items():
            category, cluster = key
            if len(merged_data) > 0:
                cluster_label = f"cluster{cluster}" if cluster else "ALL"
                detail_file = self.output_dir / f'{comparison_name}_{method_name}_{category}_{cluster_label}_genes.csv'
                merged_data.to_csv(detail_file, index=False)
                export_files[f'{category}_{cluster_label}'] = detail_file
        
        logger.info(f"Exported {len(export_files)} result files")
        return export_files
    
    def _analyze_key_transition_genes(self, de_results: pd.DataFrame) -> pd.DataFrame:
        """Analyze key transition genes in detail"""
        
        key_genes = de_results[de_results['gene_symbol'].isin(self.key_transition_genes)].copy()
        
        if len(key_genes) == 0:
            return pd.DataFrame()
        
        # Add transition relevance scores
        key_genes['is_significant'] = (
            (key_genes['padj'] < self.fdr_threshold) & 
            (key_genes['log2FoldChange'].abs() > self.logfc_threshold)
        )
        
        key_genes['transition_relevance'] = 'Low'
        
        # High relevance: significant and strong effect
        high_relevance_mask = (
            key_genes['is_significant'] & 
            (key_genes['log2FoldChange'].abs() > 2.0)
        )
        key_genes.loc[high_relevance_mask, 'transition_relevance'] = 'High'
        
        # Medium relevance: significant but moderate effect
        medium_relevance_mask = (
            key_genes['is_significant'] & 
            (key_genes['log2FoldChange'].abs() <= 2.0)
        )
        key_genes.loc[medium_relevance_mask, 'transition_relevance'] = 'Medium'
        
        # Sort by significance and effect size
        key_genes = key_genes.sort_values(
            ['is_significant', 'padj', 'log2FoldChange'], 
            ascending=[False, True, False]
        )
        
        logger.info(f"Key transition genes: {len(key_genes)} found, "
                   f"{key_genes['is_significant'].sum()} significant")
        
        return key_genes


def analyze_expression_overlap(
    de_results_file: Union[str, Path],
    loop_genes_directory: Union[str, Path],
    comparison_name: str = "prePB_vs_memB",
    method_name: str = "diffHic",
    config: Optional[Config] = None
) -> Dict[str, Any]:
    """
    Convenience function to analyze expression overlap
    
    Args:
        de_results_file: Path to differential expression results
        loop_genes_directory: Directory with loop gene files
        comparison_name: Name of comparison
        method_name: Analysis method name
        config: Configuration object
        
    Returns:
        Dictionary with analysis results
    """
    
    if config is None:
        from ..config import get_default_config
        config = get_default_config()
    
    analyzer = GeneExpressionAnalyzer(config)
    return analyzer.analyze_differential_expression(
        de_results_file, loop_genes_directory, comparison_name, method_name
    )
"""
BED file utilities for loop coordinate processing
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd

logger = logging.getLogger(__name__)


def parse_loop_coordinates(
    loop_id: str,
) -> Optional[Tuple[str, int, int, str, int, int]]:
    """
    Parse loop_id format: chr2:111430000-111440000_chr2:111500000-111510000

    Args:
        loop_id: Loop identifier string

    Returns:
        Tuple of (chr1, start1, end1, chr2, start2, end2) or None if parsing fails
    """
    try:
        # Split into two anchors
        anchor1_str, anchor2_str = loop_id.split("_")

        # Parse anchor1: chr2:111430000-111440000
        chr1 = anchor1_str.split(":")[0]
        coords1 = anchor1_str.split(":")[1]
        start1, end1 = coords1.split("-")

        # Parse anchor2: chr2:111500000-111510000
        chr2 = anchor2_str.split(":")[0]
        coords2 = anchor2_str.split(":")[1]
        start2, end2 = coords2.split("-")

        return chr1, int(start1), int(end1), chr2, int(start2), int(end2)

    except (ValueError, IndexError) as e:
        logger.warning(f"Error parsing loop_id '{loop_id}': {e}")
        return None


def create_anchor_beds(
    results_df: pd.DataFrame,
    category_name: str,
    comparison_name: str,
    output_dir: Union[str, Path],
) -> Dict[str, Path]:
    """
    Create BED files for loop anchors from results DataFrame

    Args:
        results_df: DataFrame with loop results containing 'loop_id' column
        category_name: Category name (e.g., 'up', 'down', 'common')
        comparison_name: Comparison name (e.g., 'prePB_vs_memB')
        output_dir: Output directory for BED files

    Returns:
        Dictionary mapping anchor names to BED file paths
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Prepare anchor data
    anchor1_data = []
    anchor2_data = []

    valid_loops = 0

    for idx, row in results_df.iterrows():
        # Get loop_id from different possible column names
        loop_id = None
        for col in ["loop_id", "interaction_id", "id"]:
            if col in row and pd.notna(row[col]):
                loop_id = row[col]
                break

        if loop_id is None:
            logger.warning(f"No loop_id found in row {idx}")
            continue

        # Parse coordinates
        parsed = parse_loop_coordinates(str(loop_id))
        if parsed is None:
            continue

        chr1, start1, end1, chr2, start2, end2 = parsed

        # Create unique loop names
        loop_name = f"{comparison_name}_{category_name}_loop{valid_loops + 1}"

        # Add to anchor lists (BED format: chr, start, end, name, score, strand)
        anchor1_data.append([chr1, start1, end1, f"{loop_name}_anchor1", ".", "."])

        anchor2_data.append([chr2, start2, end2, f"{loop_name}_anchor2", ".", "."])

        valid_loops += 1

    if valid_loops == 0:
        logger.warning(f"No valid loops found for {category_name} category")
        return {}

    logger.info(f"Creating BED files for {category_name}: {valid_loops} valid loops")

    # Write BED files
    bed_files = {}

    # Anchor 1 BED file
    anchor1_file = output_path / f"{comparison_name}_{category_name}_anchor1.bed"
    with open(anchor1_file, "w") as f:
        for row in anchor1_data:
            f.write("\t".join(map(str, row)) + "\n")

    bed_files[f"{category_name}_anchor1"] = anchor1_file

    # Anchor 2 BED file
    anchor2_file = output_path / f"{comparison_name}_{category_name}_anchor2.bed"
    with open(anchor2_file, "w") as f:
        for row in anchor2_data:
            f.write("\t".join(map(str, row)) + "\n")

    bed_files[f"{category_name}_anchor2"] = anchor2_file

    logger.info(f"Created BED files: {anchor1_file.name}, {anchor2_file.name}")

    return bed_files


def validate_bed_file(bed_file: Union[str, Path]) -> bool:
    """
    Validate BED file format and content

    Args:
        bed_file: Path to BED file

    Returns:
        True if valid, False otherwise
    """
    bed_path = Path(bed_file)

    if not bed_path.exists():
        logger.error(f"BED file does not exist: {bed_path}")
        return False

    try:
        with open(bed_path, "r") as f:
            lines = f.readlines()

        if len(lines) == 0:
            logger.error(f"BED file is empty: {bed_path}")
            return False

        # Check first few lines for proper format
        for i, line in enumerate(lines[:5]):
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")

            # BED files need at least 3 columns (chr, start, end)
            if len(fields) < 3:
                logger.error(f"BED file line {i+1} has < 3 fields: {bed_path}")
                return False

            # Check that start and end are integers
            try:
                start = int(fields[1])
                end = int(fields[2])

                if start >= end:
                    logger.error(
                        f"Invalid coordinates in BED file line {i+1}: start >= end"
                    )
                    return False

            except ValueError:
                logger.error(
                    f"Non-integer coordinates in BED file line {i+1}: {bed_path}"
                )
                return False

        return True

    except Exception as e:
        logger.error(f"Error validating BED file {bed_path}: {e}")
        return False


def merge_bed_files(
    bed_files: List[Union[str, Path]],
    output_file: Union[str, Path],
    sort_output: bool = True,
) -> Path:
    """
    Merge multiple BED files into one

    Args:
        bed_files: List of BED file paths to merge
        output_file: Output merged BED file path
        sort_output: Whether to sort the output by coordinates

    Returns:
        Path to merged BED file
    """
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    all_regions = []

    for bed_file in bed_files:
        bed_path = Path(bed_file)
        if not bed_path.exists():
            logger.warning(f"BED file not found, skipping: {bed_path}")
            continue

        with open(bed_path, "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    all_regions.append(line)

    # Sort by coordinates if requested
    if sort_output:

        def sort_key(line):
            fields = line.split("\t")
            try:
                return (fields[0], int(fields[1]), int(fields[2]))
            except (IndexError, ValueError):
                return (fields[0], 0, 0)

        all_regions.sort(key=sort_key)

    # Write merged file
    with open(output_path, "w") as f:
        for region in all_regions:
            f.write(region + "\n")

    logger.info(
        f"Merged {len(bed_files)} BED files into {output_path} ({len(all_regions)} regions)"
    )

    return output_path


def bed_to_dataframe(bed_file: Union[str, Path]) -> pd.DataFrame:
    """
    Load BED file into pandas DataFrame

    Args:
        bed_file: Path to BED file

    Returns:
        DataFrame with BED regions
    """
    bed_path = Path(bed_file)

    if not bed_path.exists():
        raise FileNotFoundError(f"BED file not found: {bed_path}")

    # Define column names for BED format
    bed_columns = ["chr", "start", "end", "name", "score", "strand"]

    # Read file and determine number of columns
    with open(bed_path, "r") as f:
        first_line = f.readline().strip()

    if not first_line:
        return pd.DataFrame(columns=bed_columns[:3])

    n_cols = len(first_line.split("\t"))
    column_names = bed_columns[:n_cols]

    # Read BED file
    bed_df = pd.read_csv(
        bed_path, sep="\t", header=None, names=column_names, comment="#"
    )

    return bed_df


def dataframe_to_bed(
    df: pd.DataFrame, output_file: Union[str, Path], columns: Optional[List[str]] = None
) -> Path:
    """
    Save DataFrame to BED file

    Args:
        df: DataFrame with genomic regions
        output_file: Output BED file path
        columns: Column names to include (default: first 3-6 columns)

    Returns:
        Path to output BED file
    """
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if columns is None:
        # Use first few columns that match BED format
        available_cols = df.columns.tolist()
        bed_columns = ["chr", "start", "end", "name", "score", "strand"]
        columns = []

        for col in bed_columns:
            if col in available_cols:
                columns.append(col)
            elif len(columns) < 3:
                # Need at least chr, start, end
                if len(available_cols) > len(columns):
                    columns.append(available_cols[len(columns)])

        if len(columns) < 3:
            raise ValueError("DataFrame must have at least 3 columns for BED format")

    # Save to BED file
    df[columns].to_csv(output_path, sep="\t", header=False, index=False)

    logger.info(f"Saved DataFrame to BED file: {output_path} ({len(df)} regions)")

    return output_path

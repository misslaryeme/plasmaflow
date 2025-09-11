"""
Command-line interface for PlasmaFlow
"""

import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import click
import yaml

from . import __version__, check_dependencies, get_info
from .config import Config, get_default_config, save_config
from .core import PlasmaFlowAnalysis
from .utils import setup_logging


# Global context for CLI
class CLIContext:
    def __init__(self):
        self.config_file: Optional[Path] = None
        self.config: Optional[Config] = None
        self.verbose: bool = False
        self.quiet: bool = False


@click.group()
@click.version_option(__version__)
@click.option(
    "--config", "-c", type=click.Path(exists=True), help="Configuration file path"
)
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
@click.option("--quiet", "-q", is_flag=True, help="Enable quiet mode (minimal output)")
@click.pass_context
def main(ctx, config, verbose, quiet):
    """
    PlasmaFlow: Hi-C chromatin loop analysis pipeline for plasma cell differentiation

    PlasmaFlow provides a comprehensive workflow for analyzing Hi-C chromatin
    interaction data through the lens of plasma cell differentiation, from loop
    calling to pathway enrichment analysis.
    """
    # Create CLI context
    cli_ctx = CLIContext()
    cli_ctx.verbose = verbose
    cli_ctx.quiet = quiet

    # Set up logging
    log_level = "DEBUG" if verbose else "WARNING" if quiet else "INFO"
    setup_logging(level=log_level)

    # Load config if provided
    if config:
        cli_ctx.config_file = Path(config)
        from .config import load_config

        cli_ctx.config = load_config(cli_ctx.config_file)

    ctx.obj = cli_ctx


@main.command()
def info():
    """Show PlasmaFlow package information"""

    info_data = get_info()

    click.echo("=" * 50)
    click.echo(f"PlasmaFlow v{info_data['version']}")
    click.echo("=" * 50)
    click.echo(f"Description: {info_data['description']}")
    click.echo(f"Python version: {info_data['python_version']}")
    click.echo()

    click.echo("Available modules:")
    for module in info_data["modules"]:
        click.echo(f"  - {module}")
    click.echo()

    # Check dependencies
    deps = check_dependencies()
    click.echo("Dependency status:")
    for dep, available in deps.items():
        status = "✓" if available else "✗"
        click.echo(f"  {status} {dep}")


@main.command()
@click.argument("output_file", type=click.Path())
@click.option(
    "--format",
    type=click.Choice(["yaml", "json"]),
    default="yaml",
    help="Output format for configuration file",
)
def init_config(output_file, format):
    """Initialize a new PlasmaFlow configuration file"""

    output_path = Path(output_file)

    if output_path.exists():
        if not click.confirm(f"File {output_path} already exists. Overwrite?"):
            click.echo("Configuration initialization cancelled.")
            return

    # Create default config
    config = get_default_config()

    try:
        if format == "json":
            # Convert to dict and save as JSON
            config_dict = {
                "project_name": config.project_name,
                "resolution": config.resolution,
                "random_seed": config.random_seed,
                "n_threads": config.n_threads,
                "input_dir": config.input_dir,
                "output_dir": config.output_dir,
                "temp_dir": config.temp_dir,
                "samples": config.samples,
                "cell_types": config.cell_types,
                "loop_calling": config.loop_calling,
                "quality_control": config.quality_control,
                "differential": config.differential,
                "visualization": config.visualization,
                "genomics": config.genomics,
                "enrichment": config.enrichment,
                "r_config": config.r_config,
            }

            with open(output_path, "w") as f:
                json.dump(config_dict, f, indent=2)
        else:
            # Save as YAML using the save_config function
            save_config(config, output_path)

        click.echo(f"Configuration file created: {output_path}")
        click.echo("Edit this file to customize your analysis parameters.")

    except Exception as e:
        click.echo(f"Error creating configuration file: {e}", err=True)
        sys.exit(1)


@main.command()
@click.pass_context
def run(ctx):
    """Run the complete PlasmaFlow analysis pipeline"""

    cli_ctx = ctx.obj

    if cli_ctx.config is None:
        click.echo(
            "Error: No configuration file provided. Use --config option or 'plasmaflow init-config'",
            err=True,
        )
        sys.exit(1)

    try:
        # Initialize analysis
        analysis = PlasmaFlowAnalysis(
            config=cli_ctx.config,
            log_level=(
                "DEBUG" if cli_ctx.verbose else "WARNING" if cli_ctx.quiet else "INFO"
            ),
        )

        # Run full pipeline
        click.echo("Starting PlasmaFlow analysis pipeline...")
        results = analysis.run_full_pipeline()

        # Save results
        results_file = Path(cli_ctx.config.output_dir) / "plasmaflow_results.pkl"
        analysis.save_results(results_file)

        click.echo(f"Analysis completed. Results saved to: {results_file}")

        # Print summary
        execution_times = analysis.get_execution_times()
        total_time = execution_times.get("total", 0)

        click.echo(f"Total execution time: {total_time:.2f} seconds")

        # Show step results
        success_count = 0
        for step, result in results.items():
            if isinstance(result, dict) and result.get("success", False):
                success_count += 1
                status = "✓"
            else:
                status = "✗"
            click.echo(f"  {status} {step}")

        click.echo(
            f"Successfully completed {success_count}/{len(results)} pipeline steps"
        )

    except Exception as e:
        click.echo(f"Pipeline execution failed: {e}", err=True)
        if cli_ctx.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


@main.command()
@click.option(
    "--samples", type=click.Path(exists=True), help="CSV file with sample information"
)
@click.option("--output", "-o", type=click.Path(), help="Output directory")
@click.pass_context
def loop_calling(ctx, samples, output):
    """Run only the loop calling step"""

    cli_ctx = ctx.obj

    if cli_ctx.config is None:
        click.echo("Error: Configuration required", err=True)
        sys.exit(1)

    # Override output directory if provided
    if output:
        cli_ctx.config.output_dir = str(output)

    try:
        analysis = PlasmaFlowAnalysis(cli_ctx.config)

        # Load samples if provided
        sample_data = None
        if samples:
            import pandas as pd

            samples_df = pd.read_csv(samples)
            sample_data = {}
            for _, row in samples_df.iterrows():
                sample_data[row["sample_name"]] = {
                    "cool_file": row.get("cool_file"),
                    "read_count": row.get("read_count"),
                }

        click.echo("Running loop calling analysis...")
        results = analysis.run_loop_calling(sample_data)

        if results["success"]:
            click.echo("Loop calling completed successfully")

            # Show summary
            summary = results["summary"]
            click.echo(f"Processed {len(summary)} samples:")
            for _, row in summary.iterrows():
                status = "✓" if row["success"] else "✗"
                click.echo(
                    f"  {status} {row['sample']}: {row.get('num_loops', 0)} loops"
                )
        else:
            click.echo("Loop calling failed", err=True)
            sys.exit(1)

    except Exception as e:
        click.echo(f"Loop calling failed: {e}", err=True)
        if cli_ctx.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


@main.command()
@click.option(
    "--counts",
    type=click.Path(exists=True),
    required=True,
    help="Counts matrix CSV file",
)
@click.option(
    "--metadata",
    type=click.Path(exists=True),
    required=True,
    help="Sample metadata CSV file",
)
@click.option(
    "--methods",
    multiple=True,
    type=click.Choice(["diffHic", "DESeq2", "edgeR"]),
    default=["diffHic"],
    help="Differential analysis methods to use",
)
@click.option("--output", "-o", type=click.Path(), help="Output directory")
@click.pass_context
def differential(ctx, counts, metadata, methods, output):
    """Run differential analysis"""

    cli_ctx = ctx.obj

    if cli_ctx.config is None:
        config = get_default_config()
    else:
        config = cli_ctx.config

    # Override output directory if provided
    if output:
        config.output_dir = str(output)

    try:
        analysis = PlasmaFlowAnalysis(config)

        click.echo(f"Running differential analysis with methods: {', '.join(methods)}")

        # Override methods in config
        analysis.config.differential["methods"] = list(methods)

        # Run differential analysis
        results = analysis.differential_analyzer.run_differential_analysis(
            counts_matrix=counts, metadata=metadata, methods=list(methods)
        )

        # Show results
        click.echo("Differential analysis completed:")
        for result_key, result in results.items():
            status = "✓" if result.success else "✗"
            if result.success:
                click.echo(
                    f"  {status} {result_key}: {result.n_significant} significant"
                )
            else:
                click.echo(f"  {status} {result_key}: {result.error_message}")

    except Exception as e:
        click.echo(f"Differential analysis failed: {e}", err=True)
        if cli_ctx.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


@main.command()
@click.argument("config_file", type=click.Path(exists=True))
def validate_config(config_file):
    """Validate a PlasmaFlow configuration file"""

    try:
        from .config import load_config
        from .config import validate_config as validate_config_func

        # Load config
        config = load_config(config_file)
        click.echo(f"Configuration loaded successfully: {config_file}")

        # Validate config
        issues = validate_config_func(config)

        if not issues:
            click.echo("✓ Configuration is valid")
        else:
            click.echo("Configuration issues found:")
            for issue in issues:
                click.echo(f"  ✗ {issue}")
            sys.exit(1)

    except Exception as e:
        click.echo(f"Configuration validation failed: {e}", err=True)
        sys.exit(1)


@main.command()
def check_env():
    """Check PlasmaFlow environment and dependencies"""

    click.echo("Checking PlasmaFlow environment...")
    click.echo()

    # Check Python dependencies
    deps = check_dependencies()

    click.echo("Python dependencies:")
    all_good = True
    for dep, available in deps.items():
        status = "✓" if available else "✗"
        click.echo(f"  {status} {dep}")
        if not available:
            all_good = False

    click.echo()

    # Check R and R packages
    try:
        from .utils import RInterface

        r_interface = RInterface()

        click.echo("R environment:")

        # Check R availability
        if r_interface.check_r_available():
            click.echo("  ✓ R is available")

            # Check key R packages
            r_packages = ["diffHic", "DESeq2", "edgeR", "clusterProfiler"]
            package_status = r_interface.check_packages(r_packages)

            for pkg in r_packages:
                status = "✓" if package_status.get(pkg, False) else "✗"
                click.echo(f"  {status} R package: {pkg}")
                if not package_status.get(pkg, False):
                    all_good = False
        else:
            click.echo("  ✗ R is not available")
            all_good = False

    except Exception as e:
        click.echo(f"  ✗ R environment check failed: {e}")
        all_good = False

    click.echo()

    # Check external tools
    click.echo("External tools:")

    # Check Peakachu
    import subprocess

    try:
        result = subprocess.run(
            ["peakachu", "--version"], capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            click.echo("  ✓ Peakachu is available")
        else:
            click.echo("  ✗ Peakachu not working properly")
            all_good = False
    except (subprocess.TimeoutExpired, FileNotFoundError):
        click.echo("  ✗ Peakachu not found")
        all_good = False

    # Check deepTools
    try:
        result = subprocess.run(
            ["computeMatrix", "--version"], capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            click.echo("  ✓ deepTools is available")
        else:
            click.echo("  ✗ deepTools not working properly")
            all_good = False
    except (subprocess.TimeoutExpired, FileNotFoundError):
        click.echo("  ✗ deepTools not found")
        all_good = False

    click.echo()

    if all_good:
        click.echo("✓ Environment check passed!")
    else:
        click.echo("✗ Environment check failed. Please install missing dependencies.")
        click.echo()
        click.echo("Installation suggestions:")
        click.echo("  - Python packages: pip install plasmaflow[all]")
        click.echo("  - R packages: See documentation for Bioconductor installation")
        click.echo("  - Peakachu: conda install -c bioconda peakachu")
        click.echo("  - deepTools: conda install -c bioconda deeptools")
        sys.exit(1)


# Additional utility commands
@main.group()
def utils():
    """Utility commands"""
    pass


@utils.command()
@click.argument("bedpe_file", type=click.Path(exists=True))
def validate_loops(bedpe_file):
    """Validate a BEDPE loops file"""

    try:
        from .loop_calling import validate_loops

        click.echo(f"Validating loops file: {bedpe_file}")

        quality_metrics = validate_loops([bedpe_file])

        for sample_name, metrics in quality_metrics.items():
            click.echo(f"Sample: {sample_name}")
            click.echo(f"  Total loops: {metrics.num_loops}")
            click.echo(f"  Intra-chromosomal: {metrics.intra_chromosomal}")
            click.echo(f"  Inter-chromosomal: {metrics.inter_chromosomal}")
            click.echo(f"  Mean distance: {metrics.mean_distance:.0f} bp")
            click.echo(f"  Valid coordinates: {metrics.valid_coordinates}")

    except Exception as e:
        click.echo(f"Validation failed: {e}", err=True)
        sys.exit(1)


@utils.command()
@click.argument("input_dir", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
def create_sample_sheet(input_dir, output_file):
    """Create a sample sheet from a directory of files"""

    import pandas as pd

    input_path = Path(input_dir)

    # Find cool and bedpe files
    cool_files = list(input_path.glob("*.cool"))
    bedpe_files = list(input_path.glob("*.bedpe"))

    samples = []

    for cool_file in cool_files:
        sample_name = cool_file.stem.split(".")[0]  # Remove .10000_balanced.cool

        # Look for corresponding BEDPE file
        bedpe_file = None
        for bedpe in bedpe_files:
            if sample_name in bedpe.stem:
                bedpe_file = bedpe
                break

        samples.append(
            {
                "sample_name": sample_name,
                "cool_file": str(cool_file.absolute()),
                "bedpe_file": str(bedpe_file.absolute()) if bedpe_file else "",
                "cell_type": sample_name,  # Adjust as needed
                "read_count": "",  # To be filled manually
            }
        )

    # Create DataFrame and save
    samples_df = pd.DataFrame(samples)
    samples_df.to_csv(output_file, index=False)

    click.echo(f"Sample sheet created: {output_file}")
    click.echo(f"Found {len(samples)} samples")
    click.echo("Please edit the file to add read counts and correct cell types")


if __name__ == "__main__":
    main()

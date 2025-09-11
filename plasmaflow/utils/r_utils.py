"""
R integration utilities for PlasmaFlow
"""

import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

logger = logging.getLogger(__name__)


class RInterface:
    """Interface for running R code from Python"""

    def __init__(self, r_config: Optional[Dict[str, Any]] = None):
        """
        Initialize R interface

        Args:
            r_config: R configuration dictionary
        """
        self.r_config = r_config or {}
        self.r_home = self.r_config.get("r_home")

        # Set R_HOME if specified
        if self.r_home:
            os.environ["R_HOME"] = self.r_home

        self._check_r_available()

    def _check_r_available(self) -> bool:
        """Check if R is available"""
        try:
            result = subprocess.run(
                ["R", "--version"], capture_output=True, text=True, timeout=10
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            logger.warning("R not found in PATH")
            return False

    def check_r_available(self) -> bool:
        """Public method to check R availability"""
        return self._check_r_available()

    def check_packages(self, packages: List[str]) -> Dict[str, bool]:
        """
        Check if R packages are installed

        Args:
            packages: List of package names to check

        Returns:
            Dictionary of package_name -> installed status
        """
        if not self._check_r_available():
            return {pkg: False for pkg in packages}

        # Create R script to check packages
        check_script = f"""
        packages <- c({', '.join([f'"{pkg}"' for pkg in packages])})
        installed <- sapply(packages, function(pkg) {{
            suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE))
        }})
        cat(paste(packages, installed, sep=':', collapse='\\n'))
        """

        try:
            result = subprocess.run(
                ["R", "--slave", "-e", check_script],
                capture_output=True,
                text=True,
                timeout=30,
            )

            if result.returncode == 0:
                package_status = {}
                for line in result.stdout.strip().split("\n"):
                    if ":" in line:
                        pkg, status = line.split(":", 1)
                        package_status[pkg] = status.strip().lower() == "true"
                return package_status
            else:
                logger.error(f"Error checking R packages: {result.stderr}")
                return {pkg: False for pkg in packages}

        except subprocess.TimeoutExpired:
            logger.error("Timeout checking R packages")
            return {pkg: False for pkg in packages}
        except Exception as e:
            logger.error(f"Error checking R packages: {e}")
            return {pkg: False for pkg in packages}

    def install_packages(self, packages: List[str]) -> List[str]:
        """
        Install R packages

        Args:
            packages: List of package names to install

        Returns:
            List of packages that failed to install
        """
        if not self._check_r_available():
            logger.error("R not available for package installation")
            return packages

        failed_packages = []
        bioc_packages = []
        cran_packages = []

        # Classify packages as Bioconductor or CRAN
        known_bioc = {
            "diffHic",
            "InteractionSet",
            "GenomicRanges",
            "DESeq2",
            "edgeR",
            "clusterProfiler",
            "msigdbr",
            "Biostrings",
            "rtracklayer",
        }

        for pkg in packages:
            if pkg in known_bioc:
                bioc_packages.append(pkg)
            else:
                cran_packages.append(pkg)

        # Install CRAN packages
        if cran_packages:
            cran_script = f"""
            repos <- c({', '.join([f'"{repo}"' for repo in self.r_config.get('cran_repos', ['https://cloud.r-project.org'])])})
            install.packages(c({', '.join([f'"{pkg}"' for pkg in cran_packages])}), repos=repos)
            """

            try:
                result = subprocess.run(
                    ["R", "--slave", "-e", cran_script],
                    capture_output=True,
                    text=True,
                    timeout=300,  # 5 minutes
                )

                if result.returncode != 0:
                    logger.error(f"Error installing CRAN packages: {result.stderr}")
                    failed_packages.extend(cran_packages)

            except subprocess.TimeoutExpired:
                logger.error("Timeout installing CRAN packages")
                failed_packages.extend(cran_packages)
            except Exception as e:
                logger.error(f"Error installing CRAN packages: {e}")
                failed_packages.extend(cran_packages)

        # Install Bioconductor packages
        if bioc_packages:
            bioc_script = f"""
            if (!requireNamespace("BiocManager", quietly = TRUE)) {{
                install.packages("BiocManager", repos="https://cloud.r-project.org")
            }}
            BiocManager::install(c({', '.join([f'"{pkg}"' for pkg in bioc_packages])}))
            """

            try:
                result = subprocess.run(
                    ["R", "--slave", "-e", bioc_script],
                    capture_output=True,
                    text=True,
                    timeout=600,  # 10 minutes
                )

                if result.returncode != 0:
                    logger.error(
                        f"Error installing Bioconductor packages: {result.stderr}"
                    )
                    failed_packages.extend(bioc_packages)

            except subprocess.TimeoutExpired:
                logger.error("Timeout installing Bioconductor packages")
                failed_packages.extend(bioc_packages)
            except Exception as e:
                logger.error(f"Error installing Bioconductor packages: {e}")
                failed_packages.extend(bioc_packages)

        return failed_packages

    def run_script(
        self, r_code: str, working_dir: Optional[Union[str, Path]] = None
    ) -> Dict[str, Any]:
        """
        Run R script

        Args:
            r_code: R code to execute
            working_dir: Working directory for R script

        Returns:
            Dictionary with execution results
        """
        if not self._check_r_available():
            return {"success": False, "error": "R not available", "output": None}

        # Create temporary script file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as f:
            f.write(r_code)
            script_path = f.name

        try:
            # Set working directory
            if working_dir is None:
                working_dir = tempfile.mkdtemp()

            working_dir = Path(working_dir)
            working_dir.mkdir(parents=True, exist_ok=True)

            # Run R script
            result = subprocess.run(
                ["R", "--slave", "--no-restore", "--no-save", "-f", script_path],
                cwd=str(working_dir),
                capture_output=True,
                text=True,
                timeout=3600,  # 1 hour timeout
            )

            return {
                "success": result.returncode == 0,
                "output": result.stdout,
                "error": result.stderr if result.returncode != 0 else None,
                "working_dir": str(working_dir),
                "script_path": script_path,
            }

        except subprocess.TimeoutExpired:
            return {
                "success": False,
                "error": "R script execution timed out",
                "output": None,
            }
        except Exception as e:
            return {"success": False, "error": str(e), "output": None}
        finally:
            # Clean up script file
            try:
                os.unlink(script_path)
            except:
                pass


def check_r_packages(packages: List[str]) -> Dict[str, bool]:
    """
    Check if R packages are available

    Args:
        packages: List of package names

    Returns:
        Dictionary of package availability
    """
    r_interface = RInterface()
    return r_interface.check_packages(packages)


def install_r_packages(packages: List[str]) -> List[str]:
    """
    Install R packages

    Args:
        packages: List of package names to install

    Returns:
        List of packages that failed to install
    """
    r_interface = RInterface()
    return r_interface.install_packages(packages)

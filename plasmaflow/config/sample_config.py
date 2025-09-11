"""
Sample and cell type configuration for PlasmaFlow
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


@dataclass
class CellTypeConfig:
    """Configuration for a specific cell type"""

    name: str
    display_name: str
    color: str = "#000000"
    marker_genes: List[str] = field(default_factory=list)
    description: str = ""

    # Hi-C specific parameters
    read_count: Optional[int] = None
    peakachu_model: Optional[str] = None

    # File paths
    cool_file: Optional[str] = None
    loops_file: Optional[str] = None
    expression_file: Optional[str] = None

    def validate(self) -> List[str]:
        """Validate cell type configuration"""
        issues = []

        if not self.name:
            issues.append("Cell type name cannot be empty")

        if self.cool_file and not Path(self.cool_file).exists():
            issues.append(f"Cool file does not exist: {self.cool_file}")

        if self.loops_file and not Path(self.loops_file).exists():
            issues.append(f"Loops file does not exist: {self.loops_file}")

        if self.expression_file and not Path(self.expression_file).exists():
            issues.append(f"Expression file does not exist: {self.expression_file}")

        return issues


@dataclass
class SampleConfig:
    """Configuration for sample management"""

    cell_types: Dict[str, CellTypeConfig] = field(default_factory=dict)

    def add_cell_type(
        self,
        name: str,
        display_name: str,
        read_count: Optional[int] = None,
        peakachu_model: Optional[str] = None,
        **kwargs,
    ) -> None:
        """Add a cell type to the configuration"""

        self.cell_types[name] = CellTypeConfig(
            name=name,
            display_name=display_name,
            read_count=read_count,
            peakachu_model=peakachu_model,
            **kwargs,
        )

    def get_cell_type(self, name: str) -> Optional[CellTypeConfig]:
        """Get cell type configuration by name"""
        return self.cell_types.get(name)

    def list_cell_types(self) -> List[str]:
        """Get list of all cell type names"""
        return list(self.cell_types.keys())

    def validate_all(self) -> Dict[str, List[str]]:
        """Validate all cell type configurations"""
        all_issues = {}

        for name, cell_type in self.cell_types.items():
            issues = cell_type.validate()
            if issues:
                all_issues[name] = issues

        return all_issues


def get_default_sample_config() -> SampleConfig:
    """Get default sample configuration for plasma cell analysis"""

    config = SampleConfig()

    # Memory B cells (memB)
    config.add_cell_type(
        name="memB",
        display_name="Memory B cells",
        read_count=113922292,
        peakachu_model="100million",
        color="#1f77b4",
        marker_genes=["CD19", "CD20", "CD27", "IGD"],
        description="Antigen-experienced memory B cells",
    )

    # Pre-plasmablasts (prePB)
    config.add_cell_type(
        name="prePB",
        display_name="Pre-plasmablasts",
        read_count=159514688,
        peakachu_model="150million",
        color="#ff7f0e",
        marker_genes=["CD19", "CD38", "CD138", "PRDM1"],
        description="Activated B cells transitioning to plasmablasts",
    )

    # Plasmablasts (PB)
    config.add_cell_type(
        name="PB",
        display_name="Plasmablasts",
        read_count=158184861,
        peakachu_model="150million",
        color="#2ca02c",
        marker_genes=["CD19", "CD38", "CD138", "IRF4", "XBP1"],
        description="Proliferating antibody-secreting cells",
    )

    # Plasma cells (PC)
    config.add_cell_type(
        name="PC",
        display_name="Plasma cells",
        read_count=91590221,
        peakachu_model="100million",
        color="#d62728",
        marker_genes=["CD38", "CD138", "TNFRSF17", "SDC1"],
        description="Terminally differentiated antibody-secreting cells",
    )

    return config

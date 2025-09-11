"""
Model management for Peakachu loop calling
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


class ModelManager:
    """Manager for Peakachu model files"""

    def __init__(self, model_dir: Optional[str] = None):
        self.model_dir = Path(model_dir) if model_dir else None
        self._validate_model_dir()

    def _validate_model_dir(self) -> None:
        """Validate that model directory exists"""
        if self.model_dir and not self.model_dir.exists():
            logger.warning(f"Model directory does not exist: {self.model_dir}")
            self.model_dir = None

    def get_model_path(self, model_name: str) -> Optional[Path]:
        """
        Get full path to a Peakachu model file

        Args:
            model_name: Name of the model (e.g., "100million", "150million")

        Returns:
            Path to model file or None if not found
        """
        if not self.model_dir:
            logger.error("Model directory not set")
            return None

        # Try different possible model file patterns
        model_patterns = [
            f"high-confidence.{model_name}.10kb.w6.pkl",
            f"{model_name}.10kb.w6.pkl",
            f"{model_name}.pkl",
            f"high-confidence.{model_name}.pkl",
        ]

        for pattern in model_patterns:
            model_path = self.model_dir / pattern
            if model_path.exists():
                logger.debug(f"Found model: {model_path}")
                return model_path

        logger.error(f"Model not found: {model_name} in {self.model_dir}")
        return None

    def list_available_models(self) -> List[str]:
        """List all available model names"""
        if not self.model_dir or not self.model_dir.exists():
            return []

        models = []
        for file in self.model_dir.glob("*.pkl"):
            # Extract model name from filename
            stem = file.stem
            if "high-confidence." in stem:
                # high-confidence.100million.10kb.w6.pkl -> 100million
                parts = stem.split(".")
                if len(parts) >= 2:
                    models.append(parts[1])
            else:
                # 100million.10kb.w6.pkl -> 100million
                parts = stem.split(".")
                if parts:
                    models.append(parts[0])

        return list(set(models))  # Remove duplicates

    def validate_model(self, model_name: str) -> bool:
        """Check if a model exists and is valid"""
        model_path = self.get_model_path(model_name)
        return model_path is not None and model_path.exists()


def get_recommended_model(read_count: int) -> str:
    """
    Get recommended Peakachu model based on read count

    Args:
        read_count: Number of Hi-C reads

    Returns:
        Recommended model name
    """
    # Thresholds based on Peakachu documentation
    if read_count < 75_000_000:
        return "50million"
    elif read_count < 125_000_000:
        return "100million"
    elif read_count < 175_000_000:
        return "150million"
    elif read_count < 225_000_000:
        return "200million"
    else:
        return "300million"


def get_model_recommendations() -> Dict[str, Dict[str, int]]:
    """
    Get all model recommendations with read count ranges

    Returns:
        Dictionary mapping model names to their recommended read count ranges
    """
    return {
        "50million": {"min_reads": 0, "max_reads": 75_000_000},
        "100million": {"min_reads": 75_000_000, "max_reads": 125_000_000},
        "150million": {"min_reads": 125_000_000, "max_reads": 175_000_000},
        "200million": {"min_reads": 175_000_000, "max_reads": 225_000_000},
        "300million": {"min_reads": 225_000_000, "max_reads": float("inf")},
    }


def validate_read_count_model_match(
    read_count: int, model_name: str
) -> Tuple[bool, str]:
    """
    Validate if read count matches the selected model

    Args:
        read_count: Number of Hi-C reads
        model_name: Selected model name

    Returns:
        Tuple of (is_valid, message)
    """
    recommended_model = get_recommended_model(read_count)

    if model_name == recommended_model:
        return True, f"Model {model_name} is optimal for {read_count:,} reads"

    recommendations = get_model_recommendations()
    if model_name in recommendations:
        model_range = recommendations[model_name]
        if model_range["min_reads"] <= read_count < model_range["max_reads"]:
            return True, f"Model {model_name} is acceptable for {read_count:,} reads"
        else:
            return False, (
                f"Model {model_name} is not suitable for {read_count:,} reads. "
                f"Recommended: {recommended_model}"
            )
    else:
        return (
            False,
            f"Unknown model: {model_name}. Available models: {list(recommendations.keys())}",
        )

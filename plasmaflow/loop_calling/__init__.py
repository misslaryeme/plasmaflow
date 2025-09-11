"""
Loop calling module for PlasmaFlow

This module provides functionality for calling chromatin loops from Hi-C data
using the Peakachu framework. It includes sample processing, model management,
and result validation.
"""

from .caller import LoopCaller, PeakachuRunner
from .models import ModelManager, get_recommended_model
from .utils import merge_loop_results, process_sample_batch
from .validation import LoopQualityMetrics, validate_loops

__all__ = [
    "LoopCaller",
    "PeakachuRunner",
    "ModelManager",
    "get_recommended_model",
    "validate_loops",
    "LoopQualityMetrics",
    "process_sample_batch",
    "merge_loop_results",
]

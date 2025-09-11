"""
Caching functionality for APA analysis results
"""

import hashlib
import json
import logging
import pickle
from pathlib import Path
from typing import Any, Dict, Optional, Union

logger = logging.getLogger(__name__)


class APACache:
    """Cache manager for APA results"""

    def __init__(self, cache_dir: Union[str, Path]):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _get_params_hash(self, params: Dict[str, Any]) -> str:
        """Generate hash of parameters for cache validation"""
        # Sort params for consistent hashing
        sorted_params = json.dumps(params, sort_keys=True)
        return hashlib.md5(sorted_params.encode()).hexdigest()[:8]

    def get_cache_path(self, sample_name: str, params: Dict[str, Any]) -> Path:
        """Get cache file path for given sample and parameters"""
        params_hash = self._get_params_hash(params)
        return self.cache_dir / f"{sample_name}_{params_hash}_apa.pkl"

    def save_result(self, result: "APAResult", params: Dict[str, Any]) -> None:
        """Save APA result to cache"""
        cache_path = self.get_cache_path(result.sample_name, params)

        cache_data = {"result": result, "params": params, "version": "1.0"}

        try:
            with open(cache_path, "wb") as f:
                pickle.dump(cache_data, f)
            logger.debug(f"Cached APA result: {cache_path}")
        except Exception as e:
            logger.warning(f"Failed to cache APA result: {e}")

    def load_result(
        self, sample_name: str, params: Dict[str, Any]
    ) -> Optional["APAResult"]:
        """Load APA result from cache if parameters match"""
        cache_path = self.get_cache_path(sample_name, params)

        if not cache_path.exists():
            return None

        try:
            with open(cache_path, "rb") as f:
                cache_data = pickle.load(f)

            cached_params = cache_data.get("params", {})

            # Validate parameters match
            if self._params_match(cached_params, params):
                logger.debug(f"Loaded cached APA result: {cache_path}")
                return cache_data["result"]
            else:
                logger.debug(f"Cache parameters mismatch for {sample_name}")
                return None

        except Exception as e:
            logger.warning(f"Failed to load cached result: {e}")
            return None

    def _params_match(
        self, cached_params: Dict[str, Any], current_params: Dict[str, Any]
    ) -> bool:
        """Check if cached parameters match current parameters"""
        for key, value in current_params.items():
            if cached_params.get(key) != value:
                return False
        return True

    def clear_cache(self, sample_name: Optional[str] = None) -> int:
        """Clear cache files. If sample_name provided, clear only for that sample"""
        cleared = 0

        if sample_name:
            pattern = f"{sample_name}_*_apa.pkl"
        else:
            pattern = "*_apa.pkl"

        for cache_file in self.cache_dir.glob(pattern):
            try:
                cache_file.unlink()
                cleared += 1
            except Exception as e:
                logger.warning(f"Failed to delete cache file {cache_file}: {e}")

        logger.info(f"Cleared {cleared} cache files")
        return cleared


def save_apa_cache(
    result: "APAResult", cache_file: Union[str, Path], params: Dict[str, Any]
) -> None:
    """
    Save APA result to cache file

    Args:
        result: APAResult object to cache
        cache_file: Path to cache file
        params: Parameters used for analysis
    """
    cache_data = {"result": result, "params": params, "version": "1.0"}

    try:
        cache_path = Path(cache_file)
        cache_path.parent.mkdir(parents=True, exist_ok=True)

        with open(cache_path, "wb") as f:
            pickle.dump(cache_data, f)

        logger.debug(f"Saved APA cache: {cache_path}")

    except Exception as e:
        logger.warning(f"Failed to save APA cache: {e}")
        raise


def load_cached_apa(
    cache_file: Union[str, Path], params: Dict[str, Any]
) -> Optional["APAResult"]:
    """
    Load cached APA result if parameters match

    Args:
        cache_file: Path to cache file
        params: Current analysis parameters

    Returns:
        APAResult if cache valid, None otherwise
    """
    cache_path = Path(cache_file)

    if not cache_path.exists():
        return None

    try:
        with open(cache_path, "rb") as f:
            cache_data = pickle.load(f)

        cached_params = cache_data.get("params", {})

        # Check if parameters match
        if params_match(cached_params, params):
            logger.debug(f"Using cached APA result: {cache_path}")
            return cache_data["result"]
        else:
            logger.debug("Cache parameters don't match, recomputing")
            return None

    except Exception as e:
        logger.warning(f"Failed to load APA cache: {e}")
        return None


def params_match(cached_params: Dict[str, Any], current_params: Dict[str, Any]) -> bool:
    """Check if cached parameters match current parameters"""
    for key, value in current_params.items():
        if cached_params.get(key) != value:
            logger.debug(
                f"Parameter mismatch: {key} cached={cached_params.get(key)} vs current={value}"
            )
            return False
    return True


def get_cache_info(cache_dir: Union[str, Path]) -> Dict[str, Any]:
    """Get information about cache directory"""
    cache_path = Path(cache_dir)

    if not cache_path.exists():
        return {"exists": False, "files": 0, "total_size": 0}

    cache_files = list(cache_path.glob("*_apa.pkl"))

    total_size = sum(f.stat().st_size for f in cache_files)

    return {
        "exists": True,
        "files": len(cache_files),
        "total_size": total_size,
        "total_size_mb": total_size / (1024 * 1024),
        "cache_dir": str(cache_path),
    }

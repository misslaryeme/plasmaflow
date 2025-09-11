"""
Logging utilities for PlasmaFlow
"""

import logging
import sys
from pathlib import Path
from typing import Optional, Union

import colorlog


def setup_logging(
    level: Union[str, int] = logging.INFO,
    log_file: Optional[Union[str, Path]] = None,
    format_string: Optional[str] = None,
    use_colors: bool = True,
) -> logging.Logger:
    """
    Set up logging for PlasmaFlow

    Args:
        level: Logging level (INFO, DEBUG, etc.)
        log_file: Optional file to write logs to
        format_string: Custom format string
        use_colors: Whether to use colored output for console

    Returns:
        Configured logger
    """
    # Convert string level to logging constant
    if isinstance(level, str):
        level = getattr(logging, level.upper())

    # Default format
    if format_string is None:
        format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    # Clear existing handlers
    root_logger = logging.getLogger()
    root_logger.handlers.clear()

    # Console handler with colors
    if use_colors:
        console_handler = colorlog.StreamHandler(sys.stdout)
        console_formatter = colorlog.ColoredFormatter(
            "%(log_color)s%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            log_colors={
                "DEBUG": "cyan",
                "INFO": "green",
                "WARNING": "yellow",
                "ERROR": "red",
                "CRITICAL": "red,bg_white",
            },
        )
    else:
        console_handler = logging.StreamHandler(sys.stdout)
        console_formatter = logging.Formatter(
            format_string, datefmt="%Y-%m-%d %H:%M:%S"
        )

    console_handler.setFormatter(console_formatter)
    console_handler.setLevel(level)
    root_logger.addHandler(console_handler)

    # File handler if requested
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_path)
        file_formatter = logging.Formatter(format_string, datefmt="%Y-%m-%d %H:%M:%S")
        file_handler.setFormatter(file_formatter)
        file_handler.setLevel(level)
        root_logger.addHandler(file_handler)

    root_logger.setLevel(level)

    # Set specific loggers
    plasmaflow_logger = logging.getLogger("plasmaflow")
    plasmaflow_logger.setLevel(level)

    return plasmaflow_logger


def get_logger(name: str) -> logging.Logger:
    """Get a logger for a specific module"""
    return logging.getLogger(f"plasmaflow.{name}")


class LoggerMixin:
    """Mixin class to add logging capability to any class"""

    @property
    def logger(self) -> logging.Logger:
        """Get logger for this class"""
        return get_logger(self.__class__.__name__.lower())


def log_function_call(func):
    """Decorator to log function calls"""

    def wrapper(*args, **kwargs):
        logger = get_logger(func.__module__)
        logger.debug(f"Calling {func.__name__} with args={args}, kwargs={kwargs}")
        try:
            result = func(*args, **kwargs)
            logger.debug(f"{func.__name__} completed successfully")
            return result
        except Exception as e:
            logger.error(f"{func.__name__} failed: {e}")
            raise

    return wrapper


def log_execution_time(func):
    """Decorator to log execution time of functions"""
    import time

    def wrapper(*args, **kwargs):
        logger = get_logger(func.__module__)
        start_time = time.time()

        try:
            result = func(*args, **kwargs)
            execution_time = time.time() - start_time
            logger.info(f"{func.__name__} completed in {execution_time:.2f} seconds")
            return result
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(
                f"{func.__name__} failed after {execution_time:.2f} seconds: {e}"
            )
            raise

    return wrapper

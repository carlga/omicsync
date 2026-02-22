"""
Shared utilities for Snakemake workflows.

This module provides common functions used across all osync Snakemake workflows
to ensure consistency, reduce code duplication, and improve maintainability.
"""

import os
import sys
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional, Union


def get_timestamp(format_str: str = "%Y%m%d") -> str:
    """
    Get standardized timestamp for log files and output directories.

    Args:
        format_str: strftime format string (default: "%Y%m%d")

    Returns:
        Formatted timestamp string

    Example:
        >>> timestamp = get_timestamp()
        >>> print(timestamp)
        '20260218'
    """
    return datetime.now().strftime(format_str)


def setup_workflow_paths(workflow_obj) -> tuple:
    """
    Setup and return workflow directory paths.

    Args:
        workflow_obj: Snakemake workflow object

    Returns:
        Tuple of (workflow_dir, parent_dir) as Path objects

    Example:
        >>> workflow_dir, parent_dir = setup_workflow_paths(workflow)
        >>> sys.path.insert(0, str(parent_dir))
    """
    workflow_dir = Path(workflow_obj.snakefile).parent.absolute()
    parent_dir = workflow_dir.parent
    return workflow_dir, parent_dir


def validate_config_keys(config: Dict[str, Any],
                         module: str,
                         required_keys: List[str],
                         optional_keys: Optional[List[str]] = None) -> None:
    """
    Validate required configuration keys exist for a module.

    Args:
        config: Full configuration dictionary
        module: Module name (e.g., 'runqc', 'init')
        required_keys: List of required configuration keys
        optional_keys: List of optional keys (for documentation)

    Raises:
        ValueError: If required keys are missing
        KeyError: If module not in config

    Example:
        >>> validate_config_keys(config, 'runqc', ['indir', 'outdir', 'jobs'])
    """
    if module not in config:
        raise KeyError(
            f"Module '{module}' not found in configuration. "
            f"Available modules: {list(config.keys())}"
        )

    module_config = config[module]
    missing = [k for k in required_keys if k not in module_config or module_config[k] is None]

    if missing:
        raise ValueError(
            f"[{module.upper()}] Missing required configuration keys: {missing}\n"
            f"Required keys: {required_keys}"
        )


def extract_threads_from_opt(opt_string: str, 
                             patterns: Optional[List[str]] = None) -> Optional[int]:
    """
    Extract thread count from tool option string.

    Args:
        opt_string: Tool options string (e.g., "--threads 8 --other-opt")
        patterns: List of regex patterns to match (default: common thread patterns)

    Returns:
        Extracted thread count or None if not found

    Example:
        >>> threads = extract_threads_from_opt("--threads 8 -q")
        >>> print(threads)
        8
    """
    if patterns is None:
        patterns = [
            r'(-t|--threads?)\s+(\d+)',
            r'(-p|--processors?)\s+(\d+)',
            r'(-j|--jobs?)\s+(\d+)',
            r'(-w|--workers?)\s+(\d+)',
        ]

    for pattern in patterns:
        match = re.search(pattern, opt_string)
        if match:
            return int(match.group(2))

    return None


def allocate_threads(total: int,
                    tasks: Dict[str, Union[int, float]],
                    min_per_task: int = 1) -> Dict[str, int]:
    """
    Allocate threads among multiple tasks with priorities.

    Args:
        total: Total available threads
        tasks: Dict mapping task names to priorities (higher = more threads)
               Use int for fixed allocation, float (0-1) for proportional
        min_per_task: Minimum threads per task

    Returns:
        Dict mapping task names to allocated thread counts

    Example:
        >>> tasks = {'align': 0.6, 'sort': 0.3, 'index': 0.1}
        >>> allocation = allocate_threads(10, tasks)
        >>> print(allocation)
        {'align': 6, 'sort': 3, 'index': 1}
    """
    allocation = {}
    remaining = total

    # First pass: handle fixed allocations
    fixed_tasks = {k: v for k, v in tasks.items() if isinstance(v, int)}
    for task, count in fixed_tasks.items():
        allocated = min(count, remaining)
        allocation[task] = max(allocated, min_per_task)
        remaining -= allocation[task]

    # Second pass: handle proportional allocations
    prop_tasks = {k: v for k, v in tasks.items() if isinstance(v, float)}
    for task, proportion in prop_tasks.items():
        allocated = max(min_per_task, int(total * proportion))
        allocated = min(allocated, remaining)
        allocation[task] = allocated
        remaining -= allocated

    # Distribute remaining threads to largest allocation
    if remaining > 0 and allocation:
        max_task = max(allocation, key=allocation.get)
        allocation[max_task] += remaining

    return allocation


def create_log_path(timestamp: str,
                   rule_name: str,
                   wildcards: Optional[Dict[str, str]] = None,
                   base_dir: str = "logs") -> str:
    """
    Create standardized log file path.

    Args:
        timestamp: Timestamp string (from get_timestamp())
        rule_name: Name of the Snakemake rule
        wildcards: Dictionary of wildcards (optional)
        base_dir: Base directory for logs

    Returns:
        Formatted log file path

    Example:
        >>> path = create_log_path("20260218", "align", {"sample": "S1"})
        >>> print(path)
        'logs/20260218_align_S1.log'
    """
    if wildcards:
        wildcard_str = "_".join(f"{v}" for v in wildcards.values() if v)
        return f"{base_dir}/{timestamp}_{rule_name}_{wildcard_str}.log"
    return f"{base_dir}/{timestamp}_{rule_name}.log"


def format_shell_log_msg(module: str, 
                        message: str,
                        level: str = "INFO") -> str:
    """
    Format log message for shell commands.

    Args:
        module: Module name (e.g., 'RUNQC', 'INIT')
        message: Log message
        level: Log level (INFO, WARNING, ERROR)

    Returns:
        Formatted log message

    Example:
        >>> msg = format_shell_log_msg("RUNQC", "Processing sample S1")
        >>> print(msg)
        '[RUNQC] INFO: Processing sample S1'
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return f"[{module}] {level}: {message} ({timestamp})"


def check_file_exists(filepath: Union[str, Path],
                     error_msg: Optional[str] = None) -> bool:
    """
    Check if file exists and raise informative error if not.

    Args:
        filepath: Path to file
        error_msg: Custom error message (optional)

    Returns:
        True if file exists

    Raises:
        FileNotFoundError: If file doesn't exist

    Example:
        >>> check_file_exists("reference.fa", "Reference genome not found")
    """
    filepath = Path(filepath)
    if not filepath.exists():
        if error_msg:
            raise FileNotFoundError(error_msg)
        raise FileNotFoundError(f"Required file not found: {filepath}")
    return True


def get_available_memory_gb() -> float:
    """
    Get available system memory in GB.

    Returns:
        Available memory in GB

    Example:
        >>> mem_gb = get_available_memory_gb()
        >>> print(f"Available memory: {mem_gb:.1f} GB")
    """
    try:
        import psutil
        return psutil.virtual_memory().available / (1024 ** 3)
    except ImportError:
        # Fallback if psutil not available
        return 8.0  # Conservative default


def sanitize_filename(filename: str,
                     replacement: str = "_") -> str:
    """
    Sanitize filename by replacing invalid characters.

    Args:
        filename: Original filename
        replacement: Character to replace invalid chars with

    Returns:
        Sanitized filename

    Example:
        >>> clean = sanitize_filename("sample:1_data*?.txt")
        >>> print(clean)
        'sample_1_data__.txt'
    """
    # Remove/replace invalid characters
    invalid_chars = r'[<>:"/\|?*]'
    return re.sub(invalid_chars, replacement, filename)


# Module-level constants
DEFAULT_THREAD_PATTERNS = [
    r'(-t|--threads?)\s+(\d+)',
    r'(-p|--processors?)\s+(\d+)',
    r'(-j|--jobs?)\s+(\d+)',
]

LOG_LEVELS = {
    'DEBUG': 'DEBUG',
    'INFO': 'INFO',
    'WARNING': 'WARNING',
    'ERROR': 'ERROR',
    'CRITICAL': 'CRITICAL',
}

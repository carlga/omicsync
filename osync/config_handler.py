"""Configuration file handling for osync.

This module provides functions to create, load, update, and clear
configuration files for the osync tool.
"""

import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from ruamel.yaml import YAML

from .logger import logger


# Constants
AVAILABLE_MODULES = ['init', 'getref', 'getfq', 'runqc', 'run_rnaseq', 'run_atacseq', 'all']
RECORDS_KEY = 'records'


def load_configfile(configfile: Optional[str]) -> Dict[str, Any]:
    """Load a configuration file with comments preserved.

    Args:
        configfile: Path to the configuration file

    Returns:
        Dictionary containing the configuration

    Raises:
        SystemExit: If file is missing or cannot be read
    """
    if configfile is None:
        logger.error("Missing configuration file! Please provide path to file.")
        sys.exit(1)

    config_path = Path(configfile).resolve()

    if not config_path.exists():
        logger.error(f"Configuration file not found: '{config_path}'")
        sys.exit(1)

    try:
        with open(config_path, 'r') as f:
            config = YAML().load(f)
        return config
    except Exception as e:
        logger.error(f"Failed to load configuration file: {e}")
        sys.exit(1)


def make_configfile(outfile: str) -> None:
    """Create a new configuration file from the default template.

    Args:
        outfile: Path where the new configuration file will be created

    Raises:
        SystemExit: If file already exists or cannot be created
    """
    # Load default configuration template
    default_config_path = Path(__file__).parent / "default_config.yaml"
    config = load_configfile(str(default_config_path))

    # Prepare output path
    out_path = Path(outfile).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists():
        logger.error(f"Configuration file already exists: '{out_path}'")
        logger.error("Please use a different path or delete the existing file.")
        sys.exit(1)

    # Update metadata
    config[RECORDS_KEY]["generated"] = datetime.now().strftime("%Y-%m-%d %H:%M")
    config[RECORDS_KEY]["filepath"] = str(out_path)

    try:
        with open(out_path, 'w') as f:
            YAML().dump(config, f)
        logger.info(f"Configuration file created: '{out_path}'")
    except Exception as e:
        logger.error(f"Failed to create configuration file: {e}")
        sys.exit(1)


def update_configfile(configfile: str, new_config: Dict[str, Any]) -> None:
    """Update an existing configuration file with new parameters.

    Args:
        configfile: Path to the configuration file to update
        new_config: Dictionary containing the updated configuration

    Raises:
        SystemExit: If file doesn't exist or cannot be written
    """
    config_path = Path(configfile).resolve()

    if not config_path.exists():
        logger.error(f"Cannot update non-existent file: '{config_path}'")
        sys.exit(1)

    # Update metadata
    new_config[RECORDS_KEY]["last_updated"] = datetime.now().strftime("%Y-%m-%d %H:%M")
    new_config[RECORDS_KEY]["filepath"] = str(config_path)

    try:
        with open(config_path, 'w') as f:
            YAML().dump(new_config, f)
        logger.debug(f"Configuration file updated: '{config_path}'")
    except Exception as e:
        logger.error(f"Failed to update configuration file: {e}")
        sys.exit(1)


def clear_configfile(configfile: str, module: Optional[str]) -> None:
    """Clear parameters from a configuration file for specified module(s).

    Args:
        configfile: Path to the configuration file
        module: Module name to clear ('init', 'getref', 'getfq', 'runqc', 
                'run_rnaseq', 'run_atacseq', or 'all')

    Raises:
        SystemExit: If module is invalid or file cannot be updated
    """
    # Validate module selection
    if module not in AVAILABLE_MODULES:
        if module is None:
            logger.error("Missing module information! Please specify which module to clear.")
        else:
            logger.error(
                f"Invalid module '{module}'! "
                f"Choose from: {', '.join(repr(m) for m in AVAILABLE_MODULES)}"
            )
        sys.exit(1)

    # Load and clear configuration
    config = load_configfile(configfile)

    if module == 'all':
        # Clear all modules except records
        for module_name in config.keys():
            if module_name != RECORDS_KEY:
                for param in config[module_name].keys():
                    config[module_name][param] = None
        logger.info("Cleared all module parameters")
    else:
        # Clear specific module
        for param in config[module].keys():
            config[module][param] = None
        logger.info(f"Cleared '{module}' module parameters")

    update_configfile(configfile, config)

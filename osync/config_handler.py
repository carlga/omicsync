import os
import sys
import logging
from ruamel.yaml import YAML
from datetime import datetime

from .logger import logger


def load_configfile(configfile):
    """ Loads configuration file with comments. """

    if configfile is None:
        logger.error("Missing configuration file! Please provide path to file.")
        sys.exit(1)

    configfile = os.path.abspath(configfile)

    if not os.path.exists(configfile):
        logger.error(f"Missing configuration file! '{configfile}' not found.")
        sys.exit(1)
    else:    
        with open(configfile, "r") as f:
            config = YAML().load(f)
    
    return config


def make_configfile(outfile):
    """ Makes file to store tool configuration. """

    configfile_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "default_config.yaml"
    )

    config = load_configfile(configfile_path)

    outfile = os.path.abspath(outfile)
    os.makedirs(os.path.dirname(outfile), exist_ok=True)

    # Update configuration file records
    config["records"]["generated"] = datetime.now().strftime("%Y-%m-%d %H:%M")
    config["records"]["filepath"] = outfile

    if os.path.exists(outfile):
        logger.error(f"Configuration file '{outfile}' exists! Use a different path.")
        sys.exit(1)
    else:
        with open(outfile, "w") as f:
            YAML().dump(config, f)
    

def update_configfile(configfile, new_config):
    """ Updates parameters of configuration file. """

    configfile = os.path.abspath(configfile)

    new_config["records"]["last_updated"] = datetime.now().strftime("%Y-%m-%d %H:%M")
    new_config["records"]["filepath"] = configfile

    if not os.path.exists(configfile):
        logger.error(f"Failed to update configuration! File '{configfile}' does not exist")
        sys.exit(1)
    else:
        with open(configfile, "w") as f:
            YAML().dump(new_config, f)
    

def clear_configfile(configfile, module):
    """ Clears parameters of configuration file. """

    # Check module to clear
    avail_modules = ['init', 'getref', 'getfq', 'runqc', 'run_rnaseq', 'run_atacseq', 'all']

    if module not in avail_modules:
        if module is None:
            logger.error(
                "Missing module information! Please indicate which module to clear."
            )
        else:
            logger.error(
                "Module is not valid! Choose from: 'init', 'getref', 'getfq', "
                "'runqc', 'run_rnaseq', 'run_atacseq', 'all'."
            )
        sys.exit(1)
    
    # Clear module parameters from configuration file
    config = load_configfile(configfile)

    if module == 'all':
        for m in config.keys():
            if m != 'records':
                for param in config[m].keys():
                    config[m][param] = None
    else:
        for param in config[module].keys():
            config[module][param] = None

    update_configfile(configfile, config)

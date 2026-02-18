import os
import sys
import logging

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from .logger import logger

logger = logging.getLogger(__name__)

__version__ = "0.2.1"
__author__ = "Carlos Gallardo"
__email__ = "carlos.gallardo@scilifelab.se"

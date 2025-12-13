"""
This module contains constants for the filesystem paths used by Pathways.
"""

from pathlib import Path

import platformdirs
import yaml


DATA_DIR = Path(__file__).resolve().parent / "data"

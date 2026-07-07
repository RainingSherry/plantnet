#!/usr/bin/env python3
from pathlib import Path
import sys

SCMAES_DIR = Path(__file__).resolve().parents[2]
if str(SCMAES_DIR) not in sys.path:
    sys.path.insert(0, str(SCMAES_DIR))

from common.runner import main


if __name__ == "__main__":
    main("010_longtail_prototype")


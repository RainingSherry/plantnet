#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from benchmark import MECH_DIR, write_outputs


if __name__ == "__main__":
    write_outputs(Path(MECH_DIR / "runs"))


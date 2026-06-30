#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

MECH_DIR = Path(__file__).resolve().parents[2]
if str(MECH_DIR) not in sys.path:
    sys.path.insert(0, str(MECH_DIR))

from common.trainer import VariantConfig, main
from loss import VARIANT_DEFAULTS
from model import build_model


CONFIG = VariantConfig(
    name="001_adaptive_mask",
    method_name="scMAE_plus_adaptive_mask",
    description="Original scMAE backbone with gene-statistic adaptive replacement masking.",
    defaults=VARIANT_DEFAULTS,
)


if __name__ == "__main__":
    main(CONFIG, build_model)


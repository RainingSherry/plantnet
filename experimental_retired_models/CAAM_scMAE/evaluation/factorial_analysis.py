from __future__ import annotations


def interaction_delta(y00: float, y10: float, y01: float, y11: float) -> float:
    return float(y11 - y10 - y01 + y00)


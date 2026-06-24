from pathlib import Path

import yaml

from scripts.run_formal_benchmark import DEFAULT_FORMAL_METHODS


PROJECT_ROOT = Path(__file__).resolve().parents[4]
MANIFEST_PATH = PROJECT_ROOT / "methods" / "method_manifest.yaml"


def _manifest_methods():
    data = yaml.safe_load(MANIFEST_PATH.read_text(encoding="utf-8"))
    return {entry["key"]: entry for entry in data["methods"]}


def test_caam_formal_manifest_entry_is_full_variant_only():
    methods = _manifest_methods()
    entry = methods["caam_scmae"]
    assert entry["path"] == "methods/DeepLearning/CAAM_scMAE/run.py"
    assert entry["default_in_formal"] is False
    assert entry["authenticity"] == "PENDING"
    assert entry["smoke"] == "PASS"
    assert "artifact_manifest.json" in entry["required_artifacts"]
    assert entry["extra_args"] == [
        "--variant",
        "full",
        "--benchmark_mode",
        "true",
        "--method_name",
        "caam_scmae",
    ]


def test_caam_ablation_variants_are_not_formal_methods():
    methods = _manifest_methods()
    forbidden = {"caam_scmae_control", "caam_scmae_axial", "caam_scmae_advmask"}
    assert forbidden.isdisjoint(methods)
    assert "caam_scmae" not in DEFAULT_FORMAL_METHODS

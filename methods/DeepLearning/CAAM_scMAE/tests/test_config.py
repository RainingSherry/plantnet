from methods.DeepLearning.CAAM_scMAE.registry import DEFAULT_CONFIG, VARIANTS, apply_variant


def test_variants_are_isolated():
    assert set(VARIANTS) == {"control", "axial", "advmask", "full"}
    assert apply_variant(DEFAULT_CONFIG, "control")["model"]["encoder_type"] == "mlp"
    assert apply_variant(DEFAULT_CONFIG, "axial")["model"]["mask_selector"] == "random"
    assert apply_variant(DEFAULT_CONFIG, "advmask")["model"]["encoder_type"] == "mlp"
    assert apply_variant(DEFAULT_CONFIG, "full")["model"]["mask_selector"] == "adversarial"


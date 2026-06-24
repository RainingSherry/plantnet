from methods.DeepLearning.CAAM_scMAE.mask_generator.adversarial_mask import AdversarialMaskGenerator
from methods.DeepLearning.CAAM_scMAE.models.axial_encoder import AxialEncoder
from methods.DeepLearning.CAAM_scMAE.models.caam_model import build_student
from methods.DeepLearning.CAAM_scMAE.tests.helpers import toy_assignment, toy_config


def test_model_a_has_axial_no_adversarial():
    cfg = toy_config("axial")
    model = build_student(n_genes=8, config=cfg, assignment=toy_assignment())
    assert any(isinstance(m, AxialEncoder) for m in model.modules())
    assert not any(isinstance(m, AdversarialMaskGenerator) for m in model.modules())


def test_model_b_has_no_axial():
    cfg = toy_config("advmask")
    model = build_student(n_genes=8, config=cfg)
    assert not any(isinstance(m, AxialEncoder) for m in model.modules())


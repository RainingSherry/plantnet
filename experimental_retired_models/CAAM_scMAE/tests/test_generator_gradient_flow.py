import torch

from methods.DeepLearning.CAAM_scMAE.data.dataset import CAAMExpressionDataset
from methods.DeepLearning.CAAM_scMAE.models.caam_model import build_student
from methods.DeepLearning.CAAM_scMAE.trainers.common import CAAMTrainer
from methods.DeepLearning.CAAM_scMAE.tests.helpers import toy_arrays, toy_config, toy_provider


def test_generator_step_updates_only_generator(tmp_path):
    x, batch, lib, zero = toy_arrays()
    cfg = toy_config("advmask")
    student = build_student(n_genes=x.shape[1], config=cfg)
    trainer = CAAMTrainer(
        config=cfg,
        student=student,
        train_dataset=CAAMExpressionDataset(x, batch, lib, zero),
        donor_provider=toy_provider(x, batch, lib, zero),
        full_x=torch.as_tensor(x),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )
    batch_obj = next(iter(trainer._loader()))
    stats = trainer._generator_step(batch_obj, epoch=1)
    assert stats["generator_grad_norm"] > 0
    assert stats["student_grad_norm_during_generator_step"] == 0


import numpy as np
import torch

from methods.DeepLearning.CAAM_scMAE.data.dataset import CAAMExpressionDataset
from methods.DeepLearning.CAAM_scMAE.data.context_selection import select_context_indices
from methods.DeepLearning.CAAM_scMAE.data.gene_modules import build_gene_modules
from methods.DeepLearning.CAAM_scMAE.models.caam_model import build_student
from methods.DeepLearning.CAAM_scMAE.trainers.common import CAAMTrainer
from methods.DeepLearning.CAAM_scMAE.tests.helpers import toy_arrays, toy_assignment, toy_config, toy_provider


def test_gene_modules_and_context_reproducible(tmp_path):
    x, *_ = toy_arrays()
    ids1, _ = build_gene_modules(x, 4, 2, 0, tmp_path / "a")
    ids2, _ = build_gene_modules(x, 4, 2, 0, tmp_path / "b")
    ctx1 = select_context_indices(x, 4, 2, 0, tmp_path / "a")
    ctx2 = select_context_indices(x, 4, 2, 0, tmp_path / "b")
    assert np.array_equal(ids1, ids2)
    assert np.array_equal(ctx1, ctx2)


def _run_full_short_training(save_dir):
    save_dir.mkdir(parents=True, exist_ok=True)
    np.random.seed(123)
    torch.manual_seed(123)
    x, batch, lib, zero = toy_arrays()
    cfg = toy_config("full")
    cfg["seed"] = 123
    cfg["axial"]["attention_dropout"] = 0.0
    context_indices = np.asarray([0, 1, 2, 3], dtype=np.int64)
    student = build_student(n_genes=x.shape[1], config=cfg, assignment=toy_assignment())
    trainer = CAAMTrainer(
        config=cfg,
        student=student,
        train_dataset=CAAMExpressionDataset(x, batch, lib, zero),
        donor_provider=toy_provider(x, batch, lib, zero),
        full_x=torch.as_tensor(x),
        device=torch.device("cpu"),
        save_dir=save_dir,
        context_indices=context_indices,
    )
    trainer.train()
    trainer.save_diagnostics()
    return {
        "context_indices": context_indices,
        "first_batch_indices": np.load(save_dir / "first_batch_indices.npy"),
        "first_mask_hard": np.load(save_dir / "first_mask_hard.npy"),
        "first_donor_indices": np.load(save_dir / "first_donor_indices.npy"),
        "first_loss": trainer.last_batch_debug["first_loss"],
        "context_cache_checksum": np.asarray(trainer.history["context_cache_checksum"], dtype=np.float64),
    }


def test_full_short_training_hard_reproducibility(tmp_path):
    run_a = _run_full_short_training(tmp_path / "a")
    run_b = _run_full_short_training(tmp_path / "b")
    for key in ("context_indices", "first_batch_indices", "first_mask_hard", "first_donor_indices"):
        assert np.array_equal(run_a[key], run_b[key])
    assert np.isclose(run_a["first_loss"], run_b["first_loss"])
    assert np.allclose(run_a["context_cache_checksum"], run_b["context_cache_checksum"])

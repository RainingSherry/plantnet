import numpy as np

from methods.DeepLearning.CAAM_scMAE.data.context_selection import select_context_indices
from methods.DeepLearning.CAAM_scMAE.data.gene_modules import build_gene_modules
from methods.DeepLearning.CAAM_scMAE.tests.helpers import toy_arrays


def test_gene_modules_and_context_reproducible(tmp_path):
    x, *_ = toy_arrays()
    ids1, _ = build_gene_modules(x, 4, 2, 0, tmp_path / "a")
    ids2, _ = build_gene_modules(x, 4, 2, 0, tmp_path / "b")
    ctx1 = select_context_indices(x, 4, 2, 0, tmp_path / "a")
    ctx2 = select_context_indices(x, 4, 2, 0, tmp_path / "b")
    assert np.array_equal(ids1, ids2)
    assert np.array_equal(ctx1, ctx2)


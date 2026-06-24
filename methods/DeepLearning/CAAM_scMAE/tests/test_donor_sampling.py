import numpy as np
import torch

from methods.DeepLearning.CAAM_scMAE.tests.helpers import toy_arrays, toy_provider


def test_donor_not_self_and_gene_wise():
    x, batch, lib, zero = toy_arrays()
    provider = toy_provider(x, batch, lib, zero)
    full = torch.as_tensor(x)
    idx = torch.tensor([0, 1, 2])
    out = provider.sample_batch(idx, full, torch.device("cpu"))
    donors = out["donor_indices"]
    assert donors.shape == (3, x.shape[1])
    assert torch.all(donors != idx.view(-1, 1))
    gene_ids = torch.arange(x.shape[1]).view(1, -1).expand_as(donors)
    expected = full[donors, gene_ids]
    assert torch.allclose(out["replacement"], expected)
    assert any(len(torch.unique(donors[row])) > 1 for row in range(donors.shape[0]))


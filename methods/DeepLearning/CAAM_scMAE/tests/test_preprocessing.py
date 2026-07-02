import numpy as np

from methods.DeepLearning.CAAM_scMAE.data.dataset import CAAMExpressionDataset


def test_training_dataset_returns_no_label_fields():
    ds = CAAMExpressionDataset(np.ones((3, 4), dtype=np.float32))
    sample = ds[0]
    assert "x" in sample and "index" in sample
    for forbidden in ("label", "labels", "cell_type", "true_cluster", "n_clusters"):
        assert forbidden not in sample


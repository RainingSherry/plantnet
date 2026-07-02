import numpy as np
import pytest

from methods.DeepLearning.CAAM_scMAE.data.dataset import CAAMExpressionDataset, assert_no_training_labels


def test_no_label_leakage_guard_accepts_training_batch():
    sample = CAAMExpressionDataset(np.ones((2, 3), dtype=np.float32))[0]
    assert_no_training_labels(sample)


def test_no_label_leakage_guard_rejects_labels():
    with pytest.raises(AssertionError):
        assert_no_training_labels({"x": 1, "label": 0})


import json
import subprocess
import sys

import anndata as ad
import numpy as np


def test_runner_writes_formal_artifacts(tmp_path):
    x = np.log1p(np.random.default_rng(0).poisson(2, size=(12, 8)).astype(np.float32))
    adata = ad.AnnData(X=x)
    adata.obs["cell_type"] = ["a"] * 6 + ["b"] * 6
    data_path = tmp_path / "toy.h5ad"
    out = tmp_path / "out"
    adata.write_h5ad(data_path)
    cmd = [
        sys.executable,
        "methods/DeepLearning/CAAM_scMAE/run.py",
        "--data_path",
        str(data_path),
        "--save_dir",
        str(out),
        "--n_clusters",
        "2",
        "--seed",
        "0",
        "--epochs",
        "1",
        "--variant",
        "control",
        "--no_cuda",
        "--benchmark_mode",
        "true",
    ]
    proc = subprocess.run(cmd, cwd=str(__import__("pathlib").Path.cwd()), capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr
    for name in ["metrics.json", "embedding_final.npy", "labels.npy", "args.json", "artifact_manifest.json"]:
        assert (out / name).exists()
    manifest = json.loads((out / "artifact_manifest.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["variant"] == "control"


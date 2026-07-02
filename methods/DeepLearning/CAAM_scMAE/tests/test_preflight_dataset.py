import csv
import json

import numpy as np
from anndata import AnnData

from methods.DeepLearning.CAAM_scMAE.benchmark.preflight_dataset import run_preflight


def test_preflight_dataset_writes_report_and_summary(tmp_path):
    x = np.asarray(
        [
            [1.0, 0.0, 2.0, 0.0, 3.0, 0.0],
            [0.0, 1.0, 0.0, 2.0, 0.0, 3.0],
            [2.0, 0.0, 1.0, 0.0, 4.0, 0.0],
            [0.0, 2.0, 0.0, 1.0, 0.0, 4.0],
            [3.0, 0.0, 4.0, 0.0, 1.0, 0.0],
            [0.0, 3.0, 0.0, 4.0, 0.0, 1.0],
        ],
        dtype=np.float32,
    )
    adata = AnnData(x)
    adata.obs["batch"] = ["a", "a", "a", "b", "b", "b"]
    path = tmp_path / "toy.h5ad"
    adata.write_h5ad(path)

    reports, report_path, summary_path = run_preflight(
        [str(path)],
        tmp_path / "preflight",
        max_runtime_genes=100,
        max_runtime_params=10_000_000,
        max_estimate_cells=6,
        gene_chunk_size=3,
    )

    assert report_path.exists()
    assert summary_path.exists()
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["datasets"][0]["data_path"] == str(path)
    assert reports[0]["n_cells"] == 6
    assert reports[0]["n_genes"] == 6
    assert reports[0]["quick_ablation_status"] in {
        "pass",
        "blocked_by_budget_deficit",
        "blocked_by_parameter_match",
        "blocked_by_runtime",
        "blocked_by_other",
    }
    assert "estimated_eligibility_rate" in reports[0]
    assert "mlp_parammatched_hidden_dim" in reports[0]

    with summary_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["dataset_name"] == "toy"
    assert "estimated_budget_deficit_rate" in rows[0]

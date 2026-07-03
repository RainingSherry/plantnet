import json

from methods.DeepLearning.CAAM_scMAE.benchmark.summarize_ablation import summarize


def _write_run(root, role, method, ari):
    run_dir = root / f"{role}__seed1"
    run_dir.mkdir(parents=True)
    metrics = {
        "kmeans_known_k": {"acc": 0.5 + ari, "nmi": 0.4 + ari, "ari": ari, "f1_macro": 0.3 + ari},
        "leiden_fixed": {"acc": 0.45 + ari, "nmi": 0.35 + ari, "ari": ari / 2.0, "f1_macro": 0.25 + ari},
    }
    artifact = {"method": method, "variant": role, "seed": 1}
    (run_dir / "metrics.json").write_text(json.dumps(metrics), encoding="utf-8")
    (run_dir / "artifact_manifest.json").write_text(json.dumps(artifact), encoding="utf-8")


def test_summarize_ablation_writes_summary_and_interaction_report(tmp_path):
    _write_run(tmp_path, "control", "caam_scmae_control", 0.10)
    _write_run(tmp_path, "axial", "caam_scmae_axial", 0.12)
    _write_run(tmp_path, "advmask", "caam_scmae_advmask", 0.13)
    _write_run(tmp_path, "full", "caam_scmae_full", 0.20)
    _write_run(tmp_path, "mlp_parammatched", "caam_scmae_mlp_parammatched", 0.18)

    summary_path, report_path = summarize(tmp_path)

    assert summary_path.exists()
    assert report_path.exists()
    summary_text = summary_path.read_text(encoding="utf-8")
    assert "kmeans_known_k.ari" in summary_text
    assert "leiden_fixed.f1_macro" in summary_text

    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["Y00"]["kmeans_known_k.ari"] == 0.10
    assert report["Y10"]["kmeans_known_k.ari"] == 0.12
    assert report["Y01"]["kmeans_known_k.ari"] == 0.13
    assert report["Y11"]["kmeans_known_k.ari"] == 0.20
    assert round(report["delta_AB"]["kmeans_known_k.ari"], 6) == 0.05
    assert round(report["full_minus_parammatched_mlp"]["kmeans_known_k.ari"], 6) == 0.02
    assert report["claim_status"] == "candidate_positive_interaction"
    assert report["claim_status"] != "synergy_confirmed"

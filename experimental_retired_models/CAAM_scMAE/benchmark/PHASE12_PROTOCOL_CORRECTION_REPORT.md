# Phase 12 Protocol Correction Report

## previous_phase_check

- `caam_scmae` remains `authenticity: PENDING` in `methods/method_manifest.yaml`.
- `caam_scmae` remains `default_in_formal: false`.
- No formal manifest entries were found for `caam_scmae_control`, `caam_scmae_axial`, or `caam_scmae_advmask`.
- `registry.py` no longer forces `n_top_genes=0` in `benchmark_mode`.
- `validate_formal_smoke.py` no longer requires `n_top_genes == 0`.
- `risk_and_stop_criteria.md` already treats effective budget deficit as a diagnostic unless strict mode is explicitly enabled.
- The working tree contained pre-existing non-CAAM uncommitted changes; they were not staged or modified for this phase.
- `AGENTS.md` was requested but is absent from the repository root.

## old_protocol_behavior

- `benchmark_mode` forced `input_mode=log1p`, `n_top_genes=0`, and `scale_input=false`.
- `load_caam_data(..., benchmark_mode=True)` also forced `n_top_genes=0`, so explicit benchmark configs could not use the HVG feature space.
- Budget deficit was gated through `runtime.fail_fast`, making it a default hard failure.
- `validate_formal_smoke.py` rejected corrected HVG runs because it hard-coded `preprocessing.n_top_genes == 0`.

## new_protocol_behavior

- `benchmark_mode` defaults to the HVG feature space:
  - `input_mode=log1p`
  - `n_top_genes=2000`
  - `scale_input=false`
- CLI `--n_top_genes 0` is preserved and treated as a full-gene stress-test protocol.
- `strict_effective_budget` defaults to `false`; budget deficit only fails when explicitly set to `true`.
- `zero_to_zero_rate`, `effective_corruption_rate`, `budget_deficit_rate`, `mean_abs_delta`, and `mean_abs_delta_masked` are written to `corruption_stats.json`.
- Selected genes are recorded as `selected_gene_indices.npy` and `selected_genes.txt`, with paths recorded in `resolved_config.yaml`.
- HVG selection is deterministically truncated to exactly `n_top_genes` when the dataset has more genes than the target.

## smoke commands

```bash
MPLCONFIGDIR=/tmp/matplotlib NUMBA_CACHE_DIR=/tmp/numba-cache /data/luolie/conda/envs/scssl_bench_py310/bin/python -m compileall -q methods/DeepLearning/CAAM_scMAE
```

```bash
MPLCONFIGDIR=/tmp/matplotlib NUMBA_CACHE_DIR=/tmp/numba-cache /data/luolie/conda/envs/scssl_bench_py310/bin/python - <<'PY'
# Phase 12 behavior validation: config defaults, CLI n_top_genes=0,
# strict budget gate, reproducible selected_gene_indices, validator HVG acceptance.
PY
```

```bash
CUDA_VISIBLE_DEVICES=1 MPLCONFIGDIR=/tmp/matplotlib NUMBA_CACHE_DIR=/tmp/numba-cache /data/luolie/conda/envs/scssl_bench_py310/bin/python methods/DeepLearning/CAAM_scMAE/run.py --data_path /tmp/caam_phase12_smoke/toy_smoke.h5ad --save_dir /tmp/caam_phase12_smoke/full_seed42 --dataset_name toy_smoke --method_name caam_scmae --variant full --n_clusters 2 --seed 42 --gpu 1 --benchmark_mode true --epochs 1 --batch_size 12 --latent_dim 8
```

```bash
MPLCONFIGDIR=/tmp/matplotlib NUMBA_CACHE_DIR=/tmp/numba-cache /data/luolie/conda/envs/scssl_bench_py310/bin/python methods/DeepLearning/CAAM_scMAE/benchmark/validate_formal_smoke.py /tmp/caam_phase12_smoke/full_seed42
```

## smoke status

- Compile check: PASS.
- Phase 12 behavior validation: PASS.
- Corrected-protocol `variant=full` 1-seed GPU smoke on GPU 1: PASS.
- `validate_formal_smoke.py`: PASS for `/tmp/caam_phase12_smoke/full_seed42`.
- `pytest` was not available in the server Python environments, so equivalent direct Python validation was used.

## whether smoke PASS remains valid

The old full-gene smoke should be treated as old-protocol smoke. The corrected-protocol minimal GPU smoke and validator pass, so the current smoke gate is compatible with the Phase 12 protocol. No method manifest fields were changed in this phase.

## files changed

- `methods/DeepLearning/CAAM_scMAE/registry.py`
- `methods/DeepLearning/CAAM_scMAE/data/preprocessing.py`
- `methods/DeepLearning/CAAM_scMAE/run.py`
- `methods/DeepLearning/CAAM_scMAE/trainers/common.py`
- `methods/DeepLearning/CAAM_scMAE/benchmark/validate_formal_smoke.py`
- `methods/DeepLearning/CAAM_scMAE/benchmark/preflight_dataset.py`
- `methods/DeepLearning/CAAM_scMAE/benchmark/run_ablation.py`
- `methods/DeepLearning/CAAM_scMAE/configs/benchmark_main.yaml`
- `methods/DeepLearning/CAAM_scMAE/tests/test_protocol_correction.py`
- `methods/DeepLearning/CAAM_scMAE/benchmark/PHASE12_PROTOCOL_CORRECTION_REPORT.md`

## tests run

- `pytest -q ...` attempted; failed because `pytest` is not installed.
- `python -m pytest -q ...` attempted; failed because `pytest` is not installed.
- `compileall` PASS.
- Direct Phase 12 behavior validation PASS.
- Tiny corrected-protocol full smoke PASS.
- Formal smoke validator PASS.

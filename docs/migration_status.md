# Migration Status

This document tracks the migration of scCluBench baseline methods into the isolated `methods/` directory structure.

## Summary

| Status | Count |
|---|---|
| PASS | 4 (Priority A: DEC, scDCC, scDeepCluster, scDSC) |
| ENV-BLOCKED | 4 (Priority B: scNAME, scziDesk, DESC; Priority B GNN: scGNN, scCDCG, AttentionAE_sc) |
| PLACEHOLDER | 3 (Foundation: scGPT, GeneFormer, GeneCompass) |

---

## Priority A — Smoke Test Required

| Model | Source | Target | Framework | Status | Smoke Test | Notes |
|---|---|---|---|---|---|---|
| DEC | OtherMode/.../DeepLearning/dec | methods/DeepLearning/dec | PyTorch | migrated | PASS | `--help` OK, `--no_cuda` supported |
| scDCC | OtherMode/.../DeepLearning/scDCC | methods/DeepLearning/scDCC | PyTorch | migrated | PASS | `--help` OK, `--no_cuda` supported |
| scDeepCluster | OtherMode/.../DeepLearning/scDeepCluster | methods/DeepLearning/scDeepCluster | TensorFlow | migrated, TF-gated | PASS | `--help` OK, TF required at runtime |
| scDSC | OtherMode/.../GNN/scDSC | methods/GNN/scDSC | PyTorch | migrated | PASS | `--help` OK, `--no_cuda` supported |

---

## Priority B — Env-Blocked (TensorFlow / Complex Dependencies)

| Model | Source | Target | Framework | Status | Notes |
|---|---|---|---|---|---|
| scNAME | OtherMode/.../DeepLearning/scNAME | methods/DeepLearning/scNAME | TensorFlow | env-blocked | requires TensorFlow/Keras |
| scziDesk | OtherMode/.../DeepLearning/scziDesk | methods/DeepLearning/scziDesk | TensorFlow | env-blocked | requires TensorFlow/Keras |
| DESC | OtherMode/.../DeepLearning/desc | methods/DeepLearning/desc | TensorFlow | env-blocked | requires TensorFlow/Keras |
| scGNN | OtherMode/.../GNN/scGNN | methods/GNN/scGNN | PyTorch | env-blocked | requires graph dependency check |
| scCDCG | OtherMode/.../GNN/scCDCG | methods/GNN/scCDCG | PyTorch | env-blocked | dependency check pending |
| AttentionAE_sc | OtherMode/.../GNN/AttentionAE-sc | methods/GNN/AttentionAE_sc | PyTorch | env-blocked | directory renamed (hyphen → underscore) |

---

## Priority C — Foundation Model Placeholders

| Model | Source | Target | Framework | Status | Notes |
|---|---|---|---|---|---|
| scGPT | OtherMode/.../Foundation/scGPT | methods/Foundation/scGPT | Foundation | placeholder | weights not downloaded, not auto-run |
| GeneFormer | OtherMode/.../Foundation/GeneFormer | methods/Foundation/GeneFormer | Foundation | placeholder | weights not downloaded, not auto-run |
| GeneCompass | OtherMode/.../Foundation/GeneCompass | methods/Foundation/GeneCompass | Foundation | placeholder | weights not downloaded, not auto-run |

---

## Structural Changes

- `AttentionAE-sc` renamed to `AttentionAE_sc` (合法 Python 包名)
- `methods/utils.py` 的 `save()` 增强输出：`embedding_final.npy`, `labels.npy`, `args.json`, `preprocess_config.json`
- 所有 Priority A 模型 `run.py` 支持 `--no_cuda` 参数
- 所有 Priority A 模型输出包含标准化文件：`embedding_final.npy`, `labels.npy`, `metrics.json`, `args.json`
- Priority B 模型提供 `--help` 接口，缺失依赖时给出清晰错误信息

---

## Non-Migrated (Reference Only)

| Model | Path | Notes |
|---|---|---|
| All scCluBench baselines | OtherMode/scCluBench-main/ | reference-only, not imported by runnable code |
| PlantSPADE_LGCL | methods/DeepLearning/PlantSPADE_LGCL/ | legacy, may be removed after NeighborMix_scMAE migration complete |

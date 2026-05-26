# PlantSPADE-LGCL

PlantSPADE-LGCL is the main paper-facing method in this repository. The goal is a clean, reproducible protocol for plant single-cell sparse support geometry, not a collection of unrelated model entries.

## Method Summary

PlantSPADE-LGCL decomposes expression as `X = M o A`:

- `M`: support matrix, indicating whether a cell-gene expression event is observed as non-zero.
- `A`: amplitude matrix, using shared normalization and `log1p`.
- Local view: cell-gene bipartite support graph with LightGCN-style propagation.
- Global view: TF-IDF + SVD low-rank embedding.
- Training: BPR ranking loss plus InfoNCE local-global alignment.
- Interpretation: sparse `SupportGeneAttention` over each cell support set `S_c`, with candidate explanation genes compared against DEG markers.

Full method notes are in `methods/DeepLearning/PlantSPADE_LGCL/docs/method_note.md`.

## Data

The configured 8 h5ad datasets are listed in:

```text
methods/DeepLearning/PlantSPADE_LGCL/configs/datasets_8plant.yaml
```

Expected file paths:

```text
data/CRA002977_1.h5ad
data/CRA007122.h5ad
data/SRP145013.h5ad
data/SRP171040.h5ad
data/SRP182008.h5ad
data/SRP224648.h5ad
data/SRP235541.h5ad
data/SRP309176.h5ad
```

Every run writes `dataset_profile.json`, `preprocess_config.json`, and `selected_genes.txt`.

## Main Methods and Baselines

Main controls:

- `traditional_pca`
- `phytocluster`
- `scvi`
- `scmae`
- `plantspade_lgcl_baseline`
- `plantspade_lgcl_support_attention`

Ablations:

- `plantspade_lgcl_attention_no_idf`
- `plantspade_lgcl_attention_no_amplitude`
- `plantspade_lgcl_attention_topk_64/128/256`
- `plantspade_lgcl_neg_random_zero`
- `plantspade_lgcl_neg_idf_weighted_zero`
- `plantspade_lgcl_neg_neighbor_conflict_zero`

## Run

Do not use GPU `0` or `7`. The new runners reject those IDs.

Profile datasets:

```bash
python methods/DeepLearning/PlantSPADE_LGCL/scripts/profile_datasets.py
```

Run the Arabidopsis root main group with five seeds:

```bash
python methods/DeepLearning/PlantSPADE_LGCL/scripts/run_suite.py \
  --groups arabidopsis_root_main \
  --seeds 1,2,3,4,5 \
  --methods traditional_pca,phytocluster,scvi,scmae,plantspade_lgcl_baseline,plantspade_lgcl_support_attention \
  --gpus 1,2,3,4,5,6
```

Run the full configured suite:

```bash
python methods/DeepLearning/PlantSPADE_LGCL/scripts/run_suite.py --gpus 1,2,3,4,5,6
```

Each dataset/method/seed is a separate process. By default the suite uses one concurrent run per GPU. To increase concurrency on large-memory GPUs:

```bash
python methods/DeepLearning/PlantSPADE_LGCL/scripts/run_suite.py \
  --gpus 1,2,3,4,5,6 \
  --jobs_per_gpu 2
```

Aggregate tables:

```bash
python methods/DeepLearning/PlantSPADE_LGCL/scripts/aggregate_results.py
```

## Evaluation Protocol

Main tables use fixed protocol only:

- `kmeans_known_k`
- `leiden_fixed`, default resolution `1.0`
- `louvain_fixed`, default resolution `1.0` when available

`leiden_oracle_best` is written only as a supplementary upper bound. It is selected by NMI over a resolution sweep and must not be reported as the main result.

## Outputs

Aggregation writes:

- `all_results_long.csv`
- `all_results_mean_std.csv`
- `table_main_fixed_protocol.csv`
- `table_oracle_supplement.csv`
- `table_attention_ablation.csv`
- `table_negative_sampling_ablation.csv`
- `dataset_profiles_summary.csv`

## Archive

Diffusion, maskdiffusion, DOLORIS, GraphDiffusion, unstable GNN entries, and foundation-model experiments are retained as archive/exploratory code. They are not part of the main PlantSPADE-LGCL runner or main result tables. See `methods/DeepLearning/PlantSPADE_LGCL/docs/archive_note.md`.

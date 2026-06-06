# Experiment Protocol

This protocol defines the paper-facing PlantSPADE-LGCL experiments.

## Dataset Groups

Config: `methods/DeepLearning/PlantSPADE_LGCL/configs/datasets_8plant.yaml`

- Arabidopsis root main reproducibility: `SRP182008`, `SRP235541`, `SRP171040`
- Cross-species root generalization: `SRP309176`, `SRP145013`
- Leaf/non-root generalization: `CRA002977_1`, `SRP224648`
- Special tissue robustness: `CRA007122`

Each dataset entry records `dataset_name`, `file_path`, `species`, `tissue`, `label_key`, `expected_n_clusters`, `group`, and `note`.

## Canonical Input

Every run writes:

- `dataset_profile.json`
- `preprocess_config.json`
- `selected_genes.txt`

Support matrix priority:

1. `adata.layers["counts"]`
2. `adata.raw.X`
3. raw-count-looking `adata.X`

Amplitude matrix:

1. normalize total counts to `1e4`
2. apply `log1p`
3. select shared HVGs, main value `n_top_genes=2000`

scVI is the exception among baselines: it should receive raw counts because its likelihood is count-based.

## Evaluation

Main results use fixed protocol only:

- `kmeans_known_k`: embedding probing with true number of classes `K`
- `leiden_fixed`: fixed resolution, default `1.0`
- `louvain_fixed`: fixed resolution, default `1.0`, supplementary if dependency is available

Supplementary oracle protocol:

- `leiden_oracle_best`: selected by NMI over a fixed resolution sweep; this is an upper bound and is not used in the main table.

Full sweep attachment:

- one row per Leiden resolution with `ACC`, `NMI`, `ARI`, `F1_macro`, `FMI`, `V-measure`, `homogeneity`, `completeness`, `n_pred_clusters`, and `silhouette`.

The old practice of choosing Leiden resolution by true-label NMI and reporting it as `leiden` is not part of the main protocol.

## Methods in Main Scope

Traditional:

- `traditional_pca`: PCA/SVD embedding with fixed KMeans, fixed Leiden, fixed Louvain

Deep baselines:

- `phytocluster`
- `scvi`
- `scmae`

Main method:

- `plantspade_lgcl_baseline`
- `plantspade_lgcl_support_attention`

Ablations:

- `plantspade_lgcl_attention_no_idf`
- `plantspade_lgcl_attention_no_amplitude`
- `plantspade_lgcl_attention_topk_64`
- `plantspade_lgcl_attention_topk_128`
- `plantspade_lgcl_attention_topk_256`
- `plantspade_lgcl_neg_random_zero`
- `plantspade_lgcl_neg_idf_weighted_zero`
- `plantspade_lgcl_neg_neighbor_conflict_zero`

Diffusion, maskdiffusion, DOLORIS, and GraphDiffusion experiments have been removed from the active working tree. Unstable GNNs and foundation-model experiments remain excluded from the main tables.

## Commands

Profile all configured datasets:

```bash
python methods/DeepLearning/PlantSPADE_LGCL/scripts/profile_datasets.py
```

Run the three main Arabidopsis root datasets, five seeds, main methods and baselines:

```bash
python methods/DeepLearning/PlantSPADE_LGCL/scripts/run_suite.py \
  --groups arabidopsis_root_main \
  --seeds 1,2,3,4,5 \
  --methods traditional_pca,phytocluster,scvi,scmae,plantspade_lgcl_baseline,plantspade_lgcl_support_attention \
  --gpus 1,2,3,4,5,6
```

Run all eight datasets:

```bash
python methods/DeepLearning/PlantSPADE_LGCL/scripts/run_suite.py --gpus 1,2,3,4,5,6
```

`run_suite.py` schedules one dataset/method/seed process per slot. Slots default to one per configured GPU. Use `--jobs_per_gpu 2` for two concurrent runs per GPU, or `--jobs N` to set a total number of slots round-robin across the GPU list.

Aggregate tables:

```bash
python methods/DeepLearning/PlantSPADE_LGCL/scripts/aggregate_results.py
```

Do not use GPU `0` or `7`; the new runners reject those IDs.

## Output Tables

`aggregate_results.py` writes:

- `all_results_long.csv`
- `all_results_mean_std.csv`
- `table_main_fixed_protocol.csv`
- `table_oracle_supplement.csv`
- `table_attention_ablation.csv`
- `table_negative_sampling_ablation.csv`
- `dataset_profiles_summary.csv`

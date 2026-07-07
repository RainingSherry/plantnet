# rank25_anomaly_association_boundary_full

Independent-full scMAE candidate adapted from **Anomaly Transformer: Time Series Anomaly Detection with Association Discrepancy**.

## Method Basis

The original paper proposes Anomaly-Attention, where each position has:

- a learned Gaussian **prior-association** over nearby positions;
- a self-attention **series-association** learned from the data;
- a symmetric KL association discrepancy between them.

The local report marks this paper as a cautious/background source and warns against using it directly as a main single-cell module. This implementation therefore does not build a time-series anomaly detector. It adapts the association discrepancy into a gene-module boundary-risk regularizer for scMAE.

## scMAE Gap Addressed

This candidate targets the **boundary / rare-cell / robust loss** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- genes are ordered by a covariance first-eigenvector direction computed from `log_expr`;
- ordered genes are split into module patches;
- an Anomaly-Transformer-style association block compares learned series attention against a Gaussian prior over adjacent gene modules;
- masked-view association discrepancy is encouraged to match clean-view discrepancy;
- the discrepancy is saved as `association_risk.npy` and summarized in `diagnostics.json`.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as the masked expression reconstruction target.
- No NB/ZINB/token/count objective is used.
- No generated cells are evaluated.

## NeighborMix Relationship

NeighborMix is not used in this candidate. The relationship is independent and potentially complementary. Since there is no cell mixing, `mixed_cell_fraction` is always `0.0`.

## Difference From The Original Paper

The original Anomaly Transformer is supervised by a minimax anomaly criterion for time-series windows. This implementation is a cautious structural adaptation for unordered single-cell expression matrices:

- temporal positions become ordered gene-module patches;
- association discrepancy is used for masked-view stability and diagnostics;
- no anomaly labels, time windows, or minimax training are used;
- screen results are only candidate evidence and are not formal performance claims.

## Example Smoke

```bash
CUDA_VISIBLE_DEVICES=1 python /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/rank25_anomaly_association_boundary_full/run.py --data_path /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/benchmark_data/Melanoma_5K.h5ad --save_dir /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/runs/scMAEs_independent_smoke/rank25_anomaly_association_boundary_full/Melanoma_5K/seed42 --dataset_name Melanoma_5K --label_key resolved_label --input_mode auto --n_top_genes 1000 --target_sum 10000 --scale_input true --n_clusters 9 --epochs 3 --batch_size 256 --lr 0.001 --weight_decay 0.00001 --seed 42 --gpu 1 --no_save_h5ad --smoke
```

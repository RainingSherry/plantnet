# Rank 02 MaskFeat-scMAE

This directory is the independent Rank 02 implementation. It does not use the
legacy shared `common/model.py` variant switch.

## Source

- Paper: `MaskFeat: Masked Feature Prediction for Self-Supervised Visual Pre-Training`
- Local PDF: `../参考文献/01_PDF论文_按推荐程度排序/002_很高_MaskFeat_Masked_Feature_Prediction_for_Self-Supervised_Visual_Pre-Training.pdf`
- GitHub URL from index: `https://github.com/facebookresearch/SlowFast/tree/main/projects/maskfeat`
- Local source clone: `../external_sources/rank02_slowfast_maskfeat`
- Commit: `287ec0076846560f44a9327e931a5a2360240533`

## scMAE Adaptation

The original MaskFeat predicts hand-crafted HOG features for masked visual
patches. This scRNA adaptation partitions genes into fixed-size patches,
replaces masked patch embeddings with a learned mask token, encodes the full
patch sequence with a Transformer, and predicts deterministic gene-patch
features on masked patches.

The gene-patch feature target concatenates:

- z-normalized expression values within the patch,
- local first-difference expression gradients,
- within-patch expression ranks.

This is a full masked feature prediction objective over gene patches, not a
single auxiliary statistic head.

## Files

- `model.py`: gene patch MaskFeat Transformer with mask token.
- `loss.py`: masked feature prediction loss with explicit denominator.
- `run.py`: benchmark-compatible training/evaluation entrypoint.
- `source_manifest.json`: paper/source provenance.

## Mask Semantics

`mask = 1` means the entire gene patch was replaced by the learned mask token
and is included in the masked feature-prediction denominator.


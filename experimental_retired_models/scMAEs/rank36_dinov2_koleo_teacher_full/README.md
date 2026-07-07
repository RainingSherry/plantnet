# rank36_dinov2_koleo_teacher_full

Independent-full scMAE candidate adapted from **DINOv2: Learning Robust Visual Features without Supervision**.

## Method Basis

The local report recommends DINOv2 as a teacher-student self-distillation mechanism. The paper combines DINO class-token distillation, iBOT masked patch-token distillation, EMA teacher updates, centered/sharpened teacher probabilities, and the KoLeo nearest-neighbor entropy regularizer for spreading representations. The official GitHub implementation was inspected for `DINOLoss`, `iBOTPatchLoss`, `KoLeoLoss`, `ssl_meta_arch.py`, and the default SSL config.

## scMAE Gap Addressed

This candidate targets the **teacher / semantic target / collapse-robustness** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- the student receives a strong masked gene-module view;
- the EMA teacher receives a weak expression view;
- DINO-style class-token CE aligns global cell semantics;
- iBOT-style masked module CE aligns masked local gene modules;
- KoLeo regularization spreads cell embeddings and is logged in training diagnostics.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- No NB/ZINB count likelihood is used.
- No offline token target is produced from scaled expression; soft targets come from the EMA teacher.

## NeighborMix Relationship

NeighborMix is not used. This method is independent and complementary: it stabilizes semantic targets and embedding spread but never mixes cells. `mixed_cell_fraction=0.0`.

## Difference From Original DINOv2

The original DINOv2 is a large-scale image ViT framework. This scRNA-seq implementation reconstructs only the training mechanism for gene-module tokens and uses no external DINOv2 weights.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.

# rank34_masked_siamese_proto_full

Independent-full scMAE candidate adapted from **Masked Siamese Networks for Label-Efficient Learning**.

## Method Basis

The local report recommends MSN as a teacher-student self-distillation mechanism. The paper defines an anchor view with random masking and a target view without masking, then matches their distributions over learnable prototypes. The official GitHub implementation confirms the central pieces used here: a momentum target encoder, learnable prototypes, soft-nearest-neighbor prototype probabilities, target sharpening, and mean-entropy maximization (`me-max`) to avoid prototype collapse.

## scMAE Gap Addressed

This candidate targets the **teacher / semantic target** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- the student/anchor encoder receives a strong masked expression view;
- the EMA teacher/target encoder receives a weak expression view;
- the masked anchor representation is matched to the sharpened target prototype distribution;
- `me-max` encourages use of the prototype set rather than a single collapsed assignment.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- No NB/ZINB count likelihood is used.
- No offline token or count target is produced from scaled expression; prototype targets come from the EMA teacher projection.

## NeighborMix Relationship

NeighborMix is not used. This method is independent and complementary: it aligns masked and weak target representations but never mixes cells. `mixed_cell_fraction=0.0`.

## Difference From Original MSN

The original MSN is an image ViT method that can drop masked patches for efficiency. This scRNA-seq version uses an MLP encoder over HVG expression vectors and keeps masked values as zeros because scMAE also needs explicit mask prediction and expression reconstruction.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.

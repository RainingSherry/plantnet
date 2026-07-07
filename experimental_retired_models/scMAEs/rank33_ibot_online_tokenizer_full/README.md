# rank33_ibot_online_tokenizer_full

Independent-full scMAE candidate adapted from **iBOT: Image BERT Pre-Training with Online Tokenizer**.

## Method Basis

The local report recommends iBOT as a teacher-student self-distillation mechanism for stabilizing scMAE targets. The paper describes masked prediction with an online tokenizer: the student sees a masked view, while an EMA teacher sees the clean or weakly augmented view and provides centered, sharpened probability targets for both the class token and masked patch tokens. The official GitHub implementation confirms the key mechanics used here: EMA teacher update, center buffers, teacher sharpening, student-temperature CE, and masked-token CE computed only on masked positions.

## scMAE Gap Addressed

This candidate targets the **teacher / semantic target** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- the student receives a strong gene mask view;
- the EMA teacher receives a weak expression view;
- the teacher acts as an online tokenizer for a cell-level class token and masked gene-module tokens;
- centered/sharpened teacher probabilities provide soft targets instead of fixed offline gene tokens.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as the masked expression reconstruction target.
- No NB/ZINB count likelihood is used.
- No discrete token target is produced from scaled expression; online token distributions come from the EMA teacher.

## NeighborMix Relationship

NeighborMix is not used. This method is independent and complementary: it stabilizes targets via an EMA teacher and never mixes cells. `mixed_cell_fraction=0.0`.

## Difference From Original iBOT

The original iBOT is an image ViT method with image patches and multi-crop augmentations. This implementation reconstructs the mechanism for scRNA-seq by using contiguous gene modules as tokens, a weak expression view for the teacher, and scMAE reconstruction/mask heads for the student.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.

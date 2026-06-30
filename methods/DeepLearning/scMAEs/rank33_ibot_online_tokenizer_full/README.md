# Rank 33: iBOT Online Tokenizer scMAE

This folder contains an independent full variant for **iBOT: Image BERT Pre-Training with Online Tokenizer**.

## Source Basis

- Local PDF: `methods/DeepLearning/scMAEs/参考文献/01_PDF论文_按推荐程度排序/033_中高_可用_iBOT_Image_BERT_Pre-Training_with_Online_Tokenizer.pdf`
- PDF SHA-256: `e76c32b28234434443cba27abda926ec7c69c64af240d2d7aa78427131c4273d`
- GitHub URL: `https://github.com/bytedance/ibot`
- Git commit: `da316d82636a7a7356835ef224b13d5f3ace0489`

The implementation follows the paper's online tokenizer idea: a teacher network produces centered soft targets, the student predicts masked patch token distributions, and the class token is trained with self-distillation. The source files `main_ibot.py` and `models/head.py` provide the original iBOT loss/head structure.

## scMAE Adaptation

- `model.py` implements expression patch embeddings, a student Transformer, an EMA teacher Transformer, shared iBOT projection heads for `[CLS]` and patch tokens, teacher centers, and auxiliary masked expression reconstruction.
- `loss.py` implements centered teacher soft targets, class-token cross-view self-distillation, masked patch-token distillation, masked reconstruction, and mask prediction. Mask semantics are fixed as `1 = corrupted/replaced/masked target`.
- `run.py` builds two stochastic expression views per cell by feature dropout plus small Gaussian noise, masks patches independently for both views, updates the EMA teacher and centers after each optimizer step, and writes the standard independent benchmark outputs.

## Not Reproduced

The original iBOT uses image crops, ViT patch grids, and very large output vocabularies. This adaptation keeps the core online tokenizer, shared projection head, centered teacher distributions, EMA teacher, class-token self-distillation, and masked patch-token loss, but replaces image crops with expression-vector stochastic views and uses a smaller output dimension suitable for the fixed scRNA quick-screen protocol.

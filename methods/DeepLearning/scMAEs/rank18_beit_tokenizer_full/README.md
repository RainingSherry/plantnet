# Rank 18 BEiT expression tokenizer scMAE

This folder independently adapts **BEiT: BERT Pre-Training of Image Transformers** to scRNA expression clustering.

## Source and reproduction scope

- GitHub source: `https://github.com/microsoft/unilm/tree/master/beit`.
- Local sparse source snapshot: `../external_sources/rank18_unilm_beit`, commit `833df7e7832e5064a281131ee64a481afa8e5b95`.
- PDF: `../参考文献/01_PDF论文_按推荐程度排序/018_高_BEiT_BERT_Pre-Training_of_Image_Transformers.pdf`.

## Implemented paper mechanisms

- Continuous expression patches are used as model input.
- Masked patches are replaced by a learned `[M]` mask token.
- A Transformer encoder predicts discrete target tokens only at masked patch positions.
- The loss uses masked token cross-entropy, matching BEiT's masked visual token prediction.

## Expression tokenizer adaptation

The original BEiT uses a frozen visual tokenizer/dVAE. For scRNA data, no external visual tokenizer is appropriate. This variant builds a dataset-local unsupervised tokenizer from the same HVG genes:

- Use nonnegative log-normalized expression with `scale_input=False`.
- Split genes into fixed patches.
- Compute each patch's mean expression.
- Quantize patch means by global quantile edges into a discrete vocabulary.

The encoder/evaluation input remains the same scaled HVG expression used by the shared benchmark.

## scMAE-family adaptation

An auxiliary masked expression reconstruction branch is kept for scMAE-family fairness. Mask semantics are:

`1 = expression patch replaced by learned BEiT mask token`.

No external pretrained image weights or visual tokenizers are loaded.

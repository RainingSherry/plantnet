# Rank 19 data2vec EMA teacher scMAE

This folder independently adapts **data2vec: A General Framework for Self-supervised Learning in Speech, Vision and Language** to scRNA expression clustering.

## Source and reproduction scope

- GitHub source: `https://github.com/facebookresearch/fairseq/tree/main/examples/data2vec`.
- Local sparse source snapshot: `../external_sources/rank19_fairseq_data2vec`, commit `3d262bb25690e4eb2e7d3c1309b1e9c406ca4b99`.
- PDF: `../参考文献/01_PDF论文_按推荐程度排序/019_高_data2vec_A_General_Framework_for_Self-supervised_Learning_in_Speech,_Vision_and_Language.pdf`.

## Implemented paper mechanisms

- Student encoder receives a masked expression-patch sequence with a learned mask token.
- Teacher encoder receives the full unmasked expression-patch sequence.
- Teacher parameters are an exponential moving average of the student encoder.
- Targets are the layer-normalized average of the top K teacher Transformer block outputs.
- The student predicts teacher latent targets only at masked patch positions.
- The primary objective is Smooth L1 latent regression with stop-gradient teacher targets.

## scMAE-family adaptation

An auxiliary masked expression reconstruction branch is retained for fair comparison with the scMAE-family benchmark. Mask semantics are:

`1 = expression patch replaced by learned data2vec mask token`.

No external pretrained data2vec checkpoints are loaded.

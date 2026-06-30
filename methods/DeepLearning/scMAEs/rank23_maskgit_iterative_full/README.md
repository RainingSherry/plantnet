# Rank 23: MaskGIT Iterative Expression-Token scMAE

This folder is an independent scMAE-family adaptation of:

- Paper: MaskGIT: Masked Generative Image Transformer
- Official code: https://github.com/google-research/maskgit
- Local source snapshot: `../external_sources/rank23_maskgit`
- Commit: `1db23594e1bd328ee78eadcd148a19281cd0f5b8`

## scMAE Adaptation

MaskGIT trains a bidirectional transformer to predict randomly masked discrete visual tokens, then decodes by repeatedly filling all masked positions and remasking low-confidence predictions according to a schedule. This implementation maps each cell's expression vector into contiguous gene patches and discretizes patch means with dataset-local quantile bins.

The adapted model keeps the core MaskGIT structure:

- discrete token IDs plus a learned `[MASK]` token;
- BERT-style bidirectional token Transformer;
- tied token embedding / MLM output projection;
- high-mask training with a cosine masking schedule;
- masked-token cross entropy only on masked positions;
- confidence-based iterative decode for representation extraction;
- auxiliary patch reconstruction for scMAE-compatible masked expression training.

Mask semantics are fixed as `1 = quantized expression patch replaced by [MASK]`.

## Fairness Notes

No image tokenizer checkpoint, image data, or pretrained MaskGIT weights are used. The quantizer is fit only on the benchmark training matrix for the current run. Preprocessing, KMeans evaluation, seed handling, and output files stay inside the shared scMAE-family protocol.

## Outputs

`run.py` writes the standard independent scMAEs artifacts: `args.json`, `source_manifest.json`, `quantizer.json`, `training_history.json`, `model_checkpoint.pth`, `embedding_final.npy`, `labels.npy`, `gene_names.npy`, `embedding.h5`, and `eval_fixed.csv`/`metrics.json` when evaluation is enabled.

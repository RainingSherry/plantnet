# Rank 01 ScDiVa-scMAE

This directory is the independent Rank 01 implementation. It does not use the
legacy shared `common/model.py` variant switch.

## Source

- Paper: `ScDiVa: Masked Discrete Diffusion for Joint Modeling of Single-Cell Identity and Expression`
- Local PDF: `../参考文献/01_PDF论文_按推荐程度排序/001_很高_ScDiVa_Masked_Discrete_Diffusion_for_Joint_Modeling_of_Single-Cell_Identity_and_Expression.pdf`
- GitHub URL from index: `https://github.com/SindiLab/ScDiVa`
- GitHub status: anonymous clone failed in this environment, so this implementation is based on the local PDF/report description.

## scMAE Adaptation

The model discretizes normalized expression values into per-gene expression
tokens, applies an absorbing masked discrete diffusion process, and trains a
time-conditioned Transformer denoiser. The embedding used for clustering is the
projected CLS representation at diffusion step 0.

## Files

- `model.py`: time-conditioned discrete diffusion Transformer.
- `loss.py`: token, value, and mask losses with explicit masked denominator.
- `run.py`: standard benchmark entrypoint.
- `run_rank01.py`: implementation entrypoint retained for compatibility.
- `source_manifest.json`: paper/source provenance.

## Mask Semantics

`mask = 1` means the expression token was replaced by the absorbing mask token.

# Rank 11 scMamba selective-SSM scMAE

This folder is an independent scMAE-family variant inspired by scMamba.

## Source and reproduction scope

- Paper/index item: Rank 11 `scMamba`.
- GitHub source: `https://github.com/23AIBox/scMamba`.
- Local source snapshot: `../external_sources/rank11_scmamba`, commit `4887c0a8ab060b2482384d2294fe265b633d2406`.
- Referenced source files:
  - `scmamba2/tokenizer/tokenizer.py`
  - `scmamba2/models/scmamba.py`
  - `scmamba2/models/block.py`
  - `scmamba2/models/config_scmamba.py`

The official implementation depends on `mamba_ssm` and `causal_conv1d`, which are not installed in this environment. This folder therefore does not import the official CUDA kernels and does not silently fall back to a feed-forward gate. Instead, `model.py` implements a dependency-free selective state-space scan with the same core SSM ingredients: patch tokenization, learned negative state dynamics `A`, input/output selectors `B/C`, skip `D`, positive per-token `delta`, causal depthwise convolution, residual SSM blocks, and sequence pooling.

## Combination with scMAE

- Input remains `[batch, genes]` after the shared fair preprocessing.
- Genes are split into fixed-width patches, matching scMamba's patch-based tokenization.
- Mask semantics are fixed as `1 = gene belongs to a masked/corrupted patch`.
- Masked patches are replaced by a learned mask token inside the tokenizer.
- The SSM encoder reconstructs masked expression values through a patch decoder.
- KMeans evaluation uses the pooled SSM embedding and the same known-`k` protocol as the other independent variants.

## Deliberate exclusions

- No external pretrained scMamba checkpoints are loaded.
- No multi-omics RNA/ATAC branch is used because the benchmark datasets are single-expression matrices.
- This is not a `Linear + Sigmoid` gate approximation: the recurrent state update is explicit in `SelectiveSSMLayer.selective_scan`.

## Smoke-test expectations

- `model(clean, patch_mask)["embedding"]` has shape `[batch, hidden]`.
- `loss.py` uses the masked-gene denominator `gene_mask.sum().clamp_min(1.0)`.
- A valid backward pass must produce nonzero gradients for the SSM state parameters such as `A_log`.

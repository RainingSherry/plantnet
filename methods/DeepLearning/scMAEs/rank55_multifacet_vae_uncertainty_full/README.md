# rank55_multifacet_vae_uncertainty_full

Independent-full scMAE candidate based on Multi-Facet Clustering Variational Autoencoders.

## Theory basis

MFCVAE uses multiple latent facets, each with a Mixture-of-Gaussians prior, to learn several unsupervised clusterings simultaneously. I inspected the public implementation at `https://github.com/FabianFalck/mfcvae`, especially `mfcvae.py` and `models_fc.py`.

## scMAE integration

This candidate fills scMAE's robust loss / uncertainty gap:

- The encoder branches into multiple `mu/logvar` latent facets.
- Each facet has its own learnable MoG prior.
- The concatenated latent sample feeds a scMAE decoder and mask predictor.
- Posterior variance is saved as an uncertainty diagnostic.
- An NB likelihood is enabled only when raw counts can be aligned to the selected HVG genes.

## Expression semantics

`scaled_expr` may be used only as encoder input. Masked expression reconstruction uses unscaled log-expression. NB NLL uses aligned raw counts and size factors only; if raw counts cannot be aligned, the NB term is disabled and recorded in `preprocess_config.json` and `summary.json`.

## NeighborMix

NeighborMix is not used. The relation is independent and complementary; `mixed_cell_fraction=0.0`.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `posterior_variance.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.

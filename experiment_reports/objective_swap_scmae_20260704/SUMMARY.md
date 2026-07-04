# Objective-swap scMAE experiment summary

## Aggregate

| dataset | assignment | varw | latent | n | KMeans ARI | direct ARI | PCA ARI | eff_dim | std_min |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Macosko | pca | 0.0 | 128 | 3 | 0.8806 +/- 0.0053 | NA | 0.8806 +/- 0.0053 | 44.1 +/- 0.0 | 1.076 +/- 0.000 |
| Macosko | sharpen | 0.02 | 128 | 3 | 0.1755 +/- 0.0395 | 0.4620 +/- 0.4029 | 0.8806 +/- 0.0053 | 127.3 +/- 0.7 | 1.045 +/- 0.030 |
| Macosko | sinkhorn | 0.02 | 128 | 1 | 0.1396 +/- 0.0000 | 0.1396 +/- 0.0000 | 0.8868 +/- 0.0000 | 127.0 +/- 0.0 | 1.020 +/- 0.000 |
| Melanoma_5K | pca | 0.0 | 128 | 3 | 0.6166 +/- 0.0759 | NA | 0.6166 +/- 0.0759 | 44.5 +/- 0.0 | 1.177 +/- 0.001 |
| Quake_10x_Spleen | pca | 0.0 | 128 | 3 | 0.8630 +/- 0.0248 | NA | 0.8630 +/- 0.0248 | 43.1 +/- 0.0 | 1.072 +/- 0.003 |

## Per-run details

| run | assignment | KMeans ARI | direct ARI | PCA ARI | eff_dim | std_min | dims_std>1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| macosko_sharpen_seed42 | sharpen | 0.1846 | 0.7403 | 0.8868 | 127.6 | 1.066 | 128 |
| macosko_sharpen_seed43 | sharpen | 0.2096 | 0.6458 | 0.8773 | 127.7 | 1.059 | 128 |
| macosko_sharpen_seed44 | sharpen | 0.1323 | 0.0000 | 0.8778 | 126.5 | 1.011 | 128 |
| macosko_sinkhorn_seed42 | sinkhorn | 0.1396 | 0.1396 | 0.8868 | 127.0 | 1.020 | 128 |
| pca_macosko_seed42 | pca | 0.8868 | nan | 0.8868 | 44.1 | 1.076 | 128 |
| pca_macosko_seed43 | pca | 0.8773 | nan | 0.8773 | 44.1 | 1.076 | 128 |
| pca_macosko_seed44 | pca | 0.8778 | nan | 0.8778 | 44.1 | 1.076 | 128 |
| pca_melanoma_seed42 | pca | 0.7040 | nan | 0.7040 | 44.5 | 1.178 | 128 |
| pca_melanoma_seed43 | pca | 0.5679 | nan | 0.5679 | 44.5 | 1.178 | 128 |
| pca_melanoma_seed44 | pca | 0.5778 | nan | 0.5778 | 44.5 | 1.176 | 128 |
| pca_quake_seed42 | pca | 0.8885 | nan | 0.8885 | 43.1 | 1.069 | 128 |
| pca_quake_seed43 | pca | 0.8615 | nan | 0.8615 | 43.1 | 1.075 | 128 |
| pca_quake_seed44 | pca | 0.8390 | nan | 0.8390 | 43.1 | 1.071 | 128 |

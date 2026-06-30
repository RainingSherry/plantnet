# Rank 29: Adaptive Fuzzy Clustering scMAE

This folder contains an independent full variant for **Deep Adaptive Fuzzy Clustering for Evolutionary Unsupervised Representation Learning** (Tan et al., arXiv:2103.17086).

## Source Basis

- Local PDF: `methods/DeepLearning/scMAEs/参考文献/01_PDF论文_按推荐程度排序/029_高_Deep_Adaptive_Fuzzy_Clustering_for_Evolutionary_Unsupervised_Representation_Learning.pdf`
- PDF SHA-256: `36d60ac770e78aaa81661f097c02ba71c05d3d2e40f6afc0e1b6f5ee73bfb9ae`
- GitHub URL: none listed in `02_整理索引.csv`.

The implementation follows the PDF method section: feature extraction and reconstruction are jointly trained with fuzzy clustering in the bottleneck space, fuzzy memberships represent soft cluster assignment probabilities, and a weighted adaptive entropy term is optimized with the clustering objective.

## scMAE Adaptation

- `model.py` implements a standalone MLP autoencoder for expression vectors, a trainable fuzzy membership layer, learned centroids, a learnable fuzzifier, and a fuzzy reconstruction decoder from membership-weighted centroids.
- `loss.py` implements masked expression reconstruction, fuzzy compactness, adaptive entropy, partition balance, center separation, and mask prediction. Mask semantics are fixed as `1 = corrupted/replaced/masked target`.
- `run.py` follows the independent benchmark protocol and periodically refreshes centroids from the current latent space with KMeans plus EMA, matching the DAFC idea of alternating fuzzy cluster updates and representation learning.

## Not Reproduced

The original paper is image/CNN oriented and describes pseudo-label similarity feedback and test-error based stopping. This adaptation keeps the core fuzzy bottleneck objective and adaptive entropy, but replaces image ConvNets with expression-vector MLP blocks and evaluates fixed-k cell clustering through the shared scMAE benchmark protocol.


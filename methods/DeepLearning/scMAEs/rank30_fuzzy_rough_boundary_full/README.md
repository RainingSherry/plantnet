# Rank 30: Fuzzy Rough Boundary scMAE

This folder contains an independent full variant for **Fuzzy Rough Sets Based on Fuzzy Quantification** (Theerens and Cornelis, arXiv:2212.04327 / Fuzzy Sets and Systems).

## Source Basis

- Local PDF: `methods/DeepLearning/scMAEs/参考文献/01_PDF论文_按推荐程度排序/030_高_Fuzzy_Rough_Sets_Based_on_Fuzzy_Quantification.pdf`
- PDF SHA-256: `baec986c62573f3d0fc7e622020482a73b762a67c72e9b906be66fc6549af3e8`
- GitHub URL: none listed in `02_整理索引.csv`.

The implementation follows the paper's fuzzy quantifier-based fuzzy rough set framework: lower approximations express that almost all similar objects belong to a concept, upper approximations express that some similar objects belong to a concept, and the Yager Weighted Implication-style lower quantifier supplies noise-tolerant boundary handling.

## scMAE Adaptation

- `model.py` implements a standalone expression autoencoder plus a fuzzy rough pseudo-concept head with Student-t cluster memberships and a boundary membership head.
- `loss.py` implements a batch-local fuzzy relation, Kleene-Dienes implication, RIM S-function quantifiers, YWI-style lower approximation, unary upper approximation, rough boundary width, prototype compactness, and masked expression reconstruction. Mask semantics are `1 = corrupted/replaced/masked target`.
- `run.py` follows the independent benchmark protocol and periodically refreshes cluster centers from the current latent space with KMeans plus EMA.

## Not Reproduced

The original paper is a fuzzy set theory and classification paper rather than a neural single-cell model. This adaptation does not reproduce the UCI classification experiment. It uses the paper's fuzzy rough lower/upper approximation machinery as a boundary-aware unsupervised objective for scMAE latent cell representations.


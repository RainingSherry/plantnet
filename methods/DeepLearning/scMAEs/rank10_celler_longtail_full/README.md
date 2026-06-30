# Rank 10: Celler Long-Tail scMAE

Source paper: **Celler: A Genomic Language Model for Long-Tailed Single-Cell
Annotation**.

The original method is a genomic language model for long-tailed annotation. This
independent scMAE-family variant adapts its two relevant long-tail ideas without
using external pretrained weights:

- Gaussian Inflation Loss (GInf Loss), from `hicell/loss.py`;
- hard sample mining, implemented with reconstruction difficulty and prototype
  uncertainty.

Because this benchmark is unsupervised clustering, true labels are not used for
training. Instead, the current embedding is periodically clustered into
`known_k` pseudo prototypes. Pseudo cluster counts provide the long-tail
`class_counts` for GInf Loss, and rare/hard pseudo cells receive stronger masked
reconstruction weight.

Mask semantics: `1 = expression gene replaced for masked reconstruction branch`.
The reconstruction denominator is the weighted sum of masked entries, clamped
only for degenerate smoke tests.

Fair protocol shared with other independent variants:

- repository `scMAE_family` preprocessing and fixed KMeans known-k evaluation;
- default `n_top_genes=1000`, `target_sum=10000`, `scale_input=True`;
- no dependence on old `scMAEs/common/model.py`;
- no GPU 0 or 7 should be used by the benchmark runner.

Not reproduced:

- large genomic language model pretraining and external Celler weights are not
  used, per the no-large-external-weight assumption;
- supervised annotation labels are not used during training.

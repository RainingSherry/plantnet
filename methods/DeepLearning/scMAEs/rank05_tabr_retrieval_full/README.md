# Rank 05: TabR Retrieval scMAE

Source paper: **TabR: Tabular Deep Learning Meets Nearest Neighbors**.

This independent variant adapts TabR's retrieval-augmented tabular mechanism to
masked single-cell expression reconstruction. The query is a masked cell
expression vector. Candidate cells are encoded into keys, the nearest contexts
are retrieved with self-neighbor removal, and a TabR-style context value is added
back to the query latent before decoding:

```text
value = candidate_value(candidate_expression) + T(query_key - context_key)
query_state = query_state + softmax(sim(query_key, context_key)) @ value
```

The value source differs from supervised TabR: scMAE has no target label during
self-supervised pretraining, so candidate expression/context embedding is used as
the retrieval value instead of `candidate_y`.

Mask semantics: `1 = expression gene masked before retrieval-augmented
reconstruction`. The reconstruction denominator is the number of masked gene
entries in the batch, clamped only for degenerate smoke tests.

Fair protocol shared with other independent variants:

- repository `scMAE_family` preprocessing and fixed KMeans known-k evaluation;
- default `n_top_genes=1000`, `target_sum=10000`, `scale_input=True`;
- no dependence on old `scMAEs/common/model.py`;
- no GPU 0 or 7 should be used by the benchmark runner.

Not reproduced:

- supervised label encoding and tabular regression/classification heads are not
  used because this is an unsupervised scRNA clustering task;
- FAISS is not required for the target datasets used here; exact torch distance
  search keeps the implementation explicit and avoids extra dependency drift.

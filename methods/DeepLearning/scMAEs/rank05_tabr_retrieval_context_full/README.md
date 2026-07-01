# rank05_tabr_retrieval_context_full

Independent-full scMAE candidate inspired by TabR.

The model builds a nearest-neighbor retrieval context from log-normalized
expression and uses that context during both training and inference. The encoder
receives `x`, `context`, and `x - context`, then trains with scMAE mask
prediction, masked expression reconstruction, and a small retrieval consistency
term.

NeighborMix is not used; `mixed_cell_fraction` is always `0.0`.


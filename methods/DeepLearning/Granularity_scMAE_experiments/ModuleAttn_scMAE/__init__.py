"""ModuleAttn-scMAE — gene co-expression module tokenization + module-level
self-attention encoder, on top of the DEC + per-dim variance-floor winner.

Minimal falsifiable test (2026-07): does explicit gene-module structure give an
ORTHOGONAL gain over the plain-MLP AdaptiveSwitch encoder, on Macosko (the
fine-grained dataset most likely to benefit), WITHOUT breaking Quake?
Everything except the encoder (DEC head, variance floor, force_gate=1) is held
identical for clean attribution. Genes are grouped into co-expression modules
first (solves the 'genes have no natural order' problem that breaks naive conv).
"""

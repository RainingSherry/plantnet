# 002 Prototype DEC

This variant keeps original scMAE as the main task and adds a delayed DEC-style
prototype regularizer. It warms up the scMAE embedding, initializes prototypes
with KMeans, then ramps the prototype KL weight from zero.

Default mechanism:

```text
scMAE + delayed DEC prototype KL
```

The prototype loss is high-confidence gated and capped by the warmup weight.


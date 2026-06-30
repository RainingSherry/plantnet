# 005 NeighborMix Prototype

This variant combines the two first-priority clustering mechanisms while keeping
the original scMAE task dominant: warmup mutual-KNN NeighborMix plus delayed
high-confidence DEC prototype alignment.

Default mechanism:

```text
scMAE + NeighborMix + delayed prototype DEC
```

Both auxiliary weights warm up from zero and remain below the scMAE loss weight.


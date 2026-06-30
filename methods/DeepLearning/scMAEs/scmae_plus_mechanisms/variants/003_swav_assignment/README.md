# 003 SwAV Assignment

This variant generates two independently corrupted scMAE views for the same
cell. The original reconstruction and mask-prediction objective remains primary;
a lightweight prototype head adds balanced swapped assignment as a regularizer.

Default mechanism:

```text
scMAE + two-view swapped assignment
```

The SwAV loss starts after warmup and uses Sinkhorn-balanced assignments.


# dec_with_floor

Pure DEC sanity check for the variance-floor hypothesis.

## 2026-07-02 Macosko result

Dataset: `methods/DeepLearning/scMAEs/benchmark_data/Macosko.h5ad`

Runs:

- `runs/universality_test/baseline`: `floor_weight=0.0`, default `tol=0.001`
- `runs/universality_test/with_floor`: `floor_weight=0.02`, default `tol=0.001`
- `runs/universality_test/baseline_forced_tol0`: `floor_weight=0.0`, `tol=0.0`
- `runs/universality_test/with_floor_forced_tol0`: `floor_weight=0.02`, `tol=0.0`

Findings:

- Default early-stop runs converged at epoch 1 and both produced ARI 0.239.
- Forced 300-epoch runs produced identical embeddings and identical metrics:
  ARI 0.457, NMI 0.418, ACC 0.654.
- In the floor run, `Floor: 0.000000` throughout training. Final latent std had min about 2.27 and all 32 dimensions were active.

Conclusion: this pure DEC implementation does not currently reproduce per-dimension variance collapse on Macosko, and the raw `floor_scale=1.0` constraint is inactive. Do not use this as evidence that std-floor generally rescues pure DEC; it is a negative/narrowing result for the universality claim.

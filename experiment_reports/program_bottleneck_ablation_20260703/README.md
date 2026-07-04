# Gene-program bottleneck ablation

External ablation scaffold (lives under `experiment_reports/`, not `methods/`).
Built on the current cross-dataset winner: **scMAE + DEC + per-dim std-floor**
(`AdaptiveSwitch_scMAE`, `variance_weight=0.02`, `force_gate=1.0`).

## Question

```text
std-floor keeps latent dimensions alive (necessary), but does it guarantee those
dims carry cell-type semantics (sufficient)? Does binding part of the latent to
gene-program activity add SEMANTIC ALIGNMENT on top of the std-floor's geometric
non-collapse -- and if so, is it the program SEMANTICS or just a generic aux loss?
```

This is a **minimal falsifiable test**, not the full "semantic multi-view" system.
It only asks whether a gene-program regularizer beats DEC+floor. Posterior fusion
and multi-kernel evidence are deferred until this returns positive.

## Design (Macosko, seeds 42/43/44)

| Arm | program_mode | split_mode | clustering z | program head reads | purpose |
|---|---|---|---|---|---|
| `a0_baseline` | none | none | full 128 | -- | DEC+floor winner (reference) |
| `a1_prog_w02` | nmf | none | full 128 | full 128 | program aux-regularizer, gentle (pw=0.02) |
| `a1_prog_w10` | nmf | none | full 128 | full 128 | program aux-regularizer, ~15% of base (pw=0.10) |
| `a1_shuffle_w10` | shuffled | none | full 128 | full 128 | **control**: row-permuted target, same loss form |
| `a2_fixed_w10` | nmf | fixed | first 96 | last 32 | split latent -> clustering shrinks (reverse-test) |
| `a2_extra_w10` | nmf | extra | full 128 | extra 32 | split but z_type keeps full width (confound control) |

Two `a1` weights bracket the null: if `a1 ~= a0` at both a gentle and a meaningful
weight, "no effect" cannot be dismissed as "lambda too small".

## Key implementation choices (why they matter)

- **Program target = NMF on UNSCALED log-expr** (same source as the recon target),
  then **z-scored per program column** so a few high-activity programs do not
  dominate the MSE and make `program_weight` meaningless.
- **`shuffled` permutes target ROWS** (breaks cell<->program pairing, preserves
  every program's marginal distribution). This is the clean control for
  "is the gain the program semantics or just any structured aux loss?"
  Verified: shuffled `program_r2 ~= 0` (unlearnable), nmf `program_r2 ~= 0.7`.
- **`a1` never splits the latent** -- program is a pure auxiliary head on the full
  128-dim latent, so it does NOT remove clustering capacity (compatible with the
  std-floor mechanism, whose lever is *maximizing* active clustering dims).
- **`a2` (split) is a mechanism reverse-test, not a candidate model.** Prediction:
  `a2_fixed <= a1` (and possibly `<= a0`), because moving active dims out of the
  clustering subspace should weaken the std-floor's fine-grained separation.
  `a2_extra` isolates the "split" variable from the "fewer dims" confound.
- `forward()["latent"]` is always the clustering subspace `z_type`, so the
  std-floor, DEC `student_q`, and exported KMeans embedding all follow the split
  automatically without touching the shared loss/eval code.

## Metrics

Beyond ARI/NMI/ACC and `effective_dim_pr` (std spectrum participation ratio, the
std-floor's own metric), this runner adds:

- **`cluster_aligned_eff_dim`**: participation ratio of the *between-class scatter*
  eigenspectrum (uses true labels). Measures how many dims carry discriminative
  signal, vs `effective_dim_pr` which only measures how many dims are active.
  This is the direct "semantic alignment" probe.
- **`program_r2`**: how well the head predicts the standardized NMF target
  (sanity that the program branch actually trained).

## Run

```bash
bash experiment_reports/program_bottleneck_ablation_20260703/run_all.sh 1   # GPU id
python experiment_reports/program_bottleneck_ablation_20260703/summarize.py
```

## Decision rule (see SUMMARY.md for auto-computed deltas)

```text
a1 > a0 (ARI +>=0.02, stable)   -> program aux helps; proceed to posterior fusion
a1 ~= a0                        -> program is neutral decoration; STOP
a1 < a0                         -> expression-derived program hurts; STOP
a1_shuffle ~= a1                -> gain is generic aux-loss, NOT program semantics
a2_fixed < a1                   -> split conflicts with std-floor clustering capacity
```

Reference DEC+floor numbers on Macosko (from `neighbormix_floor_ablation`, same
backbone): ARI 0.7018 +/- 0.0051, effective_dim_pr ~105.

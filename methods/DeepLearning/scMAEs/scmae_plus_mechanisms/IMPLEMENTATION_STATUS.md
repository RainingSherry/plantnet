# scMAE Plus Mechanisms Status

## Implemented

- Shared protocol code under `common/`:
  - original scMAE backbone wrapper
  - adaptive replacement masking
  - DEC prototype loss and SwAV assignments
  - global mutual-KNN NeighborMix
  - diagnostics for mask rate, embedding collapse, cluster mass, prototypes, and neighbors
- First-stage variants:
  - `001_adaptive_mask`
  - `002_prototype_dec`
  - `003_swav_assignment`
  - `004_neighbormix`
  - `005_neighbormix_prototype`
- Second-stage planned combination:
  - `006_adaptive_mask_neighbormix_prototype`
- NeighborMix stabilization pass:
  - `007_soft_reliable_neighbormix`
  - `008_pseudo_agree_neighbormix`
  - `009_ema_graph_neighbormix`
  - `010_consensus_graph_neighbormix`
  - `011_confidence_adaptive_neighbormix`
- `benchmark.py` with screen default `80` epochs and default seeds `42, 2024, 3407`.
- `collect_results.py` for mechanism result aggregation.
- Legacy independent `run_independent_benchmark.py` screen default changed from `25` to `80`.

## Smoke Verification

The following checks passed on Melanoma_5K with small smoke settings
(`n_top_genes <= 128`, `hidden_size=32`, `epochs=1`, CPU):

- `001_adaptive_mask`
- `002_prototype_dec` with `--prototype_start_epoch 1 --prototype_confidence_threshold 0.0`
- `003_swav_assignment` with `--swav_start_epoch 1`
- `004_neighbormix` with `--neighbor_start_epoch 1`
- `005_neighbormix_prototype` with `--prototype_start_epoch 1 --neighbor_start_epoch 1`
- `006_adaptive_mask_neighbormix_prototype` with `--prototype_start_epoch 1 --neighbor_start_epoch 1`
- `007_soft_reliable_neighbormix` with `--neighbor_start_epoch 1`
- `008_pseudo_agree_neighbormix` with `--neighbor_start_epoch 1`
- `009_ema_graph_neighbormix` with `--neighbor_start_epoch 1`
- `010_consensus_graph_neighbormix` with `--neighbor_start_epoch 1`
- `011_confidence_adaptive_neighbormix` with `--neighbor_start_epoch 1 --neighbor_update_interval 1`
- `benchmark.py --stage screen` dispatch smoke for `001_adaptive_mask`
- `benchmark.py --stage collect`

The current mechanism CSV files are smoke artifacts only. They are not official
80-epoch, three-seed quick-screen results.

## Melanoma_5K 80-Epoch Screen

Completed on 2026-06-30 with `seeds = 42, 2024, 3407`,
`epochs = 80`, `n_top_genes = 1000`, `batch_size = 128`, GPU 5.
All 15 first-stage ablation runs finished successfully. Raw run artifacts were
written to:

```text
/tmp/scmae_plus_screen_80_20260630
```

Reference thresholds from `全benchmark结果.csv`:

```text
scMAE Melanoma_5K ARI mean = 0.668029
scMAE + NeighborMix Melanoma_5K ARI mean = 0.710067
```

Screen summary:

| variant | ARI mean | ARI std | interpretation |
|---|---:|---:|---|
| `004_neighbormix` | 0.697623 | 0.057297 | best mean, but seed variance is high |
| `003_swav_assignment` | 0.668500 | 0.005069 | barely above scMAE; too marginal |
| `001_adaptive_mask` | 0.668053 | 0.001087 | essentially tied with scMAE |
| `002_prototype_dec` | 0.664245 | 0.003041 | below scMAE |
| `005_neighbormix_prototype` | 0.655135 | 0.017876 | below scMAE; prototype hurts NeighborMix |

No result has been appended to `全benchmark结果.csv`. `004_neighbormix`
is the only mechanism with a meaningful mean gain, but it needs stability work
before formal three-dataset promotion.

## NeighborMix Stability Follow-Up

Two follow-up checks were added under the same screen output root:

| variant | ARI mean | ARI std | interpretation |
|---|---:|---:|---|
| `004_neighbormix_nm_a095_w01` | 0.656026 | 0.017969 | weaker mixing is more stable but loses the gain |
| `004_neighbormix_nm_multiavg` | 0.662329 | 0.004329 | averaging all mutual neighbors over-smooths and drops below scMAE |

Implementation note: `common/neighbormix.py` now supports both
`--neighbor_mix_mode first` and `--neighbor_mix_mode mean`; the default is
`first` because this is the only NeighborMix form with a meaningful mean gain on
Melanoma_5K. The screen shows that averaging all mutual neighbors is stable but
over-smooths. The next refinement should therefore gate or score the first
reliable neighbor instead of averaging all reliable neighbors indiscriminately.

Current result files contain 21 successful runs:

```text
机制快筛单次结果.csv: 21 rows
机制快筛汇总结果.csv: 7 rows
机制尝试记录.csv: 21 rows
```

Additional stabilization checks were then completed, bringing the current
mechanism screen to 33 successful runs:

| variant | ARI mean | ARI std | interpretation |
|---|---:|---:|---|
| `004_neighbormix_nm_score078` | 0.666393 | 0.003846 | confidence gating stabilizes but falls just below scMAE |
| `004_neighbormix_nm_score082` | 0.631748 | 0.112195 | stronger gating is unstable and harmful |
| `004_neighbormix_nm_adaptive_mask` | 0.653446 | 0.023074 | adaptive mask hurts NeighborMix |
| `004_neighbormix_nm_swav` | 0.654409 | 0.013235 | SwAV regularization hurts NeighborMix |

Current result files:

```text
机制快筛单次结果.csv: 33 rows
机制快筛汇总结果.csv: 11 rows
机制尝试记录.csv: 33 rows
```

Conclusion so far: the only meaningful mean improvement remains
`004_neighbormix` with first reliable neighbor mixing. Its mean ARI is above
scMAE, but the std is too large for formal promotion. The stabilizers tested so
far reduce or erase the gain.

## Second-Batch Combination Check

The planned combination variant `006_adaptive_mask_neighbormix_prototype`
completed the same Melanoma_5K three-seed, 80-epoch screen:

| variant | ARI mean | ARI std | ACC mean | interpretation |
|---|---:|---:|---:|---|
| `006_adaptive_mask_neighbormix_prototype` | 0.654841 | 0.016507 | 0.724426 | below scMAE; adaptive mask plus DEC prototype suppresses the standalone NeighborMix gain |

Seed-level ARI:

```text
seed 42:   0.661599
seed 2024: 0.666897
seed 3407: 0.636028
```

This does not meet the quick-screen promotion threshold
`ARI_mean > 0.668`. It has not been appended to `全benchmark结果.csv`.
The result files now contain:

```text
机制快筛单次结果.csv: 36 rows
机制快筛汇总结果.csv: 12 rows
机制尝试记录.csv: 36 rows
```

Updated conclusion: the broad three-way combination is not the answer in its
current form. The useful signal is still reliable NeighborMix itself; adaptive
masking and DEC prototypes need either tighter confidence gating or a different
schedule before being recombined.

## Soft Reliable NeighborMix Check

`007_soft_reliable_neighbormix` tests a middle ground between fixed first-neighbor
mixing and averaging multiple mutual neighbors. It keeps only the first reliable
neighbor, but scales the neighbor contribution by edge reliability:

```text
beta_ij = (1 - mix_alpha) * reliability_ij ** neighbor_soft_power
```

Two 80-epoch, three-seed Melanoma_5K settings were completed:

| variant | ARI mean | ARI std | ACC mean | interpretation |
|---|---:|---:|---:|---|
| `007_soft_reliable_neighbormix` | 0.663174 | 0.004279 | 0.736982 | stable but below scMAE; default soft mix over-smooths |
| `007_soft_reliable_neighbormix_a088` | 0.663165 | 0.002333 | 0.737056 | even more stable but still below scMAE |

This confirms that reliability-weighted soft mixing reduces variance, but it
also removes the large positive seed seen in `004_neighbormix`. The next
NeighborMix refinement should not only soften edge strength; it should improve
edge selection, for example by pseudo-cluster agreement, boundary-aware
exclusion, or updating the graph from a more stable teacher embedding.

The result files now contain:

```text
机制快筛单次结果.csv: 42 rows
机制快筛汇总结果.csv: 14 rows
机制尝试记录.csv: 42 rows
```

## Pseudo-Cluster Edge Selection Check

`008_pseudo_agree_neighbormix` changes edge selection rather than edge strength.
At each NeighborMix graph refresh, the current scMAE embedding is clustered with
KMeans and mutual-KNN edges are kept only when the two cells agree in pseudo
cluster. A stricter `confq25` run additionally requires both cells to be above
the 25th percentile of pseudo-cluster confidence.

Three related 80-epoch, three-seed Melanoma_5K checks were completed:

| variant | ARI mean | ARI std | ACC mean | interpretation |
|---|---:|---:|---:|---|
| `008_pseudo_agree_neighbormix` | 0.655011 | 0.018839 | 0.726420 | same-cluster filtering is too weak and hurts seed 3407 |
| `008_pseudo_agree_neighbormix_confq25` | 0.662751 | 0.002924 | 0.736539 | confidence filtering stabilizes but remains below scMAE |
| `004_neighbormix_frozen_graph` | 0.664380 | 0.005916 | 0.737130 | freezing the warmup graph stabilizes but removes the high-gain seed |

These results sharpen the failure analysis: the strong mean of the original
`004_neighbormix` comes from dynamic graph updates amplifying a good trajectory
in one seed. Simple stabilizers, including frozen graph, soft edge strength, and
pseudo-cluster agreement, all reduce variance but fall back below the scMAE
threshold. The next plausible mechanism should be a teacher/EMA or consensus
graph that keeps beneficial graph refinement while avoiding noisy self-feedback.

The result files now contain:

```text
机制快筛单次结果.csv: 51 rows
机制快筛汇总结果.csv: 17 rows
机制尝试记录.csv: 51 rows
```

## EMA Graph NeighborMix Check

`009_ema_graph_neighbormix` smooths the embedding used for KNN graph refreshes:

```text
z_ema(t) = decay * z_ema(t-1) + (1 - decay) * z_current(t)
```

Two 80-epoch, three-seed Melanoma_5K settings were completed:

| variant | ARI mean | ARI std | ACC mean | interpretation |
|---|---:|---:|---:|---|
| `009_ema_graph_neighbormix` | 0.663759 | 0.003946 | 0.737056 | stable but below scMAE; high decay behaves like a smoothed/frozen graph |
| `009_ema_graph_neighbormix_d02` | 0.661184 | 0.009852 | 0.730999 | more responsive EMA is noisier and still below scMAE |

This confirms that single-trajectory embedding smoothing is not sufficient.
It reduces the variance of `004_neighbormix`, but it also removes the positive
dynamic feedback that produced the best seed. The next mechanism should move
from a single EMA trajectory to a consensus or teacher graph, for example by
requiring edges to persist across multiple recent graphs or by using a separate
teacher embedding rather than the current student embedding alone.

The result files now contain:

```text
机制快筛单次结果.csv: 57 rows
机制快筛汇总结果.csv: 19 rows
机制尝试记录.csv: 57 rows
```

## Consensus Graph NeighborMix Check

`010_consensus_graph_neighbormix` keeps only edges that persist across recent
graph refreshes. This tests whether multi-update agreement can preserve useful
dynamic graph refinement while rejecting transient noisy edges.

Three 80-epoch, three-seed Melanoma_5K settings were completed:

| variant | ARI mean | ARI std | ACC mean | interpretation |
|---|---:|---:|---:|---|
| `010_consensus_graph_neighbormix` | 0.696944 | 0.056877 | 0.771918 | preserves high mean but remains unstable; high-gain seed shifts from 3407 to 2024 |
| `010_consensus_graph_neighbormix_m3` | 0.666194 | 0.002646 | 0.738090 | stable and close, but still below scMAE |
| `010_consensus_graph_neighbormix_m3_a085` | 0.662422 | 0.001301 | 0.736539 | stronger mixing hurts the stable consensus setting |

This is the most informative NeighborMix stabilization result so far. A loose
2-of-3 consensus keeps the same kind of high-mean behavior as `004_neighbormix`,
but does not solve seed variance. A strict 3-of-3 consensus solves the variance
problem but falls just below the `ARI_mean > 0.668` screen threshold. The next
candidate should use a confidence-adaptive consensus policy, not a fixed
hit-count rule: permissive for high-confidence core cells and stricter for
boundary cells.

The result files now contain:

```text
机制快筛单次结果.csv: 66 rows
机制快筛汇总结果.csv: 22 rows
机制尝试记录.csv: 66 rows
```

## Confidence-Adaptive NeighborMix Check

`011_confidence_adaptive_neighbormix` replaces the fixed hit-count rule with an
adaptive edge policy: high-reliability 2-of-3 consensus edges are allowed, while
lower-confidence edges must satisfy strict 3-of-3 consensus. The default then
uses first-neighbor mixing with the same `mix_alpha = 0.90` strength as
`004_neighbormix`.

Three 80-epoch, three-seed Melanoma_5K checks were completed:

| variant | ARI mean | ARI std | ACC mean | interpretation |
|---|---:|---:|---:|---|
| `011_confidence_adaptive_neighbormix` | 0.665941 | 0.001270 | 0.737351 | stable but below scMAE |
| `011_confidence_adaptive_neighbormix_first_t084` | 0.666288 | 0.001824 | 0.737942 | best 011 run, still below scMAE |
| `011_confidence_adaptive_neighbormix_first_t080` | 0.665406 | 0.002917 | 0.736982 | looser reliability threshold does not recover the gain |

Conclusion: confidence-adaptive consensus improves stability, but it still
behaves like the strict consensus and soft-reliable variants: it removes the
high-variance positive seed rather than making that gain reproducible.

## NeighborMix Start-Epoch Check

Because the successful but unstable `004_neighbormix` may be caused by early
graph self-feedback, a small start-epoch sweep was run with the original
first-neighbor mixing rule:

| variant | ARI mean | ARI std | ACC mean | interpretation |
|---|---:|---:|---:|---|
| `004_neighbormix_start25` | 0.664278 | 0.007468 | 0.738681 | still unstable and below scMAE |
| `004_neighbormix_start29` | 0.659932 | 0.001483 | 0.735357 | poor; update schedule alignment matters |
| `004_neighbormix_start30` | 0.667711 | 0.003425 | 0.743408 | closest stable result, but just below 0.668029 baseline |
| `004_neighbormix_start35` | 0.665808 | 0.001381 | 0.740084 | stable but below scMAE |
| `004_neighbormix_start40` | 0.663076 | 0.005763 | 0.736465 | too late; loses the useful graph signal |

This sweep is important even though it did not produce a promoted result.
Delaying NeighborMix to epoch 30 almost solves the variance problem without new
architecture, missing the Melanoma_5K scMAE ARI mean by about `0.00032`. The
next useful search should combine delayed graph activation with a mechanism
that protects boundary cells, rather than adding DEC/SwAV/adaptive masking.

The result files now contain:

```text
机制快筛单次结果.csv: 90 rows
机制快筛汇总结果.csv: 30 rows
机制尝试记录.csv: 90 rows
```

## Next Step

Do a targeted second pass around NeighborMix, not broad architecture expansion:

```bash
python methods/DeepLearning/scMAEs/scmae_plus_mechanisms/benchmark.py \
  --stage screen \
  --variants 004_neighbormix \
  --seeds 42 2024 3407 \
  --screen_epochs 80 \
  --batch_size 128 \
  --n_top_genes 1000 \
  --gpu 1 \
  --rerun_existing
```

Recommended refinements:

- Use `004_neighbormix_start30` as the closest stable anchor, not the unstable
  default `004_neighbormix`.
- Tune first-neighbor reliability and strength around this anchor:
  `mix_alpha`, `mix_weight`, `neighbor_k`, and graph update interval.
- Add confidence gating for NeighborMix edges before mixing; do not average all
  mutual neighbors by default.
- The tested hard score gates (`0.78`, `0.82`), soft edge weighting,
  pseudo-cluster agreement, frozen graph, and single-trajectory EMA graph are
  not sufficient; multi-update consensus and confidence-adaptive consensus are
  also stable but below threshold.
- Add an explicit boundary/rare-cell protection signal next. The current edge
  filters only use graph reliability; they do not know whether a cell is near a
  KMeans boundary or belongs to a small pseudo-cluster.
- Do not add the current DEC prototype to NeighborMix as-is; it reduced mean ARI.
- Do not combine SwAV or adaptive mask with NeighborMix in their current form;
  both reduced Melanoma_5K ARI in the three-seed screen.
- Promotion still requires Melanoma_5K three-seed `ARI_mean > 0.668`, stable
  diagnostics, and no direct append to `全benchmark结果.csv`.

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

- Tune first-neighbor reliability and strength first: `mix_alpha`,
  `mix_weight`, `neighbor_k`, and `neighbor_start_epoch`.
- Add confidence gating for NeighborMix edges before mixing; do not average all
  mutual neighbors by default.
- The tested hard score gates (`0.78`, `0.82`), soft edge weighting,
  pseudo-cluster agreement, and frozen graph are not sufficient; future work
  should use a teacher/EMA or consensus graph rather than relying on a single
  current-embedding graph.
- Do not add the current DEC prototype to NeighborMix as-is; it reduced mean ARI.
- Do not combine SwAV or adaptive mask with NeighborMix in their current form;
  both reduced Melanoma_5K ARI in the three-seed screen.
- Promotion still requires Melanoma_5K three-seed `ARI_mean > 0.668`, stable
  diagnostics, and no direct append to `全benchmark结果.csv`.

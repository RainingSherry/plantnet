# scMAE structural pivot — diagnosis and roadmap (2026-07-03)

Status: **planning document**, written after 8+ negative "bolt-on mechanism"
ablations. It records (a) the upper-layer condition every negative result shares,
(b) what that does and does NOT license us to conclude, and (c) a cost-ordered
plan that attacks scMAE's own structure/task instead of hanging new pseudo-
structure on it.

This lives under `experimental_retired_models/` by user request; it is a strategy
note, not runnable benchmark code. Runners it proposes live under
`experiment_reports/` when built.

---

## 1. The shared upper-layer condition

Across GatedNeighborMix, AdaptiveGranularity, ReliableRecon, ModuleAttn,
DualAxisGated, marker-masking, and the gene-program bottleneck, five things were
held FIXED:

1. Input is the same expression matrix X.
2. scMAE main task = corrupt input -> mask prediction + reconstruction of X.
3. The encoder embedding is the sole clustering object.
4. KMeans / DEC act on that same latent; evaluation is KMeans, known-K.
5. Every new mechanism injects X-derived pseudo-structure (neighbors, programs,
   markers, reliability, modules).

So the negative results license ONLY this conclusion:

> With scMAE's structure and reconstruction task unchanged, expression-derived
> external structure mostly pollutes or deflects the embedding. The only
> intervention that helped (DEC + per-dim std-floor) acts directly on the
> geometry / capacity of the clustering space, not on added pseudo-supervision.

They do NOT license "scMAE cannot be improved." They falsify *bolt-ons*, not the
*backbone task*. What was never tested: is the corruption too coarse, is mask
prediction too easy, does the decoder hand the encoder a shortcut, should the
clustering embedding be separated from the reconstruction embedding.

## 2. The decisive new evidence (program-bottleneck ablation, Macosko)

| arm | seed42 ARI | seed43 ARI |
|---|---:|---:|
| a0 DEC+floor (reference) | 0.779 | 0.730 |
| a1 real-NMF program, lambda 0.10 | 0.458 | 0.356 |
| **a1-shuffle (row-permuted target, same loss)** | **0.632** | **0.536** |
| a2 latent split (cluster on z_type only) | 0.326 | 0.305 |

Two clean reads:
- **Real program hurts MORE than a shuffled/meaningless target.** The damage is
  the program *semantics*, not the generic aux loss. Real NMF activity pulls the
  embedding onto non-cell-type axes (library size / cell cycle / batch-like
  structure). This is the sharpest confirmation yet of "X-derived pseudo-
  structure != label structure."
- **Latent split collapses to ~0.31.** Moving active dims out of the clustering
  subspace fights the std-floor lever (whose mechanism is *maximizing* active
  clustering dims). This pre-emptively kills the "dual-subspace z_cls/z_aux" and
  "gene-program bottleneck" proposals as MAIN models.

## 3. The Phase-0 gate — RUN, and the verdict is "ceiling = embedding, NOT head"

The OTHER invariant was the head + protocol: KMeans, known-K, single-dataset
Macosko ARI. We re-clustered FROZEN embeddings with KMeans / GMM / Leiden (Leiden
resolution swept, selected by silhouette label-free). Runner:
`experiment_reports/clustering_head_reeval_20260703/`.

Macosko Leiden resolution sweep (ARI), true k = 12:

| embedding | KMeans knownK | Leiden@k~5 (res .02) | Leiden@k~11-14 (res .2-.3) | dec-vs-floor |
|---|---:|---:|---:|---|
| pure-backbone | 0.384 | 0.834 (k5) | 0.223-0.229 | — |
| dec (no floor) | 0.342 | 0.346 (k9) | 0.147-0.164 | worst everywhere |
| dec_floor | 0.696 | 0.780 (k6) | 0.265-0.276 | best everywhere |
| neighbormix_dec_floor | 0.473 | 0.592 (k8) | ~0.28 | below dec_floor |

VERDICT (careful reading overturns the naive one):
1. **The head is NOT a hidden lever.** The eye-catching backbone "0.834" is at
   k=5 for a 12-class set = a COARSE-partition + imbalance artifact (Macosko is
   ~80% rods; a 5-way split that nails the big class inflates ARI). At MATCHED
   granularity (k~11-14) backbone Leiden = 0.22, far below its own KMeans-known-K
   (0.384). At matched k, **KMeans-known-K is the BEST head** (dec_floor 0.696 vs
   Leiden-k14 0.276). Given the true K, KMeans-known-K is already near-optimal for
   these embeddings; Leiden only "wins" via the artifact.
2. This **refutes the "maybe it's the clustering head" hypothesis** cleanly, with
   data. By elimination it **strengthens "ceiling = embedding/task"** -> attacking
   scMAE corruption/decoder (Phases 1-2) is the justified direction.
3. **std-floor survives head-independently**: `dec` (no floor) is worst at EVERY
   resolution and EVERY head; `dec_floor` > `dec` everywhere. Variance collapse is
   real embedding damage, not a readout effect.
4. Real nuance to KEEP: the pure-backbone embedding carries strong HIERARCHICAL /
   coarse structure (NMI 0.74 at k5-7, the highest NMI anywhere) that known-K
   KMeans cannot exploit because forcing 12 balanced clusters splits the dominant
   class (same failure mode as Quake in the memory). DEC+floor trades some coarse
   structure for known-K-KMeans robustness. Where the structure LIVES differs by
   embedding; the readout decides which looks best. For the FINE-GRAINED 12-class
   goal at matched k, dec_floor still wins.

METHODOLOGICAL CATCH (now a hard rule, see Section 4): **silhouette-over-resolution
is a BROKEN label-free selector on imbalanced scRNA** — it degenerates to the
coarsest partition, which is confounded with imbalance-inflated ARI. Do not use it
for model selection. Use stability-based selection or fixed/known-K comparison.

## 4. Leakage discipline — LOCK before any "chase good clustering" work

The moment the goal becomes "excellent clustering numbers," the subtlest leak is
tuning on the test labels. Rules, fixed up front:

1. No hyperparameter (resolution, temperature, loss weight, epoch) chosen by test
   ARI/NMI.
2. Model/hyperparameter selection uses ONLY label-free signals — but NOT raw
   silhouette-over-resolution: Phase 0 showed it degenerates to the coarsest
   partition, confounded with imbalance-inflated ARI (backbone silhouette picked
   k=5 -> ARI 0.834, an artifact). Prefer seed/subsample **stability** at a fixed
   granularity, assignment entropy, cluster balance, eff-dim. If comparing across
   K, always report at MATCHED cluster count and report NMI beside ARI.
3. Labels are touched exactly once, for final reporting.
4. Report known-K and inferred-K separately; never silently use known K to pick K.
5. Fixed multi-dataset protocol (Macosko + Quake + Melanoma at minimum, plant sets
   next). No single-dataset tuning — Macosko's seed sd is +/-0.087; single-set
   gains inside that band are noise.

---

## 5. The plan — cost-ordered, each step falsifiable

Ordered by (information per GPU-hour), NOT by radicalness. Every step is measured
against DEC+floor on the SAME harness, multi-seed, with the Section 4 discipline.

### Phase 0 — close the diagnostic loop  [cheapest, do first]
- 0a. Finish the program ablation (seed 44 + a2_extra) and `summarize.py`. Confirm
  "real-program < shuffle < baseline" holds across 3 seeds. (in flight)
- 0b. Head re-eval on pure-backbone vs DEC vs DEC+floor vs NeighborMix+floor.
  Answers the Section 3 gate. (in flight)
- Deliverable: a one-paragraph verdict "ceiling = embedding" or "= head/protocol".
  Everything below is conditional on "= embedding".

### Phase 1 — decoder bypass  [DONE 2026-07-03, `experiment_reports/decoder_bypass_ablation_20260703/`]
Hypothesis: scMAE's `decoder(concat[latent, mask_logits])` lets the encoder off
the hook — the mask vector is a G-dim side channel, so the latent need not encode
what was corrupted. Good for reconstruction, possibly bad for a clustering latent.
- D0 `decoder(latent, mask_logits)`  (original)
- D1 `decoder(latent)`               (mask only in the BCE loss, not the decoder)
- D2 `decoder(latent, lowrank(mask_logits))`  (mask enters only as a 16-dim pool)
3 datasets × 3 seeds each. RESULT (ARI mean±sd):

| dataset | D0 concat | D1 none | D2 lowrank |
|---|---|---|---|
| Macosko | 0.702±0.004 | 0.471±0.155 | 0.627±0.174 |
| Melanoma_5K | 0.648±0.002 | **0.715±0.042** | 0.678±0.042 |
| Quake | 0.920±0.001 | 0.921±0.001 | 0.920±0.001 |

VERDICT — the mask→decoder path is primarily a **STABILIZER, not a pure shortcut**:
- D0 reproduces the winner exactly on all three, tiny variance (faithful runner).
- **Melanoma (hard tumor data): removing the bypass HELPS on average** (none +0.067,
  lowrank +0.030) — the "shortcut" hypothesis holds only here.
- Quake (easy): saturated, flat.
- Macosko (fine-grained): removing the bypass tanks the MEAN and explodes variance
  (sd 0.15-0.17) — BUT the per-seed detail is the real story: **lowrank seed43 =
  0.872, none seed43 = 0.690**. The 0.872 is the HIGHEST single Macosko ARI in the
  entire search (winner 0.70). So a much better embedding is REACHABLE without the
  full mask bypass; training is bimodal (soars or crashes by seed).

So D0's mask input regularizes/stabilizes reconstruction. Not a deployable win, but
the Melanoma mean-gain + the Macosko 0.872 outlier = the strongest constructive lead
to date. Cost was ~1 GPU-hour (6-way parallel dispatch across GPU 1-6).

### Phase 1b — stabilize the high mode  [DONE 2026-07-03]
Mapped the bimodal Macosko behavior: lowrank × mask_rank{4,8,16,32} × seeds42-49 +
none × seeds45-49. RESULT — the ~0.87 basin is reachable but the decoder-bypass
route is UNRELIABLE at every mask_rank:
- mr4 hit-0.8+ 12%, mr8 25%, mr16 25%, mr32 25% (mean 0.53-0.59, full spread 0.34-0.87)
- none: 0/5 seeds, all crash to ~0.355
So the Phase-1 lowrank seed43=0.872 was a ~1-in-8 lucky draw, NOT a stable method.
CONTRAST with Phase 2 swap_lib/swap_ndet = 100% of seeds hit ~0.867 (sd<0.003).
=> The SAME ~0.87 Macosko basin is reached RELIABLY by conditional-shuffle and only
OCCASIONALLY by decoder-bypass. Decoder-bypass is a dead end for reliability; the
corruption-donor route is the real lever into that basin.

### Phase 2 — conditional / nuisance-matched shuffle  [most fundamental; the real "fix scMAE" bet]
Hypothesis: scMAE's swap-noise draws the donor value from ALL cells for that gene,
so "is this value corrupted?" is often solvable via library size / detection /
zero-rate — technical axes, not cell-type ones. Restrict the donor pool to cells
in the same nuisance bin (library size, n_detected, zero-rate; batch if present)
so corruption detection must use finer gene-gene conditional structure.
- S0 original global per-gene shuffle
- S1 library-size-matched shuffle
- S2 detected-genes-matched shuffle
- S3 zero-rate/library-matched shuffle
DONE 2026-07-03 (`experiment_reports/conditional_shuffle_ablation_20260703/`, full
results in FINDINGS.md). RESULT (ARI mean±sd, 3 seeds; Δ vs same-ds zero-mask winner):

| arm | Macosko (win 0.702) | Melanoma (0.648) | Quake (0.920) |
|---|---|---|---|
| swap_global S0 | 0.299±0.028 (−0.40) | 0.655±0.005 (+0.01) | 0.923 (flat) |
| swap_lib S1 | **0.867±0.002 (+0.16)** | 0.602±0.076 (−0.05) | 0.915 (flat) |
| swap_ndet S2 | **0.864±0.001 (+0.16)** | 0.658±0.006 (+0.01) | 0.922 (flat) |
| swap_zerolib S3 | 0.759±0.147 (bimodal) | 0.660±0.003 (+0.01) | 0.917 (flat) |

FINDINGS:
1. Original scMAE GLOBAL swap-noise is CATASTROPHIC on fine-grained Macosko (0.30 vs
   zero-mask 0.70) — clean negative on the original corruption design; fine-grained
   data needs zero-masking, not global swap.
2. nuisance-matched swap (S1/S2) = Macosko +0.16 STABLE (all seeds ~0.867), and
   survives BOTH KMeans-knownK (0.870 vs 0.696) AND label-free Leiden (0.894 vs
   0.780) → genuine embedding improvement, not a KMeans-known-K artifact.
3. Does NOT generalize: Melanoma/Quake flat. Macosko-specific.
4. Obvious mechanism REFUTED: NMI(nuisance_bin, label) Macosko 0.11/0.13 ≈ Melanoma
   0.10/0.14, yet only Macosko benefits. "same-bin donor ≈ same-type donor" is wrong.
   Why Macosko-specific = UNKNOWN (open puzzle).
Disciplined verdict: NOT a general win (violates multi-dataset rule), NOT the paper
headline. DEC+std-floor stays the robust spine. But the reliable ~0.87 Macosko basin
(also glimpsed via decoder-bypass) is a real phenomenon worth understanding.

### Phase 2b — PLANT data conditional-shuffle  [DONE 2026-07-03, decisive NEGATIVE]
Project-relevant test. RESULT (ARI mean, 3 seeds; Δ vs zero-mask):

| arm | SRP182008 k=15 (zero 0.390) | CRA002977_1 k=7 (zero 0.683) |
|---|---|---|
| swap_global | 0.397 (+0.01) | 0.508 (−0.175) |
| swap_lib | 0.389 (−0.00) | 0.537 (−0.146) |
| swap_ndet | 0.384 (−0.01) | 0.544 (−0.139) |

**Macosko +0.16 does NOT transfer.** SRP182008 flat (like Melanoma); CRA002977_1 all
swaps HURT ~0.15 (zero-mask clearly best). Across 5 datasets: only Macosko benefits;
all others flat-or-worse. => The "attack corruption" branch does NOT help the project's
plant data. **zero-mask + DEC + std-floor is the robust choice everywhere, plant included.**

## 5b. Branch verdict (attack-corruption line CLOSED)

corruption design has large but DATASET-DEPENDENT effects; NO single corruption beats
zero-mask everywhere. "matched > global" ordering is consistent, but "swap beats
zero-mask" is Macosko-only. Macosko's uniqueness = unexplained (nuisance≈type refuted).
DEC + per-dim std-floor (zero-mask) remains the robust cross-dataset + plant spine.
Clean science kept: (a) original scMAE global swap-noise is worse than zero-mask on
fine-grained data; (b) the Macosko ~0.87 basin is real (reached by 2 perturbations,
survives label-free Leiden) but isolated and mechanism-unknown — a future-work note.

### Phase 3 — cluster-first / generative-clustering objective  [highest variance; the pivot, if Phases 1-2 stall]
Only if Phases 1-2 do not clear DEC+floor. Replace reconstruction-primary training
with a clustering-primary objective on a single latent (NO subspace split — a2
already showed split collapses):
- Option A (self-labelling): SwAV/DINO-style — prototypes + swapped/EMA-teacher
  assignment prediction, biological augmentations (dropout simulation, count
  down-sampling, gene masking). Keep std-floor as anti-collapse.
- Option B (generative-clustering): make the cluster variable EXPLAIN expression —
  `x_hat = sum_k q_k * mu_k + R(z)` (mixture-prototype decoder), so assignment is
  a first-class generative factor, not a post-hoc pull on the latent.
Risk: this is the neighborhood of scDCC / scDeepCluster / contrastive-scRNA — the
novelty must be the DIAGNOSIS (why reconstruction-first fails fine-grained
clustering), with cluster-first as the constructive proof, not a "new" method
claim on its own. Cost: high; treat as a research fork, not a quick ablation.

## 6. Explicitly de-prioritised (and why)

- Dual-subspace z_cls/z_aux, gene-program bottleneck as main model: a2 split ->
  ~0.31. Removing clustering capacity fights the std-floor lever. Dead.
- More complex attention / gene-token transformer / graph encoder / new
  NeighborMix / multi-kernel main model: all re-enter the "derive structure from
  X, assume it aligns with cell type" trap the negatives already closed.
- Contrastive/cluster-first as a standalone method-novelty claim: likely
  re-derives published work; only defensible as the sequel to the diagnosis.

## 7. Paper framing (unchanged spine, extended)

- Main line (mature, robust, multi-seed): DEC on fine-grained scRNA suffers
  PER-DIM VARIANCE COLLAPSE; a per-dim std-floor (VICReg hinge, w=0.02) is the
  exact antidote; eff-dim 60 -> 114 on Macosko; neutral on Quake/Melanoma.
- Diagnostic extension (this search): reconstruction-fidelity, dimension-liveness,
  and neighbor/program pseudo-structure are each NECESSARY-BUT-NOT-SUFFICIENT or
  actively harmful for fine-grained clustering — shown with clean same-harness
  controls incl. the shuffle control (real program < shuffled target < baseline).
- Constructive sequel (Phases 1-3, whichever clears the bar): the first
  intervention that improves scMAE from INSIDE its own task (corruption/decoder),
  or a cluster-first objective, motivated by the diagnosis.

## 8. Immediate next actions

1. Read Phase-0 outputs (`clustering_head_reeval_20260703/reeval.csv`, program
   `SUMMARY.md`). Write the Section 3 verdict.
2. If ceiling = embedding: build Phase 1 (decoder-bypass) runner under
   `experiment_reports/decoder_bypass_ablation_YYYYMMDD/`, reusing the
   AdaptiveSwitch backbone + the program-runner scaffold.
3. Only after Phase 1 reads out: decide Phase 2 vs Phase 3.
4. Keep DEC+std-floor as the frozen reference in every comparison.


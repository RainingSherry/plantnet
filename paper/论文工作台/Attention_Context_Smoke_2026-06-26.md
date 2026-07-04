# Attention / Context Smoke Decision

Date: 2026-06-26

Authoring status: development-only smoke memo for research direction. This is not a publication claim.

## 1. Why This Smoke Was Run

Phase 14 failed the AdvMask gate:

```text
gate_result = fail
recommendation = drop_or_downgrade_advmask
```

Therefore, Axial/full should not be run as an AdvMask rescue path. The remaining scientific question from the original TabPFN-inspired idea is narrower:

```text
Does the existing axial/context encoder, without AdvMask, improve over an MLP under the corrected scMAE-style protocol?
```

This smoke tests only that narrow question.

## 2. Protocol

```text
scope = development-only smoke
datasets = Limb_Muscle, Mouse_Pancreas_1, Quake_Smart-seq2_Lung
seed = 42
epochs = 3
corruption_type = scmae_shuffle
mask_selector = random
input_mode = log1p
n_top_genes = 2000
scale_input = false
validation/sealed test = not used
AdvMask/full = not run
output_root = /tmp/caam_attention_context_smoke/dev_20260626
```

Compared roles:

```text
control = MLP hidden_dim 256
axial = existing axial encoder
mlp_parammatched = MLP hidden_dim 291
```

Parameter counts:

```text
control student params = 9,080,688
axial student params = 9,154,544
parameter-matched MLP student params = 9,155,273
relative axial vs parameter-matched MLP gap ~= 0.00008
```

## 3. Results

Primary development metric:

```text
kmeans_known_k.ari
```

| dataset | role | student params | ARI | NMI | ACC | macro-F1 | Leiden ARI |
|---|---:|---:|---:|---:|---:|---:|---:|
| Limb_Muscle | control | 9,080,688 | 0.866751 | 0.859992 | 0.869788 | 0.760256 | 0.419897 |
| Limb_Muscle | axial | 9,154,544 | 0.204241 | 0.330692 | 0.418266 | 0.428134 | 0.114455 |
| Limb_Muscle | mlp_parammatched | 9,155,273 | 0.861273 | 0.858861 | 0.866462 | 0.758125 | 0.395307 |
| Mouse_Pancreas_1 | control | 9,080,688 | 0.365151 | 0.570103 | 0.498409 | 0.434805 | 0.316153 |
| Mouse_Pancreas_1 | axial | 9,154,544 | 0.303162 | 0.482826 | 0.383881 | 0.250695 | 0.163806 |
| Mouse_Pancreas_1 | mlp_parammatched | 9,155,273 | 0.370947 | 0.578567 | 0.500530 | 0.422805 | 0.450371 |
| Quake_Smart-seq2_Lung | control | 9,080,688 | 0.454009 | 0.714620 | 0.505370 | 0.378251 | 0.440553 |
| Quake_Smart-seq2_Lung | axial | 9,154,544 | 0.066389 | 0.170535 | 0.249403 | 0.222460 | 0.054413 |
| Quake_Smart-seq2_Lung | mlp_parammatched | 9,155,273 | 0.521373 | 0.699487 | 0.581742 | 0.316602 | 0.360003 |

ARI deltas:

| dataset | axial - control | axial - parameter-matched MLP |
|---|---:|---:|
| Limb_Muscle | -0.662509 | -0.657031 |
| Mouse_Pancreas_1 | -0.061989 | -0.067784 |
| Quake_Smart-seq2_Lung | -0.387620 | -0.454984 |

Mean ARI:

| role | mean ARI |
|---|---:|
| control | 0.561970 |
| axial | 0.191264 |
| mlp_parammatched | 0.584531 |

## 4. Interpretation

The existing axial/context encoder does not pass even a weak development smoke under the corrected protocol. It underperforms both the standard MLP and the parameter-matched MLP on all three development datasets for seed 42.

This result should not be overclaimed. It does not prove that all attention or TabPFN-like context modeling is useless for scRNA-seq. It does show that the current axial implementation should not be promoted as a main contribution or expanded into a full factorial experiment without a new mechanism-level redesign.

## 5. What This Supports

Supported conservative statements:

```text
1. AdvMask should remain downgraded after Phase 14.
2. The existing axial/context encoder is not currently a rescue path.
3. Parameter-matched controls are essential; the parameter-matched MLP is at least as strong as the standard MLP in this smoke.
4. The current paper route should remain protocol_analysis / diagnostic unless a new independently justified module is designed.
```

## 6. What This Does Not Support

Forbidden claims:

```text
1. Axial improves clustering.
2. TabPFN-like attention is validated by current CAAM evidence.
3. Axial + AdvMask synergy is worth testing with the current modules.
4. The smoke result is publication-level evidence.
```

## 7. Research Decision

Current decision:

```text
drop_or_downgrade_current_axial = true
do_not_run_full_factorial = true
continue_protocol_analysis_route = true
```

Recommended next action:

```text
Write the Phase 16 publication decision document.
If attention is revisited, it should be a new design question, not a continuation of the current Axial implementation.
```

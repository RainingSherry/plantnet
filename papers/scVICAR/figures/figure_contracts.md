# Figure contracts

## Figure 1 — Method

- **Core conclusion:** scVICAR trains one masked autoencoder to recover an original cell from a bounded graph-vicinal view, with either fixed or topology-adaptive perturbation.
- **Archetype:** schematic-led composite.
- **Backend:** Python/matplotlib only.
- **Final size:** 183 mm double-column, approximately 183 × 105 mm.
- **Panel map:** (a) scMAE corruption; (b) fixed graph-vicinal anchor recovery; (c) topology-informed affinity and cell-wise perturbation budget; (d) objective and bounded diffusion identity.
- **Reviewer risk:** do not depict a dynamic graph, learned uncertainty or interpolated reconstruction target.

## Figure 2 — Confirmatory clustering

- **Core conclusion:** matched-backbone fixed and adaptive vicinal recovery are evaluated on six frozen external datasets without per-dataset tuning.
- **Archetype:** quantitative grid.
- **Hero evidence:** paired dataset-level change in ARI relative to NoMix.
- **Supporting evidence:** NMI/macro-F1, win–tie–loss and runtime.
- **Statistics:** dataset is the independent unit; three model seeds are averaged first; hierarchical bootstrap interval and Holm-adjusted paired permutation test.

## Figure 3 — Components

- **Core conclusion:** edge affinity and node gating have separable and potentially interacting effects.
- **Archetype:** quantitative grid.
- **Panels:** six-variant ablation forest plot; interaction plot; embedding-geometry diagnostic.
- **Reviewer risk:** no component is credited unless the matched ablation supports it.

## Figure 4 — Graph contamination

- **Core conclusion:** topology adaptation is tested for graceful degradation as cross-class edges are injected.
- **Archetype:** asymmetric mixed-modality.
- **Hero evidence:** ARI versus contamination rate for F and T.
- **Supporting evidence:** affinity AUROC/Spearman, gate distribution and rare-class recall.
- **Reviewer risk:** labels are used only to inject/measure the stress test, not for training.

## Figure 5 — Downstream utility

- **Core conclusion:** biological utility is assessed independently of clustering ARI through marker recovery, marker-overlap annotation and low-label probes.
- **Archetype:** quantitative grid.
- **Panels:** Recovery@100; annotation macro-F1 and coverage; 10%/30% probe macro-F1.
- **Reviewer risk:** call the probe transductive, not cross-dataset transfer.

## Figure 6 — Pancreas case study

- **Core conclusion:** fixed and adaptive vicinal recovery can be inspected through cell-type markers rather than clustering metrics alone.
- **Archetype:** asymmetric mixed-modality.
- **Panels:** marker overlap heatmap, dot plot and type→cluster→annotation→gold Sankey.
- **Frozen display run:** Human_Pancreas_3, model seed 42, downstream split seed 11,
  shown for NoMix, scVICAR-F and scVICAR-T; this is not a best-seed choice.
- **Reviewer risk:** the dataset and marker panel are fixed before inspecting scVICAR results.
- **Fail-closed rule:** generation requires `verification_status=verified` in the
  frozen marker-panel JSON and the complete case-study downstream artifacts.
  The panel is now verified; incomplete downstream evidence still blocks export.

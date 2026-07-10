# Figure 1 QA — method framework

- Intended conclusion: scVICAR perturbs the training input toward a local graph
  neighborhood while recovering the original anchor; F fixes the budget and T
  separates topology-informed affinity from a cell-wise budget.
- Evidence logic: panel a establishes the scMAE anchor target; panels b/c isolate
  fixed and adaptive local corruption; panel d supplies their unified operator,
  diffusion identity, and objective.
- Visual inspection: passed at 1801×1082 PNG preview. Panel headings, equations,
  arrows, and target distinction are legible without zoom; no observed overlap
  or clipped label.
- Editable export: SVG preserves text as text (`svg.fonttype=none`); PDF uses
  TrueType font embedding (`pdf.fonttype=42`).
- Raster export: PNG 300 dpi and uncompressed TIFF 600 dpi are generated from
  the same Python source. The large TIFF is intentionally excluded from Git and
  retained in the remote figure archive.
- Color: anchor blue, neighbors green, topology affinity purple, adaptive gate
  red. Meaning is also encoded through labels and line width, not color alone.
- Review risks checked: the figure does not depict dynamic graph updates,
  learned uncertainty, cluster-guided gates, or mixed targets. The symbol
  `\widehat P` denotes the realized stochastic kernel, consistent with the
  Methods text.
- Source: `papers/scVICAR/code/figures.py::plot_method`.

# VarFloor-scMAE

Formal benchmark runner for the method previously tested as zero-mask scMAE +
DEC + per-dimension std-floor.

The method uses labels only for final fixed-K benchmark evaluation. Training is
label-free except for the benchmark-provided value of `--n_clusters`, matching
the known-K protocol used by the existing benchmark table.

Core definition:

- zero-mask scMAE reconstruction objective
- DEC sharpened target after warmup
- per-dimension standard-deviation floor with weight `0.02`
- KMeans known-K evaluation on the exported latent embedding

This runner intentionally omits the adaptive switch/soft branch used in earlier
diagnostic experiments. The promoted method is the fixed sharp path
(`force_gate=1`) with the std-floor enabled.


# CAAM-scMAE Design

This package implements the BDD-defined 2x2 internal design:

- control: MLP + random fixed-budget mask
- axial: Axial encoder + random fixed-budget mask
- advmask: MLP + adversarial mask selector
- full: Axial encoder + adversarial mask selector

The formal benchmark method is only `caam_scmae`, mapped to `--variant full`.


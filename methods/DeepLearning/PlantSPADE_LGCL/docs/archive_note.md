# Archive / Exploratory Code

The main paper protocol is under `methods/DeepLearning/PlantSPADE_LGCL`.

The following code is retained for history or negative/exploratory results, but is excluded from main runners, main tables, and the README main flow:

- `methods/DeepLearning/PlantSPADE`
- `methods/DeepLearning/_archive`
- root-level GraphDiffusion/DOLORIS helper scripts such as `infer_graphdiffusion.py`, `infer_graphdiffusion_ckpt.py`, and `launch_doloris_eval.sh`
- unstable or dependency-heavy GNN entries such as `scCDCG`, `scDSC`, `scGNN`
- foundation model experiments such as `scGPT`, `GeneFormer`, `GeneCompass`, and `scPlantLLM`

These files should not be deleted unless a separate archival cleanup task explicitly requests removal. New paper-facing experiments should be added through the PlantSPADE-LGCL config and runner stack.

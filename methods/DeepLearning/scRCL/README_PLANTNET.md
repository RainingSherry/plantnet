# scRCL PlantNet Adapter

This directory vendors the full upstream THPengL/scRCL implementation and adds
a thin PlantNet benchmark wrapper.

## Upstream provenance

- Repository: https://github.com/THPengL/scRCL
- Vendored commit: `6da37870445a391953a3c85970a5ccf53c8e4b9a`
- Vendored tree: `methods/DeepLearning/scRCL/scRCL_upstream/`
- Commit record: `methods/DeepLearning/scRCL/scRCL_upstream_commit.txt`

The upstream source is kept as a complete source snapshot without its `.git`
metadata. The benchmark wrapper does not reimplement the model; it imports and
uses upstream `load_h5_data`, `Model`, scRCL loss methods, fusion utilities and
clustering.

Two small runtime patches are applied inside the vendored snapshot and should be
treated as adapter patches rather than model replacements:

- `clustering.py`: compare `torch.device` values with equality so CPU/GPU
  routing works reliably.
- `model.py`: add optional chunked forms of the two cell-cell loss matrices via
  `loss_chunk_size`, preserving the same loss equations while avoiding H100 OOM
  on larger benchmark datasets.

## Benchmark wrapper

`methods/DeepLearning/scRCL/run.py` adapts PlantNet's standard method contract:

```bash
python methods/DeepLearning/scRCL/run.py \
  --data_path result/scmae_all_methods_20260705_full/converted_data/Pollen.h5ad \
  --save_dir /tmp/scrcl_pollen_smoke \
  --n_clusters 11 \
  --seed 42 \
  --epochs 2 \
  --gpu 1
```

The wrapper temporarily writes an AnnData file with the `cell_type1` categorical
column expected by upstream `load_h5_data()`.  Ground-truth labels are not used
for model selection during training.  The final epoch embedding is clustered and
then evaluated once through the repository-standard `methods/utils.py::save`.

## Runtime

scRCL requires a runtime with PyTorch, PyTorch Geometric, Scanpy, UMAP and POT
(`ot`).  The suite registry uses the dedicated runtime name `plantnet-scrcl` so
the existing benchmark environments are not modified.

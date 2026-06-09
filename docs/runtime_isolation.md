# Runtime Isolation

This document describes the per-method runtime isolation architecture used in this benchmark project.

---

## Why Per-Method Runtime Isolation?

### Problem: Dependency Conflicts

Different models in the scCluBench benchmark require fundamentally different library versions that **cannot coexist in a single conda environment**:

| Model | Key Dependency | Version |
|-------|---------------|---------|
| scNAME, scziDesk, scDeepCluster | TensorFlow | 2.x |
| DESC | TensorFlow | 2.x |
| scCDCG | PyTorch | 1.12 (older) |
| scGNN | PyTorch + rpy2 | (R/LTMG mode) |
| scMAE, scDCC, scDSC | PyTorch + Scanpy | (main env) |

TensorFlow 2.x and PyTorch 1.12 cannot be installed in the same environment without conflicts. TensorFlow 2.x and modern Scanpy also have complex dependency chains.

### Why Not One Big Environment?

1. **Dependency hell**: TensorFlow 2.x + PyTorch 1.12 + Scanpy + rpy2 + h5py + scipy creates unsolvable version conflicts.
2. **Floating-point reproducibility**: Mixing TF and PyTorch in the same process can cause non-deterministic behavior.
3. **Isolation for audit**: Each model's runtime is independently auditable without cross-contamination.

### Why Not Direct OtherMode Import?

1. **Provenance**: Models must run from `methods/`, not from `OtherMode/`. Using `OtherMode/` at runtime means the benchmark is not self-contained.
2. **Audit trail**: `environment.json` must record exactly which Python and packages were used.
3. **Reproducibility**: Results must be reproducible without `OtherMode/` being available.

---

## Architecture

### Runtime Registry

The mapping from logical runtime names to actual executables is stored in `envs/runtime_registry.yaml`:

```yaml
runtimes:
  plantnet-core:      # PyTorch + Scanpy (main models)
    backend: conda
    python: /path/to/scclubench-main/bin/python
    env_file: envs/plantnet-core.yml

  plantnet-tf1:      # TensorFlow 2.x for scNAME, scziDesk, scDeepCluster
    backend: conda
    python: /path/to/scssl_bench_py310/bin/python
    env_file: envs/plantnet-tf1.yml

  plantnet-desc:      # TensorFlow 2.x for DESC (separate from plantnet-tf1)
    backend: conda
    python: /path/to/scssl_bench_py310/bin/python
    env_file: envs/plantnet-desc.yml

  plantnet-sccdcg:   # PyTorch 1.12 for scCDCG
    backend: conda
    python: /path/to/scclubench-sccdcg/bin/python
    env_file: envs/plantnet-sccdcg.yml

  plantnet-scgnn:     # PyTorch for scGNN (noregu mode, no R)
    backend: conda
    python: /path/to/scclubench-main/bin/python
    env_file: envs/plantnet-scgnn.yml

  plantnet-attentionae:  # PyTorch for AttentionAE_sc
    backend: conda
    python: /path/to/scclubench-main/bin/python
```

The **example file** (`runtime_registry.example.yaml`) contains placeholders and should be committed. The **local file** (`runtime_registry.yaml`) contains absolute paths and is gitignored.

### Per-Method Configuration

Each method in `method_manifest.yaml` specifies its runtime:

```yaml
- key: scname
  runtime:
    env_name: plantnet-tf1
    backend: conda
```

### How the Runner Uses Runtime Isolation

`run_formal_benchmark.py` resolves the correct Python executable at runtime:

```python
python_bin = resolve_python_executable(method_info, runtime_registry)
# Returns the specific conda env's Python, or sys.executable as fallback

cmd = [
    python_bin,           # e.g., /path/to/scssl_bench_py310/bin/python
    str(run_py),
    "--data_path", data_path,
    # ... other args
]
```

This means each method's training runs in its own conda environment with **no cross-contamination**.

---

## TensorFlow Models: Lazy Import Pattern

TensorFlow/Keras models (scNAME, scziDesk, DESC, scDeepCluster) use **lazy import** to allow `--help` to work without TF being installed:

```python
_TF_READY = False

def _ensure_tf():
    """Lazy import TensorFlow only when actually needed (not for --help)."""
    global _TF_READY, tf, autoencoder
    if _TF_READY:
        return
    import tensorflow.compat.v1 as tf
    tf.disable_v2_behavior()
    from network import autoencoder
    _TF_READY = True

def main():
    args = parse_args()
    _ensure_tf()  # Only called when actually running
    # ... rest of training
```

This means:
- `python methods/DeepLearning/scNAME/run.py --help` → works in any Python
- `python methods/DeepLearning/scNAME/run.py ...` → runs in `plantnet-tf1` environment

---

## GNN Models: No-R Mode

scGNN has an optional LTMG mode that requires R/rpy2. For the benchmark smoke test, we use the **noregu** mode (default):

```bash
--regulized_type noregu  # No R dependency
```

This ensures the smoke test can run in a pure Python environment.

---

## Environment.json Recording

Every run outputs `environment.json` to its output directory:

```json
{
  "runtime_backend": "conda",
  "runtime_env": "plantnet-tf1",
  "python_executable": "/path/to/scssl_bench_py310/bin/python",
  "python_version": "3.10.12",
  "cuda_visible_devices": "1",
  "torch_version": "N/A",
  "tensorflow_version": "2.13.0",
  "scanpy_version": "1.9.3",
  "anndata_version": "0.9.2"
}
```

This file is automatically generated by `run_formal_benchmark.py` and is included in every benchmark run's output directory.

---

## Status.json Runtime Fields

Every `status.json` also records the runtime:

```json
{
  "status": "success",
  "method": "scname",
  "runtime_env": "plantnet-tf1",
  "python_executable": "/path/to/scssl_bench_py310/bin/python",
  "error": null
}
```

---

## How to Upgrade a Model from ENV-GATED to VERIFIED

1. **Ensure full source is migrated** to `methods/<category>/<model>/`
2. **Create/update core card** at `docs/model_core_cards/<model>.yaml`
3. **Verify smoke test**: Run in the appropriate runtime environment
4. **Update manifest**: Change `authenticity` to `VERIFIED`, `smoke` to `PASS`
5. **Do NOT** set `default_in_formal: true` unless explicitly desired

Example smoke test:
```bash
python scripts/run_formal_benchmark.py \
  --data_path data/SRP182008.h5ad \
  --dataset_name SRP182008 \
  --out_dir results/smoke/<model> \
  --n_clusters auto \
  --seeds 42 \
  --methods <model> \
  --runtime_registry envs/runtime_registry.yaml \
  --pretrain_epochs 1 \
  --epochs 1 \
  --gpu 1
```

---

## Summary Table

| Model | Runtime | Framework | Status |
|-------|---------|-----------|--------|
| scNAME | `plantnet-tf1` | TensorFlow 2.x | ENV-GATED |
| scziDesk | `plantnet-tf1` | TensorFlow 2.x | ENV-GATED |
| scDeepCluster | `plantnet-tf1` | TensorFlow 2.x | ENV-GATED |
| DESC | `plantnet-desc` | TensorFlow 2.x | ENV-GATED |
| scCDCG | `plantnet-sccdcg` | PyTorch 1.12 | PENDING_AUDITED |
| scGNN | `plantnet-scgnn` | PyTorch | PENDING_AUDITED |
| AttentionAE_sc | `plantnet-attentionae` | PyTorch | PENDING_AUDITED |
| All others | `plantnet-core` | PyTorch + Scanpy | VERIFIED |

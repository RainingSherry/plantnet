## Review guidelines

For every PR, Codex should treat the following as blocking issues:

- Unrelated refactoring.
- Changes outside the task's allowed file list.
- Broad try/except or fallback branches that hide logic errors.
- Silent compatibility code for unknown historical versions.
- Ambiguous mask semantics.
- Incorrect loss denominator.
- Tensor shape mismatch.
- Missing forward/backward smoke test for model or loss changes.
- Any preprocessing, training, or evaluation change that may affect experimental fairness.

For ML model changes, review:

- Tensor shapes.
- Mask meaning.
- Loss sign and normalization.
- Gradient path.
- Batch/gene/cell dimension consistency.
- Whether the implementation matches the task file.
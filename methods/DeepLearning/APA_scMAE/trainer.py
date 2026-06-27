from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch.utils.data import DataLoader

from methods.DeepLearning.APA_scMAE.data import APAExpressionDataset, ScMAEShuffleCorruption, assert_no_training_labels
from methods.DeepLearning.APA_scMAE.losses import generator_losses, mask_diagnostics, student_losses
from methods.DeepLearning.APA_scMAE.model import APAModel, freeze, grad_norm


def save_json(path: Path, data: Any) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, default=str)


class APATrainer:
    def __init__(
        self,
        *,
        config: dict[str, Any],
        model: APAModel,
        train_dataset: APAExpressionDataset,
        full_x: torch.Tensor,
        gene_stats: torch.Tensor,
        prototypes: torch.Tensor,
        device: torch.device,
        save_dir: Path,
    ) -> None:
        self.config = config
        self.model = model.to(device)
        self.dataset = train_dataset
        self.full_x = full_x.detach().cpu().float()
        self.gene_stats = gene_stats.to(device)
        self.prototypes = prototypes.to(device)
        self.device = device
        self.save_dir = save_dir
        self.corruption = ScMAEShuffleCorruption(
            self.full_x,
            seed=int(config["seed"]),
            atol=float(config["corruption"]["changed_tolerance_abs"]),
            rtol=float(config["corruption"]["changed_tolerance_rel"]),
        )
        self.student_optimizer = torch.optim.AdamW(
            list(self.model.shared.parameters()) + list(self.model.student.parameters()),
            lr=float(config["training"]["lr_student"]),
            weight_decay=float(config["training"]["weight_decay"]),
        )
        self.generator_optimizer = torch.optim.AdamW(
            self.model.generator.parameters(),
            lr=float(config["training"]["lr_generator"]),
            weight_decay=float(config["training"]["weight_decay"]),
        )
        self.history: dict[str, list[float]] = {
            "loss_student": [],
            "loss_generator": [],
            "loss_rec": [],
            "loss_mask": [],
            "mask_ratio": [],
            "effective_mask_ratio": [],
            "zero_to_zero_rate": [],
            "mask_entropy": [],
            "mask_gini": [],
            "top_gene_concentration": [],
            "student_grad_norm": [],
            "generator_grad_norm": [],
            "generator_update_count": [],
            "batch_seconds": [],
        }
        self.corruption_totals = {
            "selected": 0.0,
            "effective": 0.0,
            "zero_to_zero": 0.0,
            "abs_delta": 0.0,
            "positions": 0.0,
        }
        self.last_mask_stats: dict[str, float] = {}
        self.last_gradient_stats: dict[str, float] = {}

    def _loader(self) -> DataLoader:
        generator = torch.Generator()
        generator.manual_seed(int(self.config["seed"]))
        return DataLoader(
            self.dataset,
            batch_size=int(self.config["training"]["batch_size"]),
            shuffle=True,
            drop_last=False,
            generator=generator,
            num_workers=int(self.config["runtime"]["num_workers"]),
        )

    def _shared(self) -> tuple[torch.Tensor, torch.Tensor]:
        return self.model.shared_context(self.gene_stats)

    def _generator_forward(self, x: torch.Tensor, replacement: torch.Tensor, effective: torch.Tensor, *, detach_shared: bool):
        if detach_shared:
            with torch.no_grad():
                gene_vec, stat_vec = self._shared()
        else:
            gene_vec, stat_vec = self._shared()
        return self.model.generator(
            x,
            replacement,
            effective,
            gene_vec.detach() if detach_shared else gene_vec,
            stat_vec.detach() if detach_shared else stat_vec,
            self.prototypes,
            mask_ratio=float(self.config["mask"]["ratio"]),
            temperature=float(self.config["mask"]["temperature"]),
            topk_only_effective=bool(self.config["mask"].get("generator_topk_only_effective", True)),
        )

    def _student_forward(self, x_tilde: torch.Tensor):
        gene_vec, stat_vec = self._shared()
        return self.model.student(x_tilde, gene_vec, stat_vec, self.prototypes)

    def _record_corruption(self, mask_hard: torch.Tensor, effective_mask: torch.Tensor, delta: torch.Tensor) -> None:
        selected = float(mask_hard.detach().sum().cpu())
        effective = float(effective_mask.detach().sum().cpu())
        self.corruption_totals["selected"] += selected
        self.corruption_totals["effective"] += effective
        self.corruption_totals["zero_to_zero"] += selected - effective
        self.corruption_totals["abs_delta"] += float(delta.detach().sum().cpu())
        self.corruption_totals["positions"] += float(delta.numel())

    def _student_step(self, batch: dict[str, torch.Tensor]) -> dict[str, float]:
        assert_no_training_labels(batch)
        idx = batch["index"].to(self.device).long()
        x = batch["x"].to(self.device)
        donor = self.corruption.sample(idx, device=self.device)
        replacement = donor["replacement"]
        effective = donor["effective"]
        self.model.generator.zero_grad(set_to_none=True)
        freeze(self.model.generator, True)
        self.model.shared.train()
        self.model.student.train()
        self.model.generator.eval()
        with torch.no_grad():
            gen = self._generator_forward(x, replacement, effective, detach_shared=True)
        mask = gen["mask_hard"].detach()
        effective_mask = mask * effective
        x_tilde = x * (1.0 - mask) + replacement.detach() * mask
        self.student_optimizer.zero_grad(set_to_none=True)
        out = self._student_forward(x_tilde)
        losses = student_losses(
            x,
            out["x_recon"],
            out["mask_logits"],
            effective_mask,
            masked_data_weight=float(self.config["mask"]["masked_data_weight"]),
            gamma=float(self.config["training"]["gamma"]),
        )
        losses["loss_student"].backward()
        torch.nn.utils.clip_grad_norm_(
            list(self.model.shared.parameters()) + list(self.model.student.parameters()),
            float(self.config["training"]["student_grad_clip"]),
        )
        s_grad = grad_norm(self.model.student) + grad_norm(self.model.shared)
        g_grad = grad_norm(self.model.generator)
        if bool(self.config["runtime"]["fail_fast"]):
            if s_grad <= 0:
                raise RuntimeError("student/shared gradient is zero during student step")
            if g_grad != 0.0:
                raise RuntimeError("generator received gradient during student step")
        self.student_optimizer.step()
        freeze(self.model.generator, False)
        stats = mask_diagnostics(
            mask,
            effective_mask,
            gen["logits"],
            target_mask_ratio=float(self.config["mask"]["ratio"]),
            budget_deficit=gen.get("budget_deficit"),
        )
        self._record_corruption(mask, effective_mask, gen["delta"])
        self.last_mask_stats = stats
        return {
            "loss_student": float(losses["loss_student"].detach().cpu()),
            "loss_rec": float(losses["loss_rec"].detach().cpu()),
            "loss_mask": float(losses["loss_mask"].detach().cpu()),
            "student_grad_norm": float(s_grad),
            "generator_grad_norm": float(g_grad),
            **stats,
        }

    def _generator_step(self, batch: dict[str, torch.Tensor]) -> dict[str, float]:
        assert_no_training_labels(batch)
        idx = batch["index"].to(self.device).long()
        x = batch["x"].to(self.device)
        donor = self.corruption.sample(idx, device=self.device)
        replacement = donor["replacement"]
        effective = donor["effective"]
        freeze(self.model.shared, True)
        freeze(self.model.student, True)
        freeze(self.model.generator, False)
        self.model.generator.train()
        self.model.student.eval()
        self.generator_optimizer.zero_grad(set_to_none=True)
        self.model.shared.zero_grad(set_to_none=True)
        self.model.student.zero_grad(set_to_none=True)
        gen = self._generator_forward(x, replacement, effective, detach_shared=True)
        mask_st = gen["mask_st"]
        effective_mask_st = mask_st * effective
        x_tilde = x * (1.0 - mask_st) + replacement.detach() * mask_st
        with torch.enable_grad():
            gene_vec, stat_vec = self._shared()
            out = self.model.student(x_tilde, gene_vec.detach(), stat_vec.detach(), self.prototypes)
        losses = generator_losses(
            x,
            out["x_recon"],
            effective_mask_st,
            gen["logits"],
            gen["mask_soft"],
            gen["delta"],
            mask_ratio=float(self.config["mask"]["ratio"]),
            lambda_entropy=float(self.config["generator_loss"]["lambda_entropy"]),
            lambda_balance=float(self.config["generator_loss"]["lambda_balance"]),
            lambda_distortion=float(self.config["generator_loss"]["lambda_distortion"]),
            lambda_coverage=float(self.config["generator_loss"]["lambda_coverage"]),
        )
        losses["loss_generator"].backward()
        torch.nn.utils.clip_grad_norm_(self.model.generator.parameters(), float(self.config["training"]["generator_grad_clip"]))
        g_grad = grad_norm(self.model.generator)
        s_grad = grad_norm(self.model.student) + grad_norm(self.model.shared)
        if bool(self.config["runtime"]["fail_fast"]):
            if g_grad <= 0:
                raise RuntimeError("generator gradient is zero during generator step")
            if s_grad != 0.0:
                raise RuntimeError("student/shared received gradient during generator step")
        self.generator_optimizer.step()
        freeze(self.model.shared, False)
        freeze(self.model.student, False)
        self.last_gradient_stats = {
            "generator_grad_norm": float(g_grad),
            "student_grad_norm_during_generator_step": float(s_grad),
        }
        return {"loss_generator": float(losses["loss_generator"].detach().cpu()), "generator_grad_norm": float(g_grad)}

    def train(self) -> dict[str, list[float]]:
        epochs = max(1, int(self.config["training"]["epochs"]))
        loader = self._loader()
        for _epoch in range(1, epochs + 1):
            totals: dict[str, float] = {}
            n_batches = 0
            n_generator_updates = 0
            for batch_id, batch in enumerate(loader):
                start = time.perf_counter()
                step = self._student_step(batch)
                if (batch_id + 1) % int(self.config["training"]["generator_update_interval"]) == 0:
                    gen_step = self._generator_step(batch)
                    step.update(gen_step)
                    n_generator_updates += 1
                step["batch_seconds"] = time.perf_counter() - start
                for key, value in step.items():
                    totals[key] = totals.get(key, 0.0) + float(value)
                n_batches += 1
            for key in self.history:
                if key == "loss_generator":
                    self.history[key].append(totals.get(key, float("nan")) / n_generator_updates if n_generator_updates else float("nan"))
                elif key == "generator_grad_norm":
                    self.history[key].append(totals.get(key, float("nan")) / n_generator_updates if n_generator_updates else float("nan"))
                elif key == "generator_update_count":
                    self.history[key].append(float(n_generator_updates))
                else:
                    self.history[key].append(totals.get(key, 0.0) / max(1, n_batches))
        return self.history

    def corruption_stats(self) -> dict[str, float | str]:
        selected = max(1.0, self.corruption_totals["selected"])
        positions = max(1.0, self.corruption_totals["positions"])
        return {
            "corruption_type": str(self.config["corruption"]["type"]),
            "selected_positions": float(self.corruption_totals["selected"]),
            "effective_corruption_rate": float(self.corruption_totals["effective"] / selected),
            "zero_to_zero_rate": float(self.corruption_totals["zero_to_zero"] / selected),
            "mean_abs_delta": float(self.corruption_totals["abs_delta"] / positions),
            "budget_deficit_rate": float(self.last_mask_stats.get("budget_deficit_rate", 0.0)),
        }

    @torch.no_grad()
    def extract_embeddings(self, batch_size: int = 512) -> np.ndarray:
        self.model.eval()
        loader = DataLoader(self.dataset, batch_size=batch_size, shuffle=False, drop_last=False)
        gene_vec, stat_vec = self._shared()
        chunks: list[np.ndarray] = []
        for batch in loader:
            x = batch["x"].to(self.device)
            z = self.model.student.feature(x, gene_vec, stat_vec, self.prototypes)
            chunks.append(z.detach().cpu().numpy().astype(np.float32))
        return np.concatenate(chunks, axis=0)

    def save_diagnostics(self) -> None:
        save_json(self.save_dir / "training_history.json", self.history)
        save_json(self.save_dir / "mask_stats.json", self.last_mask_stats)
        save_json(self.save_dir / "gradient_stats.json", self.last_gradient_stats)
        save_json(self.save_dir / "corruption_stats.json", self.corruption_stats())

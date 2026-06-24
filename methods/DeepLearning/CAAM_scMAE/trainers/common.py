from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch.utils.data import DataLoader

from methods.DeepLearning.CAAM_scMAE.corruption.matched_donor import MatchedDonorCorruption
from methods.DeepLearning.CAAM_scMAE.data.dataset import CAAMExpressionDataset, assert_no_training_labels
from methods.DeepLearning.CAAM_scMAE.data.donor_candidates import DonorCandidateProvider
from methods.DeepLearning.CAAM_scMAE.diagnostics.attention_stats import summarize_attention
from methods.DeepLearning.CAAM_scMAE.diagnostics.gradient_stats import collect_gradient_stats
from methods.DeepLearning.CAAM_scMAE.diagnostics.mask_stats import summarize_mask
from methods.DeepLearning.CAAM_scMAE.losses.generator_objective import generator_loss
from methods.DeepLearning.CAAM_scMAE.losses.loss_bundle import student_loss_bundle
from methods.DeepLearning.CAAM_scMAE.mask_generator.adversarial_mask import AdversarialMaskGenerator
from methods.DeepLearning.CAAM_scMAE.mask_generator.random_mask import RandomFixedBudgetMask
from methods.DeepLearning.CAAM_scMAE.mask_generator.regularizers import generator_regularization
from methods.DeepLearning.CAAM_scMAE.models.common import freeze_module, grad_norm


class CAAMTrainer:
    def __init__(
        self,
        *,
        config: dict[str, Any],
        student,
        train_dataset: CAAMExpressionDataset,
        donor_provider: DonorCandidateProvider,
        full_x: torch.Tensor,
        device: torch.device,
        save_dir: Path,
        context_indices: np.ndarray | None = None,
    ) -> None:
        self.config = config
        self.student = student.to(device)
        self.dataset = train_dataset
        self.donor_provider = donor_provider
        self.full_x = full_x.to(device)
        self.device = device
        self.save_dir = save_dir
        self.context_indices = torch.as_tensor(context_indices, dtype=torch.long, device=device) if context_indices is not None else None
        self.random_mask = RandomFixedBudgetMask(float(config["mask"]["ratio"]))
        self.corruption = MatchedDonorCorruption()
        self.generator = None
        if config["model"]["mask_selector"] == "adversarial":
            self.generator = AdversarialMaskGenerator(
                n_genes=int(full_x.shape[1]),
                hidden_dim=int(config["generator"]["hidden_dim"]),
                mask_ratio=float(config["mask"]["ratio"]),
            ).to(device)
        self.student_optimizer = torch.optim.AdamW(
            self.student.parameters(),
            lr=float(config["training"]["lr_student"]),
            weight_decay=float(config["training"]["weight_decay"]),
        )
        self.generator_optimizer = None
        if self.generator is not None:
            self.generator_optimizer = torch.optim.AdamW(
                self.generator.parameters(),
                lr=float(config["training"]["lr_generator"]),
                weight_decay=float(config["training"]["weight_decay"]),
            )
        if self.generator_optimizer is not None:
            ids_a = {id(p) for group in self.student_optimizer.param_groups for p in group["params"]}
            ids_b = {id(p) for group in self.generator_optimizer.param_groups for p in group["params"]}
            if ids_a.intersection(ids_b):
                raise RuntimeError("Student and generator optimizer parameter sets overlap.")
        self.history: dict[str, list] = {
            "loss_student": [],
            "loss_generator": [],
            "loss_rec_masked": [],
            "loss_mask": [],
            "mask_ratio": [],
            "student_grad_norm": [],
            "generator_grad_norm": [],
            "context_cache_checksum": [],
        }
        self.last_batch_debug: dict[str, Any] = {}
        self.last_mask_stats: dict[str, Any] = {}
        self.last_attention_stats: dict[str, Any] = {}
        self.last_gradient_stats: dict[str, Any] = {}
        self.last_generator_stats: dict[str, Any] = {}

    def _loader(self) -> DataLoader:
        generator = torch.Generator()
        generator.manual_seed(int(self.config["seed"]))
        return DataLoader(
            self.dataset,
            batch_size=int(self.config["training"]["batch_size"]),
            shuffle=True,
            drop_last=False,
            generator=generator,
            num_workers=0,
        )

    def _temperature(self, epoch: int) -> float:
        warmup = int(self.config["training"]["student_warmup_epochs"])
        total = max(1, int(self.config["training"]["epochs"]) - warmup)
        t = min(1.0, max(0.0, (epoch - warmup) / total))
        start = float(self.config["generator"]["temperature_start"])
        end = float(self.config["generator"]["temperature_end"])
        return start + t * (end - start)

    def _refresh_context(self) -> None:
        if self.context_indices is not None:
            context_x = self.full_x[self.context_indices]
            self.student.refresh_context_cache(context_x, self.context_indices)

    def _student_mask(self, x: torch.Tensor, eligibility: torch.Tensor, epoch: int):
        if self.generator is None or epoch <= int(self.config["training"]["student_warmup_epochs"]):
            return self.random_mask(x, eligibility)
        with torch.no_grad():
            logits, hard, _soft, _st, info = self.generator(x, eligibility, self._temperature(epoch), add_gumbel=True)
        return logits.detach(), hard.detach(), info

    def _student_step(self, batch: dict, epoch: int, batch_id: int) -> dict[str, float]:
        assert_no_training_labels(batch)
        if self.generator is not None:
            self.generator.zero_grad(set_to_none=True)
        idx = batch["index"].to(self.device).long()
        x = batch["x"].to(self.device)
        donor = self.donor_provider.sample_batch(idx, self.full_x, self.device)
        logits, mask, mask_info = self._student_mask(x, donor["eligibility"], epoch)
        corrupt = self.corruption.corrupt(x, mask, donor["replacement"], donor["eligibility"], donor["donor_indices"])
        self.student.train()
        self.student_optimizer.zero_grad(set_to_none=True)
        out = self.student(corrupt.x_tilde, mask=mask, indices=idx)
        losses = student_loss_bundle(
            x,
            out["x_hat"],
            out["mask_logits"],
            mask,
            lambda_visible=float(self.config["loss"]["lambda_visible"]),
            lambda_mask=float(self.config["loss"]["lambda_mask"]),
        )
        if not torch.isfinite(losses["loss_student"]):
            raise FloatingPointError("loss_student is NaN/Inf")
        losses["loss_student"].backward()
        torch.nn.utils.clip_grad_norm_(self.student.parameters(), float(self.config["training"]["student_grad_clip"]))
        s_grad = grad_norm(self.student)
        if s_grad <= 0 and bool(self.config["runtime"]["fail_fast"]):
            raise RuntimeError("student gradient is zero")
        self.student_optimizer.step()
        g_grad = grad_norm(self.generator) if self.generator is not None else 0.0
        if self.generator is not None and g_grad != 0.0 and bool(self.config["runtime"]["fail_fast"]):
            raise RuntimeError("generator received gradient during student step")
        self.last_mask_stats = summarize_mask(mask, donor["eligibility"])
        self.last_attention_stats = summarize_attention(out.get("gene_attn"), out.get("cell_attn"))
        self.last_gradient_stats = collect_gradient_stats(self.student, self.generator)
        if batch_id == 0 and epoch == 1:
            self.last_batch_debug = {
                "first_batch_indices": idx.detach().cpu().numpy().astype(int).tolist(),
                "first_mask_hard_sum": float(mask.detach().sum().cpu()),
                "first_donor_indices_checksum": int(donor["donor_indices"].detach().sum().cpu()),
                "first_loss": float(losses["loss_student"].detach().cpu()),
            }
        return {
            "loss_student": float(losses["loss_student"].detach().cpu()),
            "loss_rec_masked": float(losses["loss_rec_masked"].detach().cpu()),
            "loss_mask": float(losses["loss_mask"].detach().cpu()),
            "mask_ratio": float(mask.detach().mean().cpu()),
            "student_grad_norm": float(s_grad),
            "generator_grad_norm": float(g_grad),
            **mask_info,
        }

    def _generator_step(self, batch: dict, epoch: int) -> dict[str, float]:
        if self.generator is None or self.generator_optimizer is None:
            return {}
        assert_no_training_labels(batch)
        idx = batch["index"].to(self.device).long()
        x = batch["x"].to(self.device)
        donor = self.donor_provider.sample_batch(idx, self.full_x, self.device)
        tau = self._temperature(epoch)
        self.generator.train()
        self.generator_optimizer.zero_grad(set_to_none=True)
        self.student.zero_grad(set_to_none=True)
        with freeze_module(self.student):
            self.student.eval()
            logits, hard, soft, st, _info = self.generator(x, donor["eligibility"], tau, add_gumbel=True)
            x_tilde = x * (1.0 - st) + donor["replacement"].detach() * st
            out = self.student(x_tilde, mask=hard, indices=idx)
            losses = student_loss_bundle(
                x,
                out["x_hat"],
                out["mask_logits"],
                hard,
                lambda_visible=float(self.config["loss"]["lambda_visible"]),
                lambda_mask=float(self.config["loss"]["lambda_mask"]),
            )
            regs = generator_regularization(
                soft,
                logits,
                donor["eligibility"],
                x,
                donor["replacement"],
                tau=tau,
                distortion_min=float(self.config["generator"]["distortion_min"]),
                distortion_max=float(self.config["generator"]["distortion_max"]),
            )
            loss_g = generator_loss(
                losses["loss_rec_masked"],
                losses["loss_mask"],
                regs,
                beta_mask_loss=float(self.config["generator"]["beta_mask_loss"]),
                lambda_coverage=float(self.config["loss"]["lambda_coverage"]),
                lambda_distortion=float(self.config["loss"]["lambda_distortion"]),
                lambda_entropy=float(self.config["loss"]["lambda_entropy"]),
            )
            if not torch.isfinite(loss_g):
                raise FloatingPointError("loss_generator is NaN/Inf")
            loss_g.backward()
        torch.nn.utils.clip_grad_norm_(self.generator.parameters(), float(self.config["training"]["generator_grad_clip"]))
        g_grad = grad_norm(self.generator)
        s_grad = grad_norm(self.student)
        if bool(self.config["runtime"]["fail_fast"]):
            if g_grad <= 0:
                raise RuntimeError("generator gradient is zero during generator step")
            if s_grad != 0.0:
                raise RuntimeError("student received gradient during generator step")
        self.generator_optimizer.step()
        self.last_generator_stats = {
            "loss_generator": float(loss_g.detach().cpu()),
            "generator_grad_norm": float(g_grad),
            "student_grad_norm_during_generator_step": float(s_grad),
            "temperature": float(tau),
            "coverage_loss": float(regs["coverage_loss"].detach().cpu()),
            "distortion_loss": float(regs["distortion_loss"].detach().cpu()),
            "entropy_loss": float(regs["entropy_loss"].detach().cpu()),
            "mask_entropy": float(regs["mask_entropy"].cpu()),
            "mask_gini": float(regs["mask_gini"].cpu()),
        }
        return self.last_generator_stats

    def train(self) -> dict[str, list]:
        epochs = max(1, int(self.config["training"]["epochs"]))
        loader = self._loader()
        for epoch in range(1, epochs + 1):
            self._refresh_context()
            totals: dict[str, float] = {}
            n_batches = 0
            for batch_id, batch in enumerate(loader):
                step = self._student_step(batch, epoch, batch_id)
                for key in ("loss_student", "loss_rec_masked", "loss_mask", "mask_ratio", "student_grad_norm", "generator_grad_norm"):
                    totals[key] = totals.get(key, 0.0) + float(step.get(key, 0.0))
                n_batches += 1
                if (
                    self.generator is not None
                    and epoch > int(self.config["training"]["student_warmup_epochs"])
                    and (batch_id + 1) % int(self.config["training"]["generator_update_interval"]) == 0
                ):
                    g_step = self._generator_step(batch, epoch)
                    if g_step:
                        totals["loss_generator"] = totals.get("loss_generator", 0.0) + float(g_step["loss_generator"])
            for key in ("loss_student", "loss_rec_masked", "loss_mask", "mask_ratio", "student_grad_norm", "generator_grad_norm"):
                self.history[key].append(totals.get(key, 0.0) / max(1, n_batches))
            self.history["loss_generator"].append(totals.get("loss_generator", 0.0) / max(1, n_batches))
            self.history["context_cache_checksum"].append(float(self.student.context_cache_checksum()))
        return self.history

    def save_diagnostics(self) -> None:
        def dump(name: str, obj: Any) -> None:
            with open(self.save_dir / name, "w", encoding="utf-8") as handle:
                json.dump(obj, handle, indent=2)

        dump("training_history.json", self.history)
        dump("mask_stats.json", self.last_mask_stats)
        dump("gradient_stats.json", self.last_gradient_stats)
        dump("attention_stats.json", self.last_attention_stats)
        if self.generator is not None:
            dump("generator_stats.json", self.last_generator_stats)
        dump("reproducibility_debug.json", self.last_batch_debug)

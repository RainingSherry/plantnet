from __future__ import annotations

import copy
import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch.utils.data import DataLoader

from methods.DeepLearning.APA_scMAE.data import APAExpressionDataset, ScMAEShuffleCorruption, assert_no_training_labels
from methods.DeepLearning.APA_scMAE.losses import generator_losses, mask_diagnostics, student_losses
from methods.DeepLearning.APA_scMAE.model import APAModel, freeze, grad_norm, straight_through_topk
from methods.DeepLearning.APA_scMAE.representation_losses import (
    balanced_assignment_loss,
    prototype_kl_loss,
    soft_assignment,
    teacher_consistency_loss,
    vicreg_losses,
)


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
        self.global_student_steps = 0
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
        self.teacher = self._build_teacher()
        self.embedding_prototypes: torch.Tensor | None = None
        self._pending_proto_update: tuple[torch.Tensor, torch.Tensor] | None = None
        self.history: dict[str, list[float]] = {
            "loss_student": [],
            "loss_generator": [],
            "loss_rec": [],
            "loss_mask": [],
            "loss_repr": [],
            "loss_repr_invariance": [],
            "loss_repr_variance": [],
            "loss_repr_covariance": [],
            "loss_teacher": [],
            "loss_proto": [],
            "loss_proto_balance": [],
            "mask_ratio": [],
            "effective_mask_ratio": [],
            "zero_to_zero_rate": [],
            "mask_entropy": [],
            "mask_gini": [],
            "top_gene_concentration": [],
            "student_grad_norm": [],
            "generator_grad_norm": [],
            "generator_update_count": [],
            "generator_update_forced": [],
            "warmup_epoch": [],
            "generator_enabled": [],
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

    def _build_teacher(self) -> APAModel | None:
        if not bool(self.config.get("teacher", {}).get("use_ema_teacher", True)):
            return None
        teacher = copy.deepcopy(self.model).to(self.device)
        teacher.eval()
        for param in teacher.parameters():
            param.requires_grad_(False)
        return teacher

    @torch.no_grad()
    def _init_embedding_prototypes(self) -> torch.Tensor | None:
        proto_cfg = self.config.get("prototype_consistency", {})
        if not bool(proto_cfg.get("use_proto_consistency", True)):
            return None
        n_cells = len(self.dataset)
        requested = int(proto_cfg.get("num_embedding_prototypes", 0))
        if requested <= 0:
            requested = int(self.config.get("prototype", {}).get("n_prototypes", 16))
        n_proto = max(1, min(int(requested), int(n_cells)))
        self.model.eval()
        loader = DataLoader(
            self.dataset,
            batch_size=max(1, int(self.config["training"]["batch_size"])),
            shuffle=False,
            drop_last=False,
        )
        gene_vec, stat_vec = self._shared()
        chunks: list[np.ndarray] = []
        for batch in loader:
            z = self.model.student.feature(batch["x"].to(self.device), gene_vec, stat_vec, self.prototypes)
            chunks.append(z.detach().cpu().numpy().astype(np.float32))
        embedding = np.concatenate(chunks, axis=0)
        if n_proto == 1:
            centers = embedding.mean(axis=0, keepdims=True)
        else:
            from sklearn.cluster import KMeans

            centers = KMeans(n_clusters=n_proto, random_state=int(self.config["seed"]), n_init=10).fit(embedding).cluster_centers_
        return torch.as_tensor(centers, dtype=torch.float32, device=self.device)

    def _ensure_embedding_prototypes_initialized(self) -> None:
        if self.embedding_prototypes is not None:
            return
        if not bool(self.config.get("prototype_consistency", {}).get("use_proto_consistency", True)):
            return
        self.embedding_prototypes = self._init_embedding_prototypes()

    def _shared(self) -> tuple[torch.Tensor, torch.Tensor]:
        return self.model.shared_context(self.gene_stats)

    def _representation_mechanisms_enabled(self) -> bool:
        return (
            bool(self.config.get("representation_loss", {}).get("use_repr_loss", True))
            or bool(self.config.get("teacher", {}).get("use_ema_teacher", True))
            or bool(self.config.get("prototype_consistency", {}).get("use_proto_consistency", True))
        )

    def _teacher_momentum(self) -> float:
        teacher_cfg = self.config.get("teacher", {})
        start = float(teacher_cfg.get("teacher_momentum_start", 0.99))
        end = float(teacher_cfg.get("teacher_momentum_end", 0.999))
        total = max(1, int(self.config["training"]["epochs"]) * max(1, len(self._loader())))
        progress = min(1.0, float(self.global_student_steps) / float(total))
        return start + (end - start) * progress

    @torch.no_grad()
    def _update_teacher(self) -> None:
        if self.teacher is None:
            return
        momentum = self._teacher_momentum()
        for teacher_param, student_param in zip(self.teacher.shared.parameters(), self.model.shared.parameters()):
            teacher_param.data.mul_(momentum).add_(student_param.data, alpha=1.0 - momentum)
        for teacher_param, student_param in zip(self.teacher.student.parameters(), self.model.student.parameters()):
            teacher_param.data.mul_(momentum).add_(student_param.data, alpha=1.0 - momentum)

    def _random_mask(self, x: torch.Tensor, replacement: torch.Tensor, effective: torch.Tensor) -> dict[str, torch.Tensor]:
        logits = torch.rand_like(x)
        if bool(self.config["mask"].get("generator_topk_only_effective", True)):
            eligibility = effective.bool()
        else:
            eligibility = torch.ones_like(effective, dtype=torch.bool)
        hard, soft, st, info = straight_through_topk(
            logits,
            float(self.config["mask"]["ratio"]),
            float(self.config["mask"]["temperature"]),
            eligibility,
        )
        return {
            "logits": logits,
            "mask_hard": hard,
            "mask_soft": soft,
            "mask_st": st,
            "delta": (replacement - x).abs(),
            **info,
        }

    def _representation_objectives(
        self,
        *,
        x: torch.Tensor,
        z_masked: torch.Tensor,
        gene_vec: torch.Tensor,
        stat_vec: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        device_zero = z_masked.sum() * 0.0
        losses: dict[str, torch.Tensor] = {
            "loss_repr": device_zero,
            "loss_repr_invariance": device_zero,
            "loss_repr_variance": device_zero,
            "loss_repr_covariance": device_zero,
            "loss_teacher": device_zero,
            "loss_proto": device_zero,
            "loss_proto_balance": device_zero,
        }
        repr_cfg = self.config.get("representation_loss", {})
        teacher_cfg = self.config.get("teacher", {})
        proto_cfg = self.config.get("prototype_consistency", {})

        use_repr = bool(repr_cfg.get("use_repr_loss", True))
        use_teacher = self.teacher is not None and bool(teacher_cfg.get("use_ema_teacher", True))
        use_proto = self.embedding_prototypes is not None and bool(proto_cfg.get("use_proto_consistency", True))
        if not (use_repr or use_teacher or use_proto):
            return losses

        z_clean: torch.Tensor | None = None
        if use_repr or (use_proto and not use_teacher):
            z_clean = self.model.student.encode_clean(x, gene_vec, stat_vec, self.prototypes)
        if use_repr:
            assert z_clean is not None
            vic = vicreg_losses(
                z_clean,
                z_masked,
                variance_margin=float(repr_cfg.get("variance_margin", 1.0)),
                eps=float(repr_cfg.get("eps", 1.0e-4)),
            )
            losses.update(vic)
            losses["loss_repr"] = (
                float(repr_cfg.get("lambda_invariance", 1.0)) * vic["loss_repr_invariance"]
                + float(repr_cfg.get("lambda_variance", 0.1)) * vic["loss_repr_variance"]
                + float(repr_cfg.get("lambda_covariance", 0.01)) * vic["loss_repr_covariance"]
            )

        if use_teacher:
            with torch.no_grad():
                teacher_gene_vec, teacher_stat_vec = self.teacher.shared_context(self.gene_stats)
                z_clean_teacher = self.teacher.student.encode_clean(x, teacher_gene_vec, teacher_stat_vec, self.prototypes)
            losses["loss_teacher"] = float(teacher_cfg.get("teacher_consistency_weight", 1.0)) * teacher_consistency_loss(
                z_masked,
                z_clean_teacher,
            )
        else:
            assert z_clean is not None
            z_clean_teacher = z_clean.detach()

        if use_proto:
            q_clean = soft_assignment(
                z_clean_teacher.detach(),
                self.embedding_prototypes.detach(),
                temperature=float(proto_cfg.get("proto_temperature", 0.2)),
            )
            q_masked = soft_assignment(
                z_masked,
                self.embedding_prototypes.detach(),
                temperature=float(proto_cfg.get("proto_temperature", 0.2)),
            )
            proto_kl = prototype_kl_loss(q_clean, q_masked)
            proto_balance = balanced_assignment_loss(q_masked)
            losses["loss_proto"] = float(proto_cfg.get("proto_loss_weight", 0.1)) * proto_kl
            losses["loss_proto_balance"] = float(proto_cfg.get("proto_balance_weight", 0.01)) * proto_balance
            self._pending_proto_update = (z_clean_teacher.detach().clone(), q_clean.detach().clone())
        return losses

    @torch.no_grad()
    def _apply_pending_embedding_prototype_update(self) -> None:
        if self._pending_proto_update is None:
            return
        z_clean_teacher, q_clean = self._pending_proto_update
        self._pending_proto_update = None
        self._update_embedding_prototypes(z_clean_teacher, q_clean)

    @torch.no_grad()
    def _update_embedding_prototypes(self, z_clean_teacher: torch.Tensor, q_clean: torch.Tensor) -> None:
        if self.embedding_prototypes is None:
            return
        weights = q_clean.sum(dim=0)
        centroids = q_clean.T.matmul(z_clean_teacher) / weights.clamp_min(1.0e-8).unsqueeze(1)
        active = weights > 1.0e-6
        if not bool(active.any()):
            return
        momentum = float(self.config.get("prototype_consistency", {}).get("proto_ema_momentum", 0.95))
        self.embedding_prototypes[active] = (
            momentum * self.embedding_prototypes[active]
            + (1.0 - momentum) * centroids[active]
        )

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

    def _student_step(self, batch: dict[str, torch.Tensor], *, warmup: bool, generator_enabled: bool) -> dict[str, float]:
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
            if warmup or not generator_enabled:
                gen = self._random_mask(x, replacement, effective)
            else:
                gen = self._generator_forward(x, replacement, effective, detach_shared=True)
        mask = gen["mask_hard"].detach()
        effective_mask = mask * effective
        x_tilde = x * (1.0 - mask) + replacement.detach() * mask
        self.student_optimizer.zero_grad(set_to_none=True)
        gene_vec, stat_vec = self._shared()
        out = self.model.student(x_tilde, gene_vec, stat_vec, self.prototypes)
        losses = student_losses(
            x,
            out["x_recon"],
            out["mask_logits"],
            effective_mask,
            masked_data_weight=float(self.config["mask"]["masked_data_weight"]),
            gamma=float(self.config["training"]["gamma"]),
        )
        repr_losses = self._representation_objectives(x=x, z_masked=out["z"], gene_vec=gene_vec, stat_vec=stat_vec)
        losses.update(repr_losses)
        losses["loss_student"] = (
            losses["loss_student"]
            + losses["loss_repr"]
            + losses["loss_teacher"]
            + losses["loss_proto"]
            + losses["loss_proto_balance"]
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
        self.global_student_steps += 1
        self._update_teacher()
        self._apply_pending_embedding_prototype_update()
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
            "loss_repr": float(losses["loss_repr"].detach().cpu()),
            "loss_repr_invariance": float(losses["loss_repr_invariance"].detach().cpu()),
            "loss_repr_variance": float(losses["loss_repr_variance"].detach().cpu()),
            "loss_repr_covariance": float(losses["loss_repr_covariance"].detach().cpu()),
            "loss_teacher": float(losses["loss_teacher"].detach().cpu()),
            "loss_proto": float(losses["loss_proto"].detach().cpu()),
            "loss_proto_balance": float(losses["loss_proto_balance"].detach().cpu()),
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
            if self._representation_mechanisms_enabled():
                z_clean = self.model.student.encode_clean(x, gene_vec.detach(), stat_vec.detach(), self.prototypes)
                representation_delta = (out["z"] - z_clean.detach()).pow(2).mean(dim=1)
            else:
                representation_delta = None
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
            lambda_reconstruction_moderate=float(self.config["generator_loss"].get("lambda_reconstruction_moderate", 0.1)),
            lambda_representation_moderate=float(self.config["generator_loss"].get("lambda_representation_moderate", 0.1)),
            representation_delta=representation_delta,
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
        warmup_epochs = max(0, int(self.config["training"].get("student_warmup_epochs", 0)))
        for _epoch in range(1, epochs + 1):
            warmup = _epoch <= warmup_epochs
            generator_enabled = (not warmup) and bool(self.config["training"].get("enable_generator_after_warmup", True))
            if not warmup:
                self._ensure_embedding_prototypes_initialized()
            totals: dict[str, float] = {}
            n_batches = 0
            n_generator_updates = 0
            generator_update_forced = 0
            last_batch: dict[str, torch.Tensor] | None = None
            for batch_id, batch in enumerate(loader):
                start = time.perf_counter()
                last_batch = batch
                step = self._student_step(batch, warmup=warmup, generator_enabled=generator_enabled)
                if generator_enabled and (batch_id + 1) % int(self.config["training"]["generator_update_interval"]) == 0:
                    gen_step = self._generator_step(batch)
                    step.update(gen_step)
                    n_generator_updates += 1
                step["batch_seconds"] = time.perf_counter() - start
                for key, value in step.items():
                    totals[key] = totals.get(key, 0.0) + float(value)
                n_batches += 1
            if generator_enabled and n_batches > 0 and n_generator_updates == 0 and last_batch is not None:
                gen_step = self._generator_step(last_batch)
                for key, value in gen_step.items():
                    totals[key] = totals.get(key, 0.0) + float(value)
                n_generator_updates += 1
                generator_update_forced = 1
            for key in self.history:
                if key == "loss_generator":
                    self.history[key].append(totals.get(key, 0.0) / n_generator_updates if n_generator_updates else None)
                elif key == "generator_grad_norm":
                    self.history[key].append(totals.get(key, 0.0) / n_generator_updates if n_generator_updates else None)
                elif key == "generator_update_count":
                    self.history[key].append(float(n_generator_updates))
                elif key == "generator_update_forced":
                    self.history[key].append(float(generator_update_forced))
                elif key == "warmup_epoch":
                    self.history[key].append(1.0 if warmup else 0.0)
                elif key == "generator_enabled":
                    self.history[key].append(1.0 if generator_enabled else 0.0)
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

import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.cluster import KMeans

from .support_mask import SupportMaskNet
from .latent_diffusion import LatentDiffusionAE


class ClusterLoss(nn.Module):
    """DEC-style clustering loss (Xie et al. 2016, Deep Embedded Clustering).

    Uses soft assignment Q and target distribution P to iteratively refine
    cluster centers. MUST be initialized with KMeans before use.

    IMPORTANT: This loss should be kept weak (weight ~0.02-0.05) and should
    NOT be active during early training. See ScSpade.train_epoch() for the
    staged schedule that governs when this loss is applied.
    """

    def __init__(self, n_clusters: int, latent_dim: int, cluster_dim: int = 16):
        super().__init__()
        self.n_clusters = n_clusters
        # Learnable cluster centers, initialized to zeros (must call init_centers() before training)
        self.cluster_centers = nn.Parameter(torch.zeros(n_clusters, cluster_dim))

    def init_centers(self, embeddings: torch.Tensor):
        """Initialize cluster centers with KMeans on pretrained embeddings.

        This MUST be called after the autoencoder is pretrained but BEFORE
        cluster loss is enabled. Random initialization would pull the latent
        space to meaningless positions.

        Args:
            embeddings: (n_cells, latent_dim) tensor from pretrained encoder.
        """
        with torch.no_grad():
            z = embeddings.detach().cpu().numpy()
        kmeans = KMeans(n_clusters=self.n_clusters, random_state=42, n_init=20)
        kmeans.fit(z)
        self.cluster_centers.data.copy_(torch.from_numpy(kmeans.cluster_centers_).float())

    @property
    def centers(self) -> torch.Tensor:
        """Return cluster centers, expanded to match embedding batch dimension."""
        return self.cluster_centers

    def get_soft_assignments(self, z: torch.Tensor) -> torch.Tensor:
        """Compute Student's t-soft assignment Q (one-hop kernel).

        Q_{ic} = (1 + ||z_i - μ_c||²)^{-1} / Σ_c' (1 + ||z_i - μ_{c'}||²)^{-1}
        """
        # Expand z to (batch, 1, latent) and centers to (1, n_clusters, latent)
        z_exp = z[:, None, :]          # (batch, 1, latent)
        c_exp = self.centers[None, :, :]  # (1, n_clusters, latent)
        # Squared distances: (batch, n_clusters)
        dist2 = ((z_exp - c_exp) ** 2).sum(dim=-1)
        q = (1.0 + dist2) ** (-1.0)
        q = q / (q.sum(dim=-1, keepdim=True) + 1e-10)
        return q

    def get_target_distribution(self, q: torch.Tensor) -> torch.Tensor:
        """Compute target distribution P from Q (DDPM-style sharpening).
        P_ic = Q_ic² / Σ Q_ic
        """
        p = q ** 2
        p = p / (p.sum(dim=-1, keepdim=True) + 1e-10)
        return p

    def forward(self, z: torch.Tensor, weight: float = 1.0) -> torch.Tensor:
        """DEC clustering loss: KL(Q || P).

        Args:
            z: Latent embeddings (batch, latent_dim).
            weight: Scalar multiplier on the loss.
        Returns:
            Scalar clustering loss.
        """
        q = self.get_soft_assignments(z)
        p = self.get_target_distribution(q).detach()  # target, no backprop through P
        loss = F.kl_div(q.log(), p, reduction="batchmean")
        return loss * weight


class ScSpade(nn.Module):
    """Support-Masked Latent Diffusion Autoencoder for single-cell clustering.

    This model integrates three components following the DOLORIS philosophy:

        1. SupportMaskNet  — predicts per-gene activation probability P(X>0).
                              The prediction is used as soft per-gene importance
                              weights for the reconstruction loss (not ground-truth mask).

        2. LatentDiffusionAE — learns a cell latent embedding z via encoder,
                              denoises z via DDPM, and reconstructs gene expression.

        3. ClusterLoss      — DEC-style soft clustering loss on the latent space.
                              Must be KMeans-initialized before activation.

    DATA FLOW:
        X
          -> SupportMaskNet -> mask_prob (P_gene_activated)
          -> soft_mask = observed_mask * (0.5 + 0.5 * mask_prob.detach())
          -> LatentDiffusionAE(
                 x=X,
                 mask=soft_mask,          <- KEY: mask_prob enters the main flow
                 return_recon=True,
                 sample_diffusion=False,
             )
          -> z (encoder output)
          -> z_denoised (single-step denoised, used in training)
          -> z_diffusion (full p_sample_loop, optional at inference)
          -> x_recon (decoder output, loss computed on active genes)

    TRAINING SCHEDULE (see train_joint.py:train_epoch):
        Epoch < warmup_epochs:
            effective_cluster_weight = 0.0
            effective_diffusion_weight = 0.0
            Only mask loss + reconstruction loss active.
        Epoch < warmup_epochs + 50:
            effective_diffusion_weight = min(diffusion_weight, 0.05)
            effective_cluster_weight = 0.0
        Epoch >= warmup_epochs + 50:
            Full loss schedule.

    This staged approach ensures the encoder first learns stable cell structure
    before diffusion and clustering losses compete for the latent space.
    """

    def __init__(
        self,
        n_genes: int,
        n_clusters: int,
        latent_dim: int = 32,
        hidden_dim: int = 256,
        diffusion_hidden_dim: int = 256,
        diffusion_steps: int = 100,
        dropout: float = 0.1,
        mask_hidden_dims: list = None,
        mask_dropout: float = 0.1,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_clusters = n_clusters
        self.latent_dim = latent_dim

        # 1. Support mask network (predicts gene activation probability)
        self.mask_net = SupportMaskNet(
            n_genes=n_genes,
            hidden_dims=mask_hidden_dims or [512, 256, 128],
            dropout=mask_dropout,
        )

        # 2. Latent diffusion autoencoder (encoder + DDPM denoiser + decoder)
        self.diffusion_ae = LatentDiffusionAE(
            n_genes=n_genes,
            latent_dim=latent_dim,
            hidden_dim=hidden_dim,
            diffusion_hidden_dim=diffusion_hidden_dim,
            diffusion_steps=diffusion_steps,
            dropout=dropout,
        )

        # 3. DEC-style clustering loss (disabled until KMeans init)
        self.cluster_loss_fn = ClusterLoss(
            n_clusters=n_clusters,
            latent_dim=latent_dim,
            cluster_dim=latent_dim,
        )

        # Track initialization state
        self._cluster_centers_initialized = False

    # ── Main forward pass ─────────────────────────────────────────────────────

    def forward(
        self,
        x: torch.Tensor,
        support: torch.Tensor = None,
        return_recon: bool = True,
        sample_diffusion: bool = False,
        mask_override: torch.Tensor = None,
    ) -> dict:
        """Full forward pass with mask-aware diffusion autoencoding.

        KEY DESIGN: mask_prob from SupportMaskNet is used as a SOFT per-gene
        weight in the reconstruction loss, via:
            soft_mask = observed_mask * (0.5 + 0.5 * mask_prob.detach())

        Using .detach() ensures that:
          - mask_net learns support structure from BCE only
          - reconstruction loss doesn't backprop through mask predictions
        This is the "conservative coupling" strategy from the design doc.

        Args:
            x: Input expression (batch, n_genes).
            return_recon: If True, decode and compute reconstruction loss.
            sample_diffusion: If True, run full diffusion loop for z_denoised.
                              If False, use single-step denoised z (training mode).
            mask_override: Optional tensor to use as soft_mask directly.
                           If None, computed from mask_net(x).
        Returns:
            dict with: z, z_denoised, x_recon, losses (mask/recon/diffusion/cluster)
        """
        # ── 1. Predict gene activation probabilities ──────────────────────
        mask_output = self.mask_net(x)
        mask_prob = torch.nan_to_num(
            mask_output["gene_activation_prob"],
            nan=0.5,
            posinf=1.0,
            neginf=0.0,
        ).clamp(1e-6, 1.0 - 1e-6)

        # Ground-truth observed support (for BCE mask loss target). Prefer raw-count
        # support from the unified dataloader; fall back to values > 0 for legacy calls.
        observed_mask = support.float() if support is not None else (x > 0).float()

        # ── 2. Build soft mask for reconstruction loss ─────────────────────
        # Key fix from design doc: mask_prob flows into the main loss pipeline
        if mask_override is not None:
            soft_mask = mask_override
        else:
            # Conservative coupling: observed * (0.5 + 0.5 * mask_prob.detach())
            # All observed non-zero positions participate; high-prob positions weighted more
            soft_mask = observed_mask * (0.5 + 0.5 * mask_prob.detach())

        # ── 3. Forward through latent diffusion AE ─────────────────────────
        ae_result = self.diffusion_ae(
            x=x,
            mask=soft_mask,
            return_recon=return_recon,
            sample_diffusion=sample_diffusion,
        )

        z = ae_result["z"]
        z_denoised = ae_result["z_denoised"]
        x_recon = ae_result["x_recon"]
        loss_diffusion = ae_result["diffusion_loss"]
        loss_recon = ae_result["recon_loss"]

        # ── 4. Compute losses ───────────────────────────────────────────────
        # Mask loss: BCE between predicted probabilities and observed support
        loss_mask = F.binary_cross_entropy(
            mask_prob, observed_mask, reduction="mean"
        )

        # Clustering loss: applied externally via train_epoch staged schedule
        loss_cluster = torch.tensor(0.0, device=x.device)
        if self.training and self._cluster_centers_initialized:
            loss_cluster = self.cluster_loss_fn(z)

        return {
            "z": z,
            "z_denoised": z_denoised,
            "x_recon": x_recon,
            "mask_prob": mask_prob,
            "soft_mask": soft_mask,
            "losses": {
                "mask": loss_mask,
                "recon": loss_recon,
                "diffusion": loss_diffusion,
                "cluster": loss_cluster,
            },
        }

    # ── Embedding extraction ───────────────────────────────────────────────────

    def get_embedding(self, x: torch.Tensor, use_diffusion: bool = False) -> torch.Tensor:
        """Extract cell embedding for downstream clustering.

        Args:
            x: Input expression (batch, n_genes).
            use_diffusion: If False (default), return raw encoder output.
                          If True, run full diffusion sampling (slower but potentially
                          smoother latent space).
        Returns:
            Cell embeddings (batch, latent_dim).
        """
        z = self.diffusion_ae.get_direct_embedding(x)
        if use_diffusion:
            with torch.no_grad():
                z = self.diffusion_ae.denoise_embedding(z)
        return torch.clamp(z, -10.0, 10.0)

    def get_direct_embedding(self, x: torch.Tensor) -> torch.Tensor:
        """Convenience: raw encoder output (no diffusion)."""
        return self.get_embedding(x, use_diffusion=False)

    @torch.no_grad()
    def get_diffusion_embedding(self, x: torch.Tensor) -> torch.Tensor:
        """Convenience: full diffusion-sampled embedding."""
        return self.get_embedding(x, use_diffusion=True)

    # ── Cluster center management ──────────────────────────────────────────────

    def init_cluster_centers(self, embeddings: torch.Tensor):
        """Initialize DEC cluster centers with KMeans (call after pretraining).

        Args:
            embeddings: (n_cells, latent_dim) tensor from pretrained encoder.
        """
        self.cluster_loss_fn.init_centers(embeddings)
        self._cluster_centers_initialized = True

    @property
    def cluster_centers_initialized(self) -> bool:
        return self._cluster_centers_initialized

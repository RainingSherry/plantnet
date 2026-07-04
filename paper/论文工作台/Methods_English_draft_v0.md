# Methods English Draft v0

Updated: 2026-06-26

> This draft follows the user's preferred algorithm-writing structure: natural-language description, LaTeX formulation, and practical code skeleton.

## 1. Problem Formulation

### Natural-Language Description

Let \(X\in\mathbb{R}^{N\times G}\) denote a preprocessed scRNA-seq expression matrix, where \(N\) is the number of cells and \(G\) is the number of genes. The goal is to learn an encoder \(f_\theta\) that maps each cell \(x_i\) to a low-dimensional embedding \(z_i\). These embeddings are then used for downstream clustering. During training, ground-truth cell-type labels are not used. Labels are used only for post hoc evaluation metrics such as ACC, NMI, ARI, and marker-overlap analysis.

The central training signal is masked reconstruction. A subset of expression entries is corrupted, the model observes the corrupted matrix, and it learns to reconstruct the original expression values while predicting which entries were corrupted. CAAM-scMAE extends this framework by adding bi-axial context modeling and a constrained mask selector.

### LaTeX Formulation

The input matrix is:

$$
X = [x_{ij}] \in \mathbb{R}^{N\times G}.
$$

The encoder maps a cell to an embedding:

$$
z_i = f_\theta(\widetilde{x}_i), \qquad z_i\in\mathbb{R}^{d}.
$$

The learned embedding matrix is:

$$
Z = f_\theta(\widetilde{X}) \in \mathbb{R}^{N\times d}.
$$

The downstream clustering function is applied after training:

$$
\widehat{y} = \mathcal{C}(Z),
$$

where \(\mathcal{C}\) may be K-means, Leiden clustering, or another fixed evaluation clustering protocol.

### Practical Code Skeleton

```python
class ClusteringEncoder(nn.Module):
    def forward(self, x):
        """Return cell embeddings used for downstream clustering."""
        raise NotImplementedError


def extract_embeddings(model, dataloader, device):
    model.eval()
    embeddings = []
    with torch.no_grad():
        for batch in dataloader:
            x = batch["x"].to(device)
            z = model.encoder(x)
            embeddings.append(z.cpu())
    return torch.cat(embeddings, dim=0).numpy()
```

## 2. scMAE-Compatible Masked Reconstruction

### Natural-Language Description

The first stage of CAAM-scMAE should reproduce the core scMAE training mechanism as closely as possible. For each gene, expression values are shuffled across cells, producing a gene-wise replacement matrix. A binary mask determines which entries are replaced by shuffled values. The model receives the corrupted expression matrix and learns both to reconstruct the original expression matrix and to predict the binary mask.

This stage is crucial because it establishes whether the local implementation can reproduce a strong masked autoencoding baseline before adding more complex context and adversarial components.

### LaTeX Formulation

For each gene \(j\), define a random permutation over cells:

$$
\pi_j:\{1,\ldots,N\}\rightarrow\{1,\ldots,N\}.
$$

The gene-wise shuffled value is:

$$
X'_{ij}=X_{\pi_j(i),j}.
$$

The binary mask is:

$$
M_{ij}\sim \mathrm{Bernoulli}(\rho),
$$

where \(\rho\) is the mask ratio.

The corrupted input is:

$$
\widetilde{X}=(1-M)\odot X + M\odot X'.
$$

The reconstruction head predicts:

$$
\widehat{X}=g_\theta(Z).
$$

The mask head predicts:

$$
\widehat{M}=h_\theta(Z).
$$

The weighted reconstruction loss is:

$$
\mathcal{L}_{rec}
=
\frac{1}{NG}
\sum_{i=1}^{N}\sum_{j=1}^{G}
\left(1+\lambda M_{ij}\right)
(\widehat{X}_{ij}-X_{ij})^2.
$$

The mask prediction loss is:

$$
\mathcal{L}_{mask}
=
\mathrm{BCEWithLogits}(\widehat{M},M).
$$

The total baseline objective is:

$$
\mathcal{L}_{base}
=
\mathcal{L}_{rec}
+\gamma\mathcal{L}_{mask}.
$$

### Practical Code Skeleton

```python
def gene_wise_shuffle(x: torch.Tensor) -> torch.Tensor:
    batch_size, n_genes = x.shape
    out = torch.empty_like(x)
    for j in range(n_genes):
        perm = torch.randperm(batch_size, device=x.device)
        out[:, j] = x[perm, j]
    return out


def random_mask(x: torch.Tensor, mask_ratio: float) -> torch.Tensor:
    return (torch.rand_like(x) < mask_ratio).float()


def corrupt_with_shuffle(x: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    x_prime = gene_wise_shuffle(x)
    return x * (1.0 - mask) + x_prime * mask


def masked_reconstruction_loss(x, x_hat, mask, mask_logits, lambda_masked=4.0, gamma=1.0):
    rec_weight = 1.0 + lambda_masked * mask
    loss_rec = (rec_weight * (x_hat - x).pow(2)).mean()
    loss_mask = F.binary_cross_entropy_with_logits(mask_logits, mask)
    return loss_rec + gamma * loss_mask
```

## 3. Bi-Axial Context Encoder

### Natural-Language Description

The bi-axial context encoder is designed to make the encoder aware of the two-dimensional structure of scRNA-seq expression matrices. The gene axis models relationships among genes or gene modules within the same cell. The cell axis injects population-level context by allowing query cells to attend to representative context cells. This is inspired by the row/column modeling idea of TabPFN, but it is adapted to unsupervised single-cell representation learning rather than supervised tabular prediction.

To keep the model scalable, the first implementation should avoid tokenizing every expression entry \(x_{ij}\). Instead, genes are compressed into a smaller number of gene modules. The model performs attention over module tokens, not over all gene-cell pairs.

### LaTeX Formulation

Let \(M_g\) be the number of gene modules, with \(M_g\ll G\). Let:

$$
A\in\mathbb{R}^{G\times M_g}
$$

denote the gene-to-module assignment matrix.

For cell \(i\), gene-level hidden states are:

$$
H_{ij} = \phi(x_{ij}) + e_j,
$$

where \(\phi(\cdot)\) projects expression values and \(e_j\) is a gene embedding.

Module tokens are obtained by:

$$
T_{im} = \sum_{j=1}^{G} A_{jm}H_{ij}.
$$

Gene-axis attention updates module tokens within each cell:

$$
T_i^{gene} = \mathrm{Attn}(T_i,T_i,T_i).
$$

Given context tokens \(T_C\) from representative cells, cell-axis attention updates the query cell by:

$$
T_i^{cell} = \mathrm{Attn}(T_i^{gene},T_C,T_C).
$$

The final cell embedding is:

$$
z_i = \mathrm{Pool}(T_i^{cell}).
$$

### Practical Code Skeleton

```python
class BiAxialContextEncoder(nn.Module):
    def __init__(self, n_genes, n_modules=64, d_model=128, n_heads=4):
        super().__init__()
        self.gene_embed = nn.Embedding(n_genes, d_model)
        self.value_proj = nn.Linear(1, d_model)
        self.module_logits = nn.Parameter(torch.randn(n_genes, n_modules))
        self.gene_axis = nn.MultiheadAttention(d_model, n_heads, batch_first=True)
        self.cell_axis = nn.MultiheadAttention(d_model, n_heads, batch_first=True)
        self.out = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Linear(d_model, d_model),
        )

    def tokenize(self, x):
        _, n_genes = x.shape
        gene_ids = torch.arange(n_genes, device=x.device)
        h = self.value_proj(x[..., None]) + self.gene_embed(gene_ids)[None, :, :]
        assignment = torch.softmax(self.module_logits, dim=1)
        return torch.einsum("bgd,gm->bmd", h, assignment)

    def forward(self, x, context_tokens=None):
        tokens = self.tokenize(x)
        tokens, _ = self.gene_axis(tokens, tokens, tokens)
        if context_tokens is not None:
            tokens, _ = self.cell_axis(tokens, context_tokens, context_tokens)
        return self.out(tokens.mean(dim=1))
```

## 4. Constrained Adversarial Mask Selector

### Natural-Language Description

The adversarial mask selector is not a free-form GAN generator. It should be treated as a constrained mask policy that selects which entries should be corrupted. Its purpose is to make the reconstruction task more informative, not to generate arbitrary unrealistic expression values. The replacement values can still come from gene-wise shuffling or matched donor cells; the selector only chooses where to apply corruption.

The selector must obey four constraints:

1. Fixed or controlled mask budget.
2. No use of cell-type labels or evaluation labels.
3. Coverage regularization to avoid always selecting the same genes.
4. Sparsity-aware diagnostics to detect whether selected positions provide real training signal.

### LaTeX Formulation

Let \(s_\phi(x_i,z_i)\in\mathbb{R}^{G}\) be the selector score for cell \(i\). A top-\(k\) mask is:

$$
M_{ij} =
\mathbf{1}
\left[
j\in \mathrm{TopK}(s_\phi(x_i,z_i), k_i)
\right].
$$

The student minimizes:

$$
\min_\theta
\mathcal{L}_{rec}(\theta,\phi)
+\gamma\mathcal{L}_{mask}(\theta,\phi).
$$

The selector maximizes the student loss under constraints:

$$
\max_{\phi\in\Phi}
\mathcal{L}_{rec}(\theta,\phi)
+\gamma\mathcal{L}_{mask}(\theta,\phi)
-\beta\mathcal{R}_{coverage}(\phi)
-\eta\mathcal{R}_{entropy}(\phi).
$$

The legal selector family is:

$$
\Phi=
\left\{
\phi:
\sum_j M_{ij}=k_i,\ 
\mathrm{Label}(x_i)\notin \mathrm{Input}_\phi,\ 
\mathrm{Coverage}(M)\ge \tau
\right\}.
$$

### Practical Code Skeleton

```python
class ConstrainedMaskSelector(nn.Module):
    def __init__(self, n_genes, d_model, hidden=512):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_genes + d_model, hidden),
            nn.ReLU(),
            nn.Linear(hidden, n_genes),
        )

    def forward(self, x, z, mask_ratio=0.3):
        scores = self.net(torch.cat([x, z.detach()], dim=1))
        k = max(1, int(round(mask_ratio * x.shape[1])))
        index = scores.topk(k, dim=1).indices
        mask = torch.zeros_like(x)
        mask.scatter_(1, index, 1.0)
        return mask, scores


def selector_regularization(scores, mask, eps=1e-8):
    probs = torch.softmax(scores, dim=1)
    entropy = -(probs * (probs + eps).log()).sum(dim=1).mean()
    gene_frequency = mask.mean(dim=0)
    coverage_penalty = gene_frequency.pow(2).mean()
    return coverage_penalty - entropy
```

## 5. Training Schedule

### Natural-Language Description

Training should be staged. First, the student model is trained with random masks to establish a stable scMAE-compatible baseline. Second, the bi-axial encoder is introduced while keeping random masks fixed. Third, the mask selector is introduced after a warm-up period so that selector training does not exploit an untrained student.

### LaTeX Formulation

Warm-up stage:

$$
\theta^\star
=
\arg\min_{\theta}
\mathcal{L}_{base}(\theta;M_{random}).
$$

Alternating stage:

$$
\phi^{t+1}
=
\arg\max_{\phi\in\Phi}
\mathcal{L}(\theta^t,\phi),
$$

$$
\theta^{t+1}
=
\arg\min_{\theta}
\mathcal{L}(\theta,\phi^{t+1}).
$$

### Practical Code Skeleton

```python
def train_step_student(batch, model, optimizer, selector=None, mask_ratio=0.3):
    x = batch["x"]
    if selector is None:
        mask = random_mask(x, mask_ratio)
    else:
        with torch.no_grad():
            z_clean = model.encoder(x)
            mask, _ = selector(x, z_clean, mask_ratio)

    x_tilde = corrupt_with_shuffle(x, mask)
    z = model.encoder(x_tilde)
    x_hat = model.decoder(z)
    mask_logits = model.mask_head(z)

    loss = masked_reconstruction_loss(x, x_hat, mask, mask_logits)
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    return float(loss.detach().cpu())
```

## 6. Evaluation Outputs

The implementation should save at least:

```text
embedding_final.npy
metrics.json
runtime.json
param_count.json
mask_diagnostics.json
training_history.json
```

`mask_diagnostics.json` should include:

```text
actual_mask_ratio_global
actual_mask_ratio_observed
gene_mask_frequency
mask_entropy
zero_to_zero_fraction
effective_changed_fraction
```

These diagnostics are necessary because a model can appear to train normally while its mask task provides weak or degenerate signal.

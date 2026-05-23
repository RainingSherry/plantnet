import torch
import torch.nn.functional as F

from ..models import GeneDecoder, LatentDiffusionPrior, SparseEncoder


def train_embedding_model(values, support, mask_probs, device, latent_dim=16, hidden_dim=256, diffusion_hidden_dim=128, diffusion_steps=100, lr=1e-3, weight_decay=1e-5, epochs=50, batch_size=256):
    x = torch.tensor(values, dtype=torch.float32)
    m = torch.tensor(mask_probs, dtype=torch.float32)
    sup = torch.tensor(support, dtype=torch.float32)
    dataset = torch.utils.data.TensorDataset(x, m, sup)
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)

    encoder = SparseEncoder(n_genes=x.shape[1], latent_dim=latent_dim, hidden_dim=hidden_dim).to(device)
    decoder = GeneDecoder(latent_dim=latent_dim, n_genes=x.shape[1], hidden_dim=hidden_dim).to(device)
    diffusion = LatentDiffusionPrior(latent_dim=latent_dim, hidden_dim=diffusion_hidden_dim, timesteps=diffusion_steps).to(device)

    params = list(encoder.parameters()) + list(decoder.parameters()) + list(diffusion.parameters())
    optim = torch.optim.AdamW(params, lr=lr, weight_decay=weight_decay)

    for _ in range(epochs):
        encoder.train(); decoder.train(); diffusion.train()
        for xb, mb, sb in loader:
            xb = xb.to(device)
            mb = mb.to(device)
            sb = sb.to(device)
            z = encoder(xb)
            recon = decoder(z)
            active = torch.clamp(sb * mb, 0.0, 1.0)
            active_weight = active.mean(dim=1, keepdim=True).clamp(min=1e-6)
            loss_recon = (((recon - xb) ** 2) * active).sum(dim=1, keepdim=True) / active_weight
            loss = loss_recon.mean() + diffusion.loss(z)
            optim.zero_grad()
            loss.backward()
            optim.step()

    encoder.eval(); decoder.eval(); diffusion.eval()
    with torch.no_grad():
        z = encoder(x.to(device))
        z_denoised = diffusion.denoise(z).cpu().numpy()
        recon = decoder(torch.tensor(z_denoised, device=device, dtype=torch.float32)).cpu().numpy()
    return {
        "encoder": encoder.cpu(),
        "decoder": decoder.cpu(),
        "diffusion": diffusion.cpu(),
        "embedding": z_denoised,
        "denoised_values": recon,
    }

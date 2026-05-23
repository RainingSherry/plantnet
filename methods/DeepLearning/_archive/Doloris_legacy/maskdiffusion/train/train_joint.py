from .train_mask import train_support_mask
from .train_embedding import train_embedding_model


def run_two_stage_training(args, bundle):
    from ..utils import get_device

    device = get_device(gpu=getattr(args, "gpu", 0), no_cuda=getattr(args, "no_cuda", False))
    mask_model, mask_probs = train_support_mask(
        values=bundle.values,
        support=bundle.support,
        device=device,
        hidden_dim=getattr(args, "hidden_dim", 256),
        lr=getattr(args, "lr", 1e-3),
        epochs=getattr(args, "mask_epochs", 20),
        batch_size=getattr(args, "batch_size", 256),
    )
    embedding_result = train_embedding_model(
        values=bundle.values,
        support=bundle.support,
        mask_probs=mask_probs,
        device=device,
        latent_dim=getattr(args, "latent_dim", 16),
        hidden_dim=getattr(args, "hidden_dim", 256),
        diffusion_hidden_dim=getattr(args, "diffusion_hidden_dim", 128),
        diffusion_steps=getattr(args, "diffusion_steps", 100),
        lr=getattr(args, "lr", 1e-3),
        weight_decay=getattr(args, "weight_decay", 1e-5),
        epochs=getattr(args, "embedding_epochs", 50),
        batch_size=getattr(args, "batch_size", 256),
    )
    return {
        "mask_model": mask_model,
        "mask_probs": mask_probs,
        "embedding": embedding_result["embedding"],
        "denoised_values": embedding_result["denoised_values"],
        "encoder": embedding_result["encoder"],
        "decoder": embedding_result["decoder"],
        "diffusion": embedding_result["diffusion"],
    }

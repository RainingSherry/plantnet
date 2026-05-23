# -*- coding: utf-8 -*-
"""
Utility functions for maskdiffusion.
"""

import os
import numpy as np
import torch


def save_embeddings(embeddings: np.ndarray, labels: np.ndarray, path: str):
    """
    Save embeddings and labels to files.

    Args:
        embeddings: Embedding matrix (n_cells, latent_dim)
        labels: Cell labels (n_cells,)
        path: Directory to save files
    """
    os.makedirs(path, exist_ok=True)
    np.save(os.path.join(path, 'embeddings.npy'), embeddings)
    np.save(os.path.join(path, 'labels.npy'), labels)


def load_embeddings(path: str) -> tuple:
    """
    Load embeddings and labels from files.

    Args:
        path: Directory containing saved files

    Returns:
        tuple: (embeddings, labels)
    """
    embeddings = np.load(os.path.join(path, 'embeddings.npy'))
    labels = np.load(os.path.join(path, 'labels.npy'))
    return embeddings, labels


def save_checkpoint(
    model: torch.nn.Module,
    optimizer: torch.optim.Optimizer,
    epoch: int,
    path: str,
    **kwargs
):
    """
    Save model checkpoint.

    Args:
        model: PyTorch model
        optimizer: Optimizer
        epoch: Current epoch
        path: Path to save checkpoint
        **kwargs: Additional items to save
    """
    checkpoint = {
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'epoch': epoch,
    }
    checkpoint.update(kwargs)
    torch.save(checkpoint, path)


def load_checkpoint(
    model: torch.nn.Module,
    optimizer: torch.optim.Optimizer,
    path: str,
) -> int:
    """
    Load model checkpoint.

    Args:
        model: PyTorch model
        optimizer: Optimizer
        path: Path to checkpoint

    Returns:
        int: Epoch number
    """
    checkpoint = torch.load(path)
    model.load_state_dict(checkpoint['model_state_dict'])
    optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
    return checkpoint.get('epoch', 0)


def set_device(gpu: int = 0, no_cuda: bool = False) -> torch.device:
    """
    Set compute device.

    Args:
        gpu: GPU device number
        no_cuda: Whether to disable CUDA

    Returns:
        torch.device
    """
    cuda = not no_cuda and torch.cuda.is_available()
    if cuda:
        return torch.device(f'cuda:{gpu}')
    return torch.device('cpu')


def count_parameters(model: torch.nn.Module) -> dict:
    """
    Count model parameters.

    Args:
        model: PyTorch model

    Returns:
        dict with parameter counts
    """
    total = sum(p.numel() for p in model.parameters())
    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    return {
        'total': total,
        'trainable': trainable,
        'frozen': total - trainable,
    }

"""cursor2_Doloris maskdiffusion — ScSpade package."""

from .data import load_sc_dataset, SCDatasetBundle
from .models import SparseEncoder, GeneDecoder, LatentDiffusionAE, SupportMaskNet, ScSpade
from .train import run_scspade_training, train_scspade_epoch, extract_embeddings
from .eval import cluster_and_evaluate, benchmark_evaluation, evaluate_support_predictions, marker_gene_enrichment

__all__ = [
    "load_sc_dataset",
    "SCDatasetBundle",
    "SparseEncoder",
    "GeneDecoder",
    "LatentDiffusionAE",
    "SupportMaskNet",
    "ScSpade",
    "run_scspade_training",
    "train_scspade_epoch",
    "extract_embeddings",
    "cluster_and_evaluate",
    "benchmark_evaluation",
    "evaluate_support_predictions",
    "marker_gene_enrichment",
]

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Dict


@dataclass(frozen=True)
class VariantConfig:
    key: str
    rank: int
    title: str
    source_paper: str
    method_name: str
    mask_strategy: str = "random"
    encoder_kind: str = "mlp"
    reconstruction: str = "mse"
    mask_ratio: float = 0.4
    masked_data_weight: float = 0.75
    mask_loss_weight: float = 0.7
    token_weight: float = 0.0
    token_bins: int = 8
    gene_feature_weight: float = 0.0
    consistency_weight: float = 0.0
    teacher_weight: float = 0.0
    graph_weight: float = 0.0
    retrieval_weight: float = 0.0
    barlow_weight: float = 0.0
    prototype_weight: float = 0.0
    fuzzy_weight: float = 0.0
    gate_l1_weight: float = 0.0
    dropout: float = 0.0
    neighbor_k: int = 15
    note: str = ""

    def with_cli_overrides(self, **kwargs) -> "VariantConfig":
        updates = {key: value for key, value in kwargs.items() if value is not None}
        return replace(self, **updates)


def _cfg(
    key: str,
    rank: int,
    title: str,
    source: str,
    **kwargs,
) -> VariantConfig:
    return VariantConfig(
        key=key,
        rank=rank,
        title=title,
        source_paper=source,
        method_name=f"scMAEs_{key}",
        **kwargs,
    )


VARIANTS: Dict[str, VariantConfig] = {
    "001_scdiva_discrete_diffusion": _cfg(
        "001_scdiva_discrete_diffusion", 1,
        "ScDiVa-inspired masked discrete denoising",
        "ScDiVa: Masked Discrete Diffusion for Joint Modeling of Single-Cell Identity and Expression",
        mask_strategy="dropout_adaptive", token_weight=0.35, reconstruction="huber",
    ),
    "002_maskfeat_gene_features": _cfg(
        "002_maskfeat_gene_features", 2,
        "MaskFeat-inspired masked gene feature prediction",
        "MaskFeat: Masked Feature Prediction for Self-Supervised Visual Pre-Training",
        mask_strategy="variance_adaptive", gene_feature_weight=0.25, reconstruction="huber",
    ),
    "003_consistency_teacher": _cfg(
        "003_consistency_teacher", 3,
        "Consistency-model-inspired EMA teacher",
        "Improved Techniques for Training Consistency Models",
        mask_strategy="random", consistency_weight=0.15, teacher_weight=0.25,
    ),
    "004_joao_graph_aug": _cfg(
        "004_joao_graph_aug", 4,
        "JOAO-inspired adaptive graph/mask augmentation",
        "JOAO: Automated Data Augmentations for Graph Contrastive Learning",
        mask_strategy="joao", graph_weight=0.18, consistency_weight=0.08,
    ),
    "005_tabr_retrieval": _cfg(
        "005_tabr_retrieval", 5,
        "TabR-inspired nearest-neighbor retrieval context",
        "TabR: Tabular Deep Learning Meets Nearest Neighbors",
        mask_strategy="random", retrieval_weight=0.30, neighbor_k=20,
    ),
    "006_scvgae_zinb_graph": _cfg(
        "006_scvgae_zinb_graph", 6,
        "scVGAE-inspired robust graph reconstruction",
        "scVGAE: ZINB-Based Variational Graph Autoencoder for Single-Cell RNA-Seq Imputation",
        mask_strategy="dropout_adaptive", reconstruction="huber", graph_weight=0.20, neighbor_k=15,
    ),
    "007_dino_self_distill": _cfg(
        "007_dino_self_distill", 7,
        "DINO-inspired self-distillation",
        "DinoBloom: A Foundation Model for Generalizable Cell Embeddings in Hematology",
        mask_strategy="random", teacher_weight=0.35, consistency_weight=0.10,
    ),
    "008_cell_hierarchy_proto": _cfg(
        "008_cell_hierarchy_proto", 8,
        "Cell-hierarchy-inspired prototype regularization",
        "scCello: Cell-ontology guided transcriptome foundation model",
        mask_strategy="variance_adaptive", prototype_weight=0.16, fuzzy_weight=0.05,
    ),
    "010_longtail_prototype": _cfg(
        "010_longtail_prototype", 10,
        "Long-tail prototype balancing",
        "Celler: A Genomic Language Model for Long-Tailed Single-Cell Annotation",
        mask_strategy="marker_safe", prototype_weight=0.22, fuzzy_weight=0.08,
    ),
    "011_mamba_gated_sequence": _cfg(
        "011_mamba_gated_sequence", 11,
        "scMamba-inspired gated sequence encoder",
        "scMamba: Scalable Foundation Model for Single-Cell Multi-Omics Integration",
        encoder_kind="gated_sequence", mask_strategy="variance_adaptive", dropout=0.05,
    ),
    "013_masked_sc_clustering": _cfg(
        "013_masked_sc_clustering", 13,
        "Cluster-aware masked single-cell modeling",
        "Masked Modeling for Single-cell Clustering of scRNA-seq Data",
        mask_strategy="variance_adaptive", prototype_weight=0.12, consistency_weight=0.08,
    ),
    "014_cicl_iter_contrast": _cfg(
        "014_cicl_iter_contrast", 14,
        "CICL-inspired iterative contrastive regularization",
        "CICL: scRNA-seq Data Clustering by Cluster-aware Iterative Contrastive Learning",
        mask_strategy="random", consistency_weight=0.12, barlow_weight=0.06,
    ),
    "015_scagc_adaptive_graph": _cfg(
        "015_scagc_adaptive_graph", 15,
        "scAGC-inspired adaptive cell graph guidance",
        "scAGC: Learning Adaptive Cell Graphs with Contrastive Guidance for Single-Cell Clustering",
        mask_strategy="dropedge", graph_weight=0.26, consistency_weight=0.08, neighbor_k=20,
    ),
    "017_interpretable_sparse_gate": _cfg(
        "017_interpretable_sparse_gate", 17,
        "Interpretable sparse gene gate",
        "Interpretable Deep Learning in Single-Cell Omics",
        mask_strategy="marker_safe", gate_l1_weight=2e-4, gene_feature_weight=0.08,
    ),
    "018_beit_tokenizer": _cfg(
        "018_beit_tokenizer", 18,
        "BEiT-inspired expression tokenizer",
        "BEiT: BERT Pre-Training of Image Transformers",
        mask_strategy="random", token_weight=0.28,
    ),
    "019_data2vec_ema": _cfg(
        "019_data2vec_ema", 19,
        "data2vec-inspired EMA latent target",
        "data2vec: A General Framework for Self-supervised Learning in Speech, Vision and Language",
        mask_strategy="variance_adaptive", teacher_weight=0.28,
    ),
    "020_multimae_targets": _cfg(
        "020_multimae_targets", 20,
        "MultiMAE-inspired multi-target reconstruction",
        "MultiMAE: Multi-modal Multi-task Masked Autoencoders",
        mask_strategy="module_block", gene_feature_weight=0.14, token_weight=0.12, reconstruction="huber",
    ),
    "022_ijepa_latent": _cfg(
        "022_ijepa_latent", 22,
        "I-JEPA-inspired latent prediction",
        "I-JEPA: Self-Supervised Learning from Images with Joint-Embedding Predictive Architecture",
        mask_strategy="module_block", teacher_weight=0.18, consistency_weight=0.12,
    ),
    "023_maskgit_iterative": _cfg(
        "023_maskgit_iterative", 23,
        "MaskGIT-inspired iterative token denoising",
        "MaskGIT: Masked Generative Image Transformer",
        mask_strategy="high_mask_curriculum", token_weight=0.32,
    ),
    "026_bgrl_graph_bootstrap": _cfg(
        "026_bgrl_graph_bootstrap", 26,
        "BGRL-inspired graph bootstrapping",
        "BGRL: Large-Scale Representation Learning on Graphs via Bootstrapping",
        mask_strategy="dropedge", graph_weight=0.16, teacher_weight=0.22, neighbor_k=20,
    ),
    "027_graph_barlow": _cfg(
        "027_graph_barlow", 27,
        "Graph Barlow Twins redundancy reduction",
        "Graph Barlow Twins",
        mask_strategy="dropedge", graph_weight=0.12, barlow_weight=0.10,
    ),
    "028_graphormer_bias": _cfg(
        "028_graphormer_bias", 28,
        "Graphormer-inspired graph-distance bias",
        "Graphormer",
        mask_strategy="dropedge", retrieval_weight=0.18, graph_weight=0.18, neighbor_k=25,
    ),
    "029_fuzzy_clustering": _cfg(
        "029_fuzzy_clustering", 29,
        "Deep adaptive fuzzy clustering head",
        "Deep Adaptive Fuzzy Clustering for Evolutionary Unsupervised Representation Learning",
        mask_strategy="variance_adaptive", fuzzy_weight=0.22, prototype_weight=0.10,
    ),
    "030_fuzzy_rough_boundary": _cfg(
        "030_fuzzy_rough_boundary", 30,
        "Fuzzy rough boundary-aware regularization",
        "Fuzzy Rough Sets Based on Fuzzy Quantification",
        mask_strategy="marker_safe", fuzzy_weight=0.18, prototype_weight=0.14, consistency_weight=0.06,
    ),
}


def get_variant(key: str) -> VariantConfig:
    if key not in VARIANTS:
        known = ", ".join(sorted(VARIANTS))
        raise KeyError(f"Unknown scMAEs variant {key!r}. Known variants: {known}")
    return VARIANTS[key]


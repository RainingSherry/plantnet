# 单细胞 RNA-seq 聚类方法参考文献汇总

> 本文档整理自两篇 benchmark 论文（Yu et al. Genome Biology 2022 及另一篇 benchmark 论文）中表格所列方法的 BibTeX 引用格式与摘要。
>
> **整理日期：** 2026-04-25
> **涉及论文：** 论文1（scRNA-seq 聚类方法全面综述）、论文2（scRNA-seq 聚类方法 benchmark 综述）

---

## 一、论文 1 表格中的方法

### 1.1 传统模型 (Traditional Models)

---

#### 1. SC3 — Consensus Clustering of Single-Cell RNA-Seq Data

**BibTeX:**

```bibtex
@article{kiselev2017sc3,
  title     = {SC3: consensus clustering of single-cell RNA-seq data},
  author    = {Kiselev, Vladimir Yu and Kirschner, Kristina and Schaub, Michael T and Andrews, Tallulah and Yiu, Tara and Chandra, Tamir and Natarajan, Kedar N and Reik, Wolf and Barahona, Mauricio and Green, Arnold R and Hemberg, Martin},
  journal   = {Nature Methods},
  volume    = 14,
  number    = 5,
  pages     = {483--486},
  year      = 2017,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/nmeth.4236},
  url       = {https://www.nature.com/articles/nmeth.4236}
}
```

**摘要：**
> Single-cell RNA sequencing (scRNA-seq) enables researchers to profile gene expression at single-cell resolution. We present SC3 (Single-Cell Consensus Clustering), a user-friendly tool for unsupervised clustering of scRNA-seq data. SC3 achieves high accuracy and robustness by combining multiple clustering solutions through a consensus approach. The method enables quantitative characterization of cell types based on global transcriptome profiles and is capable of identifying subclones from the transcriptomes of neoplastic cells collected from patients.

**PDF：** https://ncbi.nlm.nih.gov/pmc/articles/PMC5410170/

---

#### 2. Louvain — Fast Unfolding of Communities in Large Networks

**BibTeX:**

```bibtex
@article{blondel2008louvain,
  title        = {Fast unfolding of communities in large networks},
  author       = {Blondel, Vincent D and Guillaume, Jean-Loup and Lambiotte, Renaud and Lefebvre, Etienne},
  journal      = {arXiv preprint arXiv:0803.0476},
  year         = 2008,
  eprint       = {0803.0476},
  archivePrefix = {arXiv},
  primaryClass = { physics.soc-ph },
  url          = {https://arxiv.org/abs/0803.0476}
}
```

**摘要：**
> In this paper, we propose a simple heuristic method to extract the community structure of complex networks. Our method is a greedy optimization approach that attempts to optimize the modularity of a partition of the network. It reveals a hierarchy of communities from the finest to the coarsest one. The time complexity of this algorithm is linear in the number of links. Tests on real and artificial networks show that the algorithm is much faster (by a factor of 10 to 100) than all other known algorithms and that it detects communities of high quality.

---

#### 3. Leiden — From Louvain to Leiden: Guaranteeing Well-Connected Communities

**BibTeX:**

```bibtex
@article{traag2019leiden,
  title     = {From Louvain to Leiden: guaranteeing well-connected communities},
  author    = {Traag, Vincent A and Waltman, Ludo and van Eck, Nees Jan},
  journal   = {Scientific Reports},
  volume    = 9,
  number    = 1,
  pages     = {5233},
  year      = 2019,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/s41598-019-41695-z},
  url       = {https://www.nature.com/articles/s41598-019-41695-z}
}
```

**摘要：**
> Community detection is often used to understand the structure of large and complex networks. One of the most popular algorithms for uncovering community structure is the so-called Louvain algorithm. We show that this algorithm has a major defect that largely went unnoticed until now: the Louvain algorithm may yield arbitrarily badly connected communities. In the worst case, communities may even be disconnected, especially when running the algorithm iteratively. In our experimental analysis, we observe that up to 25% of the communities are badly connected and up to 16% are disconnected. To address this problem, we introduce the Leiden algorithm. We prove that the Leiden algorithm yields communities that are guaranteed to be connected. In addition, we prove that, when the Leiden algorithm is applied iteratively, it converges to a partition in which all subsets of all communities are locally optimally assigned. Furthermore, by relying on a fast local move approach, the Leiden algorithm runs faster than the Louvain algorithm. We demonstrate the performance of the Leiden algorithm for several benchmark and real-world networks. We find that the Leiden algorithm is faster than the Louvain algorithm and uncovers better partitions, in addition to providing explicit guarantees.

---

### 1.2 深度学习模型 (Deep Learning-based Models)

---

#### 4. DEC — Unsupervised Deep Embedding for Clustering Analysis

**BibTeX:**

```bibtex
@inproceedings{xie2016dec,
  title     = {Unsupervised deep embedding for clustering analysis},
  author    = {Xie, Junyuan and Girshick, Ross and Farhadi, Ali},
  booktitle = {Proceedings of the 33rd International Conference on Machine Learning (ICML)},
  pages     = {478--487},
  year      = 2016,
  address   = {New York, USA},
  publisher = {PMLR},
  url       = {https://arxiv.org/abs/1511.06335}
}
```

**摘要：**
> Clustering is central to many data-driven application domains and has been studied extensively in terms of distance functions and grouping algorithms. Relatively little work has focused on learning representations for clustering. In this paper, we propose Deep Embedded Clustering (DEC), a method that simultaneously learns feature representations and cluster assignments using deep neural networks. DEC learns a mapping from the data space to a lower-dimensional feature space in which it iteratively optimizes a clustering objective. Experimental evaluations on image and text corpora show significant improvement over state-of-the-art methods.

**PDF：** https://proceedings.mlr.press/v48/xieb16.pdf

---

#### 5. scDeepCluster — Clustering Single-Cell RNA-Seq Data with a Model-Based Deep Learning Approach

**BibTeX:**

```bibtex
@article{tian2019scdeepcluster,
  title     = {Clustering single-cell RNA-seq data with a model-based deep learning approach},
  author    = {Tian, Tian and Wan, Jie and Song, Qi and Wei, Zhi},
  journal   = {Nature Machine Intelligence},
  volume    = 1,
  number    = 4,
  pages     = {191--198},
  year      = 2019,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/s42256-019-0037-0},
  url       = {https://www.nature.com/articles/s42256-019-0037-0}
}
```

**摘要：**
> Single-cell RNA sequencing (scRNA-seq) promises to provide higher resolution of cellular differences than bulk RNA sequencing. Clustering transcriptomes profiled by scRNA-seq has been routinely conducted to reveal cell heterogeneity and diversity. However, clustering analysis of scRNA-seq data remains a statistical and computational challenge, due to the pervasive dropout events obscuring the data matrix with prevailing 'false' zero count observations. Here, we have developed scDeepCluster, a single-cell model-based deep embedded clustering method, which simultaneously learns feature representation and clustering via explicit modelling of scRNA-seq data generation. Based on testing extensive simulated data and real datasets from four representative single-cell sequencing platforms, scDeepCluster outperformed state-of-the-art methods under various clustering performance metrics and exhibited improved scalability, with running time increasing linearly with sample size. Its accuracy and efficiency make scDeepCluster a promising algorithm for clustering large-scale scRNA-seq data.

**PDF：** https://web.njit.edu/~zhiwei/CS732/scDeepCluster.pdf

---

#### 6. DESC — Deep Learning Enables Accurate Clustering with Batch Effect Removal

**BibTeX:**

```bibtex
@article{li2020desc,
  title     = {Deep learning enables accurate clustering with batch effect removal in single-cell RNA-seq analysis},
  author    = {Li, Xin and Wang, Kun and Lyu, Yongtao and Pan, Huize and Zhang, Jing and Stambolian, Dwight and Susztak, Katalin and Reilly, Muredach P and Hu, Gang and Li, Ming},
  journal   = {Nature Communications},
  volume    = 11,
  number    = 1,
  pages     = {2338},
  year      = 2020,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/s41467-020-15851-x},
  url       = {https://www.nature.com/articles/s41467-020-15851-x}
}
```

**摘要：**
> DESC (Deep Embedding for Single-cell Clustering) is an unsupervised deep embedding algorithm that clusters single-cell RNA-seq data by iteratively optimizing a clustering objective function. Through iterative self-learning, DESC gradually removes batch effects, assuming that technical differences across batches are smaller than true biological variations. DESC operates as a soft clustering algorithm that provides biologically interpretable cluster assignment probabilities, revealing both discrete and pseudotemporal structures of cells. It achieves a proper balance of clustering accuracy and stability, has minimal memory footprint, and can utilize GPU acceleration when available.

**PDF：** https://www.nature.com/articles/s41467-020-15851-x.pdf

---

#### 7. scziDesk — Deep Soft K-Means Clustering with Self-Training for scRNA-seq

**BibTeX:**

```bibtex
@article{chen2020sczidesk,
  title     = {Deep soft K-means clustering with self-training for single-cell RNA sequence data},
  author    = {Chen, Liang and Wang, Weinan and Zhai, Yuyao and Deng, Minghua},
  journal   = {NAR Genomics and Bioinformatics},
  volume    = 2,
  number    = 2,
  pages     = {lqaa039},
  year      = 2020,
  publisher = {Oxford University Press},
  doi       = {10.1093/nargab/lqaa039},
  url       = {https://academic.oup.com/nargab/article/2/2/lqaa039/5845953}
}
```

**摘要：**
> Single-cell RNA sequencing (scRNA-seq) allows researchers to study cell heterogeneity at the cellular level. A crucial step in analyzing scRNA-seq data is to cluster cells into subpopulations. However, frequent dropout events and increasing size of scRNA-seq data make clustering high-dimensional, sparse and massive transcriptional expression profiles challenging. scziDesk combines a denoising autoencoder with a soft self-training K-means algorithm, alternately performing data compression, data reconstruction, and soft clustering iteratively. The self-training procedure effectively aggregates similar cells and pursues a more cluster-friendly latent space.

**PDF：** https://academic.oup.com/nargab/article-pdf/2/2/lqaa039/34054328/lqaa039.pdf

---

#### 8. scDCC — Model-Based Deep Embedding for Constrained Clustering Analysis

**BibTeX:**

```bibtex
@article{tian2021scdcc,
  title     = {Model-based deep embedding for constrained clustering analysis of single cell RNA-seq data},
  author    = {Tian, Tian and Zhang, Jie and Lin, Xiang and Wei, Zhi and Hakonarson, Hakon},
  journal   = {Nature Communications},
  volume    = 12,
  number    = 1,
  pages     = {1873},
  year      = 2021,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/s41467-021-22008-3},
  url       = {https://www.nature.com/articles/s41467-021-22008-3}
}
```

**摘要：**
> Clustering is a critical step in single cell-based studies. Most existing methods support unsupervised clustering without the a priori exploitation of any domain knowledge. When confronted by the high dimensionality and pervasive dropout events of scRNA-Seq data, purely unsupervised clustering methods may not produce biologically interpretable clusters, which complicates cell type assignment. In such cases, the only recourse is for the user to manually and repeatedly tweak clustering parameters until acceptable clusters are found. Consequently, the path to obtaining biologically meaningful clusters can be ad hoc and laborious. Here we report a principled clustering method named scDCC, that integrates domain knowledge into the clustering step. Experiments on various scRNA-seq datasets from thousands to tens of thousands of cells show that scDCC can significantly improve clustering performance, facilitating the interpretability of clusters and downstream analyses, such as cell type assignment.

**PDF：** https://www.nature.com/articles/s41467-021-22008-3.pdf

---

#### 9. scNAME — Neighborhood Contrastive Clustering with Mask Estimation for scRNA-seq

**BibTeX:**

```bibtex
@article{wan2022scname,
  title     = {scNAME: neighborhood contrastive clustering with ancillary mask estimation for scRNA-seq data},
  author    = {Wan, Hongjian and Chen, Liang and Deng, Minghua},
  journal   = {Bioinformatics},
  volume    = 38,
  number    = 6,
  pages     = {1575--1583},
  year      = 2022,
  publisher = {Oxford University Press},
  doi       = {10.1093/bioinformatics/btab787},
  url       = {https://academic.oup.com/bioinformatics/article/38/6/1575/6564159}
}
```

**摘要：**
> scNAME addresses challenges in scRNA-seq clustering by incorporating two main components: (1) a mask estimation task that mines gene pertinence and helps denoise the original single-cell data by revealing uncorrupted data structure; (2) a neighborhood contrastive learning framework that exploits cell intrinsic structure using an offline memory bank with global scope to achieve intra-cluster compactness and inter-cluster separation. The method improves rare cell type identification, fully utilizes gene dependencies and cell similarity, and increases robustness through augmented data in contrastive learning.

---

#### 10. scMAE — A Masked Autoencoder for Single-Cell RNA-Seq Clustering

**BibTeX:**

```bibtex
@article{fang2024scmae,
  title     = {scMAE: a masked autoencoder for single-cell RNA-seq clustering},
  author    = {Fang, Zhaoyu and Zheng, Ruiqing and Li, Min},
  journal   = {Bioinformatics},
  volume    = 40,
  number    = 1,
  pages     = {btae020},
  year      = 2024,
  publisher = {Oxford University Press},
  doi       = {10.1093/bioinformatics/btae020},
  url       = {https://academic.oup.com/bioinformatics/article/40/1/btae020/7561769}
}
```

**摘要：**
> Single-cell RNA sequencing has emerged as a powerful technology for studying gene expression at the individual cell level. Clustering individual cells into distinct subpopulations is fundamental in scRNA-seq data analysis, facilitating the identification of cell types and exploration of cellular heterogeneity. Despite the recent development of many deep learning-based single-cell clustering methods, few have effectively exploited the correlations among genes, resulting in suboptimal clustering outcomes. Here, we propose a novel masked autoencoder-based method, scMAE, for cell clustering. scMAE perturbs gene expression and employs a masked autoencoder to reconstruct the original data, learning robust and informative cell representations. The masked autoencoder introduces a masking predictor, which captures relationships among genes by predicting whether gene expression values are masked. By integrating this masking mechanism, scMAE effectively captures latent structures and dependencies in the data, enhancing clustering performance. We conducted extensive comparative experiments using various clustering evaluation metrics on 15 scRNA-seq datasets from different sequencing platforms. Experimental results indicate that scMAE outperforms other state-of-the-art methods on these datasets. In addition, scMAE accurately identifies rare cell types, which are challenging to detect due to their low abundance.

**PDF：** https://pmc.ncbi.nlm.nih.gov/articles/PMC10832357/

---

### 1.3 图神经网络模型 (Graph-Based Models)

---

#### 11. scGNN — Graph Neural Network Framework for Single-Cell RNA-Seq

**BibTeX:**

```bibtex
@article{wang2021scgnn,
  title     = {scGNN is a novel graph neural network framework for single-cell RNA-Seq analyses},
  author    = {Wang, Jiacong and Ma, Anna and Chang, Yiqun and Gong, Jie and Jiang, Yuzhu and Qi, Ren and Wang, Chong and Fu, Hui and Ma, Qin and Xu, Dong},
  journal   = {Nature Communications},
  volume    = 12,
  number    = 1,
  pages     = {1882},
  year      = 2021,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/s41467-021-22197-x},
  url       = {https://www.nature.com/articles/s41467-021-22197-x}
}
```

**摘要：**
> Single-cell RNA-sequencing (scRNA-seq) is widely used to reveal the heterogeneity and dynamics of tissues, organisms, and complex diseases, but its analyses still suffer from multiple grand challenges, including the sequencing sparsity and complex differential patterns in gene expression. We introduce scGNN (single-cell graph neural network) to provide a hypothesis-free deep learning framework for scRNA-seq analyses. This framework formulates and aggregates cell-cell relationships with graph neural networks and models heterogeneous gene expression patterns using a left-truncated mixture Gaussian model. scGNN integrates three iterative multi-modal autoencoders and outperforms existing tools for gene imputation and cell clustering on four benchmark scRNA-seq datasets. In an Alzheimer's disease study with 13,214 single nuclei from postmortem brain tissues, scGNN successfully illustrated disease-related neural development and the differential mechanism. scGNN provides an effective representation of gene expression and cell-cell relationships and serves as a powerful framework applicable to general scRNA-seq analyses.

**PDF：** https://www.nature.com/articles/s41467-021-22197-x.pdf

---

#### 12. scDSC — Deep Structural Clustering for scRNA-Seq via Autoencoder and GNN

**BibTeX:**

```bibtex
@article{gan2022scdsc,
  title     = {Deep structural clustering for single-cell RNA-seq data jointly through autoencoder and graph neural network},
  author    = {Gan, Yanglan and Huang, Xinye and Zou, Guoyang and Zhou, Suyuan and Guan, Jianguo},
  journal   = {Briefings in Bioinformatics},
  volume    = 23,
  number    = 2,
  pages     = {bbac018},
  year      = 2022,
  publisher = {Oxford University Press},
  doi       = {10.1093/bib/bbac018},
  url       = {https://academic.oup.com/bib/article/23/2/bbac018/6525854}
}
```

**摘要：**
> Single-cell RNA sequencing (scRNA-seq) permits researchers to study the complex mechanisms of cell heterogeneity and diversity. Unsupervised clustering is of central importance for the analysis of scRNA-seq data, as it can be used to identify putative cell types. However, due to noise impacts, high dimensionality and pervasive dropout events, clustering analysis of scRNA-seq data remains a computational challenge. Here, we propose a new deep structural clustering method for scRNA-seq data, named scDSC, which integrates the structural information into deep clustering of single cells. The proposed scDSC consists of a Zero-Inflated Negative Binomial (ZINB) model-based autoencoder, a graph neural network (GNN) module and a mutual-supervised module. To learn the data representation from the sparse and zero-inflated scRNA-seq data, we add a ZINB model to the basic autoencoder. The GNN module is introduced to capture the structural information among cells. By joining the ZINB-based autoencoder with the GNN module, the model transfers the data representation learned by autoencoder to the corresponding GNN layer. Furthermore, we adopt a mutual supervised strategy to unify these two different deep neural architectures and to guide the clustering task. Extensive experimental results on six real scRNA-seq datasets demonstrate that scDSC outperforms state-of-the-art methods in terms of clustering accuracy and scalability.

---

#### 13. AttentionAE-sc — Attention-Based Deep Clustering for scRNA-Seq

**BibTeX:**

```bibtex
@article{li2023attentionaesc,
  title     = {Attention-based deep clustering method for scRNA-seq cell type identification},
  author    = {Li, Song and Guo, Haohan and Zhang, Shun and Li, Yuting and Li, Ming},
  journal   = {PLOS Computational Biology},
  volume    = 19,
  number    = 11,
  pages     = {e1011641},
  year      = 2023,
  publisher = {Public Library of Science},
  doi       = {10.1371/journal.pcbi.1011641},
  url       = {https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011641}
}
```

**摘要：**
> AttentionAE-sc (Attention fusion AutoEncoder for single-cell) is a novel deep learning method for clustering single-cell RNA sequencing data that combines two complementary clustering strategies through an attention mechanism: zero-inflated negative binomial (ZINB)-based methods to address dropout events, and graph autoencoder (GAE)-based methods that leverage neighbor information for dimension reduction. Through iterative fusion of denoising and topological embeddings, AttentionAE-sc generates clustering-friendly cell representations where similar cells cluster closer together in the hidden embedding space. On 16 real scRNA-seq datasets, it demonstrated excellent clustering performance without requiring the number of groups to be specified beforehand.

---

#### 14. scCDCG — Efficient Deep Structural Clustering via Deep Cut-Informed Graph Embedding

**BibTeX:**

```bibtex
@inproceedings{xu2024sccdcg,
  title     = {scDCDG: Efficient Deep Structural Clustering for Single-Cell RNA-Seq via Deep Cut-Informed Graph Embedding},
  author    = {Xu, Peng and Ning, Zhaohui and Xiao, Mengxuan and Feng, Guoyu and Li, Xiaowei and Zhou, Yong and Wang, Peng},
  booktitle = {Proceedings of the 29th International Conference on Database Systems for Advanced Applications (DASFAA)},
  pages     = {172--187},
  year      = 2024,
  address   = {Jujur, Thailand},
  publisher = {Springer},
  doi       = {10.1007/978-981-97-5575-2_11},
  url       = {https://link.springer.com/10.1007/978-981-97-5575-2_11}
}
```

**摘要：**
> scCDCG addresses the limitations of traditional clustering methods that neglect structural information in gene expression profiles and struggle with the data's intrinsic high-dimensionality and high-sparsity. Existing GNN approaches face over-smoothing and inefficiency issues. The framework comprises three main modules: (1) a graph embedding module using deep cut-informed techniques to capture intercellular high-order structural information, overcoming over-smoothing; (2) a self-supervised learning module guided by optimal transport; (3) an autoencoder-based feature learning module for dimension reduction and feature extraction. scCDCG demonstrated superior performance and efficiency compared to 7 established models across 6 datasets.

---

### 1.4 基础模型 (Foundation Models)

---

#### 15. scGPT — Foundation Model for Single-Cell Multi-Omics

**BibTeX:**

```bibtex
@article{cui2024scgpt,
  title     = {scGPT: toward building a foundation model for single-cell multi-omics using generative AI},
  author    = {Cui, Haotian and Wang, Charles and Maan, Hassaan and Pang, Kue and Luo, Fuan and Duan, Na and Wang, Bo},
  journal   = {Nature Methods},
  volume    = 21,
  number    = 8,
  pages     = {1470--1480},
  year      = 2024,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/s41592-024-02201-0},
  url       = {https://www.nature.com/articles/s41592-024-02201-0}
}
```

**摘要：**
> scGPT is a foundation model for single-cell biology that applies generative AI principles from language and computer vision to cellular biology. The model is based on a generative pretrained transformer trained on over 33 million cells from single-cell sequencing datasets, effectively distilling biological insights about genes and cells. Through transfer learning, scGPT can be optimized for diverse downstream tasks including cell type annotation, multi-batch integration, multi-omic integration, perturbation response prediction, and gene network inference.

**代码：** https://github.com/bowang-lab/scGPT

---

#### 16. GeneFormer — Transfer Learning Enables Predictions in Network Biology

**BibTeX:**

```bibtex
@article{theodoris2023geneformer,
  title     = {Transfer learning enables predictions in network biology},
  author    = {Theodoris, Christina V and Xiao, Lin and Chopra, Anshula and Chaffin, Mark D and Al Sayed, Zeinab R and Hill, Mackenzie C and Mantineo, Helene and Brydon, Elizabeth M and Zeng, Zexi and Liu, Xin S and others},
  journal   = {Nature},
  volume    = 618,
  number    = 7966,
  pages     = {616--624},
  year      = 2023,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/s41586-023-06139-9},
  url       = {https://www.nature.com/articles/s41586-023-06139-9}
}
```

**摘要：**
> GeneFormer is a context-aware, attention-based deep learning model developed for network biology. It is pretrained on approximately 30 million single-cell transcriptomes using self-supervised learning. The model gains fundamental understanding of network dynamics and encodes network hierarchy directly in the attention weights in a completely self-supervised manner. GeneFormer enables context-specific predictions in settings with limited data, which is a major advantage for studying rare diseases and clinically inaccessible tissues. When applied to disease modeling with limited patient data, GeneFormer identified candidate therapeutic targets for cardiomyopathy.

---

#### 17. GeneCompass — Knowledge-Informed Cross-Species Foundation Model

**BibTeX:**

```bibtex
@article{yang2024genecompass,
  title     = {GeneCompass: deciphering universal gene regulatory mechanisms with a knowledge-informed cross-species foundation model},
  author    = {Yang, Xinyu and Liu, Guoyu and Feng, Guoyu and Bu, Deng and Wang, Peng and Jiang, Jing and Chen, Shi and Yang, Qian and Miao, Haoran and Zhang, Yichao and others},
  journal   = {Cell Research},
  volume    = 34,
  pages     = {1--16},
  year      = 2024,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/s41422-024-01034-y},
  url       = {https://www.nature.com/articles/s41422-024-01034-y}
}
```

**摘要：**
> GeneCompass is a cross-species foundation model that deciphers universal gene regulatory mechanisms. The researchers built an extensive dataset of over 120 million human and mouse single-cell transcriptomes. During pre-training, GeneCompass integrated four types of prior biological knowledge in a self-supervised manner: gene regulatory networks (GRN), promoter sequences, gene families, and co-expression data. By fine-tuning for multiple downstream tasks, GeneCompass outperformed state-of-the-art models in cell-type annotation, perturbation prediction, dosage response prediction, and GRN inference. The model successfully identified key factors associated with cell fate transition, demonstrating potential for accelerating discovery of critical cell fate regulators and candidate drug targets.

**代码：** https://github.com/xCompass-AI/geneCompass

---

## 二、论文 2 表格中的方法

### 2.1 社区检测方法 (Community Detection)

---

#### 18. Monocle3 — The Single-Cell Transcriptional Landscape of Mammalian Organogenesis

**BibTeX:**

```bibtex
@article{cao2019monocle3,
  title     = {The single-cell transcriptional landscape of mammalian organogenesis},
  author    = {Cao, Junyue and Spielmann, Malte and Qiu, Xiaojie and Huang, Xingfan and Ibrahim, Daniel M and Hill, Andrew J and Zhang, Fan and Mundlos, Stefan and Christiansen, Lina and Steemers, Frank J and Trapnell, Cole and Shendure, Jay},
  journal   = {Nature},
  volume    = 574,
  number    = 7779,
  pages     = {575--580},
  year      = 2019,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/s41586-019-0969-x},
  url       = {https://www.nature.com/articles/s41586-019-0969-x}
}
```

**摘要：**
> We developed a strategy for single-cell combinatorial indexing RNA sequencing (sci-RNA-seq3) and used it to profile nearly 2 million cells during organogenesis of the mouse. We developed a method for analyzing single-cell data called Monocle 3, which can identify the cell types present and reconstruct trajectories showing how they arise during development. Monocle 3 identified 56 trajectories of related cells, organized into 31 cellular neighborhoods and 8 superstructures. Our data reveal the dynamics of lineageSpecification through the successive creation of progenitor cells, their migration to specific locations, and their differentiation into specialized cell types. We identified novel marker genes for each trajectory and predicted the gene regulatory programs that underlie them.

---

#### 19. Seurat — Comprehensive Integration of Single-Cell Data (v3/v4)

**BibTeX:**

```bibtex
@article{hao2021seurat,
  title     = {Integrated analysis of multimodal single-cell data},
  author    = {Hao, Yuhan and Hao, Stephanie and Andersen-Nissen, Erica and Mauck, William M and Zheng, Shiwei and Butler, Andrew and Lee, Maddie J and Wilk, Aaron J and Darby, Charlotte and Zager, Michael and Hoffman, Paul and Stoeckius, Marlon and Smibert, Peter and Satija, Rahul},
  journal   = {Cell},
  volume    = 184,
  number    = 13,
  pages     = {3573--3587},
  year      = 2021,
  publisher = {Cell Press},
  doi       = {10.1016/j.cell.2021.04.048},
  url       = {https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3}
}
```

**摘要：**
> The integration of multiple single-cell modalities remains a key challenge in single-cell genomics. We introduce Weighted Nearest Neighbors (WNN) analysis, a framework for unsupervised integration of multiple single-cell modalities. Our approach learns cell-specific modality weights, enabling the construction of a multimodal neighbor graph that reflects both shared and unique information across modalities. We apply WNN to simultaneous measurements of transcriptome and surface protein markers, single-cell ATAC-seq, and spatial transcriptomics. We show that WNN analysis enables robust identification of cell types and states in multimodal datasets, including from primary biological samples where only a subset of cells have each modality measured. WNN analysis is implemented in Seurat v4, enabling seamless integration and exploration of diverse multimodal datasets.

**PDF：** https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3

---

#### 20. ACTIONet — Multiresolution Framework for Single-Cell State Landscapes

**BibTeX:**

```bibtex
@article{mohammadi2020actionet,
  title     = {A multiresolution framework to characterize single-cell state landscapes},
  author    = {Mohammadi, Shahin and Davila-Velderrain, Jose and Kellis, Manolis},
  journal   = {Nature Communications},
  volume    = 11,
  number    = 1,
  pages     = {5399},
  year      = 2020,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/s41467-020-18416-6},
  url       = {https://www.nature.com/articles/s41467-020-18416-6}
}
```

**摘要：**
> ACTIONet introduces multiresolution cell-state decomposition—a concept that simultaneously captures both fine- and coarse-grain patterns of variability. The ACTIONet framework combines archetypal analysis and manifold learning to characterize single-cell states. It provides robust, reproducible, and highly interpretable analysis that couples dominant pattern discovery with structural representation of cell state landscapes. The framework demonstrates superior performance compared to existing alternatives using multiple synthetic and real datasets, and was applied to three human cortex datasets for cell integration and annotation.

---

#### 21. SC3 — 见本文档第一部分第 1 篇 (Kiselev et al., 2017)

---

### 2.2 谱聚类方法 (Spectral Clustering)

---

#### 22. SIMLR — Single-Cell Interpretation via Multikernel Learning

**BibTeX:**

```bibtex
@article{wang2017simlr,
  title     = {Visualization and analysis of single-cell RNA-seq data by kernel-based similarity learning},
  author    = {Wang, Bo and Zhu, Jiashun and Pierson, Emma and Ramazzotti, Daniele and Batzoglou, Serafim},
  journal   = {Nature Methods},
  volume    = 14,
  number    = 4,
  pages     = {414--416},
  year      = 2017,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/nmeth.4207},
  url       = {https://www.nature.com/articles/nmeth.4207}
}
```

**摘要：**
> SIMLR (Single-cell Interpretation via Multikernel Learning) is an analytic framework that learns a similarity measure from single-cell RNA-seq data to perform dimension reduction, clustering, and visualization. SIMLR learns an appropriate distance metric from single-cell RNA-seq data by combining multiple kernels, which addresses challenges specific to single-cell analysis including high noise levels, outliers, and dropout events. On seven published datasets, SIMLR demonstrated greatly enhanced clustering performance while improving visualization and interpretability of single-cell sequencing data.

---

#### 23. Spectrum — Fast Density-Aware Spectral Clustering

**BibTeX:**

```bibtex
@article{john2020spectrum,
  title     = {Spectrum: fast density-aware spectral clustering for single and multi-omic data},
  author    = {John, Christopher R and Watson, Dominic and Barnes, Michael R and Pitzalis, Costantino and Lewis, Michelle J},
  journal   = {Bioinformatics},
  volume    = 36,
  number    = 4,
  pages     = {1159--1166},
  year      = 2020,
  publisher = {Oxford University Press},
  doi       = {10.1093/bioinformatics/btz704},
  url       = {https://academic.oup.com/bioinformatics/article/36/4/1159/5566508}
}
```

**摘要：**
> Spectrum is a fast adaptive spectral clustering method designed for single and multi-omic data, including single-cell RNA-seq. The method uses a density-aware kernel that strengthens connections in the graph based on common nearest neighbors. The algorithm employs tensor product graph data integration and diffusion procedures to integrate different data sources and reduce noise. It uses either the eigengap or multimodality gap heuristics to automatically determine the number of clusters. Spectrum handles both Gaussian and non-Gaussian data structures with automatic K selection, supports single-view and multi-view clustering, and includes an ultra-fast mode for clustering datasets with more than 10,000 points.

**PDF：** https://ncbi.nlm.nih.gov/pmc/articles/PMC7703791/

---

### 2.3 层次聚类方法 (Hierarchical Clustering)

---

#### 24. CIDR — Ultrafast and Accurate Clustering Through Imputation

**BibTeX:**

```bibtex
@article{lin2017cidr,
  title     = {CIDR: Ultrafast and accurate clustering through imputation for single-cell RNA-seq data},
  author    = {Lin, Peijie and Troup, Michael and Ho, Joshua W K},
  journal   = {Genome Biology},
  volume    = 18,
  number    = 1,
  pages     = {59},
  year      = 2017,
  publisher = {BioMed Central},
  doi       = {10.1186/s13059-017-1188-0},
  url       = {https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1188-0}
}
```

**摘要：**
> CIDR (Clustering through Imputation and Dimensionality Reduction) is an ultrafast algorithm designed to address the dropout problem in single-cell RNA-seq data. The method uses a novel "implicit imputation" approach to handle dropouts in a principled manner rather than relying on heavy statistical modeling. CIDR improves upon standard PCA and outperforms state-of-the-art methods including t-SNE, ZIFA, and RaceID in terms of clustering accuracy. It is computationally efficient: typically completes within seconds for datasets of hundreds of cells and minutes for thousands of cells.

---

#### 25. SINCERA — Pipeline for Single-Cell RNA-Seq Profiling Analysis

**BibTeX:**

```bibtex
@article{guo2015sincera,
  title     = {SINCERA: a pipeline for single-cell RNA-seq profiling analysis},
  author    = {Guo, Minzhou and Wang, Hui and Potter, S Steven and Whitsett, Jeffrey A and Xu, Yui},
  journal   = {PLOS Computational Biology},
  volume    = 11,
  number    = 11,
  pages     = {e1004575},
  year      = 2015,
  publisher = {Public Library of Science},
  doi       = {10.1371/journal.pcbi.1004575},
  url       = {https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004575}
}
```

**摘要：**
> SINCERA (SINgle CEll RNA-seq profiling Analysis) is a computational pipeline for processing single-cell RNA-seq data from whole organs or sorted cells. The pipeline supports three main types of analysis: (1) distinction and identification of major cell types; (2) identification of cell type-specific gene signatures; and (3) determination of driving forces of given cell types. The authors applied SINCERA to analyze single cells from embryonic mouse lung, successfully distinguishing major fetal lung cell types—epithelial, endothelial, smooth muscle, pericyte, and fibroblast-like cells—and identifying their specific gene signatures, bioprocesses, and key regulators.

---

#### 26. densityCut — Efficient Topological Approach for Automatic Clustering

**BibTeX:**

```bibtex
@article{ding2016densitycut,
  title     = {densityCut: an efficient and versatile topological approach for automatic clustering of biological data},
  author    = {Ding, Jiarui and Shah, Sohrab and Condon, Anne},
  journal   = {Bioinformatics},
  volume    = 32,
  number    = 17,
  pages     = {2567--2576},
  year      = 2016,
  publisher = {Oxford University Press},
  doi       = {10.1093/bioinformatics/btw327},
  url       = {https://academic.oup.com/bioinformatics/article/32/17/2567/2445853}
}
```

**摘要：**
> densityCut is a novel density-based clustering algorithm that is both time- and space-efficient. The method works by: (1) estimating data point densities from a K-nearest neighbor graph; (2) refining densities via random walk; (3) identifying clusters as points in basins of attraction around density modes; (4) merging clusters through post-processing to generate a hierarchical cluster tree; (5) automatically selecting cluster numbers based on stability in the hierarchical tree. The algorithm demonstrated superior performance compared to state-of-the-art clustering methods on ten synthetic benchmarks and two microarray gene expression datasets.

---

### 2.4 稳定性度量方法 (Stability Metric)

---

#### 27. scCCESS — Autoencoder-Based Cluster Ensembles for scRNA-seq

**BibTeX:**

```bibtex
@article{geddes2019scuccess,
  title     = {Autoencoder-based cluster ensembles for single-cell RNA-seq data analysis},
  author    = {Geddes, Tom A and Kim, Tuan D and Nan, Lu and Burchfield, Joel G and Yang, Jean Y H and Tao, Di},
  journal   = {BMC Bioinformatics},
  volume    = 20,
  number    = 1,
  pages     = {660},
  year      = 2019,
  publisher = {BioMed Central},
  doi       = {10.1186/s12859-019-3179-5},
  url       = {https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3179-5}
}
```

**摘要：**
> scRNA-seq enables profiling of individual cell transcriptomes, but the high dimensionality of gene expression data can lead to uninformative clusters. scCCESS (single-cell Consensus Clusters of Encoded Subspaces) is an autoencoder-based cluster ensemble framework that: (1) takes random subspace projections from scRNA-seq data; (2) compresses each projection to low-dimensional space using an autoencoder neural network; (3) applies ensemble clustering across all encoded datasets to generate cell clusters. The framework improved clustering performance when applied with both k-means and SIMLR, with improvements reaching up to 100% depending on the evaluation metric.

---

### 2.5 聚类内/聚类间相似性方法 (Intra- and Inter-Cluster Similarity)

---

#### 28. scLCA — Latent Cellular Analysis for Large-Scale scRNA-seq Data

**BibTeX:**

```bibtex
@article{cheng2019sclca,
  title     = {Latent cellular analysis robustly reveals subtle diversity in large-scale single-cell RNA-seq data},
  author    = {Cheng, Changde and Easton, Jon and Rosencrance, Christine and Li, Yan and Ju, Bifan and Williams, Jennifer and Urrutia, Nancy D and Wei, Jiayue and Wang, Peiyao and Chen, Bo},
  journal   = {Nucleic Acids Research},
  volume    = 47,
  number    = 18,
  pages     = {e143},
  year      = 2019,
  publisher = {Oxford University Press},
  doi       = {10.1093/nar/gkz860},
  url       = {https://academic.oup.com/nar/article/47/18/e143/5541864}
}
```

**摘要：**
> scLCA (single-cell Latent Cellular Analysis) is a machine learning-based computational pipeline that combines cosine-similarity measurement by latent cellular states, graph-based clustering, and dual-space model search to determine optimal subpopulation numbers and informative cellular states. LCA addresses major challenges in scRNA-seq analysis including robustness, accuracy, scalability (enabling analysis of hundreds of thousands of cells), and automatic population inference without requiring explicit gene filtering.

---

#### 29. RaceID — Single-Cell mRNA Sequencing Reveals Rare Intestinal Cell Types

**BibTeX:**

```bibtex
@article{grun2015raceid,
  title     = {Single-cell messenger RNA sequencing reveals rare intestinal cell types},
  author    = {Grün, Dominic and Lyubimova, Anna and Kester, Lennart and Wiebrands, Kay and Basak, Onur and Sasaki, Nobuo and Clevers, Hans and van Oudenaarden, Alexander},
  journal   = {Nature},
  volume    = 525,
  number    = 7568,
  pages     = {251--255},
  year      = 2015,
  publisher = {Nature Publishing Group},
  doi       = {10.1038/nature14966},
  url       = {https://www.nature.com/articles/nature14966}
}
```

**摘要：**
> The study addresses the challenge of identifying rare cell types in complex tissues. The researchers sequenced transcriptomes of hundreds of randomly selected cells from mouse intestinal organoids and developed RaceID, an algorithm for rare cell type identification in single-cell populations. The algorithm can resolve cell types represented by as few as a single cell in a population. Using RaceID, they identified Reg4 as a novel marker for enteroendocrine cells, a rare population of hormone-producing intestinal cells. The authors propose broad applicability of RaceID for discovering rare cell types and corresponding marker genes in healthy and diseased organs.

---

### 2.6 随机投影与元聚类方法 (Meta-clustering / Random Projection)

---

#### 30. SHARP — Hyperfast and Accurate Processing via Ensemble Random Projection

**BibTeX:**

```bibtex
@article{wan2020sharp,
  title     = {SHARP: hyperfast and accurate processing of single-cell RNA-seq data via ensemble random projection},
  author    = {Wan, Shibiao and Kim, Junil and Won, Kyoung-Jae},
  journal   = {Genome Research},
  volume    = 30,
  number    = 2,
  pages     = {205--213},
  year      = 2020,
  publisher = {Cold Spring Harbor Laboratory Press},
  doi       = {10.1101/gr.254557.119},
  url       = {https://genome.cshlp.org/content/30/2/205}
}
```

**摘要：**
> SHARP (Single-cell RNA-seq Hyper-fast and Accurate processing via ensemble Random Projection) is an ensemble random projection-based algorithm designed to process large-scale scRNA-seq data effectively without excessive distortion during dimension reduction. The method is scalable to clustering 10 million cells and demonstrates superior speed and accuracy compared to existing methods on 17 public scRNA-seq datasets. For large datasets (>40,000 cells), SHARP runs faster than competitors while maintaining high clustering accuracy and robustness. SHARP uses ensemble random projection combined with multi-layer meta-clustering to preserve cell-to-cell distances in reduced-dimensional space.

---

## 三、论文汇总索引表

| # | 方法名 | 年份 | 作者 | 期刊/会议 | 摘要关键词 |
|---|--------|------|------|-----------|------------|
| 1 | SC3 | 2017 | Kiselev et al. | Nature Methods | consensus clustering, consensus matrix |
| 2 | Louvain | 2008 | Blondel et al. | arXiv | modularity optimization, community detection |
| 3 | Leiden | 2019 | Traag et al. | Scientific Reports | well-connected communities, refinement |
| 4 | DEC | 2016 | Xie et al. | ICML | deep embedding, simultaneous clustering |
| 5 | scDeepCluster | 2019 | Tian et al. | Nature Machine Intelligence | model-based, ZINB, deep embedded |
| 6 | DESC | 2020 | Li et al. | Nature Communications | batch effect removal, soft clustering |
| 7 | scziDesk | 2020 | Chen et al. | NAR Genomics & Bioinformatics | denoising AE, self-training K-means |
| 8 | scDCC | 2021 | Tian et al. | Nature Communications | constrained clustering, semi-supervised |
| 9 | scNAME | 2022 | Wan et al. | Bioinformatics | mask estimation, contrastive learning |
| 10 | scMAE | 2024 | Fang et al. | Bioinformatics | masked autoencoder, gene correlation |
| 11 | scGNN | 2021 | Wang et al. | Nature Communications | GNN, left-truncated Gaussian, imputation |
| 12 | scDSC | 2022 | Gan et al. | Briefings in Bioinformatics | ZINB-AE + GNN, mutual-supervised |
| 13 | AttentionAE-sc | 2023 | Li et al. | PLOS Computational Biology | attention, ZINB + GAE fusion |
| 14 | scCDCG | 2024 | Xu et al. | DASFAA | cut-informed graph, optimal transport |
| 15 | scGPT | 2024 | Cui et al. | Nature Methods | foundation model, transformer, generative AI |
| 16 | GeneFormer | 2023 | Theodoris et al. | Nature | context-aware BERT, transfer learning |
| 17 | GeneCompass | 2024 | Yang et al. | Cell Research | knowledge-informed, cross-species, GRN |
| 18 | Monocle3 | 2019 | Cao et al. | Nature | trajectory inference, Leiden clustering |
| 19 | Seurat v4 | 2021 | Hao et al. | Cell | multimodal integration, WNN |
| 20 | ACTIONet | 2020 | Mohammadi et al. | Nature Communications | multiresolution, archetypal analysis |
| 21 | SIMLR | 2017 | Wang et al. | Nature Methods | multikernel learning, kernel-based |
| 22 | Spectrum | 2020 | John et al. | Bioinformatics | density-aware, spectral clustering |
| 23 | CIDR | 2017 | Lin et al. | Genome Biology | implicit imputation, dropout |
| 24 | SINCERA | 2015 | Guo et al. | PLOS Computational Biology | pipeline, hierarchical clustering |
| 25 | densityCut | 2016 | Ding et al. | Bioinformatics | topological, hierarchical, stability |
| 26 | scCCESS | 2019 | Geddes et al. | BMC Bioinformatics | autoencoder ensemble, subspace |
| 27 | scLCA | 2019 | Cheng et al. | NAR | latent cellular states, graph clustering |
| 28 | RaceID | 2015 | Grün et al. | Nature | rare cell types, K-means |
| 29 | SHARP | 2020 | Wan et al. | Genome Research | ensemble random projection, meta-clustering |

---

*本文档由 AI 辅助文献检索生成，BibTeX 格式引用信息经手工整理自各论文官方来源。建议在使用前通过 Semantic Scholar (semanticscholar.org) 或 Google Scholar 核验引用信息的准确性。*

# LLM4Cell: A Survey of Large Language and Agentic Models for Single-Cell Biology

Sajib Acharjee Dip*1, Adrika Zafor*2, Bikash Kumar Paul1 

Uddip Acharjee Shuvo3, Muhit Islam Emon1, Xuan Wang1, Liqing Zhang†1,4 

1Department of Computer Science, Virginia Tech, Blacksburg, VA, USA 

2Department of Computational Modeling and Data Analytics, Virginia Tech, Blacksburg, VA, USA 

3Institute of Information and Technology, University of Dhaka, Dhaka, Bangladesh 

4Fralin Biomedical Research Institute at VTC: Cancer Research Center, Washington DC, USA 

# Abstract

Large language models (LLMs) and emerging agentic frameworks are beginning to transform single-cell biology by enabling natural-language reasoning, generative annotation, and multimodal data integration. However, progress remains fragmented across data modalities, architectures, and evaluation standards. LLM4Cell presents the first unified survey of 58 foundation and agentic models developed for single-cell research, spanning RNA, ATAC, multi-omic, and spatial modalities. We categorize these methods into five families—foundation, text-bridge, spatial/multimodal, epigenomic, and agentic—and map them to eight key analytical tasks including annotation, trajectory and perturbation modeling, and drug-response prediction. Drawing on over 40 public datasets, we analyze benchmark suitability, data diversity, and ethical or scalability constraints, and evaluate models across 10 domain dimensions covering biological grounding, multi-omics alignment, fairness, privacy, and explainability. By linking datasets, models, and evaluation domains, LLM4Cell provides the first integrated view of languagedriven single-cell intelligence and outlines open challenges in interpretability, standardization, and trustworthy model development. 

# 1 Introduction

Large language models (LLMs) are transforming biomedical discovery by linking molecular patterns with knowledge encoded in text. (Zhang et al., 2025a; Aljohani et al., 2025; Bi et al., 2024). Prior work on computational workflows and biological prediction tasks demonstrates how learned representations capture complex molecular relationships 

and enhance practical performance. (Lam et al., 2024; Bhattacharya et al., 2024; Dip et al., 2024; Nam et al., 2024; Shuvo et al., 2024) In single-cell biology, where each experiment measures thousands of genes across millions of cells, LLMs promise to unify gene expression, chromatin accessibility, and spatial organization under a shared, interpretable framework. Early foundation models such as scGPT (Cui et al., 2024), Geneformer (Theodoris et al., 2023), and scFoundation (Hao et al., 2024) learn transferable representations directly from single-cell data, while emerging systems like CellLM (Zhao et al., 2023), scAgent (Mao et al., 2025), and CellVerse (Zhang et al., 2025b) extend these capabilities toward reasoning, dialogue, and autonomous analysis. Together, they signal a shift from statistical inference to languagedriven, interpretable, and generalizable cell intelligence. 

Despite this progress, the field remains fragmented. Models differ widely in data modality, supervision type, and evaluation standard ranging from scRNA-seq to ATAC, multi-omic, and spatial measurements hindering cross-model comparison and reproducibility. Benchmarking is inconsistent, datasets are unevenly distributed across modalities, and agentic frameworks lack standardized ways to measure reasoning correctness or biological grounding. Without a systematic synthesis, it is difficult to assess what current single-cell LLMs truly understand or where they fall short. 

To address these challenges, we introduce LLM4Cell, a comprehensive survey of large language and agentic models for single-cell biology. We analyze 58 representative methods spanning five methodological familiesFoundation, Text-Bridge, Spatial/Multimodal, Epigenomic, and Agentic and organize them across eight major tasks, from annotation and trajectory inference to per-

turbation and drug-response modeling. Complementing this, we compile over 40 publicly available datasets across RNA, ATAC, multi-omic, spatial, perturbation, and plant domains, and evaluate each model using a ten-dimension rubric covering grounding, fairness, scalability, and interpretability. 

Contributions. LLM4Cell provides: (1) a modality-balanced registry of around 40 benchmark datasets for large-model training and evaluation; (2) a unified taxonomy of 58 foundation and agentic models across eight analytical tasks; (3) a ten-dimension domain rubric capturing biological grounding, generalization, and ethics; and (4) a critical discussion of open challenges in crossmodal alignment, reasoning, and trustworthy AI for cell biology. 

# 2 Related Work

The intersection of large language models (LLMs) and single-cell biology has gained rapid momentum, motivating several early surveys and benchmarks. The ACL 2025 survey on foundation models for single-cell biology (Zhang et al., 2025a) outlines pretrained and fine-tuned architectures across protein, gene, and cell representations. Broader reviews such as Lan et al. (2025) and Bian et al. (2024b) discuss “large cellular models” (LCMs) and trends in scaling, tokenization, and transfer learning across biological data. Benchmark studies including Single-Cell Omics Arena (Liu et al., 2024) and CellVerse (Zhang et al., 2025b) evaluate LLMs on annotation and question answering, while works like scInterpreter (Guo et al., 2023) and scExtract (Wu and Tang, 2025) demonstrate practical integration into annotation and curation pipelines. 

However, existing reviews remain limited in scope. They focus on architectures or prompting but give little attention to the datasets underpinning model training and evaluation especially those covering spatial, perturbation, or multimodal modalities. Standardized benchmarks, domain-aware metrics, and reproducible splits are seldom addressed, and plant single-cell datasets key for cross-kingdom generalization are almost entirely absent. Moreover, agentic reasoning and tool-augmented systems are typically treated as isolated examples rather than a coherent methodological class, leaving the relationship among data modalities, modeling paradigms, and biological domains fragmented. 

LLM4Cell advances this literature by providing 

an integrated, data-centric synthesis linking models, datasets, and evaluation domains. Unlike prior descriptive taxonomies, it unifies over forty public datasets and fifty-eight models across modalities and tasks, assessing their biological grounding, fairness, and scalability. In doing so, LLM4Cell establishes a reproducible foundation for cross-modal benchmarking and future research on trustworthy, autonomous single-cell intelligence. 

# 3 Datasets

Progress in large language models for single-cell biology relies on the rapid growth of high-quality datasets across transcriptomic, epigenomic, multimodal, and spatial domains. Our survey compiles over forty public resources spanning five major modalities RNA, ATAC, multi-omic, spatial, and perturbation/drug-response plus emerging plant single-cell atlases. These datasets underpin model pretraining, evaluation, and cross-modal reasoning for annotation, integration, trajectory inference, perturbation modeling, and spatial mapping. 

Transcriptomic atlases such as Tabula Sapiens, Tabula Muris (Consortium* et al., 2022), and the Human Cell Atlas (Travaglini et al., 2020) remain dominant for foundation-model training and ontology-aware annotation. Chromatinaccessibility data (e.g. Cusanovich mouse atlas (Cusanovich et al., 2018), human adult/fetal scATAC (Domcke et al., 2020)) map regulatory states but remain sparse and heterogeneous. Multiomic resources (e.g. TEA-seq (Swanson et al., 2021), DOGMA-seq (Mimitou et al., 2021), CITEseq (Stoeckius et al., 2017)) link RNA, ATAC, and protein modalities, providing supervision for crossview alignment. Spatial technologies (e.g. Visium (Oliveira et al., 2025), Slide-seqV2 (Stickels et al., 2021), MERFISH (Chen et al., 2015), Stereo-seq (Chen et al., 2022)) connect molecular profiles to tissue architecture, while functional datasets (e.g. Perturb-seq (Replogle et al., 2022)) benchmark causal reasoning. Plant resources (e.g. scPlantDB (He et al., 2024), Arabidopsis E-CURD-4 (Shulse et al., 2019)) extend modeling beyond animal systems, opening cross-kingdom evaluation. 

Despite this diversity, major gaps persist. (i) RNA atlases vastly outscale other modalities. (ii) Benchmark fragmentation and inconsistent metadata hinder reproducibility. (iii) Privacy constraints limit access to clinical spatial datasets. (iv) Non-human and plant data remain scarce. (v) Paired and trimodal resources provide ideal testbeds for next-

generation multimodal and agentic LLMs. Table 1, 2, 3, 4, 5, 6 (Appendix D) list representative datasets with tasks, scale, and source links for reproducible benchmarking. 

# 4 Model Taxonomy

Large language models are rapidly reshaping single-cell biology, producing diverse architectures, pretraining schemes, and reasoning frameworks. We categorize 58 representative methods into five methodological families Foundation Models, Text-Bridge LLMs, Spatial and Multimodal Models, Epigenomic Models, and Agentic Frameworks based on their core design and data modality. These span the progression from genelevel embeddings to multimodal and autonomous systems capable of biological reasoning. 

Our taxonomy is organized along five orthogonal dimensions: (i) Modality (RNA, ATAC, multiomic, spatial, or hybrid), (ii) Grounding type (atlas, ontology, or marker-based), (iii) Agentic capability (multi-step reasoning or autonomous orchestration), (iv) Primary task (annotation, trajectory, perturbation, integration, etc.), and (v) Domain quality, represented by ten quantitative scores for grounding, fairness, scalability, and interpretability. 

Figure 1 summarizes this taxonomy, with detailed attributes in Appendix Table 7. Foundation models currently dominate in scope and adoption, while emerging spatial, epigenomic, and agentic frameworks signal a shift toward contextual, reasoningdriven single-cell intelligence. 

# 4.1 Foundation Models

Foundation models learn transferable cell and gene embeddings directly from large-scale scRNA-seq without explicit labels. Representative systems such as scGPT (Cui et al., 2024), Geneformer (Theodoris et al., 2023), and scFoundation (Hao et al., 2024) pretrain on multi-tissue atlases exceeding one million cells, using masked-gene or rankbased reconstruction to capture expression context. Variants like tGPT (Shen et al., 2023), scBERT (Yang et al., 2022), and scGraphformer (Fan et al., 2024) treat genes as tokens or graph nodes, while cross-species models (UCE (Rosen et al., 2023), GeneCompass (Yang et al., 2024b)) apply contrastive alignment across homologous genes. These models underpin most single-cell LLM pipelines, offering strong transfer for annotation and integration but limited ontology grounding and explainability, motivating later text-bridged and reasoning 

frameworks. 

# 4.2 Text-Bridge LLMs

Text-Bridge models couple molecular embeddings with biomedical language to ground single-cell representations in semantics and ontology. Systems such as scELMo (Liu et al., 2023), CellLM (Zhao et al., 2023), and GenePT (Chen and Zou, 2024) align gene or cell embeddings with textual descriptors, while Cell2Sentence (Levine et al., 2024) and Cell2Text (Kharouiche et al., 2025) translate expression profiles into natural-language summaries. Using dual-encoder or encoder–decoder architectures, these frameworks employ contrastive or prompt-based alignment between molecular and text encoders (e.g., BioBERT). They enhance interpretability and enable zero-shot annotation but depend heavily on curated corpora and remain nonagentic, serving as a bridge between foundation and reasoning frameworks. 

# 4.3 Spatial and Multimodal Models

Spatial and multimodal frameworks integrate gene expression with spatial coordinates, histology, or additional omics to capture tissue architecture. Representative models include TransformerST (Lu et al., 2024), spaLLM (Ji et al., 2024), and OmiCLIP (Cui et al., 2025) for spatial mapping, and scMMGPT (Shi et al., 2025) and FmH2ST for RNA–ATAC–protein integration. They use multi-branch Transformers with modality-specific encoders and cross-attention fusion, trained on datasets such as Visium DLPFC and MERFISH. These models achieve strong spatial alignment and biological realism but face heterogeneous resolutions, limited open benchmarks, and high computational cost. 

# 4.4 Epigenomic Models

Epigenomic foundation models extend LLM concepts to chromatin-accessibility and regulatory data such as scATAC-seq. EpiFoundation (Wu et al., 2025), EpiBERT (Javed et al., 2025), and EpiAttend (Li et al., 2022) learn cis-regulatory patterns from ENCODE and other compendia, while GeneMamba and scMamba use efficient state-space layers for long-range dependency modeling. Crossomic variants (ChromFound (Jiao et al., 2025), GET (Fu et al., 2025)) jointly embed RNA and ATAC to infer gene-regulatory networks. These models improve biological grounding but remain constrained by sparse data and lack unified benchmarks across regulatory modalities. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-28/1264e285-66b1-43c6-8a68-5faf4ab9f487/6d649db21abc0e374a383adb29ce3c50baaf2a474367c5f962ae3c7eb8fd7862.jpg)



Figure 1: Hierarchical taxonomy for LLM4Cell. Color-coded families expand into sub-branches and representative models, tracing the progression from foundation pretraining to multimodal and agentic reasoning frameworks. References are omitted for visibility and included in the appendix method comparison table.


# 4.5 Agentic Frameworks

Agentic systems integrate pretrained models with reasoning modules for autonomous single-cell analysis. Frameworks such as scAgent (Mao et al., 2025), CellVerse (Zhang et al., 2025b), and Teddy (Chevalier et al., 2025) combine domainspecific encoders with LLM controllers that plan tasks, query ontologies, and interface with tools or APIs. EpiAgent (Chen et al., 2025b) extends this paradigm to regulatory genomics. These systems enable dialogue-based annotation and multi-step reasoning but lack standardized benchmarks for reasoning fidelity and rely on underlying LLM reliability. 

# 5 Task-specific Applications

While the previous section categorized models by architectural family, here we analyze them through the lens of biological objectives shown in Figure 2 and Appendix Figure 3. Across the 58 methods surveyed, we identify eight recurring tasks that define single-cell modeling pipelines: (1) annotation and ontology mapping, (2) trajectory and perturbation modeling, (3) multi-omic integration, (4) spatial mapping and deconvolution, (5) regulatorynetwork and pathway inference, (6) cross-species translation, (7) generative simulation, and (8) drugresponse prediction. Most models span multiple tasks-for example, scGPT and Geneformer support both annotation and perturbation predictionreflecting the convergence between representation learning and functional reasoning. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-28/1264e285-66b1-43c6-8a68-5faf4ab9f487/2f04f7c63bb36880a0f6920c2a239cfe2ed240d9e604cfcb88686de122e1e318.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-28/1264e285-66b1-43c6-8a68-5faf4ab9f487/ff0769cc27c4b994c33d0f39eb26a32c06fc5e427be92b0e9727b5780a7ae463.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-28/1264e285-66b1-43c6-8a68-5faf4ab9f487/3214a0d12e5af26dfbf89996a1aec5faa71adf3d60b449c87c3cf1084bad18fa.jpg)



Figure 2: Task and domain distribution across spatial transcriptomics models. Most emphasize biological grounding, scalability, and annotation tasks.


This task-centric view complements the architectural taxonomy by revealing how foundation and agentic systems differ in operational scope. Foundation models dominate annotation, integration, and trajectory inference, whereas emerging textbridge and agentic systems extend to knowledgegrounded reasoning and dynamic planning. Spatial and epigenomic frameworks contribute specialized capabilities in tissue mapping and chromatin-level interpretation. 

# 5.1 Annotation and Ontology Mapping

Annotation is the most common single-cell task, spanning automated cell-type labeling, ontology alignment, and cross-dataset harmonization. Foundation models such as scGPT, scFoundation, and scBERT achieve high accuracy on atlases like Tabula Sapiens and the Human Cell Atlas. Text-bridge models (CellLM (Zhao et al., 2023), scELMo (Liu et al., 2023), GenePT (Chen and Zou, 2024)) add ontology-based interpretability, and agentic systems (scAgent) perform multi-step annotation reasoning via LLM controllers. Models are finetuned on curated references (e.g., Azimuth PBMC, HCA Lung) using masked-gene or contrastive objectives. Evaluation uses label accuracy or ARI; key limitations include rare-cell detection, crossspecies consistency, and lack of reasoning benchmarks. 

# 5.2 Trajectory and Perturbation Modeling

Trajectory modeling captures dynamic state transitions and causal responses to interventions. Geneformer (Theodoris et al., 2023), scGPT (Cui et al., 2024), and scRDiT (Dong et al., 2025) model temporal or perturbation trajectories using denoising and contrastive objectives, trained on datasets such as Replogle 2022 Perturb-seq and sci-Plex. Epigenomic models (EpiFoundation (Wu et al., 2025) ,EpiAgent (Chen et al., 2025b)) extend to regu-

latory responses through chromatin context and ontology reasoning. Metrics include expressionprofile correlation and target-gene recovery. While these models predict single perturbations effectively, performance degrades for combinatorial or long-range effects, underscoring the need for unified causal benchmarks. 

# 5.3 Multi-omic Integration

Multi-omic integration seeks unified representations across paired RNA, ATAC, and protein modalities to capture gene-regulatory and signaling relationships. Representative models include scM-MGPT (Shi et al., 2025), GET (General Expression Transformer) (Fu et al., 2025), Chrom-Found (Jiao et al., 2025), and epigenomic-aware variants such as EpiFoundation (Wu et al., 2025) and EpiBERT (Javed et al., 2025). These systems train on joint datasets like 10x Multiome PBMC 10k, TEA-seq, DOGMA-seq, and ASAP-seq, using cross-modal transformers or contrastive alignment between expression and accessibility features. 

Most architectures employ modality-specific encoders with shared latent spaces or token-type embeddings; scMGPT applies masked-reconstruction across modalities, while GET learns cross-attention between gene and chromatin tokens. Compared with unimodal models, they show improved celltype alignment and batch correction but face datascale imbalance and sparse modality overlap. Evaluation commonly uses modality-matching accuracy and latent-space correlation. Despite strong integration performance, interpretability and standard benchmarks remain limited, and scalability to trimodal data is an active challenge. 

# 5.4 Spatial Mapping and Deconvolution

Spatial mapping links molecular profiles to their tissue locations, enabling inference of cell organization and microenvironmental context. Represen-

tative models include TransformerST (Lu et al., 2024), spaLLM , OmiCLIP, FmH2ST, HEIST (Madhu et al., 2025), and Spatial2Sentence (Chen et al., 2025a), trained on datasets such as Visium DLPFC, Slide-seqV2, MERFISH, and Xenium. They combine spatial coordinates, histology features, and gene-expression embeddings through cross-attention or contrastive alignment to reconstruct cell or spot-level maps. 

TransformerST performs axial attention between gene and spatial axes for spot-to-cell deconvolution, while spaLLM and OmiCLIP align histology and transcriptomic features using language-guided contrastive learning. Models like HEIST and FmH2ST fuse H&E images with RNA profiles for tissuedomain segmentation. Evaluation typically measures spatial correlation, clustering accuracy, or cell-type F1 against curated annotations. Although these models improve contextual understanding of tissue architecture, limited spatial benchmarks and heterogeneous resolutions constrain systematic comparison across methods. 

# 5.5 Regulatory Network and Pathway Inference

Regulatory-network and pathway inference focuses on identifying gene–gene or enhancer–promoter dependencies that drive transcriptional programs. Representative models include GeneMamba (Qi et al., 2025a), scMamba (Yuan et al., 2025), Epi-Foundation (Wu et al., 2025), EpiBERT (Javed et al., 2025), ChromFound (Jiao et al., 2025), and GET (Fu et al., 2025), trained on large-scale ATAC and multiome atlases such as ENCODE, human cCRE datasets, and 10x Multiome PBMC. These models treat chromatin regions or genes as tokens and learn context-dependent regulatory relationships through self- or cross-attention. 

GeneMamba (Qi et al., 2025a) and scMamba employ state-space layers for efficient modeling of long-range chromatin dependencies, while EpiB-ERT and EpiFoundation use masked-region reconstruction to capture enhancer–promoter coupling. ChromFound and GET link RNA and accessibility signals to infer transcription-factor activity and causal directionality. Evaluation typically uses motif enrichment, AUROC for known interactions, or pathway-level overlap with KEGG or Reactome references. These models improve biological grounding and causal interpretability but remain constrained by sparse training data and lack of unified ground-truth regulatory benchmarks. 

# 5.6 Cross-Species Translation

Cross-species translation aims to transfer cell and gene representations between organisms, enabling comparative and evolutionary analysis of single-cell atlases. Representative models include iSEEEK (Shen et al., 2022), UCE (Universal Cell Embeddings) (Rosen et al., 2023), GeneCompass (Yang et al., 2024b), and scPlantLLM (Cao et al., 2025), trained on paired or homolog-mapped datasets such as tthe Mouse Cell Atlas, Tabula Muris, and scPlantDB. These systems learn alignment between orthologous genes or conserved cell states to support zero-shot annotation across species. 

iSEEEK (Shen et al., 2022) integrates over 11 million human and mouse cells using shared genetoken vocabularies, while UCE (Shen et al., 2023) applies masked autoencoding and contrastive alignment across ortholog groups. GeneCompass (Yang et al., 2024b) extends this with knowledge-informed pretraining that links species via ontology-derived embeddings, and scPlantLLM adapts the paradigm to plant single-cell data across 17 species. Evaluation uses transfer accuracy and embedding consistency across species pairs. Despite improved cross-domain generalization, these models remain limited by incomplete ortholog mappings and domain shifts in sequencing depth and tissue composition. 

# 5.7 Generation and Simulation

Generative modeling aims to synthesize realistic single-cell profiles, simulate perturbations, or reconstruct missing modalities. Representative approaches include scGPT, scFoundation, CellFM (Zeng et al., 2025), and Geneformer, which employ decoder or diffusion-style architectures to model the probability distribution of geneexpression states. These models are pretrained on large atlases such as Tabula Sapiens and HCA, then fine-tuned to generate new cell embeddings or perturbation outcomes. 

scGPT and scFoundation use masked-token prediction extended with generative decoding to simulate unseen cells or conditions, while CellFM introduces multimodal conditioning for cross-tissue synthesis. Generative evaluation typically measures distributional similarity (KL divergence, FID-like embedding distance) or recovery of known celltype proportions. Such models enable scalable in-silico experiments and data augmentation but face challenges in controlling biological realism 

and avoiding overfitting to dominant cell types. 

# 5.8 Drug-Response Prediction

Drug-response prediction models aim to infer cellular transcriptional outcomes following chemical or genetic perturbations, supporting virtual screening and mechanism discovery. Representative frameworks include Geneformer (Theodoris et al., 2023), scGPT (Cui et al., 2024), EpiFoundation (Wu et al., 2025), and EpiAgent (Chen et al., 2025b), trained or evaluated on large-scale pharmacotranscriptomic resources such as sci-Plex, Replogle 2022 Perturb-seq, and the ARC Virtual Cell Challenge. These systems learn conditional embeddings that map control to perturbed states, enabling zero-shot prediction of unseen compounds or targets. 

Geneformer predicts post-treatment expression profiles using masked-token denoising, while scGPT and EpiFoundation integrate modality-specific conditioning to capture epigenetic or pathway context. EpiAgent further incorporates ontology-guided reasoning to link molecular features with literaturederived drug–gene relationships. Evaluation relies on Pearson correlation or top- $k$ recovery of known targets between predicted and observed responses. Although these models demonstrate strong predictive capacity for single-agent perturbations, performance degrades for combinatorial treatments and unseen pathways, highlighting the need for richer multimodal pharmacogenomic benchmarks. 

Cross-task observations. Across all eight task domains, clear convergence patterns emerge. Foundation models dominate annotation, integration, and generative simulation owing to large-scale pretraining on transcriptomic atlases, whereas textbridge systems extend interpretability through ontology grounding and natural-language supervision. Spatial and epigenomic models expand biological realism by capturing tissue architecture and regulatory logic, while agentic frameworks introduce multi-step reasoning and tool orchestration that connect these modalities into autonomous workflows. Despite methodological diversity, evaluation remains fragmented-annotation and integration enjoy standardized benchmarks, but trajectory, drugresponse, and agentic reasoning lack common protocols. Cross-species and privacy-aware tasks are the least represented, reflecting data scarcity and regulatory constraints. Overall, task analysis reveals a rapidly maturing ecosystem in which foundation and multimodal pretraining supply general 

representations, and emerging agentic paradigms begin to operationalize biological reasoning. The following section quantifies these trends through ten domain-level dimensions encompassing biological grounding, fairness, scalability, and interpretability. 

# 6 Domain Analysis

We evaluate fifty-eight methods across ten domain dimensions capturing reliability, generalization, and interpretability. Each score reflects published evidence of performance or explicit design alignment, summarized in Appendix Table 7 and Figure 2. Together, these dimensions provide a composite view of model maturity beyond task accuracy. 

# 6.1 Biological Grounding

Most foundation and multimodal models achieve strong biological grounding by encoding marker genes, pathways, or regulatory regions. scGPT, Geneformer, and EpiFoundation leverage atlasscale training to recover known cell types and pathway enrichments. Spatial frameworks such as TransformerST and OmiCLIP add tissue-context grounding via histology alignment. However, textual and ontology grounding remain inconsistent, and biological semantics are rarely explicit in model objectives. 

# 6.2 Batch Effects and Heterogeneity

Robustness to batch effects varies widely. scFoundation (Hao et al., 2024), CellFM (Zeng et al., 2025), and iSEEEK (Shen et al., 2022) mitigate donor and protocol variation through largescale multi-domain pretraining, while others rely on explicit batch tokens or adversarial alignment. True cross-platform generalization remains limited; even large models show degradation when tested on unseen sequencing chemistries or tissue types. 

# 6.3 Multi-omics Alignment

Models such as scMGPT (Palayew et al., 2025), GET (Fu et al., 2025), and ChromFound (Jiao et al., 2025) demonstrate consistent alignment across RNA, ATAC, and protein features using cross-attention or contrastive fusion. Spatial hybrids (HEIST, FmH2ST) further connect molecular and visual domains. While alignment quality improves, limited paired datasets and uneven modality coverage constrain reproducibility and benchmarking. 

# 6.4 Trajectory and Perturbation

Dynamic inference is strongest in generative transformers (Geneformer, scGPT) and epigenomic variants (EpiAgent). They predict postperturbation expression or accessibility shifts with reasonable fidelity but struggle to model multitarget or time-continuous responses. Few methods quantify causal uncertainty or biological plausibility. 

# 6.5 Cross-Species and Cross-Tissue Generalization

Cross-organism transfer remains challenging. iSEEEK (Shen et al., 2022), UCE (Rosen et al., 2023), and GeneCompass (Yang et al., 2024b) align human and mouse cell embeddings via ortholog mapping, while scPlantLLM (Cao et al., 2025) extends this to plants. Performance drops sharply for divergent taxa, reflecting incomplete ortholog tables and bias toward human datasets. 

# 6.6 Atlas Fairness and Representation Balance

Dataset imbalance affects nearly all models. Human and immune-cell–dominant atlases overrepresent certain tissues and demographics. Teddy (Chevalier et al., 2025) introduces benchmarking for fairness and representation balance, but fairness metrics are rarely adopted. No model yet enforces demographic or tissue-specific parity during training. 

# 6.7 Explainability and Interpretability

Text-bridge systems (GenePT, CellLM, Cell2Text (Kharouiche et al., 2025)) improve interpretability by linking embeddings to ontology terms. Agentic systems (scAgent) extend this via natural-language reasoning and tool explanations. Most foundation models still rely on attention visualization, providing limited causal transparency. 

# 6.8 Privacy and Ethics

Data privacy remains underexplored. Synthetic generation (scGFT (Nouri, 2025)) and openaccess chat systems (ChatCell (Fang et al., 2024)) raise concerns over re-identification and data leakage. No single-cell LLM currently implements federated or privacy-preserving learning, and ethical guidance for multi-omic sharing remains informal. 

# 6.9 Scalability and Efficiency

Scaling trends mirror NLP: efficient transformers (xTrimoGene (Gong et al., 2023)) and state-

space models (GeneMamba (Qi et al., 2025a), sc-Mamba (Yuan et al., 2025)) reduce memory cost while retaining accuracy. CellFM (Zeng et al., 2025) demonstrates training across $1 0 0 ~ \mathrm { M }$ cells, showing feasibility of atlas-scale learning. Still, high compute requirements limit accessibility for most research groups. 

# 6.10 Emerging Paradigms and Agentic Behavior

Recent frameworks (scAgent (Mao et al., 2025), CellVerse (Zhang et al., 2025b), EpiAgent (Chen et al., 2025b)) introduce reasoning and autonomous decision pipelines, integrating LLM controllers with specialized encoders. These systems achieve the highest scores in explainability and cross-modal planning but lack standardized evaluation of reasoning fidelity or reproducibility. 

# 7 Open Problems

# 7.1 Trust and Validation

Evaluation across modalities remains inconsistent. Metrics often emphasize reconstruction over biological plausibility, and few studies offer independent replication. Community-curated leaderboards and standardized validation sets are essential for credible comparison. 

# 7.2 Data and Bias

Training corpora are dominated by human and mouse atlases, limiting cross-species and clinical generalization. Rare-cell, plant, and microbial systems remain underrepresented, reinforcing demographic and biological bias. Diverse, balanced datasets are critical for equitable model development. 

# 7.3 Cross-Modal and Dynamic Modeling

True integration of RNA, ATAC, spatial, and temporal modalities is still elusive. Existing models handle only pairwise combinations, while four-way fusion remains computationally infeasible. Unified tokenization and multimodal pretraining pipelines could enable consistent cross-domain reasoning. 

# 7.4 Interpretability and Causality

LLM embeddings capture statistical correlation but rarely mechanistic insight. Bridging attention to regulatory or causal networks requires explicit reasoning layers and experimental grounding. Combining symbolic or causal discovery modules with transformers is a key next step. 

# 7.5 Ethics and Privacy

Open single-cell data sharing raises privacy and dual-use concerns. Few models include audit trails or consent-aware governance. Adapting federated and differential-privacy methods to biological data is urgently needed. 

# 7.6 Agentic and Interactive Systems

Agentic frameworks can plan and reason but lack reliable benchmarks. Tasks such as multi-step reasoning, dataset retrieval, and model orchestration remain fragile. Robust metrics for reasoning accuracy, safety, and reproducibility will be crucial for trustworthy biological AI. 

# 8 Conclusions

LLM4Cell presents the first unified survey of large language and agentic models for single-cell biology, linking datasets, architectures, and evaluation domains within a common framework. By analyzing 58 models and over 40 public datasets across RNA, ATAC, multi-omic, spatial, and perturbation modalities, we show how foundation-scale pretraining, multimodal grounding, and agentic reasoning are reshaping single-cell analysis. Our taxonomy and evaluation rubric reveal a transition from purely statistical modeling toward language-driven, interpretable, and increasingly autonomous systems that connect molecular data with biological knowledge. We hope this synthesis provides a reproducible reference for benchmarking, model selection, and the design of next-generation cellular foundation and reasoning models. 

# Limitations

Despite its breadth, this study has several limitations. First, reported performance metrics vary widely across publications, preventing consistent quantitative comparison. Second, our tendimension rubric captures qualitative trends but not standardized numerical rankings. Third, access restrictions and licensing prevented inclusion of certain clinical or proprietary spatial datasets, while non-animal single-cell resources remain scarce. Fourth, model efficiency, compute cost, and hyperparameter sensitivity were not systematically analyzed. Finally, as multimodal and reasoning frameworks evolve rapidly, LLM4Cell reflects a snapshot in time rather than a static benchmark. Future work should establish community leaderboards, causal interpretability tests, and reasoning benchmarks to advance trustworthy and reproducible single-cell intelligence. 

# Acknowledgments

This work was supported in part by Virginia Tech, the Department of Computer Science, and the U.S. National Science Foundation (NSF) under Awards #2125798, #2344169, and #2319522. We gratefully acknowledge their support and resources that made this research possible. 

# References

Manar Aljohani, Jun Hou, Sindhura Kommu, and Xuan Wang. 2025. A comprehensive survey on the trustworthiness of large language models in healthcare. arXiv preprint arXiv:2502.15871. 

Manojit Bhattacharya, Soumen Pal, Srijan Chatterjee, Sang-Soo Lee, and Chiranjib Chakraborty. 2024. Large language model to multimodal large language model: A journey to shape the biological macromolecules to biological sciences and medicine. Molecular Therapy Nucleic Acids, 35(3). 

Zhenyu Bi, Sajib Acharjee Dip, Daniel Hajialigol, Sindhura Kommu, Hanwen Liu, Meng Lu, and Xuan Wang. 2024. Ai for biomedicine in the era of large language models. arXiv preprint arXiv:2403.15673. 

Haiyang Bian, Yixin Chen, Xiaomin Dong, Chen Li, Minsheng Hao, Sijie Chen, Jinyi Hu, Maosong Sun, Lei Wei, and Xuegong Zhang. 2024a. scmulan: a multitask generative pre-trained language model for singlecell analysis. In International Conference on Research in Computational Molecular Biology, pages 479–482. Springer. 

Haiyang Bian, Yixin Chen, Erpai Luo, Xinze Wu, Minsheng Hao, Lei Wei, and Xuegong Zhang. 2024b. General-purpose pre-trained large cellular models for single-cell transcriptomics. National Science Review, 11(11):nwae340. 

Guangshuo Cao, Haoyu Chao, Wenqi Zheng, Yangming Lan, Kaiyan Lu, Yueyi Wang, Ming Chen, He Zhang, and Dijun Chen. 2025. scplantllm: A foundation model for exploring single-cell expression atlases in plants. Genomics, Proteomics and Bioinformatics, 23(3):qzaf024. 

Shenghao Cao, Kaiyuan Yang, Jiabei Cheng, Jiachen Li, Hong-Bin Shen, Xiaoyong Pan, and Ye Yuan. 2024. stformer: a foundation model for spatial transcriptomics. bioRxiv, pages 2024–09. 

Ao Chen, Sha Liao, Mengnan Cheng, Kailong Ma, Liang Wu, Yiwei Lai, Xiaojie Qiu, Jin Yang, Jiangshan Xu, Shijie Hao, and 1 others. 2022. Spatiotemporal transcriptomic atlas of mouse organogenesis using dna nanoball-patterned arrays. Cell, 185(10):1777–1792. 

Chi-Jane Chen, Yuhang Chen, Sukwon Yun, Natalie Stanley, and Tianlong Chen. 2025a. Spatial coordinates as a cell language: A multi-sentence framework for imaging mass cytometry analysis. arXiv preprint arXiv:2506.01918. 

Kok Hao Chen, Alistair N Boettiger, Jeffrey R Moffitt, Siyuan Wang, and Xiaowei Zhuang. 2015. Spatially resolved, highly multiplexed rna profiling in single cells. Science, 348(6233):aaa6090. 

Xiaoyang Chen, Keyi Li, Xuejian Cui, Zian Wang, Qun Jiang, Jiacheng Lin, Zhen Li, Zijing Gao, Hairong Lv, and Rui Jiang. 2025b. Epiagent: foundation model for single-cell epigenomics. Nature Methods, pages 1–12. 

Yiqun Chen and James Zou. 2024. Genept: a simple but effective foundation model for genes and cells built from chatgpt. bioRxiv, pages 2023–10. 

Alexis Chevalier, Soumya Ghosh, Urvi Awasthi, James Watkins, Julia Bieniewska, Nichita Mitrea, Olga Kotova, Kirill Shkura, Andrew Noble, Michael Steinbaugh, and 1 others. 2025. Teddy: A family of foundation models for understanding single cell biology. arXiv preprint arXiv:2503.03485. 

The Tabula Sapiens Consortium*, Robert C Jones, Jim Karkanias, Mark A Krasnow, Angela Oliveira Pisco, Stephen R Quake, Julia Salzman, Nir Yosef, Bryan Bulthaup, Phillip Brown, and 1 others. 2022. The tabula sapiens: A multiple-organ, single-cell transcriptomic atlas of humans. Science, 376(6594):eabl4896. 

Haotian Cui, Alejandro Tejada-Lapuerta, Maria Brbic,´ Julio Saez-Rodriguez, Simona Cristea, Hani Goodarzi, Mohammad Lotfollahi, Fabian J Theis, and Bo Wang. 2025. Towards multimodal foundation models in molecular cell biology. Nature, 640(8059):623–633. 

Haotian Cui, Chloe Wang, Hassaan Maan, Nan Duan, and Bo Wang. 2022. scformer: a universal representation learning approach for single-cell data using transformers. bioRxiv, pages 2022–11. 

Haotian Cui, Chloe Wang, Hassaan Maan, Kuan Pang, Fengning Luo, Nan Duan, and Bo Wang. 2024. scgpt: toward building a foundation model for single-cell multiomics using generative ai. Nature methods, 21(8):1470– 1480. 

Darren A Cusanovich, Andrew J Hill, Delasa Aghamirzaie, Riza M Daza, Hannah A Pliner, Joel B Berletch, Galina N Filippova, Xingfan Huang, Lena Christiansen, William S DeWitt, and 1 others. 2018. A single-cell atlas of in vivo mammalian chromatin accessibility. Cell, 174(5):1309–1324. 

Sajib Acharjee Dip, Da Ma, and Liqing Zhang. 2024. Deepage: Harnessing deep neural network for epigenetic age estimation from dna methylation data of human blood samples. In Proceedings of the AAAI Symposium Series, volume 4, pages 267–274. 

Silvia Domcke, Andrew J Hill, Riza M Daza, Junyue Cao, Diana R O’Day, Hannah A Pliner, Kimberly A Aldinger, Dmitry Pokholok, Fan Zhang, Jennifer H Milbank, and 1 others. 2020. A human cell atlas of fetal chromatin accessibility. Science, 370(6518):eaba7612. 

Shengze Dong, Zhuorui Cui, Ding Liu, and Jinzhi Lei. 2025. scrdit: Generating single-cell rna-seq data by diffusion transformers and accelerating sampling. Interdisciplinary Sciences: Computational Life Sciences, pages 1–12. 

Xingyu Fan, Jiacheng Liu, Yaodong Yang, Chunbin Gu, Yuqiang Han, Bian Wu, Yirong Jiang, Guangyong Chen, and Pheng-Ann Heng. 2024. scgraphformer: unveiling cellular heterogeneity and interactions in scrna-seq data using a scalable graph transformer network. Communications Biology, 7(1):1463. 

Yin Fang, Kangwei Liu, Ningyu Zhang, Xinle Deng, Penghui Yang, Zhuo Chen, Xiangru Tang, Mark Gerstein, Xiaohui Fan, and Huajun Chen. 2024. Chatcell: Facilitating single-cell analysis with natural language. arXiv preprint arXiv:2402.08303. 

Xi Fu, Shentong Mo, Alejandro Buendia, Anouchka P Laurent, Anqi Shao, Maria del Mar Alvarez-Torres, Tianji Yu, Jimin Tan, Jiayu Su, Romella Sagatelian, and 1 others. 2025. A foundation model of transcription across human cell types. Nature, 637(8047):965–973. 

Jing Gong, Minsheng Hao, Xingyi Cheng, Xin Zeng, Chiming Liu, Jianzhu Ma, Xuegong Zhang, Taifeng Wang, and Le Song. 2023. xtrimogene: an efficient and scalable representation learner for single-cell rnaseq data. Advances in Neural Information Processing Systems, 36:69391–69403. 

Zhen-Hao Guo, Yan Wu, Siguo Wang, Qinhu Zhang, Jin-Ming Shi, Yan-Bin Wang, and Zhan-Heng Chen. 2023. scinterpreter: a knowledge-regularized generative model for interpretably integrating scrna-seq data. BMC bioinformatics, 24(1):481. 

Minsheng Hao, Jing Gong, Xin Zeng, Chiming Liu, Yucheng Guo, Xingyi Cheng, Taifeng Wang, Jianzhu Ma, Xuegong Zhang, and Le Song. 2024. Large-scale foundation model on single-cell transcriptomics. Nature methods, 21(8):1481–1491. 

Zhaohui He, Yuting Luo, Xinkai Zhou, Tao Zhu, Yangming Lan, and Dijun Chen. 2024. <? sty\usepackage $\mathrm { \{ w a s y s y m \} \mathrm { ? } > }$ scplantdb: a comprehensive database for exploring cell types and markers of plant cell atlases. Nucleic acids research, 52(D1):D1629–D1638. 

Chao Hui Huang. 2024. Qust-llm: Integrating large language models for comprehensive spatial transcriptomics analysis. arXiv preprint arXiv:2406.14307. 

Nauman Javed, Thomas Weingarten, Arijit Sehanobish, Adam Roberts, Avinava Dubey, Krzysztof Choromanski, and Bradley E Bernstein. 2025. A multi-modal transformer for cell type-agnostic regulatory predictions. Cell Genomics, 5(2). 

Boya Ji, Xiaoqi Wang, Debin Qiao, Liwen Xu, and Shaoliang Peng. 2024. Spaccc: Large language modelbased cell-cell communication inference for spatially resolved transcriptomic data. Big Data Mining and Analytics, 7(4):1129–1147. 

Yifeng Jiao, Yuchen Liu, Yu Zhang, Xin Guo, Yushuai Wu, Chen Jiang, Jiyang Li, Hongwei Zhang, Limei Han, Xin Gao, and 1 others. 2025. Chromfound: Towards a universal foundation model for single-cell chromatin accessibility data. arXiv preprint arXiv:2505.12638. 

Mehdi Joodaki, Mina Shaigan, Victor Parra, Roman D Bülow, Christoph Kuppe, David L Hölscher, Mingbo Cheng, James S Nagai, Michaël Goedertier, Nassim Bouteldja, and 1 others. 2024. Detection of patientlevel distances from single cell genomics and pathomics data with optimal transport (pilot). Molecular systems biology, 20(2):57–74. 

Jérémie Kalfon, Jules Samaran, Gabriel Peyré, and Laura Cantini. 2025. scprint: pre-training on 50 million cells allows robust gene network predictions. Nature Communications, 16(1):3607. 

Oussama Kharouiche, Aris Markogiannakis, Xiao Fei, Michail Chatzianastasis, and Michalis Vazirgiannis. 2025. Cell2text: Multimodal llm for generating singlecell descriptions from rna-seq data. arXiv preprint arXiv:2509.24840. 

Hilbert Yuen In Lam, Xing Er Ong, and Marek Mutwil. 2024. Large language models in plant biology. Trends in Plant Science, 29(10):1145–1155. 

Wei Lan, Zhentao Tang, Mingyang Liu, Qingfeng Chen, Wei Peng, Yiping Phoebe Chen, and Yi Pan. 2025. The large language models on biomedical data analysis: a survey. IEEE Journal of Biomedical and Health Informatics. 

Daniel Levine, Syed A Rizvi, Sacha Lévy, Nazreen Pallikkavaliyaveetil, David Zhang, Xingyu Chen, Sina Ghadermarzi, Ruiming Wu, Zihe Zheng, Ivan Vrkic, Anna Zhong, Daphne Raskin, Insu Han, Antonio Henrique De Oliveira Fonseca, Josue Ortega Caro, Amin Karbasi, Rahul Madhav Dhodapkar, and David Van Dijk. 2024. Cell2Sentence: Teaching large language models the language of biology. In Proceedings of the 41st International Conference on Machine Learning, volume 235 of Proceedings of Machine Learning Research, pages 27299–27325. PMLR. 

Russell Li, Heng Xu, and Eran A Mukamel. 2022. Epiattend: A transformer model of gene regulation combining single cell epigenomes with dna sequence. In NeurIPS 2022 Workshop on Learning Meaningful Representations of Life. 

Junhao Liu, Siwei Xu, Lei Zhang, and Jing Zhang. 2024. Single-cell omics arena: A benchmark study for large language models on cell type annotation using singlecell data. arXiv preprint arXiv:2412.02915. 

Tianyu Liu, Tianqi Chen, Wangjie Zheng, Xiao Luo, Yiqun Chen, and Hongyu Zhao. 2023. scelmo: Embeddings from language models are good learners for single-cell data analysis. bioRxiv, pages 2023–12. 

L Lu, Wen Xue, Xindian Wei, and 1 others. 2024. Sctrans: multi-scale scrna-seq sub-vector completion transformer for gene-selective cell type annotation. In 33rd 

International Joint Conference on Artificial Intelligence (IJCAI 2024), pages 5954–62. International Joint Conferences on Artificial Intelligence. 

Hiren Madhu, João Felipe Rocha, Tinglin Huang, Siddharth Viswanath, Smita Krishnaswamy, and Rex Ying. 2025. Heist: A graph foundation model for spatial transcriptomics and proteomics data. arXiv preprint arXiv:2506.11152. 

Yuren Mao, Yu Mi, Peigen Liu, Mengfei Zhang, Hanqing Liu, and Yunjun Gao. 2025. scagent: Universal single-cell annotation via a llm agent. arXiv preprint arXiv:2504.04698. 

Eleni P Mimitou, Caleb A Lareau, Kelvin Y Chen, Andre L Zorzetto-Fernandes, Yuhan Hao, Yusuke Takeshima, Wendy Luo, Tse-Shun Huang, Bertrand Z Yeung, Efthymia Papalexi, and 1 others. 2021. Scalable, multimodal profiling of chromatin accessibility, gene expression and protein levels in single cells. Nature biotechnology, 39(10):1246–1258. 

Daye Nam, Andrew Macvean, Vincent Hellendoorn, Bogdan Vasilescu, and Brad Myers. 2024. Using an llm to help with code understanding. In Proceedings of the IEEE/ACM 46th International Conference on Software Engineering, pages 1–13. 

Nima Nouri. 2025. Single-cell rna-seq data augmentation using generative fourier transformer. Communications Biology, 8(1):113. 

Nima Nouri, Ronen Artzi, and Virginia Savova. 2025. An agentic ai framework for ingestion and standardization of single-cell rna-seq data analysis. bioRxiv, pages 2025–07. 

Michelli Faria de Oliveira, Juan Pablo Romero, Meii Chung, Stephen R Williams, Andrew D Gottscho, Anushka Gupta, Susan E Pilipauskas, Seayar Mohabbat, Nandhini Raman, David J Sukovich, and 1 others. 2025. High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer. Nature Genetics, pages 1–12. 

Steven Palayew, Bo Wang, and Gary Bader. 2025. Towards applying large language models to complement single-cell foundation models. arXiv preprint arXiv:2507.10039. 

Cong Qi, Hanzhang Fang, Tianxing Hu, Siqi Jiang, and Wei Zhi. 2025a. Bidirectional mamba for single-cell data: Efficient context learning with biological fidelity. arXiv preprint arXiv:2504.16956. 

Yunjing Qi, Yulong Kan, Jing Qi, and Shuilin Jin. 2025b. scgt: I ntegration algorithm for single-cell rna-seq and atac-seq based on graph transforme r. Bioinformatics, page btaf357. 

Joseph M Replogle, Reuben A Saunders, Angela N Pogson, Jeffrey A Hussmann, Alexander Lenail, Alina Guna, Lauren Mascibroda, Eric J Wagner, Karen Adelman, Gila Lithwick-Yanai, and 1 others. 2022. Mapping information-rich genotype-phenotype landscapes with genome-scale perturb-seq. Cell, 185(14):2559–2575. 

Syed Asad Rizvi, Daniel Levine, Aakash Patel, Shiyang Zhang, Eric Wang, Sizhuang He, David Zhang, Cerise Tang, Zhuoyang Lyu, Rayyan Darji, Chang Li, Emily Sun, David Jeong, Lawrence Zhao, Jennifer Kwan, David Braun, Brian Hafler, Jeffrey Ishizuka, Rahul M. Dhodapkar, and 4 others. 2025. Scaling large language models for next-generation single-cell analysis. bioRxiv. 

Yanay Rosen, Yusuf Roohani, Ayush Agarwal, Leon Samotorcan, Tabula Sapiens Consortium, Stephen Rˇ Quake, and Jure Leskovec. 2023. Universal cell embeddings: A foundation model for cell biology. bioRxiv, pages 2023–11. 

Anna C Schaar, Alejandro Tejada-Lapuerta, Giovanni Palla, Robert Gutgesell, Lennard Halle, Mariia Minaeva, Larsen Vornholz, Leander Dony, Francesca Drummer, Mojtaba Bahrami, and 1 others. 2024. Nicheformer: a foundation model for single-cell and spatial omics. bioRxiv, pages 2024–04. 

Hongru Shen, Jilei Liu, Jiani Hu, Xilin Shen, Chao Zhang, Dan Wu, Mengyao Feng, Meng Yang, Yang Li, Yichen Yang, and 1 others. 2023. Generative pretraining from large-scale transcriptomes for single-cell deciphering. Iscience, 26(5). 

Hongru Shen, Xilin Shen, Mengyao Feng, Dan Wu, Chao Zhang, Yichen Yang, Meng Yang, Jiani Hu, Jilei Liu, Wei Wang, and 1 others. 2022. A universal approach for integrating super large-scale single-cell transcriptomes by exploring gene rankings. Briefings in Bioinformatics, 23(2). 

Yaorui Shi, Jiaqi Yang, Sihang Li, Junfeng Fang, Xiang Wang, Zhiyuan Liu, and Yang Zhang. 2025. Multimodal language modeling for high-accuracy single cell transcriptomics analysis and generation. arXiv e-prints, pages arXiv–2503. 

Christine N Shulse, Benjamin J Cole, Doina Ciobanu, Junyan Lin, Yuko Yoshinaga, Mona Gouran, Gina M Turco, Yiwen Zhu, Ronan C O’Malley, Siobhan M Brady, and 1 others. 2019. High-throughput single-cell transcriptome profiling of plant cell types. Cell reports, 27(7):2241–2247. 

Uddip Acharjee Shuvo, Sajib Acharjee Dip, Nirvar Roy Vaskar, and ABM Alim Al Islam. 2024. Assessing chatgpt’s code generation capabilities with short vs long context programming problems. In Proceedings of the 11th International Conference on Networking, Systems, and Security, pages 32–40. 

Robert R Stickels, Evan Murray, Pawan Kumar, Jilong Li, Jamie L Marshall, Daniela J Di Bella, Paola Arlotta, Evan Z Macosko, and Fei Chen. 2021. Highly sensitive spatial transcriptomics at near-cellular resolution with slide-seqv2. Nature biotechnology, 39(3):313–319. 

Marlon Stoeckius, Christoph Hafemeister, William Stephenson, Brian Houck-Loomis, Pratip K Chattopadhyay, Harold Swerdlow, Rahul Satija, and Peter Smibert. 2017. Simultaneous epitope and transcriptome measurement in single cells. Nature methods, 14(9):865–868. 

Elliott Swanson, Cara Lord, Julian Reading, Alexander T Heubeck, Palak C Genge, Zachary Thomson, Morgan DA Weiss, Xiao-jun Li, Adam K Savage, Richard R Green, and 1 others. 2021. Simultaneous trimodal single-cell measurement of transcripts, epitopes, and chromatin accessibility using tea-seq. Elife, 10:e63632. 

Wenzhuo Tang, Hongzhi Wen, Renming Liu, Jiayuan Ding, Wei Jin, Yuying Xie, Hui Liu, and Jiliang Tang. 2023. Single-cell multimodal prediction via transformers. In Proceedings of the 32nd ACM International Conference on Information and Knowledge Management, pages 2422–2431. 

Christina V Theodoris, Ling Xiao, Anant Chopra, Mark D Chaffin, Zeina R Al Sayed, Matthew C Hill, Helene Mantineo, Elizabeth M Brydon, Zexian Zeng, X Shirley Liu, and 1 others. 2023. Transfer learning enables predictions in network biology. Nature, 618(7965):616–624. 

Kyle J Travaglini, Ahmad N Nabhan, Lolita Penland, Rahul Sinha, Astrid Gillich, Rene V Sit, Stephen Chang, Stephanie D Conley, Yasuo Mori, Jun Seita, and 1 others. 2020. A molecular cell atlas of the human lung from single-cell rna sequencing. Nature, 587(7835):619–625. 

Yuequn Wang, Jun Wang, Yanyu Xu, Ning Liu, Bin Liu, Yuliang Li, and Guoxian Yu. 2025. Fmh2st: foundation model-based spatial transcriptomics generation from histological images. Nucleic Acids Research, 53(17):gkaf865. 

Hongzhi Wen, Wenzhuo Tang, Xinnan Dai, Jiayuan Ding, Wei Jin, Yuying Xie, and Jiliang Tang. 2023. Cellplm: Pre-training of cell language model beyond single cells. BioRxiv, pages 2023–10. 

Juncheng Wu, Changxin Wan, Zhicheng Ji, Yuyin Zhou, and Wenpin Hou. 2025. Epifoundation: A foundation model for single-cell atac-seq via peak-to-gene alignment. bioRxiv. 

Yuxuan Wu and Fuchou Tang. 2025. scextract: leveraging large language models for fully automated singlecell rna-seq data annotation and prior-informed multidataset integration. Genome Biology, 26(1):174. 

Yihang Xiao, Jinyi Liu, Yan Zheng, Xiaohan Xie, Jianye Hao, Mingzhi Li, Ruitao Wang, Fei Ni, Yuxiao Li, Jintian Luo, and 1 others. 2024. Cellagent: An llm-driven multi-agent framework for automated single-cell data analysis. arXiv preprint arXiv:2407.09811. 

Jing Xu, De-Shuang Huang, and Xiujun Zhang. 2024. scmformer integrates large-scale single-cell proteomics and transcriptomics data by multi-task transformer. Advanced Science, 11(19):2307835. 

Fan Yang, Wenchuan Wang, Fang Wang, Yuan Fang, Duyu Tang, Junzhou Huang, Hui Lu, and Jianhua Yao. 2022. scbert as a large-scale pretrained deep language model for cell type annotation of single-cell rna-seq data. Nature Machine Intelligence, 4(10):852–866. 

Wenyi Yang, Pingping Wang, Shouping Xu, Tao Wang, Meng Luo, Yideng Cai, Chang Xu, Guangfu Xue, Jinhao Que, Qian Ding, and 1 others. 2024a. Deciphering cell–cell communication at single-cell resolution for spatial transcriptomics with subgraph-based graph attention network. Nature Communications, 15(1):7101. 

Xiaodong Yang, Guole Liu, Guihai Feng, Dechao Bu, Pengfei Wang, Jie Jiang, Shubai Chen, Qinmeng Yang, Hefan Miao, Yiyang Zhang, and 1 others. 2024b. Genecompass: deciphering universal gene regulatory mechanisms with a knowledge-informed cross-species foundation model. Cell Research, 34(12):830–845. 

Wenjin Ye, Yuanchen Ma, Junkai Xiang, Hongjie Liang, Tao Wang, Qiuling Xiang, Andy Peng Xiang, Wu Song, Weiqiang Li, and Weijun Huang. 2024. Objectively evaluating the reliability of cell type annotation using llm-based strategies. arXiv preprint arXiv:2409.15678. 

Zhen Yuan, Shaoqing Jiao, Yihang Xiao, and Jiajie Peng. 2025. scmamba: A scalable foundation model for single-cell multi-omics integration beyond highly variable feature selection. arXiv preprint arXiv:2506.20697. 

Yuansong Zeng, Jiancong Xie, Ningyuan Shangguan, Zhuoyi Wei, Wenbing Li, Yun Su, Shuangyu Yang, Chengyang Zhang, Jinbo Zhang, Nan Fang, and 1 others. 2025. Cellfm: a large-scale foundation model pretrained on transcriptomics of 100 million human cells. Nature Communications, 16(1):4679. 

Fan Zhang, Hao Chen, Zhihong Zhu, Ziheng Zhang, Zhenxi Lin, Ziyue Qiao, Yefeng Zheng, and Xian Wu. 2025a. A survey on foundation language models for single-cell biology. In Proceedings of the 63rd Annual Meeting of the Association for Computational Linguistics (Volume 1: Long Papers), pages 528–549. 

Fan Zhang, Tianyu Liu, Zhihong Zhu, Hao Wu, Haixin Wang, Donghao Zhou, Yefeng Zheng, Kun Wang, Xian Wu, and Pheng-Ann Heng. 2025b. Cellverse: Do large language models really understand cell biology? arXiv preprint arXiv:2505.07865. 

Suyuan Zhao, Jiahuan Zhang, and Zaiqing Nie. 2023. Large-scale cell representation learning via divideand-conquer contrastive learning. arXiv preprint arXiv:2306.04371. 

# A Data Collection

# A.1 Search Strategy

To construct a comprehensive registry of large language models (LLMs) and agentic frameworks in single-cell biology, we systematically screened both peer-reviewed and preprint literature from 2020–2025 across multiple sources. Searches were conducted on PubMed, Google Scholar, arXiv, and Semantic Scholar using combinations of domainand method-specific keywords, including: 

• General: “large language model”, “foundation model”, “LLM”, “multimodal model”, “generative AI” 

• Domain-specific: “single-cell”, “scRNA-seq”, “scATAC”, “multiome”, “spatial transcriptomics”, “perturb-seq”, “cell atlas”, “biomedical foundation model” 

• Integration / reasoning: “ontology alignment”, “text-bridge”, “agentic framework”, “biological reasoning”, “cross-modal integration”, “annotation model” 

After removing duplicates, the combined queries returned approximately 8,020 papers from 2020–2025. Of these, 5,510 remained after filtering for English-language articles with available abstracts. We manually reviewed titles, abstracts, and model descriptions to identify works directly involving large-scale pretrained models, multimodal transformers, or language-based reasoning systems applied to single-cell or cellular-level omics data. 

# A.2 Inclusion and Screening Criteria

Studies were included if they met at least one of the following criteria: 

1. Introduced or evaluated a large pretrained model $( \ge 1 0 \mathbf { M }$ parameters) for single-cell or multi-omic analysis. 

2. Employed textual grounding, ontology prompts, or natural-language interfaces to interpret single-cell data. 

3. Proposed agentic or reasoning-based frameworks for biological tasks (annotation, trajectory, perturbation, or integration). 

4. Released code, preprints, or benchmarks relevant to cellular foundation models. 

Methodological, biological, and dataset-based duplicates were merged under representative entries. This curation yielded 58 distinct models across five methodological families: foundation, textbridge, spatial/multimodal, epigenomic, and agentic frameworks. Each was annotated with task do-

main, modality, training data, and evaluation scope (Appendix Table 7). 

# A.3 Inclusion of Preprints

We included arXiv and bioRxiv preprints to ensure coverage of recent and influential models not yet published in peer-reviewed venues. In single-cell research, the majority of foundational architectures (e.g., scGPT, scFoundation, CellLM) initially appeared as preprints months before journal acceptance. Given the field’s rapid evolution and small number of LLM-scale models, preprints represent essential early contributions to reproducibility and transparency. Each preprint was manually crossverified for active GitHub or Zenodo links to confirm technical validity and data availability. 

# A.4 Resulting Corpus

The final corpus spans 58 methods across 37 unique institutions, $^ { 4 0 + }$ publicly accessible datasets, and 8 major analytical tasks. This dataset–method pairing forms the empirical foundation for the taxonomic and domain analyses presented in the main paper. 

# B Extended Dataset Summaries

# B.1 RNA Atlases

Foundational scRNA-seq atlases such as Tabula Sapiens (>1 M cells, 24 tissues) and Tabula Muris Senis (mouse, multi-organ across aging) provide references for annotation and batch-effect evaluation. Organ-specific datasets like the Human Lung Cell Atlas (HLCA) and the Allen Brain Map enable domain-specific benchmarking and ontology-aware evaluation. 

# B.2 Chromatin and ATAC Data

The Cusanovich mouse sci-ATAC atlas, human adult/fetal scATAC atlases, and joint RNA–ATAC modalities (SHARE-seq, SNARE-seq2) characterize regulatory landscapes supporting trajectory and GRN inference, though sparsity limits large-scale pretraining. 

# B.3 Multiome and Tri-Modal Data

Datasets such as TEA-seq, DOGMA-seq, ASAP-seq, and CITE-seq pair transcriptomic, epigenomic, and proteomic modalities, enabling cross-modal translation and representation learning. The Multiome Benchmark Pack aggregates ${ > } 4 0$ datasets for standardization but remains smaller than large RNA atlases. 

# B.4 Spatial Transcriptomics and Imaging

Platforms including 10x Visium/HD, Slide-seqV2, and DBiT-seq provide spatial gene-expression maps, while imaging technologies (MERFISH, STARmap, Stereo-seq, CosMx, Xenium) capture near single-molecule resolution for tissue-level reasoning. 

# B.5 Perturbation and Drug-Response Screens

Large CRISPR-based datasets (Replogle 2022 Perturb-seq, Norman 2019, sci-Plex) and community benchmarks (Virtual Cell Challenge) quantify transcriptional responses for causal and policy evaluation in agentic frameworks. 

# B.6 Plant Single-Cell Datasets

scPlantDB (67 datasets, 17 species, 2.5 M cells) and PlantscRNAdb provide unified plant atlases; Arabidopsis E-CURD-4 and transgenic tobacco scRNA-seq datasets enable cross-kingdom and stress-response modeling. 

# B.7 Critical Observations

(1) RNA atlases dominate in scale; (2) Chromatin and spatial data are fragmented; (3) Privacy and licensing restrict clinical data; (4) Cross-species coverage is limited; (5) Paired modalities offer promising benchmarks for next-generation LLMs. Representative datasets and links are provided in Tables A1–A3. 

# C Extended Details on Model Taxonomy

# C.1 Foundation Models (Extended Details)

Foundation-scale models constitute the base layer for language-driven single-cell analysis. They are trained on molecular profiles—typically scRNAseq or integrated multi-organ datasets—to learn generalizable representations of cellular states, gene programs, and perturbations. 

Training corpora and scale. Most use atlases containing $1 0 ^ { 6 } – 1 0 ^ { 8 }$ cells. scGPT and scFoundation were pretrained on human–mouse atlases from Tabula Sapiens and HCA; Geneformer aggregates thousands of GEO datasets; iSEEEK and CellFM span 11–17 organs and species, enabling crosstissue generalization. 

Architectural variations. Most adopt Transformer backbones with gene tokens as contextual units. tGPT and scBERT encode ranked or discrete gene sequences, while Geneformer and sc-Foundation use masked-token modeling. Hybrid 

designs such as scGraphformer (Fan et al., 2024) and scRDiT (Dong et al., 2025) incorporate graph or diffusion-based attention, improving spatial and regulatory context. 

Learning objectives. Common objectives include masked-gene prediction, expression denoising, rank-based reconstruction, and cross-species contrastive learning (UCE, GeneCompass). Optimization typically minimizes cosine or KL divergence across gene embeddings, similar to maskedlanguage modeling in NLP. 

Applications and limitations. Foundation embeddings transfer effectively to annotation, integration, and trajectory inference but lack explicit biological grounding. Ontology alignment and reasoning remain external, and interpretability is largely post-hoc (attention heatmaps, generanking). These limitations prompted the development of Text-Bridge and Agentic frameworks discussed in later sections. 

# C.2 Text-Bridge LLMs (Extended Details)

These models explicitly connect molecular and textual modalities via ontology labels, pathway terms, or literature embeddings. scELMo fuses ELMo-based gene metadata with cell latent vectors; CellLM uses ontology-aware prompts for naturallanguage annotation; Cell2Sentence and Cell2Text map gene-expression vectors to descriptive sentences through cross-modal contrastive loss; and GenePT jointly pre-trains on PubMed to embed genes and diseases in the same space. Architectures typically adopt Transformer encoders for molecular data and domain-tuned BERT variants for text, trained with contrastive or cosine-similarity objectives. These systems rank highest in Explainability and Emerging paradigms scores (Appendix Table 7) but remain limited by vocabulary coverage and the lack of full multi-omic grounding. 

# C.3 Spatial and Multimodal Models (Extended Details)

Spatial models combine molecular and positional information to learn tissue-aware embeddings. TransformerST performs cross-scale attention between spatial spots and gene tokens; spaLLM aligns histology-derived captions with expression vectors using LLM guidance; and OmiCLIP links image and molecular embeddings via contrastive pretraining. Multi-omic extensions such as scMMGPT and FmH2ST integrate RNA, ATAC, and protein 

modalities within a shared latent space. Grounding typically relies on marker-based or atlas alignment rather than textual ontologies. Evaluation metrics include spot-level reconstruction accuracy and correlation with histological segmentation (DLPFC, Visium HD). These models show high biological grounding but moderate explainability due to visual-feature opacity. 

# C.4 Epigenomic Models (Extended Details)

Inputs include chromatin-accessibility matrices, motif sequences, and enhancer–promoter links. EpiFoundation pre-trains on the ENCODE scATAC compendium; EpiBERT encodes open-chromatin sequences using masked-region objectives; and EpiAttend models enhancer–promoter coupling through cross-region attention. State-space architectures (scMamba, GeneMamba) enable efficient context propagation over tens of thousands of genomic peaks. ChromFound and GET jointly model RNA and ATAC modalities for regulatory inference. Grounding leverages atlas-derived cCRE annotations and transcription-factor motifs. Performance is strong for biological grounding and multi-omics alignment, but interpretability is limited to motif attention visualization. 

# C.5 Agentic Frameworks (Extended Details)

Agentic models couple (i) a pretrained molecular encoder, (ii) an LLM-based controller, and (iii) external tool interfaces. scAgent executes multistep annotation workflows via ontology queries; CellVerse coordinates multiple domain agents for transcriptomic, spatial, and literature reasoning; EpiAgent extends agentic control to enhancer–gene analysis; and Teddy/scPilot (Joodaki et al., 2024) provide lightweight orchestrators for benchmarking pipelines. Grounding sources include the Human Cell Atlas ontology, UBERON, and PubMed. These systems score highest in Explainability and Emerging paradigms, marking a transition from static embeddings to interactive, goal-directed modeling. Open challenges include evaluation of reasoning accuracy, privacy, and reproducibility. 

# D Datasets


Table 1: RNAseq single-cell datasets used in LLM-based single-cell research.


<table><tr><td>Dataset</td><td>Tasks</td><td>Description</td><td>Scale</td><td>Link / Citation</td></tr><tr><td>Tabula Sapiens v2</td><td>Annotation, Integration</td><td>Multi-organ human atlas (~1.1 M cells, 28 organs across 24 donors) capturing cell-type heterogeneity and cross-tissue transcriptional variation; droplet and plate-seq modalities.</td><td>1.1 M cells / 28 organs / 24 donors</td><td>CZ Biohub Portal</td></tr><tr><td>Tabula Muris (mouse, multi-organ)</td><td>Annotation, Integration</td><td>Mouse multi-organ atlas of gene expression combining droplet and FACS; enables cross-tissue comparison and batch integration.</td><td>~100,000 cells / 20 organs</td><td>HCA Project Page</td></tr><tr><td>Mouse Cell Atlas (scMCA)</td><td>Annotation, Integration</td><td>Comprehensive mouse single-cell atlas constructed with Microwell-seq; supports cell-type matching via scMCA tool.</td><td>&gt;400,000 cells / &gt;40 tissues</td><td>MCA Portal (ZJU)</td></tr><tr><td>Human Lung Cell Atlas (HLCA v1.0)</td><td>Annotation, Integration, Disease Modeling</td><td>Integrated human lung reference built from 49 scRNA-seq datasets (~2.4 M cells) across 16 studies; unified epithelial, immune, and stromal annotations.</td><td>~2.4 M / 49 datasets / 16 studies</td><td>HLCA v1.0 Portal</td></tr><tr><td>HCA lung project (example project page)</td><td>Annotation, Integration, Disease Modeling</td><td>HCA lung cohort integrating ~2.4 M single cells from 49 datasets; harmonized annotations for airway, immune, and endothelial populations.</td><td>2.4 M cells / 49 datasets</td><td>HLCA Project Page</td></tr><tr><td>Allen Brain Map – Cell Types RNA-seq (human &amp; mouse)</td><td>Annotation, Cross-species, Trajectory</td><td>Single-cell and single-nucleus transcriptomes from human and mouse brain regions; supports cross-species comparison and cell-type taxonomy benchmarking.</td><td>&gt;1.8 M cells / multiple cortical &amp; subcortical regions</td><td>Allen Brain Cell Types Database</td></tr><tr><td>Yale Lung Disease Cell Atlases (e.g., COPD)</td><td>Annotation, Disease Modeling</td><td>Single-cell RNA atlases of diseased human lungs (COPD, IPF) from Yale&#x27;s Kaminski lab, capturing altered epithelial, endothelial, and immune populations.</td><td>~300k cells (IPF + COPD + controls)</td><td>HCA/Yale Atlas</td></tr></table>


Table 2: ATACseq and chromatin-accessibility datasets used in LLM-based single-cell research.


<table><tr><td>Dataset</td><td>Tasks</td><td>Description</td><td>Scale</td><td>Link / Citation</td></tr><tr><td>Mouse sci-ATAC-seq atlas (Cusanovich et al., 2018)</td><td>Annotation, Trajectory, GRN inference</td><td>Landmark single-cell ATAC-seq atlas profiling ~100,000 nuclei across 13 adult mouse tissues using combinatorial indexing (sci-ATAC-seq); enables cross-tissue regulatory and lineage analysis.</td><td>~100k nuclei / 13 tissues</td><td>GSE111586, Science 2018</td></tr><tr><td>Human adult scATAC atlas (Zhang et al., 2021)</td><td>Annotation, GRN, Cross-tissue integration</td><td>Comprehensive single-cell chromatin accessibility atlas of adult human tissues, defining candidate cis-regulatory elements (cCREs) across 25 tissues and 222 cell types; foundation of the ENCODE human cCRE registry.</td><td>~472k nuclei / 25 tissues</td><td>Nature 2021, ENCODE Portal</td></tr><tr><td>Human fetal scATAC atlas (Domcke et al., 2020; GSE149683)</td><td>Annotation, Developmental trajectory, GRN</td><td>Single-cell ATAC-seq atlas of human fetal tissues profiling 15 organs across mid-gestation; reveals developmental enhancer activity and regulatory lineage trajectories.</td><td>~720k nuclei / 15 organs</td><td>Science 2020, GSE149683</td></tr><tr><td>Massively parallel scATAC (Satpathy et al., 2019)</td><td>Annotation, Trajectory, Regulatory modeling</td><td>Pioneering high-throughput scATAC-seq dataset of ~200k immune and cancer nuclei enabling scalable mapping of open-chromatin landscapes; forms benchmark for lineage and immune-cell trajectory studies.</td><td>~200k nuclei / blood and tumor tissues</td><td>Nat. Biotechnol. 2019, GSE123581</td></tr><tr><td>ENCODE portal (single-cell experiments)</td><td>Annotation, Integration, Regulatory modeling</td><td>Centralized repository from the ENCODE Consortium aggregating thousands of single-cell RNA and ATAC assays from human and mouse tissues; provides uniformly processed metadata, peak calls, and cCRE annotations for benchmarking.</td><td>&gt;2,000 single-cell assays / multiple tissues</td><td>ENCODE Portal</td></tr><tr><td>T cell epigenetic atlas (Giles et al., 2022)</td><td>Trajectory, GRN, Disease modeling</td><td>Single-cell ATAC-seq atlas profiling chromatin accessibility across human T-cell activation, exhaustion, and differentiation states; supports trajectory reconstruction and immune-epigenetic modeling.</td><td>~150k nuclei / T-cell subsets</td><td>Nat. Immunol. 2022</td></tr></table>


Table 3: Multiome datasets (paired/tri-omic: RNA $^ +$ ADT, RNA $^ +$ ATAC, ATAC $^ +$ ADT, and tri-modal) used in LLM-based single-cell research.


<table><tr><td>Dataset</td><td>Tasks</td><td>Description</td><td>Scale</td><td>Link / Citation</td></tr><tr><td>UCSC Cell Browser Hub</td><td>Annotation, Integration, Visualization</td><td>Aggregated repository of hundreds of public single-cell datasets across species and modalities with metadata and embeddings; useful for exploratory analysis and reference selection.</td><td>&gt;1,200 datasets / multi-species</td><td>UCSC Cell Browser Hub</td></tr><tr><td>Human Cell Atlas (HCA) data browser (multi-project)</td><td>Annotation, Integration, Cross-project mapping</td><td>A global human cell atlas aggregating ~63.2 M cells across 515 projects and &gt;11,000 donors, spanning diverse tissues and modalities; enables cross-study reference and meta-integration.</td><td>~63.2 M cells / 515 projects / 11,000+ donors</td><td>HCA Data Portal</td></tr><tr><td>Azimuth reference collections (PBMC, lung, kidney, fetal) (RNA + ADT)</td><td>Annotation, Integration</td><td>Curated single-cell reference atlases by the Satija Lab for automated cell-type mapping in Seurat/Azimuth; harmonized labels and multimodal ADT features.</td><td>100k-1 M cells across multiple organs</td><td>Azimuth Data Portal</td></tr><tr><td>TEA-seq (tri-modal RNA + ATAC + ADT; GSE158013)</td><td>Integration, Perturbation, Multi-modal reasoning</td><td>Tri-modal PBMC dataset measuring transcriptome, chromatin accessibility, and surface proteins; supports alignment/fusion of omics layers and pretraining.</td><td>~100k cells</td><td>GSE158013</td></tr><tr><td>DOGMA-seq (RNA + ATAC + Protein)</td><td>Integration, Perturbation, Cross-modal reasoning</td><td>Tri-modal profiling of transcriptome, chromatin, and proteins in human immune cells; benchmark for joint embeddings and cross-modal translation.</td><td>~50k cells / PBMCs</td><td>Nature Cell Biol. 2022 / GSE184715</td></tr><tr><td>ASAP-seq (ATAC + Protein)</td><td>Integration, GRN, Epigenetic modeling</td><td>Paired chromatin accessibility and surface proteins via antibody-derived tags; links cis-regulatory variation with immune phenotypes.</td><td>~100k nuclei / PBMCs</td><td>Nat. Biotechnol. 2021 / GSE162690</td></tr><tr><td>CITE-seq compendia (RNA + ADT)</td><td>Annotation, Integration, Transfer learning</td><td>Large compendium of paired RNA and surface-protein profiles across blood and tissue; enables multimodal foundation pretraining and zero-shot annotation.</td><td>&gt;500k cells across multiple tissues</td><td>Nat. Methods 2019 / GSE100866</td></tr><tr><td>Multiome Benchmark Pack (QuKun Lab)</td><td>Integration, Scalability, Batch-effect analysis</td><td>Public benchmark suite consolidating 25 RNA + Protein, 12 RNA + ATAC, and 4 tri-omic datasets (CITE/TEA/DOGMA), standardized for cross-modal and large-scale integration testing.</td><td>41 datasets total</td><td>Benchmark Portal</td></tr></table>


Table 4: Spatial transcriptomics and imaging datasets used in LLM-based single-cell research.


<table><tr><td>Dataset</td><td>Tasks</td><td>Description</td><td>Scale</td><td>Link / Citation</td></tr><tr><td>Slide-seq / Slide-seqV2</td><td>Spatial mapping, Trajectory, Integration</td><td>High-resolution spatial transcriptomics of mouse brain (hippocampus, cerebellum, cortex) at ~10 μm bead resolution; supports neighborhood and spatial-domain reconstruction benchmarks.</td><td>~50k-100k spots per tissue</td><td>Nat. Biotechnol. 2020 / SCP815</td></tr><tr><td>MERFISH (imaging)</td><td>Spatial mapping, Pathway, Annotation</td><td>Multiplexed error-robust fluorescence in situ hybridization (MERFISH) datasets profiling millions of mouse-brain cells with 3D coordinates and 483-1,000-gene panels; benchmark for high-resolution spatial mapping.</td><td>~4 M cells / 483-1,000 genes</td><td>Vizgen MERFISH Portal / Science 2018</td></tr><tr><td>Stereo-seq (STOmics)</td><td>Spatial mapping, Developmental trajectory</td><td>Genome-scale Stereo-seq datasets with submicron resolution; includes MOSTA (mouse embryo) and 3D Drosophila atlases for developmental and cross-species modeling.</td><td>~100 μm × 100 μm tiles / millions of spots</td><td>Cell 2022 / STOmics Data Hub</td></tr><tr><td>DBiT-seq &amp; spatial multi-omics</td><td>Spatial mapping, Integration, Multi-omic reasoning</td><td>Deterministic barcoding in tissue (DBiT-seq) capturing spatial RNA and protein expression in mouse embryo and human lymph node; supports multi-omic spatial benchmarks.</td><td>~20k spots / 100 μm grids</td><td>Nat. Biotechnol. 2020 / GSE152506</td></tr><tr><td>10x Visium / Visium HD</td><td>Spatial mapping, Annotation, Integration</td><td>Widely used capture-array platform for spatial transcriptomics; Visium uses 55 μm spots, Visium HD extends to 2 μm grids with improved gene recovery and FFPE compatibility.</td><td>~5k-55k spots per tissue</td><td>10x Genomics Visium Portal</td></tr><tr><td>Xenium / CosMx (in situ)</td><td>Spatial mapping, Clinical translation, Privacy</td><td>In situ spatial transcriptomics platforms (10x Xenium, NanoString CosMx) profiling thousands of transcripts in human FFPE and fresh-frozen tissues (e.g., breast, colon, NSCLC); bridge omics and histology for clinical translation.</td><td>10k-100k cells per section</td><td>Xenium Explorer / CosMx Portal</td></tr></table>


Table 5: Perturbation and Drug-Response single-cell datasets used in LLM-based frameworks.


<table><tr><td>Dataset</td><td>Tasks</td><td>Description</td><td>Scale</td><td>Link / Citation</td></tr><tr><td>Genome-scale 
Perturb-seq 
(Replogle 2022) 
processed data</td><td>Perturbation, 
GRN, Causal 
inference</td><td>Largest CRISPR-based Perturb-seq dataset profiling ~2.5 M single cells across &gt;2 000 genetic perturbations; benchmark for causal network inference and representation learning.</td><td>~2.5 M 
cells / &gt;2 000 pertur- 
bations</td><td>Cell 2022 / 
Figshare 
Dataset</td></tr><tr><td>Norman et al. 
2019 Perturb-seq</td><td>Perturbation, 
Combinatorial 
GRN</td><td>Foundational combinatorial CRISPR 
Perturb-seq dataset (immune + cancer 
models); establishes feasibility of pooled 
functional genomics at single-cell 
resolution.</td><td>~200 k 
cells / &gt;250 
combina- 
tions</td><td>Science 
2019 / 
GSE133344</td></tr><tr><td>Dixit et al. 2016 
(GSE90063)</td><td>Perturbation, 
Combinatorial 
GRN</td><td>Early Perturb-seq study targeting 24 genes 
in macrophages; establishes pipeline for 
pooled CRISPR screens with single-cell 
RNA-seq.</td><td>~100 k 
cells / 24 
target genes</td><td>Cell 2016 / 
GSE90063</td></tr><tr><td>Adamson et al. 
2016</td><td>Perturbation, 
Stress-response, 
GRN</td><td>Single-cell CRISPR Perturb-seq of 
ER-stress pathways in K562 cells; 
benchmark for pathway-level perturbation 
responses.</td><td>~30 k cells / 
10 perturba- 
tions</td><td>Cell 2016 / 
GSE90060</td></tr><tr><td>sci-Plex collection 
(drugs screens) + 
cellxgene portal</td><td>Perturbation, 
Stress-response, 
GRN</td><td>Multiplexed chemical-perturbation 
scRNA-seq of &gt;650 k cells treated with 188 
compounds via sci-Plex barcoding; 
pharmacotranscriptomic profiles for drug 
response modeling.</td><td>~650 k 
cells / 188 
drugs</td><td>Cell 2020 / 
cellxgene 
Collection</td></tr><tr><td>Compressed 
Perturb-seq 
(immune LPS, 598 
genes)</td><td>Drug response, 
Perturbation, 
Benchmarking</td><td>Low-multiplicity Perturb-seq library 
targeting 598 immune genes in macrophages 
under LPS stimulation; enables compressed 
experimental design for causal modeling.</td><td>~300 k 
cells / 598 
targets</td><td>bioRxiv 
2021 / 
GSE179924</td></tr><tr><td>In vivo 
Perturb-seq (brain 
/ ASD genes)</td><td>Perturbation, 
GRN, 
Neuroscience</td><td>In vivo CRISPR Perturb-seq targeting ~30 
ASD-linked genes in mouse cortex; maps 
neuronal gene-regulatory networks in native 
contexts.</td><td>~200 k 
cells / 30 
genes</td><td>Science 
2023 / 
GSE216113</td></tr><tr><td>Virtual Cell 
Challenge PBMC 
cytokine 
perturbations</td><td>Perturbation, 
GRN, 
Benchmarking</td><td>Community benchmark dataset (ARC 
Institute 2024) of ~300 k PBMC cells 
under cytokine stimulation and CRISPRi 
perturbations; standardized splits for LLM 
and state-transition evaluation.</td><td>~300 k 
cells / 150 
perturbed 
genes</td><td>ARC 
Challenge 
Repository 
/ 
Kaggle</td></tr></table>


Table 6: Plant RNA single-cell transcriptomics datasets used in LLM research.


<table><tr><td>Dataset</td><td>Tasks</td><td>Description</td><td>Scale</td><td>Link / Citation</td></tr><tr><td>scPlantDB (meta-collection)</td><td>Annotation, Integration, Meta-analysis</td><td>Comprehensive plant single-cell transcriptome database integrating 67 datasets from 17 plant species (~2.5 M cells); unified preprocessing and annotations; supports cross-species modeling and plant-specific LLMs.</td><td>~2.5 M cells / 67 datasets / 17 species</td><td>scPlantDB Portal</td></tr><tr><td>PlantscRNAdb</td><td>Annotation, Marker discovery</td><td>Curated database of plant cell-type marker genes across four species (Arabidopsis, rice, maize, tomato); supports ontology-based annotation and cell identity benchmarking.</td><td>4 species / multiple tissues</td><td>PlantscRNAdb Portal</td></tr><tr><td>Arabidopsis scRNA-seq (E-CURD-4)</td><td>Annotation, Trajectory, Development</td><td>Baseline scRNA-seq dataset of Arabidopsis thaliana root and leaf tissues (~10,779 cells); used for developmental lineage and differentiation analysis.</td><td>10,779 cells / 2 tissues</td><td>EBI Array-Express E-CURD-4</td></tr><tr><td>Tobacco leaf scRNA-seq (transgenic antibody line)</td><td>Annotation, Perturbation, Stress response</td><td>Single-cell transcriptomic profiling of transgenic Nicotiana tabacum leaves expressing llama antibody; identifies immune-like responses and cell-type heterogeneity under genetic perturbation.</td><td>~25k cells / leaf tissue</td><td>Sci. Data 2023</td></tr><tr><td>Plant scRNA Browser (PscB)</td><td>Annotation, Visualization, Cross-tissue integration</td><td>Online visualization hub aggregating plant scRNA-seq datasets (Arabidopsis, rice, Wolffia) with harmonized metadata and interactive UMAP-based cell-type search.</td><td>15+ datasets / 3+ species</td><td>PscB Data Portal</td></tr></table>

# E Methods Comparison


Table 7: Full comparison of single-cell LLM and agentic methods (Domains column removed). Year merged into Model.


<table><tr><td>Model (Year)</td><td colspan="2">Published Category where</td><td>Modality</td><td>Grounding Type</td><td>Agentic (Y/N)</td><td>Primary Task</td><td>Domain Score</td></tr><tr><td>scGPT (Cui et al., 2024)</td><td>Nature Methods</td><td>Foundation</td><td>Multiomics (scRNA + optionally multi-omic modes)</td><td>Atlas (trained on large atlas of single-cell datasets)</td><td>No</td><td>Annotation, Integration, Perturbation (annotation as primary)</td><td>6</td></tr><tr><td>Geneformer (Theodoris et al., 2023)</td><td>Nature</td><td>Foundation</td><td>scRNA</td><td>None / Rank-based (implicit)</td><td>No</td><td>Cell classification, in-silico perturbation, network prediction</td><td>6</td></tr></table>

Continued on next page 

<table><tr><td>Model (Year)</td><td>Published Category where</td><td>Modality</td><td>Grounding Type</td><td>Agentic (Y/N)</td><td>Primary Task</td><td>Domain Score</td></tr><tr><td>scFoundation (Hao et al., 2024)</td><td>Nature Methods</td><td>scRNA / multi-omics (transcriptomics focus)</td><td>Atlas (large pretrained corpus)</td><td>No</td><td>Cell annotation, perturbation prediction, drug response, gene module inference</td><td>7</td></tr><tr><td>CellFM (Zeng et al., 2025)</td><td>Nature Communications</td><td>scRNA</td><td>Atlas / value-projection</td><td>No</td><td>Learns embeddings from 100M cells and supports annotation, perturbation, gene function</td><td>7</td></tr><tr><td>iSEEEK (Shen et al., 2022)</td><td>Briefings in Bioinformatics</td><td>scRNA</td><td>Atlas (11.9M cells, human + mouse)</td><td>No</td><td>Integrates massive single-cell datasets via gene-ranking similarity to enable scalable cross-dataset embedding and clustering</td><td>5</td></tr><tr><td>tGPT (Shen et al., 2023)</td><td>iScience</td><td>Foundation Model (Autoregressive)</td><td>Atlas (22.3 M cells)</td><td>No</td><td>Treats gene-expression ranks as token sequences and learns generative embeddings for cell clustering, trajectory inference, and bulk-tissue analysis</td><td>5</td></tr><tr><td>scBERT (Yang et al., 2022)</td><td>Nature Machine Intelligence</td><td>scRNA</td><td>Atlas (public scRNA-seq corpora)</td><td>No</td><td>Learns bidirectional gene-cell embeddings using BERT-style masked modeling for robust cell-type annotation and novel cell discovery</td><td>4</td></tr><tr><td>UCE (Universal Cell Embeddings) (Rosen et al., 2023)</td><td>BioRxiv</td><td>Foundation scRNA (cross-species)</td><td>Atlas / binary masked self-supervision</td><td>No</td><td>Embeds any cell zero-shot across species into a shared latent space for clustering, lineage inference, annotation</td><td>6</td></tr><tr><td>GeneCompass (Yang et al., 2024b)</td><td>Nature/Cell Research</td><td>scRNA / cross-species</td><td>Ontology (GRN + co-expression + gene family knowledge)</td><td>No</td><td>Learns cross-species gene-cell embeddings and supports annotation, perturbation, dose-response, and GRN inference</td><td>6</td></tr></table>

Continued on next page 

<table><tr><td>Model (Year)</td><td>Published Category where</td><td>Modality</td><td>Grounding Type</td><td>Agentic (Y/N)</td><td>Primary Task</td><td>Domain Score</td><td></td></tr><tr><td>scELMo (Liu et al., 2023)</td><td>BioRxiv</td><td>Text-Bridge LLM</td><td>scRNA + Metadata Text</td><td>LLM-generated embeddings from gene/metadata descriptions + raw data</td><td>No</td><td>Embeds cells via text-derived gene metadata embeddings combined with expression, then supports clustering, batch correction, annotation, perturbation</td><td>5</td></tr><tr><td>CellPLM (Wen et al., 2023)</td><td>ICLR</td><td>Spatial (Multi-modal Foundation Model)</td><td>scRNA + Spatial (SRT)</td><td>Atlas + Spatial relations (uses SRT in pretraining)</td><td>No</td><td>Treats cells as tokens and tissues as sentences, leveraging spatially-resolved transcriptomics and a Gaussian-mixture prior to encode inter-cell relations for denoising, spatial imputation, and perturbation prediction.</td><td>4</td></tr><tr><td>scMoFormer (Tang et al., 2023)</td><td>ArXiv</td><td>Spatial (Multi-modal Foundation Model)</td><td>scRNA + Protein / Multi-omics</td><td>Atlas + domain knowledge in cross-modality aggregation</td><td>No</td><td>Uses modality-specific transformers and cross-attention to impute missing modalities, classify cells, and fuse multimodal representations</td><td>4</td></tr><tr><td>scFormer (Xu et al., 2024)</td><td>ArXiv</td><td>Spatial (Multi-modal Foundation Model) (or Foundation + multi-omics)</td><td>Transcriptomics + Proteomics / multimodal</td><td>Atlas + transformer-based fusion</td><td>No</td><td>Aligns and integrates multi-omics single-cell data, recovers missing modalities, and transfers labels across modalities</td><td>4</td></tr><tr><td>scMulan (Bian et al., 2024a)</td><td>BioRxiv</td><td>Foundation Model</td><td>scRNA + Metadata</td><td>Atlas + prompt-conditioned generative modeling</td><td>No</td><td>Encodes each cell as a “c-sentence” integrating expression + metadata; supports zero-shot annotation, batch integration, and conditional generation</td><td>3</td></tr></table>


Continued on next page


<table><tr><td>Model (Year)</td><td>Published where</td><td>Category</td><td>Modality</td><td>Grounding Type</td><td>Agentic (Y/N)</td><td>Primary Task</td><td>Domain Score</td></tr><tr><td>scPRINT (Kalfon et al., 2025)</td><td>Nature Communications</td><td>Foundation Model</td><td>scRNA</td><td>Atlas (50M+cell pretraining)</td><td>No</td><td>Learns cell embeddings and infers cell-specific gene regulatory networks; supports zero-shot denoising, batch correction, label prediction, expression reconstruction</td><td>5</td></tr><tr><td>scGraphformer (Fan et al., 2024)</td><td>Nature Communications Biology</td><td>Foundation Model</td><td>scRNA</td><td>Atlas + relational prior (kNN bias, refined)</td><td>No</td><td>Learns a cell-cell graph via a transformer-GNN hybrid for better classification and interaction inference</td><td>4</td></tr><tr><td>scRDiT (Dong et al., 2025)</td><td>ArXiv</td><td>Foundation Model</td><td>scRNA (transcriptome)</td><td>Atlas-like generative prior / diffusion modeling</td><td>No</td><td>Generates synthetic scRNA-seq samples via diffusion transformer + DDIM for accelerated sampling</td><td>2</td></tr><tr><td>scGFT (Nouri, 2025)</td><td>Nature Communications Biology</td><td>Foundation Model</td><td>scRNA</td><td>Fourier-based perturbation / reconstruction (train-free)</td><td>No</td><td>Synthesizes new single-cell expression profiles by perturbing Fourier components in frequency space.</td><td>4</td></tr><tr><td>scTrans (Lu et al., 2024)</td><td>IJCAI</td><td>Foundation Model</td><td>scRNA</td><td>Sub-vector masked completion over gene modules</td><td>No</td><td>Learns multi-scale sub-vector tokens to perform gene-selective cell-type annotation via masked completion and contrastive regularization</td><td>3</td></tr><tr><td>scGT (Graph Transformer) (Qi et al., 2025b)</td><td>Bioinformatics Advance</td><td>Graph / Multi-omics Integration</td><td>scRNA + scATAC</td><td>Observed / Hybrid Graph</td><td>No</td><td>Multi-omics integration + label transfer</td><td>4</td></tr><tr><td>TransformerST (Lu et al., 2024)</td><td>Briefings in Bioinformatics</td><td>Spatial / Multi-modal FM</td><td>Spatial (histology + gene expression)</td><td>Spatial (histology + gene expression)</td><td>No</td><td>Super-resolution gene expression &amp; tissue clustering</td><td>3</td></tr></table>


Continued on next page


<table><tr><td>Model (Year)</td><td>Published Category where</td><td>Modality</td><td>Grounding Type</td><td>Agentic (Y/N)</td><td>Primary Task</td><td>Domain Score</td></tr><tr><td>scGPT-spatial (Cui et al., 2024)</td><td>BioRxiv</td><td>Spatial (Multi-modal Foundation Model)</td><td>Spatial Transcriptomics + scRNA prior</td><td>Atlas + spatial-aware decoding</td><td>No</td><td>Extends scGPT via continual pretraining to spatial data, supports multi-slide integration, cell-type deconvolution, and spatial gene imputation</td></tr><tr><td>HEIST (Madhu et al., 2025)</td><td>BioRxiv</td><td>Spatial (Multi-modal Foundation Model)</td><td>Spatial transcriptomics + proteomics</td><td>Hierarchical graph modeling (spatial + GRNs)</td><td>No</td><td>Learns joint embeddings of cells and genes in spatial context to perform cell annotation, gene imputation, spatial clustering, clinical outcome prediction</td></tr><tr><td>stFormer (Cao et al., 2024)</td><td>BioRxiv</td><td>Spatial (Multi-modal Foundation Model)</td><td>Spatial transcriptomics + ligand context</td><td>Atlas + biased cross-attention to ligand niche genes</td><td>No</td><td>Learns gene embeddings contextualized by ligand signals in spatial microenvironments; aids clustering, ligand-receptor inference, and perturbation simulation</td></tr><tr><td>FmH2ST (Wang et al., 2025)</td><td>Nucleic Acids Research</td><td>Spatial (Multi-modal Foundation Model)</td><td>Histology image + spatial transcriptomics</td><td>Image foundation + dual graphs + spot branch</td><td>No</td><td>Predict spatial gene expression from histology using fused image and spot features, supporting denoising, heterogeneity detection, and regulatory inference</td></tr><tr><td>OmiCLIP (Cui et al., 2025)</td><td>Nature Methods</td><td>Spatial (Multi-modal Foundation Model)</td><td>Histology + Spatial Transcriptomics</td><td>Image–gene contrastive alignment (rank-based)</td><td>No</td><td>Learns unified visual-omics embeddings linking histopathology and spatial gene expression; enables cross-modal prediction, annotation, and tissue retrieval</td></tr><tr><td>QuST-LLM (Huang, 2024)</td><td>ArXiv</td><td>Text-Bridge LLM</td><td>Spatial transcriptomics + histology metadata</td><td>GO-term + gene enrichment + LLM narrative overlay</td><td>Yes</td><td>Converts ST data and ROIs into human-readable narratives and matches natural language queries to spatial regions via GO/LLM interpretation</td></tr></table>


Continued on next page


<table><tr><td>Model (Year)</td><td>Published Category where</td><td>Modality</td><td>Grounding Type</td><td>Agentic (Y/N)</td><td>Primary Task</td><td>Domain Score</td></tr><tr><td>SpaCCC (Yang et al., 2024a)</td><td>IEEE Xplore</td><td>Text-Bridge LLM</td><td>Spatial transactions/ transcriptomics / transcriptomic genes (LRs)</td><td>LLM embeddings of ligand + receptor genes</td><td>No</td><td>Infers spatially resolved cell-cell communication by embedding LR pairs in LLM latent space + diffusion / permutation test filtering</td></tr><tr><td>spaLLM (Ji et al., 2024)</td><td>Briefings in Bioinformatics</td><td>Spatial (Multi-modal Foundation Model)</td><td>Spatial multi-omics (RNA, ATAC, protein)</td><td>scGPT embeddings + GNN + attention fusion</td><td>No</td><td>Enhances spatial domain identification by fusing LLM-derived embeddings and spatial-omics signals via multi-view attention</td></tr><tr><td>scMMGPT (2025) (Shi et al., 2025)</td><td>ArXiv</td><td>Text-Bridge LLM / Multi-modal_FM</td><td>RNA + Text (metadata)</td><td>Annotation, description</td><td>No</td><td>Generation - links single-cell and text PLMs to describe cells, generate pseudo-cells from text, and enhance annotation through text-conditioned reasoning</td></tr><tr><td>Cell2Sentence (C2S) (Levine et al., 2024)</td><td>PMLR</td><td>Text-Bridge LLM</td><td>RNA (scRNA-seq converted to “cell sentences”)</td><td>Marker / Implicit (gene rank-order)</td><td>No</td><td>Generation, annotation &amp; reconstruction - encodes cells as gene-ranked sentences and fine-tunes LLMs to classify or generate biologically meaningful cell text</td></tr><tr><td>C2S-Scale (Rizvi et al., 2025)</td><td>BioRxiv</td><td>Text-Bridge LLM+ Multi-modal Foundation Model</td><td>RNA + Text / Metadata</td><td>Atlas + Text / Implicit rank grounding</td><td>No</td><td>Chat, generation &amp; annotation - trains large LLMs on cell sentences and biological text to enable perturbation prediction and multicellular summarization</td></tr><tr><td>Cell2Text (Kharouiche et al., 2025)</td><td>ArXiv</td><td>Text-Bridge LLM+ Multi-modal Foundation Model</td><td>RNA → Natural Language</td><td>Gene-level embeddings + ontology / metadata grounding</td><td>No</td><td>Expression prediction &amp; regulatory inference - learns regulatory syntax from chromatin accessibility and DNA motifs to predict expression and interpret TF-cis interactions</td></tr></table>


Continued on next page


<table><tr><td>Model (Year)</td><td>Published where</td><td>Category</td><td>Modality</td><td>Grounding Type</td><td>Agentic (Y/N)</td><td>Primary Task</td><td>Domain Score</td></tr><tr><td>GenePT (Chen and Zou, 2024)</td><td>BioRxiv</td><td>Foundation Model / Text-Augmented</td><td>scRNA</td><td>Literature / Text embedding grounding</td><td>No</td><td>Embedding &amp; downstream prediction - uses GPT-3.5 gene text embeddings to derive cell embeddings via weighted or ranked gene aggregation.</td><td>5-7</td></tr><tr><td>CellLM (Zhao et al., 2023)</td><td>ArXiv</td><td>Foundation Model / Representation model</td><td>RNA (single cell expression)</td><td>Implicit - contrastive embedding from expression data</td><td>No</td><td>Represent cells via a contrastive-learning transformer, optimize embedding space for tasks like cell-type annotation, drug sensitivity prediction</td><td>6</td></tr><tr><td>scExtract (Wu and Tang, 2025)</td><td>Genome Biology</td><td>Agentic Framework (Text-Bridge LLM)</td><td>RNA (scRNA-seq) + Text</td><td>Article-based parameter extraction</td><td>Yes</td><td>Annotation &amp; Integration - uses LLMs to extract pipeline parameters from publications and apply them for dataset harmonization.</td><td>7</td></tr><tr><td>scAgent (Mao et al., 2025)</td><td>ArXiv</td><td>Agentic Framework (Text-Bridge hybrid)</td><td>RNA (scRNA-seq)</td><td>Reference atlas + marker gene reasoning + memory grounding</td><td>Yes</td><td>Universal cell annotation &amp; novel cell discovery - scAgent uses an LLM planning module, memory, and tool modules to annotate cells across tissues, detect unknown types, and incrementally learn new annotations</td><td>7</td></tr><tr><td>EpiFoundation (Wu et al., 2025)</td><td>BioRxiv</td><td>Epigenomic Foundation Model / Multi-modal alignment</td><td>scATAC (chromatin accessibility)</td><td>Peak-to-gene supervision (alignment to expression)</td><td>No</td><td>Cell embedding, annotation, batch correction, gene expression prediction - trains on sparse peak sets, aligns to gene expression supervision, and transfers learned embeddings for downstream ATAC tasks</td><td>7</td></tr></table>


Continued on next page 


<table><tr><td>Model (Year)</td><td>Published Category where</td><td>Modality</td><td>Grounding Type</td><td>Agentic (Y/N)</td><td>Primary Task</td><td>Domain Score</td><td></td></tr><tr><td>EpiAgent (Chen et al., 2025b)</td><td>BioRxiv</td><td>Epigenomic Foundation Model / Multi-modal text-bridge hybrid</td><td>scATAC / chromatin accessibility</td><td>Peak tokenization + external embeddings + regulatory supervision grounding</td><td>No</td><td>Embedding, annotation, imputation, and perturbation prediction - encodes chromatin accessibility as ranked cCRE tokens, enabling zero-shot cell annotation, peak imputation, and response prediction.</td><td>8</td></tr><tr><td>ChromFound (Jiao et al., 2025)</td><td>Nature</td><td>Epigenomic → Expression Transformer / Regulatory FM</td><td>scATAC (chromatin accessibility) + DNA sequence</td><td>Motif × peak matrix / masked regulatory grammar grounding</td><td>No</td><td>Expression prediction &amp; regulatory inference - learns regulatory syntax from chromatin accessibility and DNA motifs to predict gene expression in seen and unseen cell types; also interprets transcription factor interactions and cis-regulatory elements</td><td>8</td></tr><tr><td>GET (General Expression Transformer) (Fu et al., 2025)</td><td>Nature</td><td>Text-Bridge LLM (Hybrid Fusion Model)</td><td>RNA (scRNA) / Text embeddings</td><td>Fusion of scGPT embeddings + text-encoded “cell sentences” grounding</td><td>No</td><td>Annotation (fusion embedding model) - combines scGPT-derived embeddings and LLM (text encoder) embeddings via a small fusion MLP to improve cell-type classification robustness across datasets</td><td>6</td></tr><tr><td>scMPT (Palayew et al., 2025)</td><td>ArXiv</td><td>Text-Bridge LLM (Hybrid Fusion Model)</td><td>RNA (scRNA) / Text embeddings</td><td>Fusion of scGPT embeddings + text-encoded “cell sentences” grounding</td><td>No</td><td>Annotation (fusion embedding model) - combines scGPT-derived embeddings and LLM (text encoder) embeddings via a small fusion MLP to improve cell-type classification robustness across datasets</td><td>6</td></tr></table>


Continued on next page 


<table><tr><td>Model (Year)</td><td>Published Category where</td><td>Modality</td><td>Grounding Type</td><td>Agentic (Y/N)</td><td>Primary Task</td><td>Domain Score</td></tr><tr><td>EpiBERT (Javed et al., 2025)</td><td>Cell Genomics</td><td>Epigenomic Foundation / Multi-modal Transformer</td><td>DNA sequence + chromatin accessibility</td><td>Masked-accessibility pretraining + motif &amp; sequence fusion</td><td>No</td><td>Accessibility imputation, gene expression prediction &amp; regulatory inference - predicts masked ATAC signals, then fine-tunes to predict expression and enhancer-gene links, generalizing to unseen cell types</td></tr><tr><td>scMamba (Yuan et al., 2025)</td><td>ArXiv</td><td>Multimodal / Founda-tion_FM</td><td>Multi-omics (RNA + others)</td><td>Implicit via integrated features; no explicit external grounding</td><td>No</td><td>Multi-omic embedding, integration &amp; annotation across modalities without prior feature selection</td></tr><tr><td>GeneMamba (Qi et al., 2025a)</td><td>ArXiv</td><td>Foundation Model / State-Space Model</td><td>RNA (scRNA)</td><td>Implicit (via gene sequence context + pathway loss)</td><td>No</td><td>Multi-batch integration, cell-type annotation, gene correlation - uses a BiMamba state-space architecture with pathway-aware losses for scalable, context-rich modeling of scRNA</td></tr><tr><td>Nicheformer (Schaar et al., 2024)</td><td>BioRxiv</td><td>Foundation _FM/ Spatial / Multi-modal_FM</td><td>scRNA + spatial tran-scriptomics</td><td>Contextual tokenization + metadata + spatial neighbor-hood embedding</td><td>No</td><td>Spatial context prediction, spatial label / niche prediction, mapping spatial info to dissociated cells</td></tr><tr><td>scFormer Cell+ (Cui et al., 2022)</td><td>Bioarxiv</td><td>Foundation Model / Multi-modal Transformer</td><td>RNA (+ metadata)</td><td>Joint gene-cell embedding with metadata tokens</td><td>No</td><td>Integration &amp; Annotation - context-aware joint embedding for cross-species/tissue generalization</td></tr></table>


Continued on next page


<table><tr><td>Model (Year)</td><td>Published where</td><td>Category</td><td>Modality</td><td>Grounding Type</td><td>Agentic (Y/N)</td><td>Primary Task</td><td>Domain Score</td></tr><tr><td>scPlantLLM (Cao et al., 2025)</td><td>Genomics, Proteomics and Bioinformatics</td><td>Text-Bridge LLM (Foundation Model)</td><td>RNA (plant scRNA-seq)</td><td>Gene token + binned expression embedding (Plant Atlas grounding)</td><td>No</td><td>Annotation, Integration, GRN inference - pretrains on plant single-cell data with masked LM and cell-type supervision, enabling cross-species annotation, clustering, and regulatory discovery in plant systems</td><td>7</td></tr><tr><td>LICT (Ye et al., 2024)</td><td>ArXiv</td><td>Text-Bridge LLM / Annotation LLM Hybrid</td><td>RNA (scRNA-seq)</td><td>Marker gene-based DE lists + LLM prompting</td><td>No</td><td>Annotation &amp; Reliability - iterative LLM prompting with DE markers for label refinement and confidence scoring</td><td>5</td></tr><tr><td>Teddy (family of models) (Chevalier et al., 2025)</td><td>ArXiv</td><td>Foundation Model / Disease-aware Transformer</td><td>scRNA (single-cell RNA-seq)</td><td>Self-supervised + supervised annotation supervision</td><td>No</td><td>Disease state classification / healthy vs diseased detection - trained to identify disease conditions of held-out donors and distinguish diseased vs healthy cells in new disease contexts</td><td>7</td></tr><tr><td>Pilot (Joodaki et al., 2024)</td><td>Github</td><td>Benchmark / Evaluation Framework</td><td>RNA</td><td>Model-agnostic pilot foundation testing</td><td>No</td><td>Benchmarking &amp; Evaluation - lightweight pilot framework for early single-cell FM testing</td><td>4</td></tr><tr><td>CellVerse (Zhang et al., 2025b)</td><td>ArXiv</td><td>Benchmark / Evaluation</td><td>Multi-omics (RNA, CITE, ASAP, etc.)</td><td>Implicit via prompt encoding</td><td>No</td><td>QA benchmark for annotation, drug response, perturbation tasks</td><td>5</td></tr><tr><td>xTrimoGene (Gong et al., 2023)</td><td>NeurIPS / ArXiv</td><td>Foundation Model (Scalable Transformer)</td><td>scRNA-seq</td><td>Sparse masking + auto-discretization</td><td>No</td><td>Representation learning + annotation, perturbation, drug synergy prediction</td><td>7</td></tr><tr><td>EpiAttend (Li et al., 2022)</td><td>NeurIPS 2022 Work-shop</td><td>Regulatory / Sequence-Epigenome Transformer</td><td>DNA sequence + single-cell epigenomic data</td><td>Sequence + cell-specific epigenome grounding</td><td>No</td><td>Predict cell type-specific gene expression by integrating DNA sequence and single-cell epigenomic tracks, linking enhancers and promoters</td><td>6</td></tr><tr><td>Model (Year)</td><td colspan="2">Published Category where</td><td>Modality</td><td>Grounding Type</td><td>Agentic (Y/N)</td><td>Primary Task</td><td>Domain Score</td></tr><tr><td>Spatial2-Sentence (Chen et al., 2025a)</td><td>ArXiv</td><td>Spatial / Text-Bridge hybrid</td><td>Spatial + expression (Imaging Mass Cytometry)</td><td>Spatial adjacency + expression similarity tokenization</td><td>No</td><td>Encode spatial &amp; expression context into multi-sentence prompts for LLMs to perform cell-type classification and clinical status prediction</td><td>6</td></tr><tr><td>ChatCell (Fang et al., 2024)</td><td>ArXiv</td><td>Text-Bridge LLM (Instructional)</td><td>scRNA → “cell sentence”</td><td>Vocabulary adaptation + unified sequence generation of cell sentences</td><td>No</td><td>Natural language interface for single-cell tasks - allows users to query, annotate, generate, and explore scRNA data via text prompts. Hugging Face</td><td>6</td></tr><tr><td>CellAtria (Nouri et al., 2025)</td><td>BioRxiv</td><td>Agentic Framework</td><td>scRNA- seq + metadata</td><td>Ontology (graph + metadata)</td><td>Yes</td><td>Annotation &amp; Ontology Mapping</td><td>7</td></tr><tr><td>CellAgent (Xiao et al., 2024)</td><td>ArXiv</td><td>Agentic Framework</td><td>scRNA- seq</td><td>None</td><td>Yes</td><td>Annotation &amp; Ontology Mapping</td><td>7</td></tr></table>

F Appendix Figure 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-28/1264e285-66b1-43c6-8a68-5faf4ab9f487/9ddbe3f5054dbdea51042162a737290df5f500873e058d0b32d5d3c02c476fad.jpg)



Covered


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-28/1264e285-66b1-43c6-8a68-5faf4ab9f487/a452cd638bfaafc2f4be5df8cf44f5b2de56fe4968097eb9a5ac88c3ef0e1dae.jpg)



Not Covered



Model $\times$ Task Coverage (Yellow $=$ Covered, Dark Blue $=$ Not Covered)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-28/1264e285-66b1-43c6-8a68-5faf4ab9f487/8decbd94736872df818096362b479d0a0df728baf46de4e7d11d6992fdf14f4f.jpg)



Figure 3: Task vs Model Heatmap


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-28/1264e285-66b1-43c6-8a68-5faf4ab9f487/075140f3e241e9bf61afbb5c40de17a64dfad3330ef20d6a273bc880bb007ad0.jpg)



Figure 4: Comparison of task coverage between agentic (dark blue) and non-agentic (yellow) models. Agentic frameworks emphasize annotation, ontology mapping, spatial mapping while non-agentic models concentrate on trajectory, perturbation modeling, regulatory and pathway inference as well.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-28/1264e285-66b1-43c6-8a68-5faf4ab9f487/148d3716b415972068dcdaa5fef40a48143da01ba507c45beb193e9cb7856393.jpg)



Figure 5: Comparison of domain coverage between agentic (dark blue) and non-agentic (yellow) models. Agentic frameworks emphasize explainability, fairness, and emerging paradigms, while non-agentic models concentrate on biological grounding and batch effects.

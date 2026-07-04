# Protocol-Analysis Literature Matrix for CAAM/scMAE

Date: 2026-06-26

Status: writing and research-planning memo. This is not experimental evidence.

Current manuscript route:

```text
protocol_analysis / diagnostic paper
```

Current forbidden route:

```text
positive CAAM-scMAE method paper centered on AdvMask, Axial, or their synergy
```

This matrix is organized by how each paper helps or hurts the current manuscript claim. It should be used to write Introduction, Related Work, Discussion, and reviewer-response material.

## 1. Working Claim

Safe current claim:

> In scRNA-seq masked autoencoding for clustering, feature-space choice, corruption semantics, effective corruption diagnostics, and mask-selection policy can materially change development conclusions. More complex corruption or learned masking mechanisms can be active by diagnostics while still failing to improve clustering representations beyond seed-level variation.

Unsafe current claims:

```text
CAAM-scMAE improves over scMAE.
AdvMask improves clustering.
Current Axial improves clustering.
Axial + AdvMask has positive synergy.
Known-K ARI is fully unsupervised evidence.
Development triage validates the final method.
```

## 2. Direct Single-Cell Masked/Clustering Evidence

| Paper | Problem | Core method | Relation to CAAM/scMAE | Supports or weakens our hypothesis | Borrowable experiment design | Main risk exposed | Manuscript use | Gap sentence |
|---|---|---|---|---|---|---|---|---|
| scMAE: a masked autoencoder for single-cell RNA-seq clustering | Learn clustering-oriented scRNA-seq representations with masked reconstruction | Gene-wise shuffle corruption, reconstruction, and mask prediction | Direct base model and strongest internal baseline | Supports masked autoencoding as a useful starting point; weakens any claim that a larger module is automatically needed | Keep scMAE-style shuffle and mask-prediction diagnostics; compare all variants against random-mask MLP/scMAE-compatible control | If we do not beat or clarify scMAE protocol behavior, the paper becomes a weak variant study | Introduction, Method, Results | Although scMAE shows that masked corruption can support clustering, it remains unclear which protocol details make the pretext task informative under sparse HVG inputs. |
| scCluBench | Benchmark clustering algorithms across many scRNA-seq datasets | Unified benchmark with multiple datasets, algorithms, metrics, and failure modes | Motivates protocol fairness and broad validation | Supports our shift from single-dataset method claim to benchmark/protocol analysis | Separate development, validation, and sealed test; report known-K and non-oracle metrics distinctly | Reviewers will ask for broader validation before publication-level claims | Introduction, Experimental Setup, Limitations | Existing benchmarks rank algorithms, but do not isolate how masked-autoencoder corruption semantics and effective corruption diagnostics affect representation quality. |
| scNAME | Single-cell clustering with neighborhood contrastive learning and ancillary mask estimation | Neighborhood contrastive signal plus ancillary mask-related objective | Closest neighbor/mask-adjacent baseline | Supports the idea that mask-related objectives can matter; weakens a naive claim that our mask predictor is novel by itself | Include mask objectives and neighborhood baselines in Related Work; do not claim mask estimation alone as novelty | If we ignore scNAME, reviewers may say the mask idea is already covered | Related Work | Prior masked or ancillary objectives motivate mask-aware learning, but they do not resolve whether learned mask selection improves clustering under controlled corruption budgets. |
| scDeepCluster | Deep clustering for scRNA-seq | Autoencoder plus model-based clustering | Standard deep clustering baseline | Supports need for strong baselines; weakens paper if we only compare internal variants | Report against representative deep baselines after validation approval | Full formal baseline grid is not yet run | Related Work, Validation Plan | Deep clustering methods demonstrate learned embeddings are useful, but they often entangle representation learning with clustering assumptions rather than isolating corruption protocol effects. |
| DESC | Deep clustering with batch-effect removal | Autoencoder-based clustering and iterative self-training/batch correction | Baseline and batch-sensitivity reference | Supports discussing batch/fairness risks | Include batch-aware evaluation only after final protocol is frozen | If CAAM changes preprocessing, fairness becomes suspect | Related Work, Discussion | Batch-aware deep clustering highlights that preprocessing and objective design can change downstream conclusions, motivating protocol-level controls. |
| scDCC | Constrained clustering for scRNA-seq | Deep embedding with pairwise constraints | Baseline showing label/constraint leakage risk | Weakens any overly broad "unsupervised" claim if known labels or constraints leak | Use as contrast: we do not use labels/constraints in training | Known-K evaluation still needs careful wording | Related Work, Reproducibility | Constrained methods can improve clustering with extra supervision, so a label-free masked protocol must state exactly where labels are excluded. |
| scGNN / scDSC / scCDCG | Graph-based single-cell clustering | Cell graph, GNN, structural or cut-informed embedding | Baselines stronger than simple MLP on some benchmarks | Weakens a pure MLP-method novelty claim; supports diagnostic route instead | Compare final protocol to graph baselines only under validation/formal approval | Graph construction can dominate performance; not fair to hand-wave | Related Work, Limitations | Graph clustering methods show that population structure matters, but also make it harder to separate representation gains from graph-construction choices. |
| SC3, Louvain, Leiden, Seurat-style graph pipelines | Classical and graph-community clustering | Consensus or graph-community clustering | Required non-deep reference points | Supports reporting fixed-Leiden non-oracle summaries | Always distinguish known-K K-means from fixed-resolution Leiden | Known-K ARI alone is oracle-like | Experimental Setup, Results | A protocol-analysis paper should report both oracle-like development diagnostics and non-oracle graph clustering views. |

## 3. Masked Modeling and Corruption Background

| Paper | Problem | Core method | Relation to CAAM/scMAE | Supports or weakens our hypothesis | Borrowable experiment design | Main risk exposed | Manuscript use | Gap sentence |
|---|---|---|---|---|---|---|---|---|
| BERT | Self-supervised representation learning from masked tokens | Predict masked tokens from context | General masked modeling background | Supports masked prediction as representation learning; weakens any claim that masking itself is novel | Explain masked pretext-task design without overclaiming novelty | Token masking analogy can be misleading for sparse continuous gene expression | Related Work | Unlike language tokens, scRNA-seq entries are sparse continuous measurements, so zero-to-zero and effective corruption must be diagnosed explicitly. |
| Masked Autoencoders Are Scalable Vision Learners | Masked reconstruction for visual representation | High-ratio masking and reconstruction from visible patches | General MAE reference | Supports reconstruction-based pretraining; weakens architecture novelty unless adapted to scRNA-seq sparsity | Use reconstruction pretext-task framing and ablation discipline | Vision MAE assumptions do not transfer directly to zero-inflated matrices | Related Work | The success of MAE-style objectives in dense modalities does not answer which corruption semantics are informative for sparse gene-expression matrices. |
| Denoising / dropout-aware single-cell autoencoders such as DCA | Denoise sparse scRNA-seq counts | ZINB/dropout-aware reconstruction | Background for sparsity and zero inflation | Supports our emphasis on zero semantics | Discuss why zero-to-zero corruption is diagnostic, not necessarily training failure | Reviewers may ask why not use ZINB loss | Discussion, Future Work | Denoising models address count noise, but clustering-oriented masked corruption requires separate diagnostics for whether a masked entry actually changes information. |
| scVI / scvi-tools ecosystem | Probabilistic representation for single-cell data | VAE and likelihood-based modeling | Strong alternative representation-learning family | Weakens a narrow autoencoder-only framing | Position our paper as protocol analysis, not universal representation replacement | Need not beat scVI unless final method claim is reopened | Related Work | Probabilistic single-cell models handle uncertainty, while our current analysis asks a narrower question about masked corruption protocol effects. |

## 4. Foundation Models and Context Modeling

| Paper | Problem | Core method | Relation to CAAM/scMAE | Supports or weakens our hypothesis | Borrowable experiment design | Main risk exposed | Manuscript use | Gap sentence |
|---|---|---|---|---|---|---|---|---|
| TabPFN | Fast prediction on small tabular datasets | Foundation model over tabular tasks with row/feature context | Original inspiration for context/attention idea | Supports exploring row/feature context; weakens direct transfer claims because TabPFN is supervised tabular prediction, not scRNA-seq clustering | Parameter-matched attention-vs-MLP comparisons; task-size and feature-count constraints | Current Axial smoke failed; do not claim TabPFN-like modeling is validated | Introduction, Future Work | Tabular foundation models motivate context-aware designs, but scRNA-seq clustering needs controlled tests showing context modeling improves representation rather than parameter count. |
| Geneformer | Transfer learning for network biology | Gene-rank/token foundation model | Large-scale single-cell representation reference | Supports broader foundation-model context; weakens small-model novelty if we claim general foundation status | Treat as background, not direct baseline unless resources allow | Foundation models may not improve clustering without task-specific evaluation | Related Work | Large single-cell foundation models provide transferable representations, but they do not replace task-specific evidence for clustering protocols. |
| scGPT | Single-cell multi-omics foundation model | Generative pretraining on cell/gene tokens | Foundation model reference | Supports the importance of gene/cell tokenization and context | Discuss as future baseline or representation comparator | We cannot imply we are building comparable-scale model | Related Work, Future Work | Foundation models shift the field toward contextual representation learning, but smaller masked autoencoders still need rigorous protocol controls. |
| scFoundation | Large-scale foundation model on transcriptomics | Large model trained on broad transcriptomic data | Foundation model reference | Supports evaluating foundation baselines later; weakens "top-tier model innovation" unless our mechanism is sharper | Add as validation/future baseline, not current development evidence | Without validation, comparisons would be speculative | Related Work | Transferable transcriptomic foundation models raise the bar for new methods, making protocol-sensitive evidence more important before proposing new modules. |
| Single-cell foundation model evaluations/surveys | Evaluate utility of foundation models across tasks | Survey or benchmark foundation representations | Helps position route as evaluation-sensitive | Supports our argument that representation claims must be task-tested | Use as Discussion support for not assuming attention/foundation gains | Need up-to-date citation cleanup before submission | Discussion | Even foundation-model papers require task-specific evaluation; our negative CAAM evidence follows the same discipline at smaller scale. |

## 5. Learned/Adversarial Masking and Negative Evidence

| Paper | Problem | Core method | Relation to CAAM/scMAE | Supports or weakens our hypothesis | Borrowable experiment design | Main risk exposed | Manuscript use | Gap sentence |
|---|---|---|---|---|---|---|---|---|
| Generative Adversarial Nets | Learn generative models by adversarial training | Generator-discriminator minimax | Historical background only | Supports adversarial framing in broad sense; weakens calling AdvMask a full GAN because our selector does not generate expression values | Use only as high-level background if needed | "GAN" wording invites mode-collapse and training-stability criticism | Related Work or omit | Our mask selector is better described as constrained adversarial mask selection, not a full generative model. |
| WGAN-GP / adversarial stabilization | Stabilize adversarial objectives | Gradient penalty and training constraints | Background for why adversarial modules need constraints | Supports our constrained-budget design | Report generator gradients, entropy, Gini, top-gene concentration | Nonzero gradients do not imply useful representation gains | Methods, Discussion | Adversarial training diagnostics are necessary for implementation validity, but not sufficient for downstream clustering utility. |
| Hard example mining / learned masking / Gumbel top-k family | Choose difficult or informative training examples/features | Learned selection under budget or relaxation | Conceptual neighbor of AdvMask | Supports the hypothesis that mask choice can matter; current results weaken it under our protocol | Require random-mask control, equal budget, and no label leakage | Learned selectors can find shortcuts or overfit difficulty rather than biology | Future Work | A learned mask can optimize difficulty without improving clustering, so downstream representation gates are essential. |

## 6. Current Evidence-to-Literature Alignment

| Our evidence | Literature alignment | Interpretation |
|---|---|---|
| Phase12 corrected HVG/log1p/no-scale protocol | scCluBench, Seurat/scRNA-seq preprocessing norms | Feature-space choice is not a cosmetic setting; it changes benchmark fairness. |
| Phase13 corruption triad chose scMAE shuffle over matched/nonzero-aware donor | scMAE plus masked modeling background | The original scMAE corruption remains a strong baseline; "more biologically plausible" replacement is not automatically better. |
| nonzero-aware donor improved diagnostics but not stable clustering | Denoising/zero-inflation literature plus our diagnostics | Effective corruption is necessary to measure, but insufficient as a quality proxy. |
| Phase14 AdvMask generator gradients positive but effect-size gate failed | adversarial training / learned masking background | Active adversarial training does not equal useful representation learning. |
| Axial smoke lost to parameter-matched MLP | TabPFN and foundation-model context literature | Context modeling remains a plausible future direction, but current Axial implementation is not evidence for it. |
| post-hoc label F1/flow diagnostics match ARI direction | single-cell biological interpretation expectations | Label recovery diagnostics can strengthen analysis, but no marker/cell-type claim is yet supported. |

## 7. Related Work Paragraph Skeleton

Masked autoencoding has become a promising self-supervised objective for single-cell clustering, with scMAE showing that gene-wise corruption and mask prediction can yield useful clustering representations. However, broader clustering benchmarks such as scCluBench and long-standing concerns about single-cell clustering evaluation show that performance claims depend heavily on preprocessing, metric choice, and benchmark protocol. Existing deep clustering methods, including scDeepCluster, DESC, scDCC, scGNN, scDSC, scCDCG, and scNAME, demonstrate that learned embeddings, graph structure, neighborhood information, and mask-related objectives can all affect clustering. In contrast to proposing another larger architecture by default, this study isolates protocol choices inside masked autoencoding: feature space, corruption semantics, effective corruption diagnostics, and learned mask selection. Inspired by masked modeling and tabular/context foundation models, we test whether more complex corruption and learned masking provide clustering gains, but our development evidence shows that diagnostic activity can diverge from representation quality.

## 8. Reviewer-Facing Positioning

If a reviewer says "this is not a new method":

```text
Correct. The current manuscript should not be sold as a new state-of-the-art model. Its contribution is a controlled protocol analysis showing that masked-autoencoder improvements can fail at the corruption/mask-selection level even when implementation diagnostics are active.
```

If a reviewer says "why not tune AdvMask/Axial more":

```text
Tuning after seeing development failures would contaminate the validation route. The disciplined response is to freeze the protocol-analysis route, validate it, and only reopen attention or learned masking as a new mechanism with its own development split and parameter-matched controls.
```

If a reviewer says "known-K is oracle":

```text
Known-K metrics are explicitly development diagnostics. The manuscript must report fixed-Leiden/non-oracle summaries and post-hoc label diagnostics separately, and it must not describe known-K ARI as fully unsupervised performance.
```

## 9. Source Links to Verify Before Submission

Primary or publisher/DOI links already represented in the current BibTeX or local paper workspace:

- scMAE: https://doi.org/10.1093/bioinformatics/btae020
- scCluBench: https://doi.org/10.1609/aaai.v40i2.37110
- scNAME: https://doi.org/10.1093/bioinformatics/btac011
- scDeepCluster: https://doi.org/10.1038/s42256-019-0037-0
- scDCC: https://doi.org/10.1038/s41467-021-22008-3
- DESC: https://doi.org/10.1038/s41467-020-15851-3
- scGNN: https://doi.org/10.1038/s41467-021-22197-x
- scDSC: https://doi.org/10.1093/bib/bbac018
- TabPFN: https://www.nature.com/articles/s41586-024-08328-6
- Geneformer: https://doi.org/10.1038/s41586-023-06139-9
- scGPT: https://doi.org/10.1038/s41592-024-02201-0
- scFoundation: https://doi.org/10.1038/s41592-024-02305-7

Zotero local API status on this run:

```text
api_running = false
base_url = http://127.0.0.1:23119
```

No Zotero library writes or imports were performed.

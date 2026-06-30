# scMAEs phase 1 status

This branch is moving from the legacy shared lightweight variant framework to
independent method folders. The legacy `common/model.py + variant_defs.py`
implementation must not be reported as a valid new result.

Implemented independent folders:

- Rank 01: `rank01_scdiva_discrete_diffusion_full`.
  - Time-conditioned masked discrete diffusion model over expression tokens.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Synthetic forward/backward smoke: PASS.
  - Small Melanoma real-data smoke: PASS (`n_top_genes=64`, `depth=1`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 32: ACC `0.657212`, NMI `0.604449`, ARI `0.530250`, F1 macro `0.503932`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.
  - GitHub clone status: `https://github.com/SindiLab/ScDiVa` anonymous clone failed; implementation uses local PDF/report.

- Rank 02: `rank02_maskfeat_gene_features_full`.
  - MaskFeat-style gene patch masking with learned mask token and masked patch feature prediction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Synthetic forward/backward smoke: PASS.
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `depth=1`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, `n_top_genes=1000`: ACC `0.403723`, NMI `0.367150`, ARI `0.251842`, F1 macro `0.341047`. This is far below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.
  - Source clone: `external_sources/rank02_slowfast_maskfeat`, commit `287ec0076846560f44a9327e931a5a2360240533`.

- Rank 03: `rank03_consistency_models_full`.
  - Consistency-model adaptation with Karras sigma schedule, boundary-condition denoising scalings, adjacent-time EMA target consistency, Pseudo-Huber robust regression, and mask prediction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Synthetic forward/backward smoke: PASS.
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `depth=2`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.618879`, NMI `0.653005`, ARI `0.525432`, F1 macro `0.506649`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.
  - Source clone: `external_sources/rank03_openai_consistency_models`, commit `e32b69ee436d518377db86fb2127a3972d0d8716`.

- Rank 04: `rank04_joao_graphcl_full`.
  - JOAO GraphCL adaptation with batch-local KNN cell graph, GCN encoder, five graph/expression augmentations, learned augmentation probability simplex update, NT-Xent loss, and masked expression reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Synthetic forward/backward smoke: PASS.
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `depth=2`, `epochs=1`, `skip_eval=true`), with `aug_P` updated from uniform to a non-uniform policy.
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.612896`, NMI `0.581485`, ARI `0.473741`, F1 macro `0.481978`. Final `aug_P` concentrated on cell feature drop/subgraph keep, but the result is below the current scMAE Melanoma reference ARI of about `0.668`.
  - Source clone: `external_sources/rank04_graphcl_automated`, commit `8f3c2ac7831b88693e932c924428d0c3fe065894`.

- Rank 05: `rank05_tabr_retrieval_full`.
  - TabR retrieval adaptation with query/candidate key encoder, exact nearest-neighbor context search, self-neighbor removal, TabR-style context value, attention residual injection, and masked expression reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Synthetic forward/backward smoke: PASS.
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `context_size=8`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.699756`, NMI `0.676540`, ARI `0.569282`, F1 macro `0.633085`. This is the strongest quick-screen result so far, but still below the current scMAE Melanoma reference ARI of about `0.668`.
  - Source clone: `external_sources/rank05_tabr`, commit `17baa9082506f8e7a0f8d11bb1e08212926a1507`.

- Rank 06: `rank06_scvgae_zinb_graph_full`.
  - scVGAE adaptation with batch-local KNN cell graph, GCN encoder, variational latent mean/log-variance, ZINB mean/dropout/dispersion graph heads, KL loss, and masked expression reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Synthetic forward/backward smoke: PASS.
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `latent_size=16`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.582539`, NMI `0.580844`, ARI `0.425145`, F1 macro `0.456258`. This is below the current scMAE Melanoma reference ARI of about `0.668`.
  - GitHub clone status: `https://github.com/STOmics/scVGAE` anonymous clone failed; implementation uses local PDF/report.

Skipped before next implementation:

- Rank 07 `DinoBloom`: hematology image foundation model; no fair no-large-external-weight scRNA reproduction.
- Rank 08 `scCello`: ontology/foundation supervision would alter unsupervised clustering fairness.
- Rank 09 `LangCell`: language-cell foundation pretraining requires external setup/weights outside this benchmark scope.
- Rank 12 `A Survey on Foundation Language Models for Single-cell Biology`: survey/related-work source, not an implementable model variant.

- Rank 10: `rank10_celler_longtail_full`.
  - Celler long-tail adaptation with Gaussian Inflation Loss over unsupervised pseudo prototypes, KMeans pseudo-cluster counts, rare pseudo-cluster reweighting, hard sample mining, and masked expression reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Synthetic forward/backward smoke: PASS.
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `depth=1`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.628850`, NMI `0.688649`, ARI `0.564727`, F1 macro `0.515831`. This is close to Rank 05 but still below the current scMAE Melanoma reference ARI of about `0.668`.
  - Source clone: `external_sources/rank10_hiceller`, commit `6b27dee9cf55c3a16a84a7a4f649647a2ab8dca2`.

- Rank 11: `rank11_scmamba_ssm_full`.
  - scMamba adaptation with patch tokenization, learned mask token, explicit selective state-space scan over gene patches, learned negative state dynamics `A`, input/output selectors `B/C`, skip `D`, positive per-token `delta`, causal depthwise convolution, residual SSM blocks, and masked patch reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: `external_sources/rank11_scmamba`, commit `4887c0a8ab060b2482384d2294fe265b633d2406`.
  - Official `mamba_ssm`/`causal_conv1d` dependencies are unavailable in the environment; this implementation uses an explicit dependency-free selective scan and does not use a `Linear + Sigmoid` substitute.
  - Synthetic forward/backward smoke: PASS (`embedding=(8,32)`, finite loss, nonzero `A_log` gradient).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `state_size=8`, `depth=1`, `patch_size=16`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.481055`, NMI `0.414006`, ARI `0.283059`, F1 macro `0.389079`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 13: `rank13_masked_sc_clustering_full`.
  - Mask-sc adaptation with 2D expression matrix patch embeddings, fixed-count random patch masking, visible-only Transformer encoder, sequence-guided Transformer decoder, learned mask tokens, and reconstruction of sequence-level target features.
  - The frozen sequence-level target encoder is trained inside `run.py` with a contrastive-sc style self-supervised objective before mask-sc training, because no GitHub URL or checkpoint is listed in the local index.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Synthetic forward/backward smoke: PASS (`embedding=(8,32)`, masked prediction shape equals target feature shape, finite loss, nonzero 2D patch projection gradient).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `decoder_size=32`, `target_dim=16`, `encoder_depth=1`, `decoder_depth=1`, `target_pretrain_epochs=1`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, target sequence pretraining 5 epochs plus mask-sc 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.596277`, NMI `0.622791`, ARI `0.496678`, F1 macro `0.453728`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 14: `rank14_cicl_iter_contrast_full`.
  - CICL adaptation with gene-patch Transformer encoder, two Gaussian/dropout augmented views, iterative KMeans pseudo labels, Student-t clustering probabilities, projection head, instance InfoNCE, cluster-aware contrastive loss, DEC-style target-distribution KL stabilization, and auxiliary scMAE masked reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: `external_sources/rank14_circle`, commit `aba15bc81dc0b4999f56c7a82dd5f13bf109f27c`.
  - Synthetic forward/backward smoke: PASS (`embedding=(8,32)`, finite combined CICL loss, nonzero patch Transformer projection gradient).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `projection_size=16`, `depth=1`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.510082`, NMI `0.502244`, ARI `0.370997`, F1 macro `0.373784`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 15: `rank15_scagc_adaptive_graph_full`.
  - scAGC adaptation with mini-batch KNN initial graph, RBF similarity, Gumbel-TopK straight-through adaptive graph sampling, symmetric adjacency, TAGCN 0..K-hop graph encoder, inner-product graph decoder, ZINB decoder, temporal graph contrastive guidance, Student-t KL clustering objective, and auxiliary scMAE masked reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - GitHub URL: none listed in `02_整理索引.csv`; implementation follows the local PDF/report.
  - Synthetic forward/backward smoke: PASS (`embedding=(10,16)`, adaptive adjacency `(10,10)`, finite combined scAGC loss, nonzero TAGCN projection gradient).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `latent_size=16`, `hop_order=1`, `graph_top_k=8`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.615333`, NMI `0.646286`, ARI `0.527687`, F1 macro `0.480268`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

Additional skipped items:

- Rank 16 `Computational Methods for Single-Cell Multi-Omics Integration and Alignment`: review/background source.
- Rank 17 `Interpretable Deep Learning in Single-Cell Omics`: interpretability survey/direction source, not a standalone model variant.

- Rank 18: `rank18_beit_tokenizer_full`.
  - BEiT adaptation with expression patch embeddings, learned `[M]` mask token, Transformer encoder, dataset-local quantile expression tokenizer, masked patch token cross-entropy, and auxiliary masked expression reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: sparse `external_sources/rank18_unilm_beit`, commit `833df7e7832e5064a281131ee64a481afa8e5b95`.
  - Synthetic forward/backward smoke: PASS (`embedding=(8,32)`, token logits `(8,8,16)`, finite masked token loss, nonzero patch embedding gradient).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `depth=1`, `patch_size=8`, `vocab_size=16`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.551518`, NMI `0.545871`, ARI `0.408118`, F1 macro `0.443030`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 19: `rank19_data2vec_ema_full`.
  - data2vec adaptation with EMA teacher/student encoders, full-input teacher targets, masked-input student prediction, top-K layer-normalized teacher target averaging, SmoothL1 latent regression only at masked patch positions, EMA decay schedule, and auxiliary masked expression reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: sparse `external_sources/rank19_fairseq_data2vec`, commit `3d262bb25690e4eb2e7d3c1309b1e9c406ca4b99`.
  - Synthetic forward/backward smoke: PASS (`embedding=(8,32)`, target `(8,8,32)`, finite latent loss, nonzero student patch projection gradient, no teacher gradients, EMA update executes).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `depth=2`, `patch_size=8`, `average_top_k_layers=2`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.412807`, NMI `0.256247`, ARI `0.161320`, F1 macro `0.374823`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 20: `rank20_multimae_targets_full`.
  - MultiMAE adaptation with RNA-internal `expr`, `rank`, and `stat` tasks, separate task input adapters, Dirichlet visible-token sampling, a shared visible-token Transformer encoder, task-specific cross-attention decoders, and masked-token-only SmoothL1 losses.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: sparse `external_sources/rank20_multimae`, commit `66910f5b5ba236f5e731883db85fe4f24ee01106`.
  - Synthetic forward/backward smoke: PASS (`embedding=(8,32)`, finite loss, nonzero input-adapter, decoder, and encoder gradients).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `depth=1`, `decoder_size=32`, `patch_size=8`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.379127`, NMI `0.334000`, ARI `0.215508`, F1 macro `0.268009`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 21: `rank21_audiomae_local_window_full`.
  - AudioMAE adaptation with expression vectors padded to a 2D gene grid, Conv2d patch embedding, AudioMAE-style 2D row/column masking, visible-patch MAE encoder, decoder mask-token unshuffle, shifted local-window attention decoder, and masked-patch-only reconstruction loss that excludes padded gene positions from the denominator.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: sparse `external_sources/rank21_audiomae`, commit `bd60e29651285f80d32a6405082835ad26e6f19f`.
  - Synthetic forward/backward smoke: PASS (`embedding=(8,32)`, prediction shape valid, finite loss, nonzero patch-embed, local-window decoder, and encoder gradients).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `depth=1`, `decoder_size=32`, `decoder_depth=2`, `patch_size=4`, `window_size=2`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 25 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.357855`, NMI `0.325007`, ARI `0.195674`, F1 macro `0.348075`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 22: `rank22_ijepa_gene_context_full`.
  - I-JEPA adaptation with a gene-patch context encoder, EMA target encoder, non-overlapping context/target masks, predictor target tokens with positional embeddings, target-block SmoothL1 latent prediction, and auxiliary scMAE masked reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: `external_sources/rank22_ijepa`, commit `52c1ae95d05f743e000e8f10a1f3a79b10cff048`.
  - Synthetic forward/backward smoke: PASS (`embedding=(8,32)`, prediction `(8,8,32)`, finite loss, nonzero context encoder and predictor gradients).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `depth=1`, `predictor_size=32`, `predictor_depth=1`, `patch_size=8`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 80 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.587857`, NMI `0.596559`, ARI `0.454131`, F1 macro `0.493923`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 23: `rank23_maskgit_iterative_full`.
  - MaskGIT adaptation with dataset-local quantile expression tokens, learned `[MASK]` token, BERT-style bidirectional token Transformer, tied MLM token logits, cosine high-mask training schedule, confidence-based iterative decode, and auxiliary masked patch reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: `external_sources/rank23_maskgit`, commit `1db23594e1bd328ee78eadcd148a19281cd0f5b8`.
  - Synthetic forward/backward smoke: PASS (`embedding=(16,32)`, logits `(16,8,16)`, finite loss, nonzero token embedding gradient, iterative decode valid).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `depth=1`, `patch_size=8`, `vocab_size=16`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 80 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.243297`, NMI `0.148566`, ARI `0.062727`, F1 macro `0.211855`. This is far below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

Additional skipped items:

- Rank 24 `Flow Matching for Generative Modeling`: full flow-matching training would add a heavy continuous generative objective with weak direct comparability to fixed-k scRNA clustering; skipped as a primary independent scMAE benchmark variant.
- Rank 25 `Anomaly Transformer`: anomaly-association priors are aimed at anomaly detection rather than representation learning for fixed-k cell clustering; skipped as a primary independent scMAE benchmark variant.

- Rank 26: `rank26_bgrl_graph_bootstrap_full`.
  - BGRL adaptation with batch-local KNN cell graphs, two graph/expression views, edge dropout, feature masking, online GraphSAGE-style encoder, independently initialized EMA target encoder, MLP predictor, symmetric stopped-gradient cosine bootstrap loss, and auxiliary masked-expression reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: `external_sources/rank26_bgrl`, commit `60f9f19ad0598f9163ad70ebbde3e7297760901e`.
  - Synthetic forward/backward smoke: PASS (`embedding=(12,16)`, finite loss, nonzero online encoder and predictor gradients, no target encoder gradients).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `latent_size=16`, `predictor_hidden=32`, `graph_top_k=8`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 80 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.608243`, NMI `0.637522`, ARI `0.518232`, F1 macro `0.464919`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 27: `rank27_graph_barlow_twins_full`.
  - Graph Barlow Twins adaptation with batch-local KNN cell graphs, two edge/feature dropout views, shared dense graph encoder, projection head, batch-normalized cross-correlation matrix, diagonal invariance loss, off-diagonal redundancy-reduction loss, and auxiliary masked-expression reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: `external_sources/rank27_graph_barlow_twins`, commit `ec62580aa89bf3f0d20c92e7549031deedc105ab`.
  - Synthetic forward/backward smoke: PASS (`embedding=(12,16)`, projection `(12,16)`, finite loss, nonzero encoder and projector gradients).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `latent_size=16`, `projection_size=16`, `projector_hidden=32`, `graph_top_k=8`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 80 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.625526`, NMI `0.558700`, ARI `0.487346`, F1 macro `0.440827`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 28: `rank28_graphormer_bias_full`.
  - Graphormer adaptation with batch-local KNN cell graphs, graph token, degree centrality embeddings, shortest-path spatial attention bias, edge-type attention bias, graph-token virtual-distance bias, per-head attention bias added before softmax, graph reconstruction, and auxiliary masked-expression reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: `external_sources/rank28_graphormer`, commit `a04573c40705fb174db261bb746a8258d00992f5`.
  - Synthetic forward/backward smoke: PASS (`embedding=(12,32)`, graph embedding `(32,)`, finite loss, nonzero attention projection and spatial-bias embedding gradients).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `depth=1`, `num_heads=4`, `graph_top_k=8`, `max_distance=4`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 80 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.643474`, NMI `0.650870`, ARI `0.514985`, F1 macro `0.496722`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 29: `rank29_adaptive_fuzzy_clustering_full`.
  - Deep Adaptive Fuzzy Clustering adaptation with MLP expression autoencoder, trainable fuzzy membership layer, learnable fuzzifier, fuzzy centroid reconstruction, weighted adaptive entropy, partition-balance and center-separation regularizers, and periodic KMeans/EMA centroid evolution.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - GitHub URL: none listed in `02_整理索引.csv`; implementation follows the local PDF/report. PDF SHA-256 `36D60AC770E78AAA81661F097C02BA71C05D3D2E40F6AFC0E1B6F5EE73BFB9AE`.
  - Synthetic forward/backward smoke: PASS (`embedding=(12,16)`, membership `(12,4)`, finite loss, nonzero encoder and fuzzy-center gradients).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `latent_size=16`, `depth=1`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 80 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.580767`, NMI `0.642314`, ARI `0.450376`, F1 macro `0.518267`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

- Rank 30: `rank30_fuzzy_rough_boundary_full`.
  - Fuzzy Rough Sets Based on Fuzzy Quantification adaptation with MLP expression autoencoder, Student-t pseudo concepts, batch-local latent fuzzy relation, Kleene-Dienes implication, RIM S-function quantifiers, YWI-style lower approximation, unary upper approximation, rough boundary width loss, and periodic KMeans/EMA centroid updates.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - GitHub URL: none listed in `02_整理索引.csv`; implementation follows the local PDF/report. PDF SHA-256 `BAEC986C62573F3D0FC7E622020482A73B762A67C72E9B906BE66FC6549AF3E8`.
  - Synthetic forward/backward smoke: PASS (`embedding=(12,16)`, lower/upper approximations `(12,4)`, relation `(12,12)`, finite loss, nonzero encoder, fuzzy-center, and boundary-head gradients).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `latent_size=16`, `depth=1`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 80 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.582539`, NMI `0.617691`, ARI `0.486931`, F1 macro `0.447473`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

Additional skipped items:

- Rank 31 `scPilot`: LLM/API-backed omics-native reasoning and tool-use framework for cell-type annotation, trajectory inference, and GRN prediction rather than a trainable masked-expression representation model. Source clone `external_sources/rank31_scpilot`, commit `402b64ee029bd7bfb2d8b3d6296ec25d3e374a4f`. The repo requires API keys/providers and external data downloads; using it would alter unsupervised clustering fairness, so it is skipped as a primary scMAE benchmark variant and can only be considered later as an interpretability/failure-analysis adjunct.
- Rank 32 `ChatCell`: supervised T5/LLM natural-language interface over Cell2Sentence-style gene sequences, with Hugging Face model weights and instruction data for random/pseudo-cell generation, cell-type annotation, and drug sensitivity prediction. Source clone `external_sources/rank32_chatcell`, commit `f7203340709c31a36fda0a350f9b8c7eac636258`. It is not an unsupervised masked-expression clustering backbone and would require supervised instruction labels/pretrained LLM weights outside the fair scMAE protocol, so it is skipped as a primary benchmark variant.

- Rank 33: `rank33_ibot_online_tokenizer_full`.
  - iBOT adaptation with expression patch embeddings, a student Transformer, EMA teacher Transformer, shared iBOT projection head for `[CLS]` and patch tokens, centered teacher distributions, class-token cross-view self-distillation, masked patch-token online-tokenizer distillation, and auxiliary masked-expression reconstruction.
  - Required files present: `model.py`, `loss.py`, `run.py`, `README.md`, `source_manifest.json`.
  - Source clone: `external_sources/rank33_ibot`, commit `da316d82636a7a7356835ef224b13d5f3ace0489`.
  - Local PDF SHA-256 `e76c32b28234434443cba27abda926ec7c69c64af240d2d7aa78427131c4273d`.
  - Synthetic forward/backward smoke: PASS (`embedding=(12,32)`, finite loss, nonzero student and projection-head gradients, no teacher gradients, center update executed).
  - Small Melanoma real-data smoke: PASS (`n_top_genes=128`, `hidden_size=32`, `out_dim=32`, `bottleneck_size=16`, `patch_size=8`, `depth=1`, `epochs=1`, `skip_eval=true`).
  - Melanoma_5K quick screen: FAIL for formal promotion. Seed 42, 80 epochs, GPU 5, `n_top_genes=1000`, batch size 128: ACC `0.424108`, NMI `0.355382`, ARI `0.243218`, F1 macro `0.305138`. This is below the current scMAE Melanoma reference ARI of about `0.668`, so it should not enter the formal three-dataset benchmark without redesign.

Not yet implemented:

- Ranks 34+.

Next required step:

- Implement Rank 34 (`Masked Siamese Networks for Label-Efficient Learning`) with an independent MSN-style masked Siamese/self-distillation objective, then run smoke and Melanoma_5K quick screen.

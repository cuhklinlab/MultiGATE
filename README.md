# MultiGATE


![Overview](https://github.com/aqlkzf/typoraimg/raw/main//imgmac/MultiGATE_framework.png)

MultiGATE is a two-level graph-attention auto-encoder designed for spatial multi-omics analysis. It extracts the latent embeddings of the pixels/spots in spatial multi-omics data, while simultaneously incorporating the regulatory relationship of the cross-modality features through the cross-modality attention mechanism and the spatial relationship of the pixels/spots through the within-modality attention mechanism. In addition to reconstruction loss, a CLIP contrastive loss aligns embeddings across modalities.  MultiGATE yields (i) latent representations of pixels for clustering and visualization and (ii) cross-modality attention scores for cross-modality regulatory inference.

## Reproduce

The `/Reproduce` directory contains scripts and resources to reproduce the main figures and results of MultiGATE, including spatial clustering results, cis-regulation, trans-regulation analysis, long-range interaction detection, and protein-gene interaction findings. Please refer to the instructions and scripts in this directory to replicate the key analyses and visualizations presented in the paper.

## Usage && installation

Please follow the [Tutorials](https://multigate.readthedocs.io/en/latest) for installation and Usage.


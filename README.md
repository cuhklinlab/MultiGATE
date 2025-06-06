# MultiGATE

> MultiGATE is a novel method that utilizes a two-level graph attention auto-encoder to integrate multi-modality and spatial information effectively.

![Overview](./fig/MultiGATE_framework.png)

MultiGATE is a two-level graph attention auto-encoder designed for spatial multi-omics analysis. 
It integrates chromatin accessibility or protein information with gene expression data using an attention mechanism. 
The model employs cross-modality attention at the first level to capture inter-modality relationships, and withinmodality attention at the second level to incorporate spatial information. 
In addition to reconstruction loss, CLIP loss is employed to enhance embedding mapping. This comprehensive approach enables the aggregation of spatial information and feature expression, facilitating joint clustering for precise spatial analysis using both modalities.

## Usage && installation

Please follow the [Tutorials](https://multigate.readthedocs.io/en/latest) for installation and Usage.


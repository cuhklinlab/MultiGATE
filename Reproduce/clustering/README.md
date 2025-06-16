# Clustering Results Reproduction

This directory contains scripts and notebooks to reproduce the clustering results for 7 datasets presented in the MultiGATE paper. These datasets span five different spatial multi-omics platforms/technologies and cover both stereotypical tissues (brain and spleen) and non-stereotypical samples (breast cancer).

## Datasets Overview

The following 7 datasets are included in this reproduction suite:

### Original 4 Benchmark Datasets:
1. **Human Hippocampus Dataset** (Spatial ATAC–RNA) - `human_hippData/`
2. **P22 Mouse Brain Dataset** (Spatial ATAC–RNA) - `P22mouseData/`
3. **Murine Spleen Dataset** (SPOTS) - `SPOTS_spleenData/`
4. **Metastatic Melanoma Dataset** (Slide-tags) - `slidetags_data/`

### Additional Datasets (Revised Version):
5. **Spatial RNA + Metabolomics Dataset** (SMA technology) - `spatial_metabolic_RNA_data/`
6. **Breast Cancer Dataset** (Visium-CytAssist RNA–protein) - `breast_cancer_10xdata/`
7. **Simulation Dataset** (breast cancer-patterned spatial ATAC + RNA) - `Simulated_data/`

## Platforms and Technologies

The datasets cover five different spatial multi-omics platforms:
- **Spatial ATAC–RNA**: Human hippocampus, P22 mouse brain, simulation datasets
- **Slide-tags**: Metastatic melanoma dataset
- **Visium-CytAssist RNA–protein**: Breast cancer dataset
- **SMA**: Spatial RNA + metabolomics dataset
- **SPOTS**: Murine spleen dataset

## Molecular Combinations

Different molecular combinations are represented:
- Spatial ATAC + RNA
- Spatial RNA + protein/ADT
- Spatial RNA + metabolite

## Directory Structure

```
clustering/
├── MultiGATE/                    # Core MultiGATE implementation
│   ├── MultiGATE.py             # Main MultiGATE class
│   ├── Train_MultiGATE.py       # Training module
│   ├── model_MultiGATE.py       # Model architecture
│   ├── utils.py                 # Utility functions
│   └── genomics.py              # Genomics-specific functions
├── human_hippData/              # Human hippocampus (Spatial ATAC–RNA)
│   ├── TrainMultiGATE.ipynb     # Training notebook
│   ├── spatialplot.ipynb        # Spatial visualization
│   └── human_hippResults/       # Output results
├── P22mouseData/                # P22 mouse brain (Spatial ATAC–RNA)
│   ├── RunMultiGATE.ipynb       # Training and analysis
│   ├── spatialplot.ipynb        # Spatial visualization
│   └── P22Mousebrain_Results/   # Output results
├── SPOTS_spleenData/            # Murine spleen (SPOTS)
│   ├── RunMultiGATE.ipynb       # Training and analysis
│   ├── spatialplot.ipynb        # Spatial visualization
│   └── SPOTS_spleenResults/     # Output results
├── slidetags_data/              # Metastatic melanoma (Slide-tags)
│   ├── RunMultiGATE.ipynb       # Training and analysis
│   ├── spatialplot.ipynb        # Spatial visualization
│   └── slidetags_Results/       # Output results
├── spatial_metabolic_RNA_data/  # Spatial RNA + metabolomics (SMA)
│   ├── trainSTSM_MultiGATE.py   # Training script
│   ├── spatialplot.ipynb        # Spatial visualization
│   ├── InputData/               # Input data files
│   └── STSMresults/             # Output results
├── breast_cancer_10xdata/       # Breast cancer (Visium-CytAssist RNA–protein)
│   ├── train_MultiGATE.py       # Training script
│   ├── spatialplot.ipynb        # Spatial visualization
│   ├── Inputdata/               # Input data files
│   └── breastcancerResults/     # Output results
└── Simulated_data/              # Simulation (breast cancer-patterned ATAC + RNA)
    ├── train_simulation.py      # Training script
    ├── spatialplot.ipynb        # Spatial visualization
    └── simulation_dataResults/  # Output results
```





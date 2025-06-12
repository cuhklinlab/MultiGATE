# Gene-Protein Interaction Analysis Pipeline

Analysis pipeline for gene-protein interactions using attention mechanisms from spatial multi-omics data.

## Scripts

1. **`gene_pathway_processor.py`** 
   - Processes raw gene pathway data
   - Maps genes to pathway annotations 
   - Integrates attention scores with gene information
   - Outputs standardized gene-pathway datasets

2. **`attention_analyzer.py`** 
   - Analyzes attention matrices from gene-protein networks
   - Generates pathway comparison visualizations
   - Performs statistical analysis across pathway types
   - Creates publication-ready plots

## Quick Start

```bash
# 1. Process pathway data
python gene_pathway_processor.py

# 2. Analyze and visualize
python attention_analyzer.py
```

## Data Files

### Input Files
**Small Input Files (included in repository):**
- `input_*.csv` - Gene pathway data and attention scores
- `*.npz` - Gene-peak attention matrices

**Large Data Files (download required):**
- `rna_data0320.h5ad` (142.8 MB) - RNA single-cell data
- `protein_data0320.h5ad` - Protein single-cell data

### Download Large Data Files

Due to GitHub's file size limitations, large data files are stored externally:

**Download Link:** [https://drive.google.com/file/d/15fi3c_Eg69Suf3JfhnsipMdzUFaV2FLM/view?usp=drive_link](https://drive.google.com/file/d/15fi3c_Eg69Suf3JfhnsipMdzUFaV2FLM/view?usp=drive_link)

**Instructions:**
1. Click the download link above
2. Download the data files to the `protein-gene_interaction/reproduce/` directory
3. Ensure files are named exactly as `rna_data0320.h5ad` and `protein_data0320.h5ad`
4. Run the analysis scripts as described in Quick Start

### Output Files
**Output:** `output_*.csv`, visualizations in `results_single_analysis_3types/`

## Output

```
results_single_analysis_3types/
├── output_analysis_results.csv
├── output_pathway_comparison.png
└── output_pathway_comparison.pdf
```

 
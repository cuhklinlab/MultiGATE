import os
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Train MultiGATE model with customizable hyperparameters')
parser.add_argument('--gpu', type=str, default='1', help='GPU device ID to use')
parser.add_argument('--base_path', type=str, default='./output0311/', help='Base path for output files')
parser.add_argument('--rad_cutoff', type=float, default=503, help='Radius cutoff for spatial network calculation')
parser.add_argument('--num_epoch', type=int, default=400, help='Number of training epochs')

parser.add_argument('--temp', type=float, default=6, help='Temperature parameter for training')
parser.add_argument('--bp_width', type=int, default=450, help='Base pair width parameter')

args = parser.parse_args()

# Set GPU device
os.environ["CUDA_VISIBLE_DEVICES"] = args.gpu

import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, '/lustre/project/Stat/s1155202250/fastfolder/code/st/MultiGATE/MultiGATEgithub0607/MultiGATE/Reproduce/clustering')

import MultiGATE
import scanpy as sc;import matplotlib.pyplot as plt
import os
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np


base_path = args.base_path 

os.makedirs(base_path, exist_ok=True)

rna_file = "./output0311/simulated_data_adata1.h5ad" 
adata1 = sc.read_h5ad(rna_file)



# Load ATAC data
atac_file="./output0311/simulated_data_adata2.h5ad"
adata2 = sc.read_h5ad(atac_file)


adata1 = adata1[:, adata1.var['highly_variable']]
# keep all ADTs in adata2, therefore we do not select highly variable features here.
adata2 = adata2[:, adata2.var['highly_variable']]

MultiGATE.Cal_gene_peak_Net_new(adata1, adata2, 150000, file='/lustre/project/Stat/s1155202250/fastfolder/code/st/MultiGATE/tutorial/data_tutorial/human/gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz')
adata1.uns['gene_peak_Net'] = adata2.uns['gene_peak_Net']


adata1, adata2 = MultiGATE.train_MultiGATE(
    adata1, adata2, 
    n_epochs=args.num_epoch,
    verbose=False,
    type='ATAC_RNA', 
    temp=args.temp,
    bp_width=args.bp_width,
    save_attention=False,
)

if 'MultiGATE_attention' in adata1.uns:
    del adata1.uns['MultiGATE_attention']
if 'MultiGATE_gene_peak_attention' in adata1.uns:
    del adata1.uns['MultiGATE_gene_peak_attention']

# Clean up adata1 object before saving
# Remove attention matrices if present
if 'MultiGATE_attention' in adata1.uns:
    del adata1.uns['MultiGATE_attention']
if 'MultiGATE_gene_peak_attention' in adata1.uns:
    del adata1.uns['MultiGATE_gene_peak_attention']

if 'MultiGATE_attention' in adata2.uns:
    del adata2.uns['MultiGATE_attention']
if 'MultiGATE_gene_peak_attention' in adata2.uns:
    del adata2.uns['MultiGATE_gene_peak_attention']


# Convert categorical categories to strings
for col in adata1.obs.columns:
    if pd.api.types.is_categorical_dtype(adata1.obs[col]):
        # Convert categories to string type
        adata1.obs[col] = adata1.obs[col].cat.rename_categories(
            adata1.obs[col].cat.categories.astype(str)
        )

for col in adata2.obs.columns:
    if pd.api.types.is_categorical_dtype(adata2.obs[col]):
        # Convert categories to string type
        adata2.obs[col] = adata2.obs[col].cat.rename_categories(
            adata2.obs[col].cat.categories.astype(str)
        )

print(adata1)
# save adata1 and adata2
adata1.write(os.path.join(base_path, 'simulated_data_adata1.h5ad'))
adata2.write(os.path.join(base_path, 'simulated_data_adata2.h5ad'))
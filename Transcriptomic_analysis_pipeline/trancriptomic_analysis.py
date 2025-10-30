#!/usr/bin/env python3
"""
Transcriptomic Data Analysis Pipeline (Python version)
Author: [Your Name]
Description:
    End-to-end transcriptomic analysis using Python.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from statsmodels.stats.multitest import multipletests
from scipy import stats
import os

# --- Setup ---
DATA_DIR = "data/"
OUTPUT_DIR = "results/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Load Data ---
adata = sc.read_loom(os.path.join(DATA_DIR, "expression_data.loom"))  # Example file format
print(adata)

# --- Quality Control ---
sc.pl.violin(adata, ['n_genes', 'n_counts'], jitter=0.4)
sc.pl.highest_expr_genes(adata, n_top=20)

# --- Normalization & Log transform ---
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata)

# --- Experimental Design ---
# Example: 'group' column in adata.obs contains 'control' and 'treatment'
groups = adata.obs['group'].unique()

control = adata[adata.obs['group'] == 'control']
treatment = adata[adata.obs['group'] == 'treatment']

# --- Differential Expression ---
expr_control = np.mean(control.X, axis=0)
expr_treatment = np.mean(treatment.X, axis=0)

# T-test
t_stat, p_val = stats.ttest_ind(control.X, treatment.X, axis=0, equal_var=False)
logFC = np.log2((expr_treatment + 1e-9) / (expr_control + 1e-9))

# Multiple testing correction
adj_p = multipletests(p_val, method='fdr_bh')[1]

# Create DataFrame
results = pd.DataFrame({
    'gene': adata.var_names,
    'logFC': logFC,
    'pval': p_val,
    'adj_pval': adj_p
})

# Filter significant genes
deg = results[(results['adj_pval'] <= 0.05) & (abs(results['logFC']) >= np.log2(2))]
deg.to_csv(os.path.join(OUTPUT_DIR, "differentially_expressed_genes.csv"), index=False)

print(f"{len(deg)} differentially expressed genes identified.")

# --- Visualization ---
# Volcano plot
plt.figure(figsize=(8,6))
sns.scatterplot(x='logFC', y=-np.log10(deg['adj_pval']), data=deg)
plt.title("Volcano Plot of Differential Expression")
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 Adjusted P-value")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "volcano_plot.png"))

# Heatmap
top_genes = deg.head(50)['gene']
sc.pl.heatmap(adata, var_names=top_genes, groupby='group', cmap='viridis', show=True)

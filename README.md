# Bioinformatics-projects

# 1) 🧬 Transcriptomic Data Analysis Pipeline

This repository provides **R and Python implementations** of an end-to-end transcriptomic data analysis workflow.  
It supports the analysis of gene expression data from microarrays or RNA-seq to identify **differentially expressed genes (DEGs)** and visualize expression patterns.
The goal of the script is to analyze gene expression data from microarray experiments (specifically, Affymetrix .CEL files) to:

- Preprocess and normalize the data

- Identify genes that are differentially expressed between conditions (e.g., control vs disease)

- Annotate those genes with biological identifiers

- Visualize expression patterns with heatmaps and clustering
---

- Install dependencies:  pip install scanpy pandas numpy statsmodels matplotlib seaborn
#
- Run the script:  python transcriptomic_analysis.py

## 🚀 Features
- Load and preprocess raw gene expression data  
- Perform background correction and normalization  
- Conduct differential expression analysis (control vs treatment)  
- Annotate significant genes using Ensembl  
- Generate visualizations:
  - Quality control plots  
  - Volcano plots  
  - Heatmaps with clustering  
---
## 🧠 R Version
Uses the **Bioconductor** ecosystem:
- [`affy`](https://bioconductor.org/packages/affy/)
- [`limma`](https://bioconductor.org/packages/limma/)
- [`biomaRt`](https://bioconductor.org/packages/biomaRt/)
- [`gplots`](https://cran.r-project.org/package=gplots`)
- [`dendextend`](https://cran.r-project.org/package=dendextend`)

To run:
```bash
Rscript transcriptomic_analysis.R
------------------------------------------------------------------------------------------------------------------------------------------------

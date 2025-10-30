# Bioinformatics Projects

# 1)Transcriptomic Data Analysis Pipeline

This repository presents a reproducible and modular pipeline for transcriptomic data analysis, implemented in both R and Python.
The project aims to demonstrate practical expertise in computational biology, data preprocessing, differential gene expression analysis, and biological data visualization.

The workflow is designed for microarray or RNA-seq datasets (e.g., Affymetrix .CEL files) and enables the identification of differentially expressed genes (DEGs) between experimental conditions such as control vs disease.

## Objectives

This project illustrates the core analytical steps of a bioinformatics transcriptomic workflow, including:
- Preprocessing and normalization of high-dimensional gene expression data
- Statistical modeling for the detection of DEGs using robust linear models
- Annotation of DEGs with biological identifiers and gene metadata
- Exploratory visualization through clustering and heatmaps
- Cross-validation of results between R and Python implementations

This pipeline emphasizes data integrity, reproducibility, and scientific rigor, mirroring best practices used in research environments.

 ## Key Features

- Quality control and exploratory data visualization (histograms, boxplots, MA plots)

- Background correction and multiple normalization strategies (RMA, quantile, LOESS)

- Differential expression analysis using limma and empirical Bayes statistics

- Gene annotation using Ensembl via biomaRt (exported as DEG.txt)

- Heatmap visualization and hierarchical clustering of top DEGs

- Parallel Python implementation leveraging Scanpy and Statsmodels

## R Implementation

Developed using the Bioconductor ecosystem, the R version focuses on transparency, interpretability, and statistical robustness.

Dependencies:
- affy – preprocessing and background correction
- limma – differential expression analysis
- biomaRt – gene annotation
- gplots, dendextend – visualization and clustering

--

Run the pipeline:

Rscript transcriptomic_analysis.R

--

Output files include:
- Normalized expression matrices
- DEG tables (DEG.txt)
- Heatmaps and cluster dendrograms

## Python Implementation

The Python version mirrors the analytical flow using modern data science frameworks.
It provides flexibility for integration with RNA-seq data and machine learning-based analyses.

Dependencies:

pip install scanpy pandas numpy statsmodels matplotlib seaborn

Run the pipeline:

python transcriptomic_analysis.py

The Python workflow performs:
- Data preprocessing and normalization
- DEG identification using linear modeling
- Visualization with heatmaps and PCA


## Repository Structure
Bioinformatics-projects/
│
├── transcriptomic_analysis.R        # R implementation (Bioconductor)
├── transcriptomic_analysis.py       # Python implementation (Scanpy)
├── data/                            # Input data (e.g., .CEL or count matrices)
├── results/                         # Output DEG tables and visualizations
└── README.md                        # Project overview

## References

Smyth, G. K. (2004). Linear Models and Empirical Bayes Methods for Assessing Differential Expression in Microarray Experiments. Statistical Applications in Genetics and Molecular Biology.

Bioconductor Project: https://bioconductor.org
Scanpy Documentation: https://scanpy.readthedocs.io

## Scientific Context

This project was developed as part of an academic initiative focused on bioinformatics research and data-driven biology.
It integrates techniques from transcriptomics, statistics, and machine learning to extract biological insights from high-throughput data.

The pipeline can be extended to incorporate:
- RNA-seq differential expression (DESeq2 or EdgeR)
- Gene Ontology (GO) and pathway enrichment analysis
- Machine learning-based classification of gene expression profiles

# 2) RNA Structure :  RNA Molecular Graph Data Exploration and Generation


This repository contains tools for **loading, inspecting, and preparing RNA molecular graph data** stored in serialized (`.pickle`) format using the **NetworkX** library.

The goal is to:
- Understand the structure and attributes of RNA graphs derived from Protein Data Bank (PDB) data,
- Detect chemically modified nucleotides and specialized interactions,
- Prepare the foundation for **dataset generation** for machine learning or structural analysis.

---

##  Project Overview

RNA molecular structures can be represented as **graphs**, where:
- **Nodes** correspond to nucleotides (`A`, `C`, `G`, `U`, or modified bases),
- **Edges** represent chemical or structural interactions between nucleotides.

Each `.pickle` file in the `NetworkxGraph/` directory stores a single molecular graph, serialized as a **NetworkX** object.

---

##  Notebook: `GenerateData_working_fulldata.ipynb`

This notebook is primarily a **data exploration and verification tool**.  
It helps confirm that the RNA graph data is structured correctly and consistent across all files.

# Author

This project reflects the my technical expertise and research experience in:
- Computational biology and bioinformatics
- Data analytics and machine learning
- Scientific programming in R and Python
- Biological data visualization and interpretation

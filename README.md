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

# 2) RNA Structure 

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

Each `.pickle` file in the `NetworkxGraph/` directory (the one provided is far from complete, you can reach out to me obtain the all dataset) stores a single molecular graph, serialized as a **NetworkX** object.

---
## > RNA Molecular Graph Analysis

This repository contains a Python script that loads and analyzes a molecular graph (e.g., RNA structure) stored in a pickle file using **NetworkX**.

---

##  Overview

The script 'algo.py' performs the following operations:

1. **Loads** a serialized molecular graph from a `.pickle` file.  
2. **Iterates through all nodes** to identify **chemically modified nucleotides**.  
3. **Traverses all edges** to find specific **chemical interactions** labeled `'CHH'`.

---

##  Code Explanation
### 1. Imports
- pickle is used to load serialized Python objects from files.
- networkx is a library for creating and manipulating graphs/networks.

### 2. Load the graph from a pickle file
- This opens a file called 1S72.pickle in binary read mode.
- pickle.load(f) deserializes the file and loads it as a Python object.
- Here, it’s expected to be a NetworkX graph object, likely representing a molecular structure (e.g., RNA).

### 3. Iterate over all nodes and access node attributes

- graph.nodes(data=True) iterates over each node and its associated attributes.

- Each node is assumed to be a tuple (position, chain):
 - position → probably the nucleotide’s position in the sequence.
 - chain → the chain identifier in the molecule.
- data is a dictionary with attributes of the node, like nucleotide.

### 4. Check for chemically modified nucleotides
- Checks if the nucleotide is not one of the standard RNA bases (A, C, G, U).
- If it’s something else (chemically modified), it prints information about that node: its position, chain, and nucleotide type.
So this part identifies modified nucleotides in the RNA structure.

### 5. Iterate over edges and check for “CHH” labels
- graph.edges(data=True) iterates over all edges and their attributes.
- Checks if the edge has a label CHH.
- If yes, it prints the source node, target node, and the edge label.
  
This part identifies specific types of interactions/edges labeled 'CHH' in the molecular graph.

### Summary
- Loads a molecular structure (RNA) graph from a pickle file.
- Iterates through all nodes to find chemically modified nucleotides.
- Iterates through all edges to find edges labeled 'CHH' (possibly a specific hydrogen bond or chemical interaction).


## Example Output
node=(12, 'A'), data['nucleotide']='m3C', chain='A', position=12

source=(5, 'B'), target=(8, 'B'), data['label']='CHH'


This output indicates:
- A modified nucleotide m3C at position 12 on chain A.
- A chemical interaction (CHH) between nucleotides at positions 5 and 8 on chain B.
  
<img width="1536" height="1024" alt="graph_example" src="https://github.com/user-attachments/assets/5ffc80e4-d338-411c-8254-2d67e19823bd" />

---
##  > Notebook: `GenerateData_working_fulldata.ipynb`



This notebook doesn’t perform analysis yet, is primarily a **data exploration and verification tool**. 
It helps you:
- Understand the structure and attributes inside each RNA graph,
- Check for modified nucleotides,
- Prepare for downstream analysis (e.g., detecting motifs, analyzing interactions, or generating datasets).
  
It is designed to inspect and explore RNA molecular graph data stored in NetworkX graph objects, which have been serialized as .pickle files.

---
## Code overview
### 1. Imports libraries

It imports standard scientific Python libraries (NumPy, Pandas, Seaborn, Matplotlib) and networkx for working with graphs.

###  2. Lists files in the “NetworkxGraph” directory

It lists all .pickle files inside the NetworkxGraph/ folder — each of these likely represents a molecular graph derived from a PDB (Protein Data Bank) structure.

### 3. Loads one example graph
Loads a specific file (e.g., 2KH9.pickle) to examine its contents.
The loaded object is a NetworkX Graph that represents a molecule (probably RNA or a similar structure).

### 4. Inspects graph nodes
i = 0

for node, data in graph.nodes(data=True): [...]

This loops through every node in the graph and prints:
- The node identifier (usually a tuple like (position, chain)),
- The node’s associated attributes (stored in a dictionary, e.g. {'nucleotide': 'A'}).

This helps understand what information is stored per nucleotide.

### 5. Loops through several files to preview graph contents

For the first 10 pickle files:
- Loads each graph.
- Prints each node’s index and nucleotide type (e.g., A, C, G, U, or modified bases).

This part is purely exploratory — the goal is to check consistency and understand the data structure across multiple molecular graphs.

### Summary of what the notebook does
1. Load libraries	Prepare environment for data inspection.
2. List .pickle files	Identify available molecular graph datasets.
3. Load example graph	Open and visualize one molecular graph.
4. Explore nodes	Print node details (position, chain, nucleotide).
5. Batch check multiple graphs	Verify data consistency across several files.

---
# Author

This project reflects the my technical expertise and research experience in:
- Computational biology and bioinformatics
- Data analytics and machine learning
- Scientific programming in R and Python
- Biological data visualization and interpretation
  
CISSE Khaly Bécaye Ba, MSc in Computer Science, Bioinformatics

khalybecayecisse@gmail.com

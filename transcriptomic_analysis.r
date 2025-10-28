##########################################################
# Transcriptomic Data Analysis Pipeline
# Author: [Your Name]
# Description: End-to-end microarray analysis workflow
# using Bioconductor (Affymetrix, limma, biomaRt)
##########################################################

# --- Setup -------------------------------------------------------------

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Core Bioconductor packages
BiocManager::install(c("affy", "limma", "biomaRt"), ask = FALSE)

# CRAN visualization packages
install.packages(c("gplots", "dendextend"), dependencies = TRUE)

# Load libraries
library(affy)
library(limma)
library(gplots)
library(dendextend)
library(biomaRt)

# --- Define directories ------------------------------------------------
data_dir <- "data/"
results_dir <- "results/"
dir.create(results_dir, showWarnings = FALSE)

setwd(data_dir)

# --- Load Affymetrix data ---------------------------------------------
cat("Reading .CEL files...\n")
cel <- ReadAffy()
summary(cel)

# --- Quality Control ---------------------------------------------------
par(mfrow = c(1, 2))
hist(cel, main = "Raw intensity distribution")
boxplot(cel, horizontal = TRUE, main = "Boxplot of raw intensities")

# --- Background correction and normalization ---------------------------
cat("Performing background correction and normalization...\n")
bg_corrected <- bg.correct(cel, method = "rma")
normalized <- normalize(bg_corrected, method = "quantiles")

# --- Visualization after normalization --------------------------------
par(mfrow = c(1, 2))
hist(normalized, main = "Normalized intensities")
boxplot(normalized, horizontal = TRUE, main = "Boxplot of normalized data")

# --- Experimental design ----------------------------------------------
cat("Creating design and contrast matrices...\n")

# Example: 7 controls, 7 treatments
design <- model.matrix(~ 0 + factor(c(rep("control", 7), rep("treatment", 7))))
colnames(design) <- c("control", "treatment")
contrast.matrix <- makeContrasts(treatment - control, levels = design)

print(design)
print(contrast.matrix)

# --- Expression summarization -----------------------------------------
cat("Computing expression set...\n")
eset <- computeExprSet(normalized, pmcorrect.method = "pmonly", summary.method = "avgdiff")
expr <- log2(exprs(eset))

# --- Differential Expression Analysis ---------------------------------
cat("Running differential expression analysis with limma...\n")
fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
efit <- eBayes(fit2)

deg <- topTable(efit,
                adjust = "BH",
                p.value = 0.05,
                number = Inf,
                sort.by = "logFC")

cat("Number of significant genes:", nrow(deg), "\n")

# --- Gene annotation ---------------------------------------------------
cat("Annotating genes using biomaRt...\n")
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

annotation <- getBM(
  attributes = c("affy_hg_u133_plus_2",
                 "hgnc_symbol",
                 "chromosome_name",
                 "ensembl_gene_id",
                 "external_gene_name"),
  filters = "affy_hg_u133_plus_2",
  values = rownames(deg),
  mart = ensembl
)

# Merge expression data with annotation
merged <- merge(annotation, expr[rownames(deg), ],
                by.x = "affy_hg_u133_plus_2",
                by.y = "row.names")

write.table(merged,
            file = file.path(results_dir, "differentially_expressed_genes.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Annotated results saved in 'results/differentially_expressed_genes.txt'\n")

# --- Heatmap Visualization --------------------------------------------
cat("Generating heatmap of top genes...\n")

# Subset for top 50 DE genes
top_genes <- head(rownames(deg), 50)
heatmap_data <- expr[top_genes, ]

# Scale rows (z-score)
scaled_data <- t(scale(t(heatmap_data)))

# Plot heatmap
heatmap.2(as.matrix(scaled_data),
          trace = "none",
          density.info = "none",
          col = colorRampPalette(c("yellow", "blue"))(64),
          main = "Heatmap of Top Differentially Expressed Genes",
          margins = c(8, 10),
          cexRow = 0.6)

# --- Clustering --------------------------------------------------------
cat("Performing hierarchical clustering...\n")
hc <- hclust(dist(scaled_data))
plot(hc, main = "Hierarchical Clustering of DE Genes", cex = 0.5)
abline(h = 7, col = "red", lty = 2)

cat("Analysis complete.\n")
##########################################################

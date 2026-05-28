# pan-myocarditis-snRNA-seq
> Single-cell RNA-seq analysis of mouse heart immune cells in myocarditis

---

## 1. System & Hardware Requirements

### System Requirements
* **Operating System:** Linux (Ubuntu 20.04.6 LTS)
* **Software Environment:** R version 4.3.1
* **Core Dependencies (R packages):** * `Seurat` (v4.3.0.1), `ROGUE` (v1.0), `lisi` (v1.0)
    * `clusterProfiler` (v4.8.3), `org.Mm.eg.db` (v3.17.0), `msigdbr` (v7.5.1)
    * `GSVA` (v1.48.3), `AUCell` (v1.22.0), `pySCENIC` (v0.12.1)
    * `CellChat` (v1.6.1), `Monocle 2` (v2.28.0), `Monocle 3` (v1.0.0)

### Hardware Requirements
* **Operating Environment:** Tested on a high-performance Linux server equipped with 104 CPU cores and 755 GB RAM.
* **Recommended Minimum:** To process the full pan-myocarditis single-cell dataset and run intensive tasks (such as large-scale integration or cell-cell communication inference), a Linux workstation/server with at least **16 CPU cores and 64 GB RAM** is highly recommended.

---

## 2. Installation Guide

### Instructions
Most dependencies are standard R/Bioconductor packages and can be installed via CRAN or Bioconductor. Please run the following scripts in your R console to set up the environment:

```R
# 1. Install CRAN packages
cran_packages <- c("Seurat", "lisi", "msigdbr", "remotes", "devtools")
install.packages(cran_packages)

# 2. Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_packages <- c("clusterProfiler", "org.Mm.eg.db", "GSVA", "AUCell", "monocle")
BiocManager::install(bioc_packages)

# 3. Install specific versions or GitHub repositories (ROGUE, CellChat, Monocle3)
remotes::install_github("PaulingLiu/ROGUE")
remotes::install_github("sqjin/CellChat")

# For Monocle 3 installation
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'coolwarm'))
remotes::install_github('cole-trapnell-lab/monocle3')

```

### Typical Installation Time

- **Standard Desktop (for Demo dependencies)**: ~20–30 minutes to install standard CRAN/Bioconductor dependencies.

---

## 3. Demo Instructions

This demo uses an optimized subset of the immune cells (`demo_immune_cells.rds`) to quickly reproduce key figures and verify the pipeline execution.

### Instructions to run
1. **Download the demo data:** Place `demo_immune_cells.rds` in your working directory (or adjust the path in the script).
2. **Run the demo script:**
   * *Option A (Within R/RStudio console):*
     ```R
     source("demo_main.R")
     ```
   * *Option B (From the terminal/command line):*
     ```bash
     Rscript demo_main.R
     ```

### Expected Output
Upon successful execution, the script will generate the following figures:
* **UMAP plot:** A 2D scatter plot demonstrating 6 major immune cell types (*B cells, T cells, NK cells, DCs, Macrophages/Monocytes, and Neutrophils*), with each type represented by 100 downsampled cells.
* **Violin plot:** Three QC quality panels showing the clear distribution of gene counts, UMI counts, and mitochondrial gene percentage.
* **FeaturePlot:** Four spatially-resolved expression maps of canonical markers (`Cd79a`, `Cd3d`, `Cd68`, `S100a8`) to highlight cell-type specificity.
* **DotPlot (Optional):** Summarized expression level and percentage of the four key markers across all identified cell lineages.

### Expected Runtime
- **~ Less than 20 seconds** on a typical desktop computer.


# 🧬 AMR_TZ Project — Scripts

This folder contains R scripts and notebooks used in the **AMR_TZ** project, which investigates antimicrobial resistance (AMR) in Tanzanian populations using metagenomic data and phylogenetic analysis.

## 📁 Contents

| File | Description |
|------|-------------|
| `create_trees.pbs` | align gene variants from flankophile and create phylogenetic trees with IQ-TREE. |
| `gene_identity.R` | Script for plotting gene identity across populations. |
| `permanova_trees.R` | Performs PERMANOVA analysis on phylogenetic trees to assess host population effects. |
| `plot_selected_gene_trees.R` | Generates visualizations of gene trees, highlighting host-specific clustering. |
| `variants_vs_samples.Rmd` | Plots number of variants vs number of samples per population. |
| `variants_vs_samples.nb.html` | Rendered HTML output of the above notebook. |

## 🧪 Requirements

- R (≥ 4.0)
- Packages:
  - `ggplot2`
  - `ggtree`
  - `vegan`
  - `ape`
  - `dplyr`
  - `tidyr`
  - `readr`
  - `phytools`
  - `ggpubr`

You can install missing packages using:
```r
install.packages(c("ggplot2", "vegan", "ape", "dplyr", "tidyr", "readr", "phytools", "ggpubr"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")
```

## 🚀 Usage

In each script the working directory must be set in the top of the script. This directory should contain the "1_hits_all.tsv" from running flankophile and the treefiles created from create_trees.pbs

## 🌍 Project Context

These scripts support the analysis of host-specific AMR gene variants and phylogenetic clustering in samples from:
- Adult village populations
- Hadza hunter-gatherers
- School children
- Chimpanzees

## 📌 Repository

This folder is part of the [AMR_TZ GitHub repository](https://github.com/marieljespersen/AMR_TZ).

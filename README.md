# PCA & Volcano Plot Shiny App

A Shiny web application for interactive bioinformatics visualization — PCA plots and Volcano plots for differential expression analyses.

**Live app:** https://jwkim0.shinyapps.io/PCA-Volcano-shinyapp/

## Features

### Volcano Plot
- Upload CSV or TSV output from DESeq2, edgeR, limma, etc.
- Auto-detects log2FC, p-value, and gene name columns
- Adjustable p-value and log2FC thresholds
- Points colored by significance (Up / Down / NS)
- Label top N most significant genes (via ggrepel)
- Custom plot title
- Download plot as PNG or SVG

### PCA Plot
- Upload an expression matrix (features × samples) and a metadata CSV
- Auto-runs PCA with optional log2 transformation (configurable pseudocount) and feature scaling
- Color and shape points by any metadata column (e.g. condition, batch)
- Toggle sample labels on/off
- Select which PCs to plot on each axis (% variance explained shown)
- Custom plot title
- Download plot as PNG or SVG

### General
- "Load toy data" button on each tab for instant demo without uploading files

## Running the App Locally

```r
shiny::runApp()
```

## Input File Formats

### Volcano Plot
A CSV or TSV with at least one log2 fold change column and one p-value column.

| gene | log2FoldChange | pvalue | padj |
|------|---------------|--------|------|
| GeneA | 2.3 | 0.001 | 0.01 |
| GeneB | -1.5 | 0.03 | 0.08 |

### PCA Plot — Expression Matrix
Rows = features (genes), columns = samples. First column = feature IDs.

| GeneID | Sample1 | Sample2 | Sample3 |
|--------|---------|---------|---------|
| Gene1  | 1020    | 980     | 1050    |
| Gene2  | 130     | 115     | 140     |

### PCA Plot — Metadata
One row per sample. Must include a column of sample IDs matching the expression matrix column names.

| SampleID | Condition | Batch |
|----------|-----------|-------|
| Sample1  | Control   | A     |
| Sample2  | Treatment | A     |

Toy datasets (`toy_expression.csv`, `toy_metadata.csv`, `toy_de_results.csv`) are included for testing.

## Requirements

- R (≥ 4.0)
- shiny
- ggplot2
- ggrepel
- readr

Install dependencies:

```r
install.packages(c("shiny", "ggplot2", "ggrepel", "readr"))
```

## License

MIT

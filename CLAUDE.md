# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

A Shiny web application for bioinformatics visualization — PCA (Principal Component Analysis) plots and Volcano plots for differential expression analyses. Compatible with DESeq2, edgeR, limma output.

## Running the App

```r
shiny::runApp()
```

```bash
Rscript -e "shiny::runApp('.')"
```

Install dependencies if needed:
```r
install.packages(c("shiny", "ggplot2", "readr"))
```

## Architecture

All code lives in a single `app.R` (UI + server). Split into `ui.R` / `server.R` only if the file grows large.

**Volcano tab** — user uploads a CSV/TSV (e.g. DESeq2 results). `guess_col()` auto-detects log2FC and p-value columns by name pattern. `col_selectors` UI is rendered dynamically after upload. The reactive chain is: `data` → `volcano_data` (adds `neg_log10_p` and `category`) → `make_volcano` (ggplot) → `output$volcano` / `output$download_volcano`.

**PCA tab** — two file uploads: expression matrix (features × samples, first column = feature IDs) and metadata CSV. Dynamic controls (`pca_controls`) render only after both files are loaded. `pca_validation` checks that sample IDs in the expression matrix columns match those in the chosen metadata ID column, and surfaces mismatches as a red error above the plot. The reactive chain is: `pca_expr` + `pca_meta` → `pca_validation` → `pca_result` (runs `prcomp`, merges metadata scores) → `make_pca` → `output$pca_plot` / `output$download_pca`.

## Code Style

- 2-space indentation
- UTF-8 encoding
- Leave comments on non-obvious logic

## Test Data

`toy_expression.csv` and `toy_metadata.csv` are included for manual testing. Expression matrix has 9 samples (Ctrl/Treat/KO × 3 replicates) and the metadata has matching `SampleID`, `Condition`, and `Batch` columns.

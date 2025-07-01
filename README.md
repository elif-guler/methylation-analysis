#  Analysis of CpG Methylation Differences Between Healthy and Diseased States

This repository contains an R script to analyze differential DNA methylation using Illumina HumanMethylation450k array data, comparing healthy and diseased samples.

## Files

- `analysis.R` – Main R script that runs the complete methylation analysis pipeline.

## Project Overview

This project aims to identify CpG sites that are differentially methylated between healthy and diseased individuals. The steps include:

- Sample and probe quality control  
- Detection p-value filtering (threshold: 0.05)  
- Quantile normalization (`preprocessQuantile`)  
- Statistical testing using **t-tests**  
- Visualization of differential methylation results

## Tools & Packages

Developed in **R** using **Bioconductor** packages:

- `minfi`  
- `gplots`  
- `qqman`  
- `IlluminaHumanMethylation450kmanifest`

## How to Run the Analysis

1. Clone the repository:
2. Place your raw IDAT files or methylation matrix into the working directory.
3. Open `analysis.R` in RStudio.
4. Install missing packages (if any).
5. Modify file paths and sample group labels as needed.
6. Run the script to execute the full pipeline.

## Visualizations Generated

The script will generate the following plots:

- Density plots of Beta and M values (raw vs normalized)  
- 6-panel quality control plot  
- PCA plots colored by group, sex, and Sentrix ID  
- Boxplots of p-values (raw, uncorrected, corrected)  
- Volcano plot of significant CpG sites  
- Manhattan plot of genome-wide results  
- Heatmap with hierarchical clustering of top differentially methylated probes


#Analysis of CpG Methylation Differences Between Healthy and Diseased States

This repository contains the R script developed for analyzing differential methylation between healthy and diseased samples using Illumina DNA methylation array data.

##Files:
analysis.R: Main R script performing the full methylation analysis pipeline

##Project Summary
DNA methylation analysis to detect CpG sites that show significant differences between healthy and diseased individuals. The steps in the pipeline include:

- Quality control of samples and probes

- Normalization using preprocessQuantile

- Filtering based on detection p-values

- Differential methylation analysis using t-tests

- Visualization of results (e.g., heatmaps, volcano plots)

##Tools & Packages:
R / RStudio

###Bioconductor packages:

- minfi

- gplots

- qqman

- IlluminaHumanMethylation450kmanifest

##Pipeline:

- Clone the repository

- Place your raw data

- Install required R packages and input files

- Run the analysis

##Visualizations

The script generates several plots, which will appear in your R graphics device:

- Density plots of raw Beta and M values for different groups.

- 6-panel QC plot comparing raw and normalized data.

- PCA plots colored by group, sex, and Sentrix ID.

- Boxplots of raw, uncorrected, and corrected p-values.

- Volcano plots highlighting significant probes.

- Manhattan plots showing genomic significance.

- Heatmaps with hierarchical clustering.

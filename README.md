# Bulk RNA-seq Data Analysis
This repository contains the input files and analysis code (Analysis.R) used for demonstrating the key steps involved in Bulk RNA-seq data analysis.
The aim is to provide an overview of the typical RNA-seq analysis workflow, starting from processed count data and proceeding to exploratory analysis, differential expression and pathway analysis.

## Repository Structure
1. script1_preprocessing_pipeline.sh: Code for downloading and processing of the raw FASTQ files.
2. script2_feature_counts.sh: Code for generating the feature counts matrix from STAR BAM files.
3. script2_analysis.R:
   Contains the code required for performing the following:
   1. Principal Component Analysis (PCA)
   2. Differential Expression Analysis (DEG)
   3. Pathway Enrichment Analysis: Over Representation Analysis (ORA), Gene Set Variation Analysis (GSVA), Gene Set Enrichment       Analysis (GSEA)
4. metadata.csv: This is the file containing metadata information for the samples used in the analysis.
5. raw_counts.csv: This contains the raw counts matrix for all the samples used in the analysis.
6. hallmark.gmt: This is the Hallmark gene set file obtained from the MSigDB database for pathway analysis.

## References:
1. https://hbctraining.github.io/rnaseq-cb321/lessons/analysis_methods.html
2. https://silicogene.com/blog/bulk-rna-seq-data-analysis/

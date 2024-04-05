# He_Larval_scRNA-seq

## Overview

This repository contains all the code used to run the analyses presented in ____.

## Required packages


## Analysis pipeline

The code for these analyses is separated across several files and repositories, many of which depend on the others. The following is the suggested order for performing the analyses.

1) Perform the initial scRNA-seq analysis steps to generate the Seurat object for the HE dataset.
    a) Run "timepoint-integration_and_scTransform.R" to integrate the single cell gene expression counts tables from each of the sample time points.
    b) Run "umap_and_clustering.Rmd" to obtain the UMAP projection of the dataset and assign cell type identities.
2) Perform the broader analyses needed for assessing specific cell types in downstream steps.
    a) Run "run-interproscan.sh" (on a compute cluster) to assign GO annotations to protein coding genes in the He genome.
    b) Run "interproscan-output-processing.Rmd" to convert the InterProScan results into useable output for later analyses that use GO terms.
    c) Run "tf-list-curation.Rmd" to generate a candidate list of transcription factor genes in the He genome.
    d) Run the Waddington-OT pipeline in the "waddington-ot" folder (see the separate README.md file). The results of this are used for analyzing the undifferentiated cell cluster.
3) Perform more detailed analyses of specific larval cell types
    a) Run "skeletogenic-cell-analysis.Rmd" to analyze the skeletogenic cell lineage.
    b) Run "undifferentiated-cluster-analysis.Rmd" to analyze the undifferentiated cell cluster.
    c) Run "neural-analysis.Rmd" to analyze the neural cluster.
4) Perform the analyses used for analyzing GRN and transcription factor gene expression profiles across the scRNA-seq time course.
    a) Run "mfuzz-pipeline.Rmd" to perform the Mfuzz analysis of the pseudobulked scRNA-seq data.
    b) Run "mfuzz-go-analysis.Rmd" to compare enriched biological processes between different Mfuzz expression profiles.
    c) Run "mfuzz-grn-analysis.Rmd" to analyze the expression of sea urchin embryonic GRN genes across different stages of development.
    d) Run "candidate-tf-analysis.Rmd" to identify candidate transcription factors that may regulate adult rudiment development in He.







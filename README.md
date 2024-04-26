# He_Larval_scRNA-seq

## Overview

This repository contains all the code used to run the analyses presented in "Contrasting the development of larval and adult body plans during the evolution of biphasic lifecycles in sea urchins" (2024).


## Required packages

The following packages and versions were used for analyses in R versions 4.2.0 and 4.3.0, depending on the package availability in each version of R:
- Seurat v4.3.0 (R v4.2.0)
- SeuratObject v4.1.3 (R v4.2.0)
- sctransform v0.3.5 (R v4.2.0)
- AnnotationForge v1.40.2 (R v4.2.0)
- Mfuzz v2.58.0 (R v4.2.0)
- DESeq2 v1.38.0 (R v4.2.0)
- SingleCellExperiment v1.20.1 (R v4.2.0)
- tidyverse v2.0.0 (R v4.2.0 and v4.3.0)
- scCustomize v2.0.1 (R v4.3.0)
- clusterProfiler v4.10.0 (R v4.3.0)
- enrichplot v1.22.0 (R v4.3.0)
- ComplexHeatmap v2.18.0 (R v4.3.0)
- ape v5.7-1 (R v4.3.0)

For the Waddington-OT pipeline, wot v1.0.8 was run in Python v3.8.

For the gene function analysis, InterProScan v5.64-96.0 was run on the Duke Compute Cluster.


## Analysis pipeline

The code for these analyses is separated across several files and repositories, many of which depend on the others. The following is the suggested order for performing the analyses.

1) Perform the initial scRNA-seq analysis steps to generate the Seurat object for the HE dataset.

    a) Run "1.1_timepoint-integration_and_scTransform.R" to integrate the single cell gene expression counts tables from each of the sample time points. This script creates the Seurat object used for all downstream analyses. The counts tables for all timepoints are combined into the same object, and are filtered and standardized using the scTransform pipeline. The final output will be a Seurat object ready for clustering and UMAP projection.
   
    b) Run "1.2_umap_and_clustering.Rmd" to obtain the UMAP projection of the dataset and assign cell type identities. This script contains the code used to run the UMAP and clustering algorithms on the Seurat object containing the normalized data for all 12 timepoints in this study. It also contains code for identifying the marker genes for each cluster, assigning cell type identities to clusters, and generating figures used in the main paper (portions of Figure 2, S1, and S2, as well as Table S1).

2) Perform the broader analyses needed for assessing specific cell types in downstream steps.

    a) Run InterProScan on a compute cluster to assign GO annotations to protein coding genes in the He genome. The folder "/data-files/" contains the file with the peptide models for He ("Hery_peptide_models.fasta"). The following command was run on the Duke Compute Cluster to perform this analysis:

   ./interproscan.sh -i ./Hery_peptide_models.fasta -goterms -cpu 16
   
    b) Run "2.1_interproscan-output-processing.Rmd" to convert the InterProScan results into useable output for later analyses that use GO terms. This script provides the code needed to convert the gene function annotation output (of the protein-coding genes in the He genome) from InterProScan into an `orgDb` package, using the tools provided by the `AnnotationForge` package. This will make it easy to filter the gene function annotation output for GO terms of interest and to perform GO over representation analyses.
   
    c) Run "2.2_tf-list-curation.Rmd". This script provides the code needed to generate the list of candidate transcription factors found in the He genome, based on GO term annotations and a previously used list.
   
    d) Run the Waddington-OT pipeline in the "waddington-ot" folder (see the separate README.md file). The results of this are used for analyzing the undifferentiated cell cluster.

4) Perform more detailed analyses of specific larval cell types

    a) Run "3.1_skeletogenic-cell-analysis.Rmd". This script provides the ability to run the analyses of the skeletogenic cell lineage in the scRNA-seq dataset. The script contains the code needed to generate all the plots used in Figures 3 and 4.
   
    b) Run "3.2_undifferentiated-cluster-analysis.Rmd". This script provides the ability to run the analyses of the undifferentiated cell cluster in the He scRNA-seq dataset. The script contains the code needed to generate all the plots used in Figure 5.
   
    c) Run "3.3_neural-analysis.Rmd". This script provides the ability to run the analyses of the neural cell lineage in the He scRNA-seq dataset. The script contains the code needed to generate all the plots used in Figure 6.
   
5) Perform the analyses used for analyzing GRN and transcription factor gene expression profiles across the scRNA-seq time course.
   
    a) Run "4.1_mfuzz-pipeline.Rmd". This script provides the code needed to perform fuzzy c-means clustering of temporal gene expression profiles for genes in the He scRNA-seq dataset using the package `Mfuzz`.
   
    b) Run "4.2_mfuzz-go-analysis.Rmd" to compare enriched biological processes between different Mfuzz expression profiles. This script provides the code needed to perform the Gene Ontology (GO) over representation analysis for the genes expressed in the different Mfuzz clusters (which requires running "mfuzz-pipeline.Rmd" first). The plots generated in this script appear in Figure 7.
   
    c) Run "4.3_mfuzz-grn-analysis.Rmd". This script provides the code needed to analyze the temporal expression profiles of sea urchin embryonic GRN genes, using the results from the Mfuzz fuzzy c-means clustering pipeline. The code needed to generate the plots for Figure 8A-B.
   
    d) Run "4.4_candidate-tf-analysis.Rmd" to identify candidate transcription factors that may regulate adult rudiment development in He. This script provides the code needed to analyze the expression profiles (using the outputs from the Mfuzz pipeline) for the transcription factors in the list generated by the "2.2_tf-list-curation.Rmd" script. This data will be used to generate a candidate list of TFs that may regulate adult rudiment development in He. The plots generated by this script are used in Figure 8C and 8E.







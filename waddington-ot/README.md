# Waddington-OT analysis

This folder contains the code needed to run the Waddington-OT pipeline on the He scRNA-seq data.

In order to run the pipeline, the He scRNA-seq Seurat object will firt have to be created and assigned cell type annotations (performed in "timepoint-integration_and_scTransform.R" and "umap_and_clutering.Rmd").

Once the Seurat object is created, run "data-preparation/seurat-export-for-wot.R" to generate the input files for the Waddington-OT pipeline. These should all appear in the "waddington-ot/data" folder.

Then the Waddington-OT analysis can be performed by running the .ipynb files in the "waddington-ot/" folder in numerical order (01-04). The final output of this pipeline will be the transition table used in Figure 5 during the analysis of the undifferentiated cell cluster.

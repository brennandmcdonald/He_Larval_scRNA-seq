# Combining timepoints into a single object for downstream clustering and UMAP projection
# Filters out ribosomal genes from variable genes used for clustering

# load packages
library(tidyverse)
library(Seurat)
library(sctransform)
Sys.setenv(VROOM_CONNECTION_SIZE=500072)


# load the raw data
he.6.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He6hpf_anno.csv", row.names = 1)
he.9.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He9hpf_anno.csv", row.names = 1)
he.12.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He12hpf_anno.csv", row.names = 1)
he.16.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He16hpf_anno.csv", row.names = 1)
he.20.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He20hpf_anno.csv", row.names = 1)
he.24.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He24hpf_anno.csv", row.names = 1)
he.30.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He30hpf_anno.csv", row.names = 1)
he.36.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He36hpf_anno.csv", row.names = 1)
he.42.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He42hpf_anno.csv", row.names = 1)
he.48.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He48hpf_anno.csv", row.names = 1)
he.54.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He54hpf_anno.csv", row.names = 1)
he.60.csv <- read.csv("/work/bdm50/scRNA-seq/csv_anno/He60hpf_anno.csv", row.names = 1)


he.6 <- CreateSeuratObject(counts = he.6.csv, project = "6hpf", 
                           row.names = gene.names.6, min.cells = 3, min.features = 200)
he.9 <- CreateSeuratObject(counts = he.9.csv, project = "9hpf", 
                           row.names = gene.names.9, min.cells = 3, min.features = 200)
he.12 <- CreateSeuratObject(counts = he.12.csv, project = "12hpf", 
                            row.names = gene.names.12, min.cells = 3, min.features = 200)
he.16 <- CreateSeuratObject(counts = he.16.csv, project = "16hpf", 
                            row.names = gene.names.16, min.cells = 3, min.features = 200)
he.20 <- CreateSeuratObject(counts = he.20.csv, project = "20hpf", 
                            row.names = gene.names.20, min.cells = 3, min.features = 200)
he.24 <- CreateSeuratObject(counts = he.24.csv, project = "24hpf", 
                            row.names = gene.names.24, min.cells = 3, min.features = 200)
he.30 <- CreateSeuratObject(counts = he.30.csv, project = "30hpf", 
                            row.names = gene.names.30, min.cells = 3, min.features = 200)
he.36 <- CreateSeuratObject(counts = he.36.csv, project = "36hpf", 
                            row.names = gene.names.36, min.cells = 3, min.features = 200)
he.42 <- CreateSeuratObject(counts = he.42.csv, project = "42hpf", 
                            row.names = gene.names.42, min.cells = 3, min.features = 200)
he.48 <- CreateSeuratObject(counts = he.48.csv, project = "48hpf", 
                            row.names = gene.names.48, min.cells = 3, min.features = 200)
he.54 <- CreateSeuratObject(counts = he.54.csv, project = "54hpf", 
                            row.names = gene.names.54, min.cells = 3, min.features = 200)
he.60 <- CreateSeuratObject(counts = he.60.csv, project = "60hpf", 
                            row.names = gene.names.60, min.cells = 3, min.features = 200)


he.6$group <- "6hpf"
he.9$group <- "9hpf"
he.12$group <- "12hpf"
he.16$group <- "16hpf"
he.20$group <- "20hpf"
he.24$group <- "24hpf"
he.30$group <- "30hpf"
he.36$group <- "36hpf"
he.42$group <- "42hpf"
he.48$group <- "48hpf"
he.54$group <- "54hpf"
he.60$group <- "60hpf"


# merge the datasets
he.integrated <- merge(he.6, y = c(he.9, he.12, he.16, he.20, he.24, he.30, he.36, he.42, he.48, he.54, he.60), 
                       add.cell.ids = c("6hpf", "9hpf", "12hpf", "16hpf", "20hpf", "24hpf", "30hpf", "36hpf", "42hpf", "48hpf", "54hpf", "60hpf"), 
                       project = "HE_6-60hpf_Integrated")

he.integrated <- PercentageFeatureSet(he.integrated, pattern = "\\b\\w*Rp[sl]\\w*\\b", col.name = "percent.Rb")

he.integrated$group <- factor(he.integrated$group, levels = c("6hpf", "9hpf", "12hpf", "16hpf", "20hpf", "24hpf", "30hpf", "36hpf", "42hpf", "48hpf", "54hpf", "60hpf"))

Idents(he.integrated) <- "Stage"


# QC and selecting cells for further analysis
RNA_nFeature_Count_vlnplot <- VlnPlot(he.integrated, features = c("nFeature_RNA", "nCount_RNA"))
RNA_nFeature_Count_scatter <- FeatureScatter(he.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

save(RNA_nFeature_Count_vlnplot, file = "/work/bdm50/scRNA-seq/plots/RNA_nFeature_Count_vlnplot.Rda")
save(RNA_nFeature_Count_scatter, file = "/work/bdm50/scRNA-seq/plots/RNA_nFeature_Count_scatter.Rda")

he.integrated <- subset(he.integrated, subset = nFeature_RNA > 200 & nCount_RNA < 10000 & nFeature_RNA < 4000)


# Normalize and scale the data
he.integrated <- SCTransform(he.integrated, method = "glmGamPoi", verbose = T, variable.features.n = 6000 , vars.to.regress = "percent.Rb" )



# Run PCA
he.integrated <- RunPCA(object = he.integrated, npcs = 200, features = VariableFeatures(object = he.integrated))

pca_VizDim <- VizDimLoadings(he.integrated, dims = 1:4, reduction = "pca")
pca_DimPlot <- DimPlot(he.integrated, reduction = "pca")
pca_DimHeatmap <- DimHeatmap(he.integrated, dims = 1:4, cells = 500, balanced = T)

save(he.integrated, file = "/work/bdm50/scRNA-seq/HE_6-60hpf_integrated-SCT6k.Rda")


save(pca_VizDim, file = "/work/bdm50/scRNA-seq/plots/pca_VizDim.Rda")
save(pca_DimPlot, file = "/work/bdm50/scRNA-seq/plots/pca_DimPlot.Rda")
save(pca_DimHeatmap, file = "/work/bdm50/scRNA-seq/plots/pca_DimHeatmap.Rda")



# Determine the dimensionality of the dataset
elbow_plot <- ElbowPlot(he.integrated, ndims = 200)

save(elbow_plot, file = "/work/bdm50/scRNA-seq/plots/elbow_plot.Rda")






save(he.integrated, file = "/work/bdm50/scRNA-seq/HE_6-60hpf_integrated-SCT6k.Rda")


# This script creates the Seurat object used for all downstream analyses. 
# The counts tables for all timepoints are combined into the same object,
# and are filtered and standardized using the scTransform pipeline.
# The final output will be a Seurat object ready for clustering and UMAP projection.

# load packages
library(tidyverse)
library(Seurat)
library(sctransform)
Sys.setenv(VROOM_CONNECTION_SIZE=500072)


# load the raw data
he.6.csv <- read.csv("./counts-tables/He6hpf_anno.csv", row.names = 1)
he.9.csv <- read.csv("./counts-tables/He9hpf_anno.csv", row.names = 1)
he.12.csv <- read.csv("./counts-tables/He12hpf_anno.csv", row.names = 1)
he.16.csv <- read.csv("./counts-tables/He16hpf_anno.csv", row.names = 1)
he.20.csv <- read.csv("./counts-tables/He20hpf_anno.csv", row.names = 1)
he.24.csv <- read.csv("./counts-tables/He24hpf_anno.csv", row.names = 1)
he.30.csv <- read.csv("./counts-tables/He30hpf_anno.csv", row.names = 1)
he.36.csv <- read.csv("./counts-tables/He36hpf_anno.csv", row.names = 1)
he.42.csv <- read.csv("./counts-tables/He42hpf_anno.csv", row.names = 1)
he.48.csv <- read.csv("./counts-tables/He48hpf_anno.csv", row.names = 1)
he.54.csv <- read.csv("./counts-tables/He54hpf_anno.csv", row.names = 1)
he.60.csv <- read.csv("./counts-tables/He60hpf_anno.csv", row.names = 1)

# create separate Seurat objects for each timepoint
he.6 <- CreateSeuratObject(counts = he.6.csv, project = "6hpf", min.cells = 3, min.features = 200)
he.9 <- CreateSeuratObject(counts = he.9.csv, project = "9hpf", min.cells = 3, min.features = 200)
he.12 <- CreateSeuratObject(counts = he.12.csv, project = "12hpf", min.cells = 3, min.features = 200)
he.16 <- CreateSeuratObject(counts = he.16.csv, project = "16hpf", min.cells = 3, min.features = 200)
he.20 <- CreateSeuratObject(counts = he.20.csv, project = "20hpf", min.cells = 3, min.features = 200)
he.24 <- CreateSeuratObject(counts = he.24.csv, project = "24hpf", min.cells = 3, min.features = 200)
he.30 <- CreateSeuratObject(counts = he.30.csv, project = "30hpf", min.cells = 3, min.features = 200)
he.36 <- CreateSeuratObject(counts = he.36.csv, project = "36hpf", min.cells = 3, min.features = 200)
he.42 <- CreateSeuratObject(counts = he.42.csv, project = "42hpf", min.cells = 3, min.features = 200)
he.48 <- CreateSeuratObject(counts = he.48.csv, project = "48hpf", min.cells = 3, min.features = 200)
he.54 <- CreateSeuratObject(counts = he.54.csv, project = "54hpf", min.cells = 3, min.features = 200)
he.60 <- CreateSeuratObject(counts = he.60.csv, project = "60hpf", min.cells = 3, min.features = 200)


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

# Create a new feature for ribosomal genes
he.integrated <- PercentageFeatureSet(he.integrated, pattern = "\\b\\w*Rp[sl]\\w*\\b", col.name = "percent.Rb")

he.integrated$group <- factor(he.integrated$group, levels = c("6hpf", "9hpf", "12hpf", "16hpf", "20hpf", "24hpf", "30hpf", "36hpf", "42hpf", "48hpf", "54hpf", "60hpf"))

Idents(he.integrated) <- "Stage"


# QC and selecting cells for further analysis
VlnPlot(he.integrated, features = c("nFeature_RNA", "nCount_RNA"))
FeatureScatter(he.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

he.integrated <- subset(he.integrated, subset = nFeature_RNA > 200 & nCount_RNA < 10000 & nFeature_RNA < 4000)


# Normalize and scale the data using SCTransform and the "glmGamPoi" method
# Regress out the effect of the ribosomal genes
he.integrated <- SCTransform(he.integrated, method = "glmGamPoi", verbose = T, variable.features.n = 6000 , vars.to.regress = "percent.Rb" )


# Run PCA
he.integrated <- RunPCA(object = he.integrated, npcs = 200, features = VariableFeatures(object = he.integrated))

VizDimLoadings(he.integrated, dims = 1:4, reduction = "pca")
DimPlot(he.integrated, reduction = "pca")
DimHeatmap(he.integrated, dims = 1:4, cells = 500, balanced = T)


# Determine the dimensionality of the dataset
ElbowPlot(he.integrated, ndims = 200)


# Save the Seurat object for the next steps
save(he.integrated, file = "HE_6-60hpf_integrated-SCT6k.Rda")


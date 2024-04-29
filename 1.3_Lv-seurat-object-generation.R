# This script creates the Lytechinus variegatus Seurat object used for all 
# downstream analyses. 

# The scRNA-seq data for this was retrieved from Massri et al. (2021).

# The counts tables for all timepoints are combined into the same object,
# and are filtered and standardized using the SCTransform pipeline.
# The Seurat object is then clustered and run through the UMAP algorithm.


# Load packages
library(Seurat)
library(tidyverse)
library(sctransform)
Sys.setenv(VROOM_CONNECTION_SIZE=500072)

# Load the raw data
Lv2hpf <- read.csv("/counts_tables/Lv-2hpf_anno.csv", row.names=1)
Lv3hpf <- read.csv("/counts_tables/Lv-3hpf_anno.csv", row.names=1)
Lv4hpf <- read.csv("/counts_tables/Lv-4hpf_anno.csv", row.names=1)
Lv5hpf <- read.csv("/counts_tables/Lv-5hpf_anno.csv", row.names=1)
Lv6hpf <- read.csv("/counts_tables/Lv-6hpf_anno.csv", row.names=1)
Lv7hpf <- read.csv("/counts_tables/Lv-7hpf_anno.csv", row.names=1)
Lv8hpf <- read.csv("/counts_tables/Lv-8hpf_anno.csv", row.names=1)
Lv9hpf <- read.csv("/counts_tables/Lv-9hpf_anno.csv", row.names=1)
Lv10hpf <- read.csv("/counts_tables/Lv-10hpf_anno.csv", row.names=1)
Lv11hpf <- read.csv("/counts_tables/Lv-11hpf_anno.csv", row.names=1)
Lv12hpf <- read.csv("/counts_tables/Lv-12hpf_anno.csv", row.names=1)
Lv13hpf <- read.csv("/counts_tables/Lv-13hpf_anno.csv", row.names=1)
Lv14hpf <- read.csv("/counts_tables/Lv-14hpf_anno.csv", row.names=1)
Lv15hpf <- read.csv("/counts_tables/Lv-15hpf_anno.csv", row.names=1)
Lv16hpf <- read.csv("/counts_tables/Lv-16hpf_anno.csv", row.names=1)
Lv18hpf <- read.csv("/counts_tables/Lv-18hpf_anno.csv", row.names=1)
Lv20hpf <- read.csv("/counts_tables/Lv-20hpf_anno.csv", row.names=1)
Lv24hpf <- read.csv("/counts_tables/Lv-24hpf_anno.csv", row.names=1)

# create separate Seurat objects for each timepoint
Lv2hpf_seurat <- CreateSeuratObject(counts = Lv2hpf, project = "Lv2hpf", min.cells = 3, min.features = 200)
Lv3hpf_seurat <- CreateSeuratObject(counts = Lv3hpf, project = "Lv3hpf", min.cells = 3, min.features = 200)
Lv4hpf_seurat <- CreateSeuratObject(counts = Lv4hpf, project = "Lv4hpf", min.cells = 3, min.features = 200)
Lv5hpf_seurat <- CreateSeuratObject(counts = Lv5hpf, project = "Lv5hpf", min.cells = 3, min.features = 200)
Lv6hpf_seurat <- CreateSeuratObject(counts = Lv6hpf, project = "Lv6hpf", min.cells = 3, min.features = 200)
Lv7hpf_seurat <- CreateSeuratObject(counts = Lv7hpf, project = "Lv7hpf", min.cells = 3, min.features = 200)
Lv8hpf_seurat <- CreateSeuratObject(counts = Lv8hpf, project = "Lv8hpf", min.cells = 3, min.features = 200)
Lv9hpf_seurat <- CreateSeuratObject(counts = Lv9hpf, project = "Lv9hpf", min.cells = 3, min.features = 200)
Lv10hpf_seurat <- CreateSeuratObject(counts = Lv10hpf, project = "Lv10hpf", min.cells = 3, min.features = 200)
Lv11hpf_seurat <- CreateSeuratObject(counts = Lv11hpf, project = "Lv11hpf", min.cells = 3, min.features = 200)
Lv12hpf_seurat <- CreateSeuratObject(counts = Lv12hpf, project = "Lv12hpf", min.cells = 3, min.features = 200)
Lv13hpf_seurat <- CreateSeuratObject(counts = Lv13hpf, project = "Lv13hpf", min.cells = 3, min.features = 200)
Lv14hpf_seurat <- CreateSeuratObject(counts = Lv14hpf, project = "Lv14hpf", min.cells = 3, min.features = 200)
Lv15hpf_seurat <- CreateSeuratObject(counts = Lv15hpf, project = "Lv15hpf", min.cells = 3, min.features = 200)
Lv16hpf_seurat <- CreateSeuratObject(counts = Lv16hpf, project = "Lv16hpf", min.cells = 3, min.features = 200)
Lv18hpf_seurat <- CreateSeuratObject(counts = Lv18hpf, project = "Lv18hpf", min.cells = 3, min.features = 200)
Lv20hpf_seurat <- CreateSeuratObject(counts = Lv20hpf, project = "Lv20hpf", min.cells = 3, min.features = 200)
Lv24hpf_seurat <- CreateSeuratObject(counts = Lv24hpf, project = "Lv24hpf", min.cells = 3, min.features = 200)

# ADDING STAGE ANNOTATIONS FOR EACH TIME SERIES LIBRARY
Lv2hpf_seurat$Stage <- "2_HPF"
Lv3hpf_seurat$Stage <- "3_HPF"
Lv4hpf_seurat$Stage <- "4_HPF"
Lv5hpf_seurat$Stage <- "5_HPF"
Lv6hpf_seurat$Stage <- "6_HPF"
Lv7hpf_seurat$Stage <- "7_HPF"
Lv8hpf_seurat$Stage <- "8_HPF"
Lv9hpf_seurat$Stage <- "9_HPF"
Lv10hpf_seurat$Stage <- "10_HPF"
Lv11hpf_seurat$Stage <- "11_HPF"
Lv12hpf_seurat$Stage <- "12_HPF"
Lv13hpf_seurat$Stage <- "13_HPF"
Lv14hpf_seurat$Stage <- "14_HPF"
Lv15hpf_seurat$Stage <- "15_HPF"
Lv16hpf_seurat$Stage <- "16_HPF"
Lv18hpf_seurat$Stage <- "18_HPF"
Lv20hpf_seurat$Stage <- "20_HPF"
Lv24hpf_seurat$Stage <- "24_HPF"

# merge the datasets
control_2to24_raw <-merge(Lv2hpf_seurat, y= c(Lv3hpf_seurat, Lv4hpf_seurat, 
                                              Lv5hpf_seurat, Lv6hpf_seurat, 
                                              Lv7hpf_seurat, Lv8hpf_seurat, 
                                              Lv9hpf_seurat, Lv10hpf_seurat,
                                              Lv11hpf_seurat, Lv12hpf_seurat,
                                              Lv13hpf_seurat, Lv14hpf_seurat,
                                              Lv15hpf_seurat, Lv16hpf_seurat,
                                              Lv18hpf_seurat, Lv20hpf_seurat,
                                              Lv24hpf_seurat),
                                add.cell.ids = c("2_HPF", "3_HPF", "4_HPF" ,
                                                 "5_HPF", "6_HPF", "7_HPF",
                                                 "8_HPF", "9_HPF", "10_HPF",
                                                 "11_HPF", "12_HPF", "13_HPF",
                                                 "14_HPF", "15_HPF", "16_HPF",
                                                 "18_HPF", "20_HPF", "24_HPF"),
                                project = "Control")

# Put custom order to the groups
control_2to24_raw$Stage <- factor(control_2to24_raw$Stage,  
                                  levels = c('2_HPF', '3_HPF','4_HPF','5_HPF','6_HPF',
								                  '7_HPF','8_HPF','9_HPF', '10_HPF','11_HPF', '12_HPF',
								                  '13_HPF','14_HPF','15_HPF','16_HPF','18_HPF','20_HPF',
								                  '24_HPF'))

control_2to24_raw$orig.ident <- factor(control_2to24_raw$orig.ident,  
                                       levels = c('Lv2hpf', 'Lv3hpf','Lv4hpf','Lv5hpf','Lv6hpf',
                                       'Lv7hpf','Lv8hpf','Lv9hpf', 'Lv10hpf','Lv11hpf', 'Lv12hpf',
                                       'Lv13hpf','Lv14hpf','Lv15hpf','Lv16hpf','Lv18hpf','Lv20hpf',
                                       'Lv24hpf'))

control_2to24_raw$experiment <- "CONTROL"


# Create a new feature for ribosomal genes
control_2to24_raw <- PercentageFeatureSet(control_2to24_raw, pattern = "\\b\\w*Rp[sl]\\w*\\b", col.name = "percent.Rb")

# QC and selecting cells for further analysis
VlnPlot(control_2to24_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.Rb"), ncol = 3)

control_2to24_raw <- subset(control_2to24_raw, subset = nFeature_RNA > 200 & nCount_RNA < 10000 )

VlnPlot(control_2to24_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.Rb"), ncol = 3)


# Normalize and scale the data using SCTransform
# Regress out the effect of the ribosomal genes
lv.integrated <- SCTransform(control_2to24_raw, vst.flavor = "v2", verbose = T, variable.features.n = 2000,  vars.to.regress = "percent.Rb")

# Run PCA
lv.integrated <- RunPCA(lv.integrated,  npcs = 200, verbose = T)

# Non-linear dimensional reduction (UMAP)
lv.integrated <- RunUMAP(lv.integrated,  dims = 1:185)

# Perform clustering with several resolution values
lv.integrated <- FindNeighbors(lv.integrated, dims = 1:185)
lv.integrated <- FindClusters(lv.integrated, resolution = 3)
lv.integrated <- FindClusters(lv.integrated, resolution = 2)
lv.integrated <- FindClusters(lv.integrated, resolution = 1)


# Save the Seurat object for downstream analyses
save(lv.integrated, file = "LV_2-24hpf_integrated-SCT2k_185dim.Rda")




# SAVING CONTROL DATA FOR WOT

library(Matrix)

getwd()

setwd("./WaddingtonOT")

dir.create("data")

DATA_PATH <- "WaddingtonOT_Full/data/"

list.files("data")

#R_FILE_PATH <- paste(DATA_PATH, "Micromereless_dataset_single_cell.Rda", sep="")

MATRIX_SAVE_PATH <- paste(DATA_PATH, "HERY_SCT.mtx", sep="/")
ANNO_SAVE_PATH <- paste(DATA_PATH, "HERY_seurat_anno.csv", sep="/")
VAR_SAVE_PATH <- paste(DATA_PATH, "HERY_var.csv", sep="/")
UMAP_SAVE_PATH <- paste(DATA_PATH, "HERY_umap.csv", sep="/")

data <- he.integrated

getwd()



writeMM(data@assays$SCT@data, file=MATRIX_SAVE_PATH)
write.csv(data@meta.data, ANNO_SAVE_PATH)
write.csv(data@assays$SCT@meta.features, VAR_SAVE_PATH)
write.csv(data[["umap"]]@cell.embeddings, UMAP_SAVE_PATH)


head(he.integrated@meta.data$celltype)

write.table(data.frame( rownames(he.integrated@meta.data)  , he.integrated@meta.data$celltype), file ="cell_annotations_hery.txt", row.names=F, col.names=T, quote=F,sep='\t')

Idents(he.integrated) <- "celltype"

tibble(cellID = colnames(he.integrated), clusterID = Idents(he.integrated)) %>%
  # write out to a file, with today's date
  write_csv(file = sprintf("WaddingtonOT_Full/cell_annotations_hery_%s.csv", Sys.Date()))



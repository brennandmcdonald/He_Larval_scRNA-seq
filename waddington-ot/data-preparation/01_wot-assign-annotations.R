# Annotate cluster identities in the 6-60hr integrated Seurat object.
# Then export those identities in a format for WaddingtonOT


CreateSeuratObject(load(
  file = "6-60hr\ Integrated\ UMAP/HE_6-60hpf_integrated-SCT6k_195dim.Rda", 
  verbose = T))

head(he.integrated@meta.data)

Idents(he.integrated) <- "seurat_clusters"


#### Annotations updated 01/12/2024
new.cluster.ids <- c("Animal Pole Domain",
                     "Animal Pole Domain",
                     "Larval Ectoderm",
                     "Larval Ectoderm",
                     "Larval Ectoderm",
                     "Pigment",
                     "Larval Ectoderm",
                     "Ciliary Band",
                     "Early Ectoderm",
                     "Larval Ectoderm",
                     "Animal Pole Domain",
                     "Ciliary Band",
                     "Vegetal Ectoderm",
                     "Larval Ectoderm",
                     "Larval Ectoderm",
                     "Ciliary Band",
                     "Larval Ectoderm",
                     "Early Endomesoderm",
                     "Rudiment Endomesoderm",
                     "Gut",
                     "Left Coelomic Pouch",
                     "Animal Pole Domain",
                     "Blastocoelar",
                     "Vegetal Ectoderm",
                     "Vestibular Ectoderm",
                     "Vestibular Ectoderm",
                     "Blastocoelar",
                     "Left Coelomic Pouch",
                     "Early Endomesoderm",
                     "Pigment",
                     "Animal Pole Domain",
                     "Neural",
                     "Rudiment Endomesoderm",
                     "Vegetal Ectoderm",
                     "Larval Ectoderm",
                     "Unknown",
                     "Vegetal Ectoderm",
                     "Ciliary Band",
                     "Early Mesoderm",
                     "Neural",
                     "Undifferentiated",
                     "Right Coelomic Pouch",
                     "Animal Pole Domain",
                     "Left Coelomic Pouch",
                     "Pigment",
                     "Left Coelomic Pouch",
                     "Larval Ectoderm",
                     "Rudiment Endomesoderm",
                     "Early Endomesoderm",
                     "Skeletal",
                     "Vestibular Ectoderm",
                     "Neural",
                     "Unknown",
                     "Neural",
                     "Neural",
                     "Unknown",
                     "Neural")

names(new.cluster.ids) <- levels(he.integrated)
he.integrated <- RenameIdents(he.integrated, new.cluster.ids)

Idents(he.integrated)

#save(he.integrated, file = "6-60hr\ Integrated\ UMAP/HE_6-60hpf_integrated-SCT6k_195dim.Rda")

# FOR Hery DATASET CELLS
head(he.integrated@meta.data)
Idents(he.integrated)
cellids <-  as.data.frame(Idents(he.integrated))
head(cellids)
write.table(cellids, file ="WaddingtonOT_Full/cell_annotations_He_6_60hpf_v2.txt", row.names=T, col.names=T, quote=F,sep='\t')


he.integrated <- AddMetaData(
  object = he.integrated,
  metadata = Idents(he.integrated),
  col.name = "celltype"
)

head(he.integrated@meta.data)
getwd()
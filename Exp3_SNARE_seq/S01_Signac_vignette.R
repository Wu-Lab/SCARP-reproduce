library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
set.seed(1234)

# ===========================================================================================
# This vignette was downloaded from https://satijalab.org/signac/articles/snareseq.html
# ===========================================================================================

# ================================================================
#               data  loading
# ================================================================

# load processed data matrices for each assay
rna <- Read10X("./Raw data/GSE126074_AdBrainCortex_rna/", gene.column = 1)
atac <- Read10X("./Raw data/GSE126074_AdBrainCortex_atac/", gene.column = 1)
fragments <- "./Raw data/GSE126074_AdBrainCortex_atac/fragments.sort.bed.gz"

# create a Seurat object and add the assays
snare <- CreateSeuratObject(counts = rna)

snare[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = fragments
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(snare[["ATAC"]]) <- annotations




# ================================================================
#              Quality control
# ================================================================

DefaultAssay(snare) <- "ATAC"
snare <- TSSEnrichment(snare)
snare <- NucleosomeSignal(snare)
snare$blacklist_fraction <- FractionCountsInRegion(
  object = snare,
  assay = 'ATAC',
  regions = blacklist_mm10
)

Idents(snare) <- "all"  # group all cells together, rather than by replicate

VlnPlot(
  snare,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 5
)



snare <- subset(
  x = snare,
  subset = blacklist_fraction < 0.03 &
    TSS.enrichment < 20 &
    nCount_RNA > 800 &
    nCount_ATAC > 500
)

VlnPlot(
  snare,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 5
)



# ================================================================
#         Gene expression data processing
# ================================================================

DefaultAssay(snare) <- "RNA"

snare <- FindVariableFeatures(snare, nfeatures = 3000)
snare <- NormalizeData(snare)
snare <- ScaleData(snare)
snare <- RunPCA(snare, npcs = 30)
snare <- RunUMAP(snare, dims = 1:30, 
                 reduction.name = "umap.rna")
snare <- FindNeighbors(snare, dims = 1:30)
snare <- FindClusters(snare, resolution = 0.5, algorithm = 3)



# ================================================================
#               Integration with scRNA-seq data
# ================================================================
# label transfer from Allen brain
allen <- readRDS("./Raw data/GSE126074_AdBrainCortex_rna/allen_brain.rds")
allen <- UpdateSeuratObject(object =  allen)

# use the RNA assay in the SNARE-seq data for integration with scRNA-seq
DefaultAssay(snare) <- 'RNA'

transfer.anchors <- FindTransferAnchors(
  reference = allen,
  query = snare,
  dims = 1:30,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen$subclass,
  weight.reduction = snare[['pca']],
  dims = 1:30
)

snare <- AddMetaData(object = snare, metadata = predicted.labels)


# label clusters based on predicted ID
new.cluster.ids <- c(
  "L2/3 IT",
  "L4",
  "L6 IT",
  "L5 CT",
  "L4",
  "L5 PT",
  "Pvalb",
  "Sst",
  "Astro",
  "Oligo",
  "Vip/Lamp5",
  "L6 IT.2",
  "L6b",
  "NP"
)
names(x = new.cluster.ids) <- levels(x = snare)
snare <- RenameIdents(object = snare, new.cluster.ids)
snare$celltype <- Idents(snare)

DimPlot(snare, group.by = 'celltype', label = TRUE, reduction = 'umap.rna') + 
  ggtitle('RNA(Signac)')

# 600*400


# ===================================================
#     computing silhouette
# ===================================================
library(factoextra)
library(cluster)


SNARE_Celltype <- data.frame(snare@meta.data[["celltype"]])
rownames(SNARE_Celltype) <- rownames(snare@meta.data)
colnames(SNARE_Celltype) <- 'celltype'

help_id <-  data.frame(1:length(unique(SNARE_Celltype$celltype)))
rownames(help_id) <- unique(SNARE_Celltype$celltype)
temp_name <- c()
for (i in 1:dim(SNARE_Celltype)[1]) {
  temp_name <-
    c(temp_name, help_id[SNARE_Celltype$celltype[i], ])
}
SNARE_Celltype$id <- temp_name
rm(temp_name, i, help_id)

x <- SNARE_Celltype$id
names(x) <- rownames(SNARE_Celltype)

sil <- silhouette(x, dist(snare@reductions[['umap.rna']]@cell.embeddings))
fviz_silhouette(sil)

# 600*400



# ==============================================
# save for SCARP in python

write.csv(SNARE_Celltype, './Processed data/SNARE_Celltype.csv')






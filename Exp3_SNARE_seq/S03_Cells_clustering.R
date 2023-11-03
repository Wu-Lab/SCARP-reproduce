library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
set.seed(1234)



# =====================================================================
#                     import  my  data and filter data
# =====================================================================
SCARP_ATAC_Cells_df <-
  read.table(
    './Results/SCARP_ATAC_Cells_df.txt',
    sep = ',',
    row.names = 1,
    header = T
  )

# load processed data matrices for each assay
rna <- Read10X("./Raw data/GSE126074_AdBrainCortex_rna/", gene.column = 1)
atac <- Read10X("./Raw data/GSE126074_AdBrainCortex_atac/", gene.column = 1)
fragments <- "./Raw data/GSE126074_AdBrainCortex_atac/fragments.sort.bed.gz"

# create a Seurat object and add the assays
snare <- CreateSeuratObject(counts = rna)

# create a Seurat object and add the assays
snare[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = fragments
)

snare <- snare[,rownames(SCARP_ATAC_Cells_df)]

# add cell type annotation
SNARE_cell_type <- read.csv('./Processed data/SNARE_Celltype.csv', row.names = 1)

snare$standard_Seurat_lablled <-
  factor(SNARE_cell_type[rownames(snare@meta.data),'celltype'],
            levels = c(
              "L2/3 IT",
              "L4",
              "L6 IT",
              "L5 CT",
              "L5 PT",
              "Pvalb",
              "Sst",
              "Astro",
              "Oligo",
              "Vip/Lamp5",
              "L6 IT.2",
              "L6b",
              "NP"))





# ===========================================================================
#                                 SCARP
# ============================================================================
snare[["SCARP_feat"]] <- CreateDimReducObject(
  embeddings = as.matrix(SCARP_ATAC_Cells_df),
  assay = 'ATAC',
  key = 'scarp_',
)

choose_dim <- dim(SCARP_ATAC_Cells_df)[2]

snare <- FindNeighbors(snare,
                       reduction = "SCARP_feat",
                       dims = 1:choose_dim)

snare <- RunUMAP(
  snare,
  reduction = 'SCARP_feat',
  dims = 1:choose_dim,
  reduction.name = "umap.atac.scarp"
)


DimPlot(snare,
        reduction = 'umap.atac.scarp',
        label = F,
        group.by = 'standard_Seurat_lablled') +
  ggtitle('ATAC(SCARP)')


library(loomR)
sdata.loom <- as.loom(x = snare, filename = "./sdata_SCARP.loom", verbose = FALSE)
# Always remember to close loom files when done
sdata.loom$close_all()



# ===================================================
#     computing silhouette
# ===================================================
library(factoextra)
library(cluster)

x <- SNARE_cell_type$id
names(x) <- rownames(snare@meta.data)

sil <- silhouette(x, dist(snare@reductions[['umap.atac.scarp']]@cell.embeddings))
fviz_silhouette(sil)

# 600*400

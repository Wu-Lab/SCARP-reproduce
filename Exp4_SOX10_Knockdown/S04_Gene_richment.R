library(clusterProfiler)
library(org.Hs.eg.db)
# library(topGO)
library(DOSE)
# library(doseplot)
library(patchwork)
library(enrichplot)
# library(ReactomePA)
library(DO.db)
library(forcats)
library(ggstance)
library(ggnewscale)
library(ggplot2)

# ===============================================================================
diff_genes <- read.csv('./Results/cytoscape/SCARP default node.csv')
diff_genes <- bitr(diff_genes$commonName, "SYMBOL", "ENSEMBL", "org.Hs.eg.db")


# ===============================================================================
back_genes <- read.csv("./Processed data/annotated_genes_1000bp.csv")
back_genes <- union(back_genes$ENSEMBL, diff_genes$ENSEMBL)
  

# ===============================================================================
CC <- enrichGO(gene = diff_genes$ENSEMBL,  #Gene list
               keyType = "ENSEMBL",  #Gene ID type
               OrgDb=org.Hs.eg.db,  #Species
               ont = "CC",
               pvalueCutoff = 1,  #pvalue threshold
               pAdjustMethod = "fdr",  #Multiple hypothesis testing correction method
               universe = back_genes,
               minGSSize = 1,   #Minimum set of genes to annotate
               maxGSSize = 500,  #Maximun set of genes to annotate
               qvalueCutoff = 1,
               readable = TRUE)  #Gene ID to gene name conversion

MF <- enrichGO(gene = diff_genes$ENSEMBL,  #Gene list
               keyType = "ENSEMBL",  #Gene ID type
               OrgDb=org.Hs.eg.db,  #Species
               ont = "MF",
               pvalueCutoff = 1,  #pvalue threshold
               pAdjustMethod = "fdr",  #Multiple hypothesis testing correction method
               universe = back_genes,
               minGSSize = 1,   #Minimum set of genes to annotate
               maxGSSize = 500,  #Maximun set of genes to annotate
               qvalueCutoff = 1,
               readable = TRUE)  #Gene ID to gene name conversion


BP <- enrichGO(gene = diff_genes$ENSEMBL,  #Gene list
               keyType = "ENSEMBL",  #Gene ID type
               OrgDb=org.Hs.eg.db,  #Species
               ont = "BP",
               pvalueCutoff = 1,  #pvalue threshold
               pAdjustMethod = "fdr",  #Multiple hypothesis testing correction method
               universe = back_genes,
               minGSSize = 1,   #Minimum set of genes to annotate
               maxGSSize = 500,  #Maximun set of genes to annotate
               qvalueCutoff = 1,
               readable = TRUE)  #Gene ID to gene name conversion

# keytypes(org.Hs.eg.db)
back_genes <- bitr(
  gene = back_genes,
  fromType = "ENSEMBL",
  toType = c('ENTREZID', 'SYMBOL'),
  OrgDb = 'org.Hs.eg.db'
)

diff_genes <- bitr(
  gene = diff_genes$ENSEMBL,
  fromType = "ENSEMBL",
  toType = c('ENTREZID', 'SYMBOL'),
  OrgDb = 'org.Hs.eg.db'
)


R.utils::setOption("clusterProfiler.download.method",'auto') 

KEGG <- enrichKEGG(gene = diff_genes$ENTREZID,
                 keyType = "kegg",
                 organism = 'hsa',  # homo sapiens
                 pvalueCutoff = 1,
                 pAdjustMethod = "fdr",
                 universe = back_genes$ENTREZID,
                 minGSSize = 1,
                 qvalueCutoff = 1)


# ===============================================================================
dotplot(CC,  #
        x = "GeneRatio",  # x stick
        color = "p.adjust",  #right y stick
        showCategory = 15,  # show top 20 dots
        title = "Cellular Components"
) + scale_color_continuous(low='purple', high='green') + aes(shape=Count < 10)

dotplot(MF,  #
        x = "GeneRatio",  # x stick
        color = "p.adjust",  #right y stick
        showCategory = 15,  # show top 20 dots
        title = "Molecular Function"
) + scale_color_continuous(low='purple', high='green') + aes(shape=Count < 10)

dotplot(BP,  #
        x = "GeneRatio",  # x stick
        color = "p.adjust",  #right y stick
        showCategory = 15,  # show top 20 dots
        title = "Biological Process"
) + scale_color_continuous(low='purple', high='green') + aes(shape=Count < 10)

dotplot(KEGG,  #
        x = "GeneRatio",  # x stick
        color = "p.adjust",  #right y stick
        showCategory = 15,  # show top 20 dots
        title = "KEGG"
) + scale_color_continuous(low='purple', high='green') + aes(shape=Count < 10) 
# 550*550


# =====================
for (ii in 1:length(KEGG@result[["geneID"]])){
  set_ <- strsplit(KEGG@result[["geneID"]][ii],'/')
  
  for(i in 1:length(set_[[1]])){
    set_[[1]][i] <- diff_genes[diff_genes$ENTREZID==set_[[1]][i],]$SYMBOL
  }
  
  KEGG@result[["geneID"]][ii] <- paste(set_[[1]],collapse='/')
}
rm(set_,i, ii)


cnetplot(KEGG,
         showCategory = 5, 
         colorEdge = T,
         node_label="all",
         circular = T)

# 1000*800


write.csv(KEGG@result, './Results/KEGG_SCARP.csv')

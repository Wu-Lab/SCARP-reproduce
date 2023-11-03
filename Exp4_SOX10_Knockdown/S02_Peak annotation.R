library(ChIPseeker)
library(ggplot2)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggupset)
library(ggimage)
require(clusterProfiler)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene



#=========================================================================
#             Peaks annotation (bed file + txdb object)
#=========================================================================
back_peaks <- readPeakFile('./Processed data/Sox10KD_peaks.bed')
back_x_annogene <- annotatePeak(back_peaks,
                                tssRegion = c(-3000, 3000),
                                TxDb = txdb,
                                annoDb = 'org.Hs.eg.db')
back_genes <- as.data.frame(back_x_annogene)
write.csv(back_genes, './Processed data/annotated_genes_3000bp_alltype.csv')
# write.csv(back_genes, './Processed data/annotated_genes_1000bp_alltype.csv')

back_genes <- back_genes[grepl('Promoter', back_genes$annotation),]
write.csv(back_genes, './Processed data/annotated_genes_1000bp.csv')




#=========================================================================
#             Data Visualization 
#=========================================================================
plotAnnoPie(back_x_annogene)
upsetplot(back_x_annogene, vennpie=TRUE)
plotDistToTSS(back_x_annogene,
              title = 'Distribution of ATAC-seq loci relative to TSS',
              ylab = 'Peaks (%) (5\'-> 3\')')

# only keep promoter sites
promoter_anno <- subset(back_x_annogene, abs(distanceToTSS)<1000)
promoter_anno <- subset(promoter_anno, grepl('Promoter', annotation))
y_promoter_anno <- as.data.frame(promoter_anno)

plotAnnoPie(promoter_anno) # 450*300
upsetplot(promoter_anno, vennpie=TRUE) # 700*500
plotDistToTSS(promoter_anno,
              title = 'Distribution of ATAC-seq loci relative to TSS (only promoters)',
              ylab = 'Peaks (%) (5\'-> 3\')') # 550*250





library(TCGAbiolinks)
library(plyr)
library(limma)
library(biomaRt)
library(SummarizedExperiment)
library(survminer)
require(survival)
library(clusterProfiler)




# =============================================================================
# clinical data
clinical_data <- read.table('./TCGA/survival_SKCM_survival.txt', sep = '\t', header = TRUE)
rownames(clinical_data) <- clinical_data[, 1]


# =============================================================================
# expression data
Expr_matrix <- read.table('./TCGA/TCGA-SKCM.htseq_fpkm.tsv',
                          header = TRUE, row.names = 1)
new_name <- c()
old_name <- strsplit(rownames(Expr_matrix), ".", fixed = TRUE)
for (i in 1:length(old_name)) {
  new_name <- c(new_name, old_name[[i]][1])
}
rownames(Expr_matrix) <- new_name
rm(new_name, old_name)

new_name <- c()
old_name <- strsplit(colnames(Expr_matrix), '.', fixed = TRUE)
for (i in 1:length(old_name)) {
  temp_chr <- paste(old_name[[i]], collapse = '-')
  new_name <- c(new_name, substr(temp_chr, 1, nchar(temp_chr)-1))
}
colnames(Expr_matrix) <- new_name
rm(new_name, old_name, i, temp_chr)


# =============================================================================
# gene id transformation, only keep the genes that can be transformed into geneid
gene <- bitr(rownames(Expr_matrix), "ENSEMBL", "SYMBOL", "org.Hs.eg.db")
Expr_matrix <- cbind(rownames(Expr_matrix), Expr_matrix)
colnames(Expr_matrix)[1] <- "ENSEMBL"
Gene_expr_df <- merge(gene, Expr_matrix, by = "ENSEMBL")
rm(gene)


# delete duplicate genes
Gene_expr_df <- Gene_expr_df[-which(duplicated(Gene_expr_df[, 2])), ]
rownames(Gene_expr_df) <- Gene_expr_df[, 2]
Gene_expr_df <- Gene_expr_df[, c(-1, -2)]

# filter
clinical_data <- clinical_data[intersect(clinical_data$sample,
                                         colnames(Gene_expr_df)),]
Gene_expr_df <- Gene_expr_df[, intersect(clinical_data$sample,
                                         colnames(Gene_expr_df))]

# the name of TCGA sample splited by '-'，
# The first three columns are patient numbers，
# The fourth column is the type: 11: normal; 01: tumor，
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
group <- strsplit(colnames(Gene_expr_df), "-")
class <- sapply(group, function(I) {I[4]})
class[which(grepl("11",class)==TRUE)]<-"Solid Tissue Normal"
class[which(grepl("06",class)==TRUE)]<-"Metastatic"
class[which(grepl("07",class)==TRUE)]<-"Additional Metastatic"
class[which(grepl("01",class)==TRUE)]<-"Primary Solid Tumor"

clinical_data <- cbind(clinical_data, class)

# merge survival data to expression data according to sample information
Gene_expr_df <- rbind(Gene_expr_df, clinical_data$X_PATIENT)
rownames(Gene_expr_df)[dim(Gene_expr_df)[1]] <- 'X_PATIENT'

df_OS <- merge(clinical_data, t(Gene_expr_df), by = "X_PATIENT")

df_OS_dropdu <- df_OS[-which(duplicated(df_OS[, 1])), ]
rownames(df_OS_dropdu) <- df_OS_dropdu[, 1]
df_OS_dropdu <- df_OS_dropdu[, -1]


# Metastatic: 367
# Additional Metastatic:1
# Primary Solid Tumor: 103
# Solid Tissue Normal: 1






#########################################################################
#########################################################################

# Dividing the sample according to genes
# gene <- 'MITF'
# gene <- 'TRIM69' 
# gene <- 'STK17B'
# gene <- 'SLC14A2'
# gene <- 'P4HTM'
# gene <- 'AKAP13'
# gene <- 'ZCCHC24'
# gene <- 'CTBP2'
# gene <- 'ZC3H6'
# gene <- 'COQ10B'
# gene <- 'XPNPEP3'
gene <- 'FBXO38'


# sample screen
df_tumor <- df_OS_dropdu[df_OS_dropdu$class!='Solid Tissue Normal', ]

gexp <- as.numeric(df_tumor[, gene])
sample_div <- quantile(gexp, 0.5)
survival_var_choose <- gexp
survival_var_choose[gexp <= sample_div] <- 'low expr'
survival_var_choose[gexp > sample_div] <- 'high expr'
print(paste0(
  'propotion:',
  sum(survival_var_choose == 'high expr') / length(gexp)
))
colnames(df_tumor)[1] <- gene
df_tumor[, 1] <- survival_var_choose



# ==============================================================================

# tow group
fit <- survfit(Surv(OS.time, OS) ~ FBXO38, data = df_tumor)
# save: 650, 510
ggsurvplot(
  fit,
  size = 1,
  linetype = "strata",
  break.time.by = 2000,
  risk.table.col = "strata",
  legend = c(0.85, 0.95),
  legend.title = "",
  palette = "Dark2",
  conf.int = TRUE,
  pval = TRUE,
  pval.coord = c(7000, 0.7),
  pval.size = 4,
  pval.method = TRUE,
  pval.method.size = 4,
  pval.method.coord = c(7000, 0.8),
  risk.table = TRUE,
  surv.median.line = "hv",
  legend.labs = c("high expression", "low expression"),
  risk.table.y.text.col = TRUE,
  xlab = "Follow up time(d)",
  title = gene
)
# 700*600


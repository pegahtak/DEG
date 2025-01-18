# PCA for CAR constructs:

# single run
rm (list=ls())
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(purrr)
library(dplyr)
library(biomaRt)
library (org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(purrr)
library(readxl)
library(tibble)
library(xlsx)
library(purrr)
library(dplyr)
rm(list = ls())
source("create_conditions.R")

# setwd("/Volumes/My_Passport/paper2_result/new_pipeline/")
getwd()
#import data
data<- read.table("data/countMatrix_all.txt", header = T)
geneIDs<- data$Geneid
gene_IDs <- read.table("/Volumes/My_Passport/paper2_result/new_pipeline/Ensemble_IDs.txt", header = T , fill = T)
gene_IDs[gene_IDs==""]<-NA
colnames(gene_IDs) <- c("ensembl_gene_id", "hgnc_symbol")

data$Geneid<- NULL

conditions<-creat_conditions(data)

rownames(data)<- geneIDs
colnames(data)<- unlist(map(strsplit(colnames(data), split = "_CKD"), 1))
data<- data[complete.cases(data),]
geneIDs<- rownames(data)

#convert ensemble IDs to gene names
data <- data %>%
  rownames_to_column(var = "ensembl_gene_id")
data <- data %>%
  left_join(gene_IDs, by = "ensembl_gene_id")
data$hgnc_symbol[data$hgnc_symbol==""]<-NA
data<- data[complete.cases(data), ]

# remove repeated genes
ind.dup<- which(!duplicated(data$hgnc_symbol) & duplicated(data$hgnc_symbol, fromLast=T))
if (length(ind.dup)>0)
{
  data<- data[-ind.dup,]
}

Genes<- data$hgnc_symbol
ENSEMBLE<-data$ensembl_gene_id
data$ensembl_gene_id<-NULL
rownames(data)<- Genes
data$hgnc_symbol<-NULL




group1<- "dualCAR_CD8_monomer"
group2<- "Full_CD8_Mb_wtIL2"
group3<- "Full_CD8_control"
# group4<- "fulldel_CD8_V7"
Trx_group<-c(group1)

group5<- "dualCAR_CD8_control"
group6<- "dualCAR_CD8_V9_Mb"
group7<- "dualCAR_CD8_V7"
group8<- "Full_CD4_control"
Ctrl_group<- c(group5, group6, group7)


data_ind<- unlist(conditions[c(Trx_group , Ctrl_group)])
# filter data
data_files<- data[,data_ind]
colnames(data_files)
data_files$dualCAR_CD8_V7_2<- NULL
# make direction
dir_name <- "dualCAR"#change back
dir.create(dir_name, showWarnings = FALSE)


# create meta data for desq2
meta<- matrix("", nrow = ncol(data_files) , ncol=1)
colnames(meta)<- c("sample_type" )
meta<- as.data.frame(meta)
sample_names<- colnames(data_files)

rownames(meta)<- sample_names
meta<- as.data.frame(meta)
meta[,1]<- substr(rownames(meta) , 1 , nchar(rownames(meta))-2 )
# meta$samplenames<- c(paste0(group1, "_" , 1:length(Trx)) , paste0(group2, "_" , 1:length(Ctrl))) ## add to code
meta$samplenames<- rownames(meta)

meta

data_files<-data_files%>%rownames_to_column(var="SYMBOL")
rownames(data_files)<- data_files$SYMBOL

data_files$SYMBOL<-NULL
# deseq with rlog 
dds <- DESeqDataSetFromMatrix(data_files, meta, design = ~ sample_type)
dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts , file.path(dir_name , paste0("normalized_counts.txt") ) 
            , sep="\t" ,quote = F, row.names = T)
rld <- rlog(dds, blind=TRUE)# regularized log transformed
log_counts<- assay(rld)
pca_result <- prcomp(log_counts, scale. = F)
pca_data<- pca_result$rotation %>% as.data.frame()%>%dplyr::select(PC1 , PC2)
pca_data$Group<- substr(rownames(pca_data), 1, nchar(rownames(pca_data))-2 )
write.table(log_counts , file.path(dir_name , paste0("log_norm_counts.txt")) ,
            sep = "\t", quote = F , row.names = T) 
pca_data <- prcomp(t(assay(rld)))
loadings <- pca_data$rotation
top20_genes<- rownames(loadings[order(abs(loadings[, "PC1"]) , decreasing = T),])[1:20]
# top20_gene_symbols <- gene_symbols %>%
#   filter(ensemble_ID %in% top20_genes) %>%
#   pull(symbol)
write.table(top20_genes, file.path(dir_name, paste0("top20_PC1_",dir_name, ".txt"))
            , col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
pdf(file.path(dir_name, paste0("PCA_",dir_name, ".pdf")), height = 5, width = 5)
print(plotPCA(rld, intgroup="sample_type" , pcsToUse=1:2 , ntop= 10000) +
        ggtitle(dir_name) )
 # ggplot(pca_data, aes(x = PC1, y = PC2,colour = Group)) + geom_point()
dev.off()

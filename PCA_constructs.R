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
# data <- data %>%
#   rownames_to_column(var = "ensembl_gene_id")
# data <- data %>%
#   left_join(gene_IDs, by = "ensembl_gene_id")
# data$hgnc_symbol[data$hgnc_symbol==""]<-NA
# data<- data[complete.cases(data), ]
# 
# # remove repeated genes
# ind.dup<- which(!duplicated(data$hgnc_symbol) & duplicated(data$hgnc_symbol, fromLast=T))
# if (length(ind.dup)>0)
# {
#   data<- data[-ind.dup,]
# }
# 
# Genes<- data$hgnc_symbol
# ENSEMBLE<-data$ensembl_gene_id
# data$ensembl_gene_id<-NULL
# rownames(data)<- Genes
# data$hgnc_symbol<-NULL




group1<- "Full_CD8_control"
group2<- "Costimonly_CD8_control"
group3<- "CD3only_CD8_control"
group4<- "fulldel_CD8_control"
Trx_group<-c(group1, group2, group3, group4)

group5<- "Full_CD8_V7"
group6<- "Costimonly_CD8_V7"
group7<- "CD3only_CD8_V7"
group8<- "fulldel_CD8_V7"
Ctrl_group<- c(group5, group6, group7, group8)


data_ind<- unlist(conditions[c(Trx_group , Ctrl_group)])
# filter data
data_files<- data[,data_ind]
colnames(data_files)
# data_files$dualCAR_CD8_V7_2<- NULL
# make direction
dir_name <- "Constructs"#change back
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
pca_result <- prcomp(t(log_counts), scale. = F)
pca_data <- as.data.frame(pca_result$x)
pca_data$Group <- colData(rld)$sample_type
# pca_data<- pca_result$rotation %>% as.data.frame()%>%dplyr::select(PC1 , PC2)
# pca_data$Group<- substr(rownames(pca_data), 1, nchar(rownames(pca_data))-2 )
write.table(log_counts , file.path(dir_name , paste0("log_norm_counts.txt")) ,
            sep = "\t", quote = F , row.names = T) 
pca_log <- prcomp(t(assay(rld)))
loadings <- pca_log$rotation
top20_genes<- rownames(loadings[order(abs(loadings[, "PC1"]) , decreasing = T),])[1:20]
# top20_gene_symbols <- gene_symbols %>%
#   filter(ensemble_ID %in% top20_genes) %>%
#   pull(symbol)
write.table(top20_genes, file.path(dir_name, paste0("top20_PC1_",dir_name, ".txt"))
            , col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


sdev<- pca_result$sdev
variance_explained <- sdev^2
proportion_variance_explained <- 100*variance_explained / sum(variance_explained)
proportion_variance_explained<- proportion_variance_explained%>%round(digits = 2)%>%t()%>%as.data.frame()
names(proportion_variance_explained)<- paste0("PC", 1:length(proportion_variance_explained))


# pdf(file.path(dir_name, paste0("PCA_",dir_name, ".pdf")), height = 5, width = 5)
print(plotPCA(rld, intgroup="sample_type" , pcsToUse=1:2 , ntop= 10000) +
        ggtitle(dir_name) )
 # ggplot(pca_data, aes(x = PC1, y = PC2,colour = Group)) + geom_point()
# dev.off()

my_cols<- c(Rested="darkgrey", wtIL2="pink", Mb_wtIL2="purple", V7="green", monomer="darkgreen",Mb_V9="cyan")
my_shapes<- c(Full_CD8=22, Costimonly_CD8=21,CD3only_CD8=24, fulldel_CD8=25, dualCAR_CD8=7, Full_CD4=23)

pca_data$Construct<- factor(c(rep("Full_CD8", 4), rep("Costimonly_CD8", 4), rep("CD3only_CD8", 4), rep("fulldel_CD8", 4)
                            ,rep("Full_CD8", 4), rep("Costimonly_CD8", 4), rep("CD3only_CD8", 4), rep("fulldel_CD8", 4)))
pca_data$Treatment<- factor(c(rep("Rested", 16), rep("V7", 16)))
pdf(file.path(dir_name, paste0("PCA.1_2:3_",dir_name, ".pdf")), height = 5, width = 5)

ggplot(pca_data, aes(PC2, PC3)) + 
  geom_point(aes(fill = Treatment, shape = Construct), size = 3, stroke = 0.75, colour = "black", show.legend = T ) +  # Black borders
  xlab(paste0("PC2 ", proportion_variance_explained$PC2, "%")) +
  ylab(paste0("PC3 ", proportion_variance_explained$PC3, "%")) +
  scale_fill_manual(name = "Treatment", values = my_cols) +  # Ensure correct fill colors
  scale_shape_manual(values = my_shapes) +  # Custom shapes for Construct
  theme(legend.key = element_rect(colour = "transparent"))+
  guides(
    fill = guide_legend(override.aes = list(shape = 21, stroke = 0)),  # Fixes the Treatment legend
    shape = guide_legend()  # Keeps Construct legend as is
  ) 


dev.off()
graphics.off()







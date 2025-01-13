
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

# convert ensemble IDs to gene names
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

for ( i in 1:(length(conditions) - 1))
{
  for( j in (i+1):(length(conditions)))
  {
    # get indices
    if ( i == j)
      next
    
    Trx<- conditions[[i]]
    Ctrl<- conditions[[j]]
  
    group1<- names(conditions)[i]
    group2<- names(conditions)[j]
  
    
    # filter data
    data_files<- data[,c(Trx , Ctrl)]
   
    colnames(data_files)
    
    # make direction
    dir_name <- paste0(group1, "_vs_", group2) #change back
   
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
    #meta$Batch<- c( rep(batches$Full_CD8_V7 , 4) , rep(batches$dualCAR_CD8_V7 , 4) )
    meta
    
    data_files<-data_files%>%rownames_to_column(var="SYMBOL")
    rownames(data_files)<- data_files$SYMBOL
    # data_files<- data_files%>%left_join(gene_IDs , by="ensembl_gene_id")
    # data_files$hgnc_symbol[data_files$hgnc_symbol==""]<-NA
    # data_files<- data_files[complete.cases(data_files),]
    # dup.id<- which(duplicated(data_files$ensembl_gene_id))
    # if(length(dup.id) >0)
    #   data_files<- data_files[-dup.id, ]
    # rownames(data_files)<- data_files$ensembl_gene_id
    #data_files$ensembl_gene_id<-NULL
    data_files$SYMBOL<-NULL
    # deseq with rlog 
    dds <- DESeqDataSetFromMatrix(data_files, meta, design = ~ sample_type)
    dds <- estimateSizeFactors(dds)
    
    normalized_counts <- counts(dds, normalized=TRUE)
    write.table(normalized_counts , file.path(dir_name , paste0("normalized_counts.txt") ) 
                , sep="\t" ,quote = F, row.names = T)
    rld <- rlog(dds, blind=TRUE)# regularized log transformed
    log_counts<- assay(rld)
    write.table(log_counts , file.path(dir_name , paste0("log_norm_counts.txt")) ,
                sep = "\t", quote = F , row.names = T) 
    pdf(file.path(dir_name, paste0("PCA_", group1, "_vs_", group2, ".pdf")), height = 5, width = 5)
    print(plotPCA(rld, intgroup="sample_type" , pcsToUse=1:2 , ntop= 500) +
            ggtitle(paste0(group1 , " vs " , group2)) )
    dev.off()
    #graphics.off()
    
    # get top 20 genes of PC1 
    pca_data <- prcomp(t(assay(rld)))
    loadings <- pca_data$rotation
    
    top20_genes<- rownames(loadings[order(abs(loadings[, "PC1"]) , decreasing = T),])[1:20]
    # top20_gene_symbols <- gene_symbols %>%
    #   filter(ensemble_ID %in% top20_genes) %>%
    #   pull(symbol)
    write.table(top20_genes, file.path(dir_name, paste0("top20_PC1_", group1, "_vs_", group2, ".txt"))
                , col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    
    # correlation heatmap
    rld_mat <- assay(rld) 
    rld_cor <- cor(rld_mat)
    pdf(file.path(dir_name, paste0("heatmap_", group1, "_vs_", group2, ".pdf"))
        , height = 5, width = 5)
    pheatmap(rld_cor)
    dev.off()
    #graphics.off()
    
    
    dds <- DESeqDataSetFromMatrix(data_files, meta, design = ~ sample_type)
    dds <- DESeq(dds, )
    dds <- DESeq(dds, test="Wald")
    
    lfc_cutoff <- 1.5
    padj_cutoff <- 0.05
    res <- results(dds, contrast = c("sample_type", group1, group2)) 
    res_tab<- res %>%
      as.data.frame()
    # 
    # SYM<- gene_symbols %>%
    #   filter(ensemble_ID %in% rownames(res_tab))%>% pull (symbol)
    res_tab<- res_tab%>%rownames_to_column(var="Symbol")
    
    sig_results <- res_tab  %>%
      filter(padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff) 
    
    write.table(res_tab ,file.path(dir_name, paste0("res_lfc_",group1,"_vs_", group2,  ".txt") ) , quote = F , sep = "\t")
    #merge with gene symbols
    # sig_results <- sig_results %>%
    #   left_join(gene_symbols, by = c("gene" = "ensemble_ID"))
    
    #save
    write.xlsx(sig_results, file.path(dir_name, paste0("DEGs_",group1 , "_vs_", group2 ,".xlsx") ),
               )
    
    
    #lab = gene_symbols$symbol
    # lables_sig<- unique(c("IL2", "CSF2","IL3","IFNG","CXCL8", "XCL1", "XCL2", "CCL4" , sig_results$Symbol))
    # lab_sig_2<-c("IL2", "CSF2","IL3","IFNG","CXCL8", "XCL1", "XCL2", "CCL4", "IL31")
    # Volcano plot
    
    #Slabs<- sig_results%>% arrange(abs(log2FoldChange))%>%pull(Symbol)%>% head(20)
    pdf(file.path(dir_name, paste0("Volcano_", group1, "_vs_", group2, ".pdf")), height =10, width = 15)
    res_tableOE_tb <- res %>% data.frame() 
    res_tableOE_tb<- res_tableOE_tb[complete.cases(res_tableOE_tb),]
    res_tableOE_tb <- res_tableOE_tb %>% 
      plyr::mutate(threshold_OE = padj < padj_cutoff & abs(log2FoldChange) >= lfc_cutoff)
    res_tableOE_tb<- res_tableOE_tb%>% rownames_to_column(var="symbol")
    x<-res_tableOE_tb$pvalue
    res_tableOE_tb$pvalue<- res_tableOE_tb$padj
    res_tableOE_tb$padj<-NULL
    
    volcano_plot <- EnhancedVolcano(res_tableOE_tb, lab = res_tableOE_tb$symbol,
                                    x = 'log2FoldChange', xlim = c(-9 , 9),
                                    y = 'pvalue', title = paste0(group1, " vs ", group2),
                                    pCutoff = padj_cutoff, FCcutoff = lfc_cutoff, pointSize = 2,
                                    labSize = 4, colAlpha = 4/5 , legendLabSize = 25 , axisLabSize = 25)
    
    # Modify the font size of the axis numbers directly using the theme() function
    volcano_plot <- volcano_plot + theme(axis.text.x = element_text(size = 25),
                                         axis.text.y = element_text(size = 25))
    print(volcano_plot)
    dev.off()
    #graphics.off()
    
  }
}

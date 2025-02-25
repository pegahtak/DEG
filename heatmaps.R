# heatmap of Full CD8 + with averages!!! with fold change
rm(list=ls())
setwd("/Volumes/My_Passport/paper2_result/new_pipeline/DEG/")
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(tximportData)
#library(readr)
#library(tximport)
#library(purrr)
library(dplyr)
#library(biomaRt)
#library (org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(purrr)
library(readxl)
library(tibble)
library(xlsx)
library(RColorBrewer)
library(tibble)
rm(list = ls())

# read gene set  
geneDat<- read.xlsx("/Volumes/My_Passport/paper2_result/new_pipeline/DEG/Full_CD8_V7_vs_fulldel_CD8_V7/DEGs_Full_CD8_V7_vs_fulldel_CD8_V7.xlsx"
                    , sheetIndex=1)
geneDat$NA.<-NULL


# select top 250  upregulated genes
# change desc(log2FoldChange to desc(abs(log2FoldChange))
# topGenes<- geneDat %>% arrange(desc(log2FoldChange))%>%pull(Symbol)%>% head(100)
topGenes<- geneDat %>% arrange(desc(log2FoldChange) )%>%pull(Symbol)
topGenes50 <- geneDat %>% arrange(desc(log2FoldChange) )%>%pull(Symbol)%>% head(50)

control<- "fulldel_CD8_V7"
sample_names<- rev(c("Full_CD8_V7","Full_CD8_Mb_wtIL2", "Full_CD8_wtIL2"))  # using rev to put the first column on right
fig_name<- "wtIL2_bridging"

files_names<- paste0(sample_names, "_vs_",control )
list_files<-as.list(paste0("/Volumes/My_Passport/paper2_result/new_pipeline/DEG/" , files_names, "/res_lfc_", files_names , ".txt"))
x_list <- lapply(list_files, read.delim)

rm(x)# for code testing
x <- lapply(x_list, function(file) {
  subset(file, Symbol %in% topGenes)[, c("log2FoldChange", "Symbol"), drop = FALSE]
})
x<- lapply(x, function(file) file[!duplicated(file$Symbol), ])

final_table <- do.call(cbind, x)
typeof(final_table)
rownames(final_table)<- NULL
final_table<- column_to_rownames(final_table,var="Symbol")
final_table<- final_table[, grep("log2FoldChange", colnames(final_table)) ]


final_table<- final_table %>%
  mutate(sd = rowSds(as.matrix(across(everything()))))%>%
 arrange(desc(sd))%>% as.data.frame()%>%dplyr::select(-sd)

topGenes250<- rownames(final_table)%>%head(250)
#topGenes50<- final_table%>%rownames()%>%head(50)

#rownames(final_table)<- final_table$Symbol
#Symbols<-  final_table$Symbol
#final_table<- final_table[,-grep("Symbol" , colnames(final_table))]
colnames(final_table)<- files_names

# remove rows with NA
#final_table <- final_table[complete.cases(final_table), ]
#write.csv(final_table, paste0(fig_name,".csv")  , quote = F)

#final_table<- final_table[ order(final_table[,ncol(final_table)] , decreasing = T), ]
colnames(final_table)<- rev(c("V7","Mb_wtIL2", "wtIL2"))
final_table_50<- final_table[which(rownames(final_table)%in%topGenes50),]
final_table<- final_table[complete.cases(final_table),]

t<- final_table
dir.create("fig6")
setwd("fig6/")
write.csv(final_table ,paste0(fig_name, ".csv") , quote = F)
# removing outliers
t[t > 20] <- NA
#t[t< -10]<- NA
# ind<- which(rownames(t)=="RNA5SP474")
# t<- t[-ind, ]
# ind<- which(rownames(t)=="RNA5SP295")
# t<- t[-ind, ]

# r_ind<- which(rownames(final_table)%in%remove)
# t<- final_table[-r_ind, ]

mycols<- colorRampPalette(c("blue", "white", "red"))(100)
# myBreaks<- c(-4, -2, 4, 10)
#ph.breaks <- seq(from = -10, to = 10, length.out = 201)
ph.breaks <- c(seq(from = min(t, na.rm=T), to = -0.1, length.out = 50), seq(from=0, to=max(t, na.rm = T), length.out=51) )
pdf(paste0(fig_name , "all.pdf"), width = 4 , height = 20)
#par(mar=c(20,20,20,20))
print(pheatmap(t , cluster_cols = T , show_colnames = T,
               fontsize_row = 5 , fontsize_col = 15 , cluster_rows = T, fontsize = 15 , legend = T
               , annotation_legend = T , angle_col = 315 , na_col="purple" , clustering_method = "ward.D", breaks = ph.breaks ))
dev.off()
graphics.off()


t50<- final_table_50
t50[t50 > 19] <- NA
#t50<- t50[complete.cases(t50),]
ph.breaks50<- c(seq(from = -6, to = -0.1, length.out = 50), seq(from=0, to=10, length.out=51) )
pdf(paste0(fig_name , "top50.pdf"), width = 4 , height = 20)
#par(mar=c(20,20,20,20))
print(pheatmap(t50 , cluster_cols = T , show_colnames = T,
               fontsize_row = 10 , fontsize_col = 10 , cluster_rows = T, fontsize = 15 , legend = T
               , annotation_legend = T , angle_col = 315 , na_col="purple" , clustering_method = "ward.D" ))
dev.off()
graphics.off()

t <- t[complete.cases(t), ]
vars<-rowVars(as.matrix(t))
ind<- which(vars> median(vars))
t<- t[ind, ]
ph.cols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
# pdf(paste0(fig_name , "_withLFC_main.pdf"), width = 4.5 , height = 15.5)
# par(mar=c(1,1,1,1))
pdf(paste0(fig_name , "_withLFC_main.pdf"), width = 50 , height = 155)
print(pheatmap(t, cluster_cols = T , show_colnames = T,
               fontsize_row = 70, fontsize_col = 100 , cluster_rows = T, fontsize = 20 , legend = T
              , annotation_legend = T , angle_col = 315 , na_col="purple" , clustering_method = "average", breaks=ph.breaks  ) )
dev.off()
graphics.off()

#, breaks = ph.breaks , color =ph.cols 
# paletteLength <- 50
# myBreaks <- c(seq(min(t, na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1), 



#               seq(max(t, na.rm = T)/paletteLength, max(t, na.rm = T), length.out=floor(paletteLength/2)) ,20
#breaks<- c()



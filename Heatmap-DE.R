### Heatmap-DE.R -- Andrew R Gross -- 2018-01-29
### A basic script for generating heatmaps of expression data

### Header   #############################################################################################
library(gplots)
library(RColorBrewer)
library(biomaRt)
library(DESeq2)
library(Hmisc)
library(reshape)

### Functions   #############################################################################################
stats.calculator <- function(dataframe) {
  min <- apply(expression.filt,1,min)
  ctr.med <- apply(expression.filt[which(Treatment == 'CTR')],1,median)
  dis.med <- apply(expression.filt[which(Treatment == 'DIS')],1,median)
  stats.df <- data.frame(min, ctr.med, dis.med)
}

### Input #############################################################################################
# Input a table of expression data
expression.df <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)

# Input de results
de.results <- read.table("z:/Data/RNAseq HT neurons and tissue/2nd rerun/DE_genes_ANDREW/Obs_iHT_vs_Ctr_iHT--1184.txt", header=TRUE, row.names = 1)

### Format #############################################################################################
expression.df <- expression.df[c(2,3,19,20,1,11,12,16,17,18)]
expression.df <- round(expression.df*100)

### Generate DE results #############################################################################################
Treatment <- c('DIS','DIS','DIS','DIS','CTR','CTR','CTR','CTR','CTR','CTR')
column.metadata <- data.frame(Treatment, row.names = names(expression.df))
expression.se <- DESeqDataSetFromMatrix(countData = as.matrix(expression.df), colData = column.metadata, design = ~ Treatment)
de.results <- results(DESeq(expression.se)); sum(de.results$padj < 0.1, na.rm=TRUE)
de.results <- de.results[order(de.results$log2FoldChange),]
de.results <- as.data.frame(subset(de.results, padj < 0.1))

### Format, stage 2 #############################################################################################
### Filter rows by p value
length(which(de.results$padj <= 0.03))
de.results.filt <- de.results[de.results$padj <= 0.03,]
expression.filt <- expression.df[match(row.names(de.results.filt),row.names(expression.df)),]
expression.stats <- stats.calculator(expression.filt)

### Filter rows by min values
expression.filt <- expression.filt[which(expression.stats$min >=10),]
expression.stats <- stats.calculator(expression.filt)

### Convert to log fold changes
expr.logfc <- log(expression.filt/expression.stats$ctr.med)
## Optional: normalize
expr.logfc <- expr.logfc-apply(expr.logfc,1,median)

### Plot #############################################################################################
### Reshape
expr.logfc$gene <- row.names(expr.logfc)
expr.logfc <- expr.logfc[rev(row.names(expr.logfc)),]
expr.logfc$gene <- factor(expr.logfc$gene, levels = expr.logfc$gene, ordered = TRUE)
expr.logfc.m <- melt(expr.logfc)

### Plot
ggplot(expr.logfc.m, aes(x = variable, y = gene)) +
  geom_tile(aes(fill = value)) +
  #scale_fill_gradientn(colors = c('blue', 'white', 'red'), limits = c(-3.3,3.3)) +
  #scale_fill_gradientn(colors = c('green', 'black', 'red'), limits = c(-3.1,3.1)) +
  scale_fill_gradientn(colors = c('green', 'black', 'red')) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(position = "top") +
  labs(title = 'Differential Expression') +
  theme(#axis.text.x = element_text(size = 12),
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 11),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    panel.border = element_rect(color = "black", fill = NA)) 

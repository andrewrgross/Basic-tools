### Ensembl to Entrez ID converter -- Andrew R Gross -- 2018-08-21
### This script converts Ensembl IDs to Entrez IDs

####################################################################################################################################################
### Header

library(gage)
#library(gageData)
library(biomaRt)
library(org.Hs.eg.db)
#library(AnnotationDbi)

ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
listMarts(host="www.ensembl.org")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

####################################################################################################################################################
### Main Function

convert.to.entrez <- function(df.to.convert) {
  data("egSymb")
  x <- org.Hs.egENSEMBL                          # Get the entrez gene IDs that are mapped to an Ensembl ID
  mapped_genes <- mappedkeys(x)                  # Convert to a list
  EN.to.EZ.dictionary <- as.list(x[mapped_genes])
  ##########################################################################
  ### Batch conversion, homebrew
  entrez.ids <- c()
  for (ensembl.id in row.names(df.to.convert)) {
    #print(ensembl.id)
    entrez.id.pos <- grep(ensembl.id, EN.to.EZ.dictionary)[1]
    #print(entrez.id.pos)
    entrez.id <- names(EN.to.EZ.dictionary)[entrez.id.pos]
    if (length(entrez.id.pos) == 0) {entrez.id <- NA}
    #print(entrez.id)
    entrez.ids <- c(entrez.ids, entrez.id)
  }
  
  df.to.convert$loop <- entrez.ids
  
  ##########################################################################
  ### Batch conversion via sym2eg
  df.to.convert$sym2eg <- sym2eg(df.to.convert$Gene)
  #df.to.convert$sym2eg <- sym2eg(df.to.convert.w.sym$Symbol)
  
  ##########################################################################
  ### Join two lists
  entrez.join <- c()
  for (row.num in 1:nrow(df.to.convert)) {
    id <- df.to.convert$sym2eg[row.num]
    if (is.na(id) == TRUE) {id = df.to.convert$loop[row.num]}
    entrez.join <- c(entrez.join, id)
  }
  df.to.convert$join <- entrez.join
  ##########################################################################
  ### Remove NAs
  unconverted.pos <- which(is.na(df.to.convert$join))
  if (length(unconverted.pos) > 0) {
    df.to.convert <- df.to.convert[-unconverted.pos,]
    print(paste(length(unconverted.pos), 'NAs removed.'))
  }
  
  ### Examine non-matches
  match <- df.to.convert$loop == df.to.convert$sym2eg
  non.match.count = length(grep(FALSE, match))
  summary(match)
  print(paste('Conversion of', nrow(df.to.convert), 'IDs produced', non.match.count, 'conflicts.'))
  
  ### Examine non-unique IDs
  non.unique.count = length(df.to.convert$join)-length(unique(df.to.convert$join))
  print(paste(non.unique.count, 'Non unique IDs generated'))
  if (non.unique.count == 0) {
    row.names(df.to.convert) <- df.to.convert$join
    df.to.convert <- df.to.convert[1:(ncol(df.to.convert)-3)]
    print('ID conversion complete')
  }
  if (non.unique.count >0) {
    id.count.table <- table(df.to.convert$join)
    id.count.table.order <- order(id.count.table, decreasing = TRUE)
    id.count.table <- id.count.table[id.count.table.order]
    non.unique.pos <- which(id.count.table>1)
    non.unique.ids <- names(non.unique.pos)
    print('The following non unique IDs were generated.  Please resolve manually:')
    print(head(non.unique.ids,20))
  }
  return(df.to.convert)
}

####################################################################################################################################################
### Convert IDs
### Input Data
setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Input files/")

results.als <- read.csv('DEGs - ALS 0.05 pv.csv', row.names= 1)
results.ko <- read.csv('DEGs - KO 0.05 pv.csv', row.names = 1)

##########################################################################
### Convert
#df.to.convert <- convert.to.entrez(head(results.w.samples,100))
results.als.ez <- convert.to.entrez(results.als)
results.ko.ez <- convert.to.entrez(results.ko)

##########################################################################
### Assess and cleanup conversion

grep('9782', results.ko.ez$join)
results.ko.ez[grep('9782', results.ko.ez$join),]
results.ko.ez2 <- results.ko.ez[-5097,]
row.names(results.ko.ez2) <- results.ko.ez2$join
results.ko.ez <- results.ko.ez2[1:11]

##########################################################################
### Output

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Input files/')

#write.csv(results.als.ez, 'DEGs - ALS 0.05 padj ENTREZ.csv')
#write.csv(results.ko.ez, 'DEGs - KO 0.05 padj ENTREZ.csv')

write.csv(results.als.ez, 'DEGs - ALS 0.05 pv ENTREZ.csv')
write.csv(results.ko.ez, 'DEGs - KO 0.05 pv ENTREZ.csv')

####################################################################################################################################################
### Manual Troubleshooting tools

##########################################################################
### Single conversion
#ensembl.id <- row.names(df.to.convert)[1]
#entrez.id.pos <- grep(ensembl.id, EN.to.EZ.dictionary)
#entrez.id <- names(EN.to.EZ.dictionary)[entrez.id.pos]

##########################################################################
### Assess conversion results
summary(df.to.convert)
summary(is.na(df.to.convert$sym2eg))
summary(is.na(df.to.convert$loop))
summary(is.na(df.to.convert$join))

### Examine unconverted
unconverted.pos <- which(is.na(df.to.convert$join))
unconverted <- df.to.convert[unconverted.pos,]

### Drop unconverted
df.to.convert <- df.to.convert[-unconverted.pos,]

### Examine non-matches
df.to.convert$match <- df.to.convert$loop == df.to.convert$sym2eg
summary(df.to.convert$match)
non.match.pos <- which(df.to.convert$match==FALSE)
non.matches <- df.to.convert[non.match.pos,]#[c(1,4,5,6,7,8,9,10,11,12,14,15,16)]
#test[-grep('Description', names(test))]
non.matches$loop.back <- eg2sym(non.matches$loop)
non.matches$sym.back <- eg2sym(non.matches$sym2eg)
head(non.matches,30)

### Examine non-unique IDs
length(df.to.convert$join)-length(unique(df.to.convert$join))
id.count.table <- table(df.to.convert$join)
id.count.table.order <- order(id.count.table, decreasing = TRUE)
id.count.table <- id.count.table[id.count.table.order]
non.unique.pos <- which(id.count.table>1)
non.unique.ids <- names(non.unique.pos)

##########################################################################
### Isolate non-unique IDs
non.unique.ids.df <- df.to.convert[0,]
for (entrez.id in non.unique.ids) {
  temp.df <- df.to.convert[grep(entrez.id, df.to.convert$join),]
  non.unique.ids.df <- rbind(non.unique.ids.df, temp.df)
}
head(non.unique.ids.df,50)

### Correct non-unique IDs manually
REPLACEMENT.ID <- 8209
NON.UNIQUE.ID.TO.REPLACE <- 
grep(REPLACEMENT.ID, df.to.convert$join)
df.to.convert[row.names(non.unique.ids.df)[1],]$join <- REPLACEMENT.ID
df.to.convert[row.names(non.unique.ids.df)[1],]

### Drop unwanted non-unique rows
NON.UNIQUE.ID.TO.DROP <- "ENSG00000280987"
row.to.drop <- grep(NON.UNIQUE.ID.TO.DROP, row.names(df.to.convert))
df.to.convert[row.to.drop,]
df.to.convert <- df.to.convert[-row.to.drop,]




#### SCRATCHWORK



df.to.convert <- results.KO
# df.convert <-results.ALS

#addGene <- function(dataframe) {
#  genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
#  genes <- genes[match(row.names(dataframe),genes[,1]),]
#  Gene <- c()
#  for (rowNumber in 1:length(genes[,1])) {
#    newGene <- genes[rowNumber,][,2]
#    Gene <- c(Gene, newGene)
#  }
#  dataframe[length(dataframe)+1] <- Gene
#  names(dataframe)[ncol(dataframe)] <- "Gene"
#  return(dataframe)
#}
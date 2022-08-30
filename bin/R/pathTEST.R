### pathTEST.R ###
# This script will take the output from seurat differential analysis / findMarker analysis to generate pathway analysis  
###################

## Loading requires packages
library(gage)
require(GOstats)
require(dplyr)
require(stringr)
require(reshape)

## Db for pathway
library("AnnotationDbi")
## if you are dealing with human
# library("org.Hs.eg.db")
library("org.Mm.eg.db")
require(GO.db)

## If you are using human, replace with "org.Hs.egPATH"
## ALso for GO analysis, replace with "org.Mm.egGODB")
mapped_genes <- mappedkeys(org.Mm.egPATH) ## KEGG


## If running on RStudio:
## you can skip this step
args = commandArgs(TRUE)
group1 = args[2]
group2 = args[3]
subtype = args[1]

## If running on RStudio
# group1  = "Control"
# group = "Disease"
# subtype = "CM" / celltype

## differential expression analysis output
## If you are unsure, replace with the following
# de = read.table(file=defile,header=TRUE)
## replace defile with the path to the seurat differential output
de = read.table(paste("DE_",group1,"_",group2,"_",subtype,".txt",sep=""),header=TRUE)

## add a column on entrez id
## Note: if you are using human, just replace org.Mm.eg.db with org.Hs.eg.db
de$entrez = mapIds(org.Mm.eg.db,keys = rownames(de),column = "ENTREZID", keytype = "SYMBOL", multiVals="first")
## remove any rows with no entrez id assigned
res = de[ !is.na(de$entrez),]

## preparing input for the test.
foldchanges = res$avg_log2FC
names(foldchanges) = res$entrez
head(foldchanges)

## subset genes (upregulated and significant 0.05)
up.entrez = res$entrez[ res$p_val_adj<0.05 & res$avg_log2FC>0.5 ]
down.entrez = res$entrez[ res$p_val_adj<0.05 & res$avg_log2FC<0.5 ]

## run KEGG analysis on up-regulated data
## Note: replace annotation with org.Hs.eg.db if dealing with human
params <- new('KEGGHyperGParams',geneIds = up.entrez, universeGeneIds = mapped_genes, pvalueCutoff = 0.1, testDirection = 'over',annotation = "org.Mm.eg.db")
hgOver <- hyperGTest(params)
result <- summary(hgOver)
returnGENE <- function(x){ return(paste(unlist(mget(up.entrez[up.entrez %in% unlist(mget(x,revmap(org.Mm.egPATH),ifnotfound=NA))],org.Mm.egSYMBOL)),collapse=","))}
result$gene_list = unlist(lapply(result$KEGGID,returnGENE))
write.table(result,file=paste(subtype,group1,group2,"GOstats_up_KEGG.txt",sep="_"),quote=FALSE,sep="\t",row.names=F)

## run KEGG analysis on up-regulated data
## Note: replace annotation with org.Hs.eg.db if dealing with human
params <- new('KEGGHyperGParams',geneIds = down.entrez, universeGeneIds = mapped_genes, pvalueCutoff = 0.1,testDirection = 'over',annotation = "org.Mm.eg.db")
hgOver <- hyperGTest(params)
result <- summary(hgOver)
returnGENE <- function(x){ return(paste(unlist(mget(down.entrez[down.entrez %in% unlist(mget(x,revmap(org.Mm.egPATH),ifnotfound=NA))],org.Mm.egSYMBOL)),collapse=","))}
result$gene_list = unlist(lapply(result$KEGGID,returnGENE))
write.table(result,file=paste(subtype,group1,group2,"GOstats_down_KEGG.txt",sep="_"),quote=FALSE,sep="\t",row.names=F)


### Repeat the analysis but with GO
mapped_genes <- mappedkeys(org.Mm.egGO)

params <- new('GOHyperGParams',geneIds = up.entrez,universeGeneIds = mapped_genes,ontology = 'BP',pvalueCutoff = 0.001,conditional = FALSE,testDirection = 'over',annotation = "org.Mm.eg.db")
hgOver <- hyperGTest(params)
result <- summary(hgOver)
returnGENE <- function(x){ return(paste(unlist(mget(up.entrez[up.entrez %in% unlist(mget(x,revmap(org.Mm.egGO),ifnotfound=NA))],org.Mm.egSYMBOL)),collapse=","))}
result$gene_list = unlist(lapply(result$GOBPID,returnGENE))
write.table(result,file=paste(subtype,group1,group2,"GOstats_up_GOBP.txt",sep="_"),quote=FALSE,sep="\t",row.names=F)

params <- new('GOHyperGParams',geneIds = down.entrez,universeGeneIds = mapped_genes,ontology = 'BP',pvalueCutoff = 0.001,conditional = FALSE,testDirection = 'over',annotation = "org.Mm.eg.db")
hgOver <- hyperGTest(params)
result <- summary(hgOver)
returnGENE2 <- function(x){ return(paste(unlist(mget(down.entrez[down.entrez %in% unlist(mget(x,revmap(org.Mm.egGO),ifnotfound=NA))],org.Mm.egSYMBOL)),collapse=","))}
result$gene_list = unlist(lapply(result$GOBPID,returnGENE2))
write.table(result,file=paste(subtype,group1,group2,"GOstats_down_GOBP.txt",sep="_"),quote=FALSE,sep="\t",row.names=F)

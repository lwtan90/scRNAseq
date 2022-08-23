require(stringr)
require(reshape)
library(gage)
#library(gageData)
library("AnnotationDbi")
library("org.Mm.eg.db")
require(GOstats)
require(GO.db)
require(dplyr)

args = commandArgs(TRUE)
group1 = args[2]
group2 = args[3]
subtype = args[1]
de = read.table(paste("DE_",group1,"_",group2,"_",subtype,".txt",sep=""),header=TRUE)

de$entrez = mapIds(org.Mm.eg.db,keys = rownames(de),column = "ENTREZID", keytype = "SYMBOL", multiVals="first")
res = de[ !is.na(de$entrez),]

foldchanges = res$avg_log2FC
names(foldchanges) = res$entrez
head(foldchanges)

mapped_genes <- mappedkeys(org.Mm.egPATH)

up.entrez = res$entrez[ res$p_val_adj<0.05 & res$avg_log2FC>0.5 ]
params <- new('KEGGHyperGParams',geneIds = up.entrez, universeGeneIds = mapped_genes, pvalueCutoff = 0.1, testDirection = 'over',annotation = "org.Mm.eg.db")
hgOver <- hyperGTest(params)
result <- summary(hgOver)
returnGENE <- function(x){ return(paste(unlist(mget(up.entrez[up.entrez %in% unlist(mget(x,revmap(org.Mm.egPATH),ifnotfound=NA))],org.Mm.egSYMBOL)),collapse=","))}
result$gene_list = unlist(lapply(result$KEGGID,returnGENE))
write.table(result,file=paste(subtype,group1,group2,"GOstats_up_KEGG.txt",sep="_"),quote=FALSE,sep="\t",row.names=F)

down.entrez = res$entrez[ res$p_val_adj<0.05 & res$avg_log2FC<(-0.5) ]
params <- new('KEGGHyperGParams',geneIds = down.entrez, universeGeneIds = mapped_genes, pvalueCutoff = 0.1,testDirection = 'over',annotation = "org.Mm.eg.db")
hgOver <- hyperGTest(params)
result <- summary(hgOver)
returnGENE <- function(x){ return(paste(unlist(mget(down.entrez[down.entrez %in% unlist(mget(x,revmap(org.Mm.egPATH),ifnotfound=NA))],org.Mm.egSYMBOL)),collapse=","))}
result$gene_list = unlist(lapply(result$KEGGID,returnGENE))
write.table(result,file=paste(subtype,group1,group2,"GOstats_down_KEGG.txt",sep="_"),quote=FALSE,sep="\t",row.names=F)


### for KEGG
#data(kegg.sets.hs)
#keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
#up.pathways = data.frame(id=rownames(keggres$greater), keggres$greater)
#write.table(up.pathways, file=paste(group1,group2,subtype,"gage_up_pathway_KEGG.txt",sep="_"),sep="\t",quote=FALSE)
#down.pathways = data.frame(id=rownames(keggres$less), keggres$less)
#write.table(down.pathways, file=paste(group1,group2,subtype,"gage_down_pathway_KEGG.txt",sep="_"),sep="\t",quote=FALSE)
#data(go.sets.hs)
#data(go.subs.hs)
#gobpsets = go.sets.hs[go.subs.hs$BP]
#gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
#up.pathways = data.frame(id=rownames(gobpres$greater), gobpres$greater)
#write.table(up.pathways, file=paste(group1,group2,subtype,"gage_up_pathway_gobp.txt",sep="_"),sep="\t",quote=FALSE)
#down.pathways = data.frame(id=rownames(gobpres$less), gobpres$less)
#write.table(down.pathways, file=paste(group1,group2,subtype,"gage_down_pathway_gobp.txt",sep="_"),sep="\t",quote=FALSE)
## GOStats

mapped_genes <- mappedkeys(org.Mm.egGO)

up.entrez = res$entrez[ res$p_val_adj<0.05 & res$avg_log2FC>0.5 ]
params <- new('GOHyperGParams',geneIds = up.entrez,universeGeneIds = mapped_genes,ontology = 'BP',pvalueCutoff = 0.001,conditional = FALSE,testDirection = 'over',annotation = "org.Mm.eg.db")
hgOver <- hyperGTest(params)
result <- summary(hgOver)
returnGENE <- function(x){ return(paste(unlist(mget(up.entrez[up.entrez %in% unlist(mget(x,revmap(org.Mm.egGO),ifnotfound=NA))],org.Mm.egSYMBOL)),collapse=","))}
result$gene_list = unlist(lapply(result$GOBPID,returnGENE))
write.table(result,file=paste(subtype,group1,group2,"GOstats_up_GOBP.txt",sep="_"),quote=FALSE,sep="\t",row.names=F)

down.entrez = res$entrez[ res$p_val_adj<0.05 & res$avg_log2FC<(-0.5) ]
params <- new('GOHyperGParams',geneIds = down.entrez,universeGeneIds = mapped_genes,ontology = 'BP',pvalueCutoff = 0.001,conditional = FALSE,testDirection = 'over',annotation = "org.Mm.eg.db")
hgOver <- hyperGTest(params)
result <- summary(hgOver)
returnGENE2 <- function(x){ return(paste(unlist(mget(down.entrez[down.entrez %in% unlist(mget(x,revmap(org.Mm.egGO),ifnotfound=NA))],org.Mm.egSYMBOL)),collapse=","))}
result$gene_list = unlist(lapply(result$GOBPID,returnGENE2))
write.table(result,file=paste(subtype,group1,group2,"GOstats_down_GOBP.txt",sep="_"),quote=FALSE,sep="\t",row.names=F)

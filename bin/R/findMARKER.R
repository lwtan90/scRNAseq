require(Seurat)
require(Signac)

##load("step4_regrouped.RData")


heart.combined = readRDS("SCT.prepped.RDS")

args = commandArgs(TRUE)
table(Idents(heart.combined))
##Idents(heart.combined) = heart.combined$integrated_snn_res.0.5
##table(Idents(heart.combined))
##heart.combined <- PrepSCTFindMarkers(heart.combined)
##saveRDS(heart.combined,file="SCT.prepped.RDS")

heart.combined.markers <- FindMarkers(heart.combined, assay="SCT", ident.1 = args[1], min.pct = 0.2, logfc.threshold = 0.2, only.pos=TRUE)

write.table(heart.combined.markers,file=paste("marker_",args[1],".txt",sep=""),sep="\t",quote=FALSE)


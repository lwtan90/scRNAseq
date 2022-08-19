require(Seurat)
require(Signac)

heart.combined = readRDS("SCT.prepped.RDS")


args = commandArgs(TRUE)
table(Idents(heart.combined))

heart.combined = subset(heart.combined, idents = args[1])
heart.combined = subset(heart.combined, sample %in% c("Control","LMNA1"))
Idents(heart.combined) = heart.combined$sample
table(Idents(heart.combined))
DefaultAssay(heart.combined) = "SCT"

heart.combined <- PrepSCTFindMarkers(heart.combined)
heart.combined.markers <- FindMarkers(heart.combined, assay="SCT", ident.1 = "LMNA1", ident.2 = "Control", min.pct = 0.2, logfc.threshold = 0.2)

write.table(heart.combined.markers,file=paste("DE_control_LMNA1_",args[1],".txt",sep=""),sep="\t",quote=FALSE)


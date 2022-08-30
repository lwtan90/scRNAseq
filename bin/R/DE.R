### DE.R ###
# This script is created to perform differential gene expression analysis for single-cell data
# Required packages: Seurat, Signac (if running multiome / atac)
# The analysis also assumed that you have normalized the data using SCTransform.
# Output will be the typical differential output from Seurat.
# Note: ident.1 will be the "case" and ident.2 will be "reference/control", so the log2FC will be ident.1-ident.2.
# Feel free to change the logfc.threshold and min.pct. Read Seurat manual is needed
# Places to edit:
# 1. P1
# 2. P2
############

require(Seurat)
require(Signac)

# p1
# the R object has to be editted
heart.combined = readRDS("SCT.prepped.RDS")

## Optional if running using RStudio because you dont need user argumnt. This is to run on cluster
args = commandArgs(TRUE)

## just to see if the current identity of cells is labelled corrected.
## If incorrect,please use the following command (assumeing that the "celltype" column contains the identity that you are interested in:
## Idents(heart.combined) = heart.combined$celltype
table(Idents(heart.combined))

## This is optional
heart.combined = subset(heart.combined, idents = args[1])
heart.combined = subset(heart.combined, sample %in% c("Control","LMNA1"))
Idents(heart.combined) = heart.combined$sample
table(Idents(heart.combined))

## It is good to make sure that your active assay is SCT / RNA but not integrated
DefaultAssay(heart.combined) = "SCT"

## I sugegst running PrepSCTFindMarkers once, and store it as SCT.prepped.RDS because this takes a long time.
## Next you just read the "prepped" R object as in p1 above.
##heart.combined <- PrepSCTFindMarkers(heart.combined)

## The actual step for differential analysis. They use wilcoxon rank sum test.
heart.combined.markers <- FindMarkers(heart.combined, assay="SCT", ident.1 = "LMNA1", ident.2 = "Control", min.pct = 0.2, logfc.threshold = 0.2)

## P2
## please reformat the outputfile
write.table(heart.combined.markers,file=paste("DE_control_LMNA1_",args[1],".txt",sep=""),sep="\t",quote=FALSE)


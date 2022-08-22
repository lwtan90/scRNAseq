require(Seurat)
require(Signac)
require(Seurat)
require(ggplot2)
require(sctransform)

options(bitmapType="cairo")

args = commandArgs(TRUE)
celltype = args[1]

### read the RData, feel free to change
heart = readRDS("../SCT.prepped.RDS")

### This script will subset the subtype listed in the mc
mc = subset(heart, idents = celltype)

### Typical Integration Protocol
mc_list = SplitObject(mc,split.by="sample")
mc_list = lapply(mc_list,SCTransform, vars.to.regress = "percent.mt")
features  <- SelectIntegrationFeatures(mc_list, nfeatures = 2000)
mc_list = PrepSCTIntegration(object.list=mc_list, anchor.features=features)
anchors <- FindIntegrationAnchors(mc_list, normalization.method = "SCT", anchor.features = features)
mc_combined = IntegrateData(anchorset= anchors, normalization.method="SCT")

### Intermediary file. DOnt worry
saveRDS(mc_combined,file=paste(celltype,".RDS",sep=""))

### cluster at resolution of 0.5, PC30
mc_combined <- RunPCA(mc_combined)
mc_combined <- RunUMAP(mc_combined, reduction = "pca", dims = 1:30)
mc_combined <- FindNeighbors(mc_combined, dims = 1:30)
mc_combined <- FindClusters(mc_combined, resolution=0.5)

#### Final output for all the subsequent analysis
saveRDS(mc_combined,file=paste(celltype,".RDS",sep=""))


### Visualize the UMAP
png(paste(celltype,"_UMAP.png",sep=""),width=2000,height=1800,res=300)
p2 <- DimPlot(mc_combined, reduction = "umap", label = TRUE, repel = TRUE)
print(p2)
dev.off()



### Integrate.R ###
# This script performs integration of single-cell datasets prior to differential gene expression analysis
# Need to tailor this script depending on the number of samples, and also whether it is multiome or just scRNA-seq
# Here I am running integration on multiome datasets for DCM with 3 samples (three condition)
# If you are running scRNA-seq, just run Section 1.
# If you are running 10x Multiome, run Section 1 and Section 2.
###################

### Loaing required packages
## for scRNA-seq
library(Seurat)
library(ggplot2)
library(patchwork)
require(data.table)
require(sctransform)
library(glmGamPoi)

## for Multiome
library(Signac)

## change this accoarding to organism
# load this just to get one of the function required
library(EnsDb.Rnorvegicus.v79)
library(BSgenome.Rnorvegicus.UCSC.rn6)

## For PNG file
options(bitmapType="cairo")
set.seed(1234)

### If dealing with mouse, change GTF path
gtf = rtracklayer::import("/labs/joewu/wlwtan/annotation/hg38/arc/gencode.v38.filtered.annotation.biotype.gtf")
gene.coords <- gtf[gtf$type == 'gene']
seqlevelsStyle(gene.coords) <- 'UCSC'
annotation <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

### Load individual libraries (must be sctransformed)
control = readRDS("../SS5/JXZ43/seurat/LSI_final_heart_ATAC.RDS")
control

LMNA1 = readRDS("../SS6/JXZ43/seurat/LSI_final_heart_ATAC.RDS")
LMNA1

LMNA2 = readRDS("../SS6/JXZ43/seurat/LSI_final_heart_ATAC.RDS")
LMNA2

## change this accordingly
control$dataset <- "control"
LMNA1$dataset <- "LMNA1"
LMNA2$dataset <- "LMNA2"


heart.list <- list(control=control,LMNA1=LMNA1, LMNA2=LMNA2)
features <- SelectIntegrationFeatures(object.list=heart.list,nfeatures=3000)
heart.list = PrepSCTIntegration(object.list=heart.list, anchor.features=features)
heart.anchors <- FindIntegrationAnchors(object.list=heart.list, normalization.method="SCT", anchor.features=features)
heart.combined <- IntegrateData(anchorset=heart.anchors,normalization.method="SCT")
save(heart.combined,features,file="step2.merged.RData")
####
heart.combined <- RunPCA(heart.combined)
heart.combined <- RunUMAP(heart.combined,reduction="pca",dims=1:30)
heart.combined <- FindNeighbors(heart.combined,reduction="pca",dims=1:30)
heart.combined <- FindClusters(heart.combined,resolution=0.5)
heart.combined$sample = rep("Control")
heart.combined$sample[grep("_2",rownames(heart.combined@meta.data))] = "LMNA1"
heart.combined$sample[grep("_3",rownames(heart.combined@meta.data))] = "LMNA2"
save(heart.combined,file="step3.merged.RData")
}

if(0){
######load("step3.merged.RData")
load("step3.merged.RData")
png("UMAP_RNA_sample.png",width=4000,height=1200,res=300)
p2 <- DimPlot(heart.combined, reduction = "umap", label = TRUE, repel = TRUE, split.by="sample")
print(p2)
dev.off()
quit()
}

if(1){
##load("step3.merged.RData")
heart.combined = readRDS("SCT.prepped.RDS")
#heart.combined$sample = rep("control")
#heart.combined$sample[grep("_2",rownames(heart.combined@meta.data))] = "LMNA1"
### reassigned
#heart.combined <- FindClusters(heart.combined,resolution=0.5)
##heart.combined$regroup = heart.combined$seurat_clusters
####heart.combined$regroup = rep(-1)
####heart.combined$regroup[ heart.combined$seurat_clusters %in% c(1,4,6,11,13,16)] = "CM"
####heart.combined$regroup[ heart.combined$seurat_clusters %in% c(22,14,15,3,5)] = "MC"
####heart.combined$regroup[ heart.combined$seurat_clusters %in% c(7,12,21)] = "Muscle-like"
####heart.combined$regroup[ heart.combined$seurat_clusters %in% c(9,23)] = "SMC"
####heart.combined$regroup[ heart.combined$seurat_clusters %in% c(0,17)] = "FB"
####heart.combined$regroup[ heart.combined$seurat_clusters %in% c(10,2,8)] = "EC"
####heart.combined$regroup[ heart.combined$seurat_clusters %in% c(18)] = "UNK1"
####heart.combined$regroup[ heart.combined$seurat_clusters %in% c(19)] = "UNK2"
####heart.combined$regroup[ heart.combined$seurat_clusters %in% c(20)] = "UNK3"
####Idents(heart.combined) = heart.combined$regroup
heart.combined$regroup[ heart.combined$regroup %in% c("UNK3")] = "Neuron"
heart.combined$regroup[ heart.combined$regroup %in% c("Muscle-like")] = "Pericyte"
heart.combined$regroup[ heart.combined$regroup %in% c("SMC")] = "Tcell"
heart.combined$regroup[ heart.combined$regroup %in% c("UNK1")] = "Endocardial"
heart.combined$regroup[ heart.combined$regroup %in% c("UNK2")] = "Lymphatic_EC"
Idents(heart.combined) = heart.combined$regroup
table(Idents(heart.combined))
#png("UMAP_RNA_regrouped.png",width=1800,height=1500,res=300)
#p1 <- DimPlot(heart.combined, reduction = "umap", label = TRUE, repel = TRUE)
#p1
#dev.off()
save(heart.combined,file="step4_regrouped.RData")
quit()}



###load("step3.merged.RData")
###
###celltype = read.table("celltype.marker.txt",sep="\t")
###names(celltype) = c("gene","celltype","hgnc")
###png("MARKER_featureplot.png",width=6000,height=6000,res=300)
###FeaturePlot(heart.combined, features=celltype$gene,raster=FALSE)
###dev.off()
###
###
###
###### Assign cell types
##### classify cell type
###heart.rna = heart.combined
###Idents(heart.rna) = paste("RNA",heart.rna$integrated_snn_res.1,sep="_") ## please check
###avg.heart.rna = as.data.frame(log1p(AverageExpression(heart.rna)$integrated))
###head(avg.heart.rna)
##### read in marker
###ct = data.frame(celltype=unique(celltype$celltype))
###cf = c("red","orange","yellow","green","blue","magenta","purple","black","grey","brown","turquoise")
###ct$color = cf[as.numeric(as.factor(ct$celltype))]
###celltype$color = ct$color[match(celltype$celltype,ct$celltype)]
###
##### extract DE
###avg.heart.rna = avg.heart.rna[ rownames(avg.heart.rna)%in%celltype$gene, ]
###avg.heart.rna = avg.heart.rna[ rowSums(avg.heart.rna>0)>0, ]
###
###hc = hclust(as.dist(1-cor(avg.heart.rna)),method="ward.D2")
###hr = hclust(as.dist(1-cor(t(avg.heart.rna))),method="ward.D2")
###
###celltype = celltype[match(rownames(avg.heart.rna),celltype$gene),]
###rownames(avg.heart.rna) = celltype$hgnc[match(rownames(avg.heart.rna),celltype$gene)]
###
###require(gplots)
###png("marker_heatmap_subtype.png",width=2000,height=3500,res=300)
###heatmap.2(as.matrix(avg.heart.rna),Rowv=as.dendrogram(hr),Colv=as.dendrogram(hc),col=redgreen(100),scale="row",trace="none",density.info="none", RowSideColors=as.character(celltype$color))
###dev.off()
###
###
###
##### name it heart because you can override later in case u make a LMNA1stake
##### Assign the cell types based on RNA integration
##### Find DE between condition for MI
###heart <- RenameIdents(
###        object = heart.combined,
###        "0"="CM",
###        "1"="CM",
###        "2"="CM",
###        "3"="EC",
###        "4"="CM",
###        "5"="CM",
###        "6"="PC",
###        "7"="FB",
###        "8"="CM",
###        "9"="CM",
###        "10"="CM",
###        "11"="EC",
###        "12"="EC",
###        "13"="MC",
###        "14"="CM",
###        "15"="EC"
###)
###heart$celltype = Idents(heart)
###png("UMAP_RNA_sample.png",width=3000,height=1500,res=300)
###p1 <- DimPlot(heart, reduction = "umap", group.by = "celltype")
###p2 <- DimPlot(heart, reduction = "umap", label = TRUE, repel = TRUE)
###p1 + p2
###dev.off()
###save(heart,file="step4.named.RData")
###

###load("step4.named.RData")
###
##### Transfer the labels from RNA to each ATAC
#### Match the celltype identity to control / LMNA1
#### apparently the order matters
#### if you merge control first, the cell name will have a "_1"
##### check if this is true, some ATAC barcode LMNA1ght not be in control
###load(file="step1.P1.Sham.RData")
###load(file="step1.P1.MI.RData")
###
###control$label = paste(rownames(control@meta.data),"_1",sep="")
###control$celltype = Idents(heart)[match(control$label,names(Idents(heart)))]
###
###LMNA1$label = paste(rownames(mi@meta.data),"_2",sep="")
###LMNA1$celltype = Idents(heart)[match(mi$label,names(Idents(heart)))]
###
###
##### Integrate Sham and MI ATAC
##### Merge the peaks list from Sham and MI first
###DefaultAssay(control) <- "peaks"
###control.peaks = granges(control)
###DefaultAssay(LMNA1) <- "peaks"
###LMNA1.peaks = granges(mi)
###peaks = reduce(c(control.peaks,LMNA1.peaks))
###
##### Recount the peaks in both datasets
##### Make [["peaks"]] entry
###reCountATAC <- function(obj,fragpath)
###{
###	counts <- FeatureMatrix(
###		fragments = Fragments(obj),
###		features = peaks,
###		cells = colnames(obj)
###	)
###
###	# create a new assay using the MACS2 peak set and add it to the Seurat object
###	obj[["peaks"]] <- CreateChromatinAssay(
###		counts = counts,
###		fragments = fragpath,
###		annotation = annotation
###	)
###	obj <- FindTopFeatures(obj, LMNA1n.cutoff = 10)
###	obj <- RunTFIDF(obj)
###	obj <- RunSVD(obj)
###	return(obj)
###	
###}
###control <- reCountATAC(control,"../JXZ43/JXZ43/outs/atac_fragments.tsv.gz")
###LMNA1   <- reCountATAC(mi,"../JXZ44/JXZ44/outs/atac_fragments.tsv.gz")
###
###save(control,file="step5.P1.Sham.RData")
###save(LMNA1,file="step5.P1.MI.RData")
###
load("step5.P1.Sham.RData")
load("step5.P1.MI.RData")
###
#### first add dataset-identifying metadata
control$dataset <- "P1_CT"
LMNA1$dataset <- "P1_MI"
###
###
#### merge
###atac.combined <- merge(control,LMNA1)
###
#### process the combined dataset
###atac.combined <- FindTopFeatures(atac.combined, LMNA1n.cutoff = 10)
###atac.combined <- RunTFIDF(atac.combined)
###atac.combined <- RunSVD(atac.combined)
###atac.combined <- RunUMAP(atac.combined, reduction = "lsi", dims = 2:30)
###save(atac.combined,file="step6.ataccombined.RData")
###
###png("UMAP_merged.png",width=1500,height=2000,res=300)
###p1 <- DimPlot(atac.combined, group.by = "dataset")
###p1
###dev.off()
###
###
#### find integration anchors
###integration.anchors <- FindIntegrationAnchors(
###	object.list = list(control, LMNA1),
###	anchor.features = rownames(control),
###	reduction = "rlsi",
###	dims = 2:30
###)
###save(integration.anchors, file="step7.integration.RData")
###
#### integrate LSI embeddings
###integrated <- IntegrateEmbeddings(
###	anchorset = integration.anchors,
###	reductions = atac.combined[["lsi"]],
###	new.reduction.name = "integrated_lsi",
###	dims.to.integrate = 1:30
###)
###save(integrated, file="step8.integrated.RData")
###
#### create a new UMAP using the integrated embeddings
###integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
###p1 <- DimPlot(integrated, group.by = "dataset")
###p2 <- DimPlot(integrated, group.by = "celltype")
###png("UMAP_ATAC_sample.png",width=3000,height=1500,res=300)
###p1 + p2
###dev.off()
###save(integrated, file="step8.integrated.RData")
###
###
#### perform label transfer
#### quantify gene activity
###load("step8.integrated.RData")
###load("step4.named.RData")
###
###gene.activities <- GeneActivity(integrated, features = VariableFeatures(heart))
###
#### add gene activities as a new assay
###integrated[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
###
#### normalize gene activities
###DefaultAssay(integrated) <- "ACTIVITY"
###integrated <- NormalizeData(integrated)
###integrated <- ScaleData(integrated, features = rownames(integrated))
###
###transfer.anchors <- FindTransferAnchors(reference = heart, query = integrated, features = VariableFeatures(object = heart),
###    reference.assay = "integrated", query.assay = "ACTIVITY", reduction = "cca")
###
###
###celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = heart$celltype,
###    weight.reduction = integrated[["integrated_lsi"]], dims = 2:30)
###
###integrated <- AddMetaData(integrated, metadata = celltype.predictions)
###
##### Perform Differential RNA analysis
###p1 <- DimPlot(integrated, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("label transfer")
###png("UMAP_ATAC_label_transfer.png",width=3000,height=1500,res=300)
###p1
###dev.off()
###
###save(integrated, file = "step9_final_labelled.object.RData")

## Perform Differential RNA analysis
load("step8.integrated.RData")
load("step4.named.RData")
Idents(heart) = heart$celltype
DefaultAssay(heart) = "RNA"
heart$celltype.dataset <- paste(Idents(heart), integrated$samplename, sep = "_")
Idents(heart) <- "celltype.dataset"

require(ggplot2)
celltype.list = unique(heart$celltype)

for(i in 1:length(celltype.list))
{
	for(j in (i+1):length(celltype.list))
	{
		print(paste(celltype.list[i],celltype.list[j]))
		group1 = paste(celltype.list[i],"_P1_MI",sep="")
		group2 = paste(celltype.list[j],"_P1_MI",sep="")
		filename = paste("RNA_marker_",celltype.list[i],"_",celltype.list[j],".txt",sep="")
		DA <- FindMarkers(heart, ident.1 = group1, ident.2 = group2, verbose = TRUE)
		write.table(DA,file=filename,sep="\t",quote=FALSE)
	}
}
quit()

for(x in celltype.list)
{
	print(x)
	group1 = paste(x,"_P1_MI",sep="")
	group2 = paste(x,"_P1_CT",sep="")
	filename = paste("RNA_DA_",group2,group1,".txt",sep="")
	volfile = paste("RNA_volcano_",group2,group1,".png",sep="")
	pctfile = paste("RNA_PCT_",group2,group1,".png",sep="")
	volfile2 = paste("RNA_volcano_",group2,group1,"_diffproportion.png",sep="")
	DA <- FindMarkers(heart, ident.1 = group1, ident.2 = group2, verbose = TRUE, LMNA1n.pct=0, logfc.threshold=0.1)
	write.table(DA,file=filename,sep="\t",quote=FALSE)
	DA$sig = DA$p_val_adj < 0.05
	p1 <- ggplot(DA,aes(x=avg_log2FC,y=-log10(p_val)))+geom_point(aes(color=sig))+theme_bw()+theme(panel.grid=element_blank())+xlab(paste("logFC",group1,"-",group2))+ylab("-log10(P)")+scale_color_manual(values=c("grey30","red"))
	p1 <- p1 + geom_hline(yintercept=0) + geom_vline(xintercept=0)
	png(volfile,width=1800,height=1500,res=300)
	print(p1)
	dev.off()
	p1 <- ggplot(DA,aes(x=avg_log2FC,y=-log10(p_val)))+geom_point(aes(color=sig,size=abs(pct.1-pct.2)))+theme_bw()+theme(panel.grid=element_blank())+xlab(paste("logFC",group1,"-",group2))+ylab("-log10(P)")+scale_color_manual(values=c("grey30","red"))
	p1 <- p1 + geom_hline(yintercept=0) + geom_vline(xintercept=0)
	png(volfile2,width=1800,height=1500,res=300)
	print(p1)
	dev.off()
	p1 <- ggplot(DA,aes(x=pct.1,y=pct.2))+geom_point(aes(color=sig,size=-log10(p_val)))+theme_bw()+theme(panel.grid=element_blank())+xlab(paste("Proportion ",group1))+ylab(paste("Proportion ",group2))+scale_color_manual(values=c("grey30","red"))
	p1 <- p1 + geom_abline(slope=1,intercept=0,linetype="dashed")
	png(pctfile,width=1800,height=1500,res=300)
	print(p1)
	dev.off()
	
}


## Perform Differential ATAC analysis
load("step8.integrated.RData")
Idents(integrated) = integrated$celltype
DefaultAssay(integrated) = "peaks"
integrated$celltype.dataset <- paste(Idents(integrated), integrated$dataset, sep = "_")
Idents(integrated) <- "celltype.dataset"

require(ggplot2)
celltype.list = unique(integrated$celltype)
for(x in celltype.list)
{
	print(x)
	filename = paste(x,"_DA.txt",sep="")
	group1 = paste(x,"_P1_MI",sep="")
	group2 = paste(x,"_P1_CT",sep="")
	filename = paste("ATAC_DA_",group2,group1,".txt",sep="")
	volfile = paste("ATAC_volcano_",group2,group1,".png",sep="")
	pctfile = paste("PCT_",group2,group1,".png",sep="")
	volfile2 = paste("ATAC_volcano_",group2,group1,"_diffproportion.png",sep="")
	DA <- FindMarkers(integrated, ident.1 = group1, ident.2 = group2, verbose = TRUE, LMNA1n.pct=0, logfc.threshold=0.1)
	write.table(DA,file=filename,sep="\t",quote=FALSE)
	DA$sig = DA$p_val_adj < 0.05
	p1 <- ggplot(DA,aes(x=avg_log2FC,y=-log10(p_val)))+geom_point(aes(color=sig))+theme_bw()+theme(panel.grid=element_blank())+xlab(paste("logFC",group1,"-",group2))+ylab("-log10(P)")+scale_color_manual(values=c("grey30","red"))
	p1 <- p1 + geom_hline(yintercept=0) + geom_vline(xintercept=0)
	png(volfile,width=1800,height=1500,res=300)
	print(p1)
	dev.off()
	p1 <- ggplot(DA,aes(x=avg_log2FC,y=-log10(p_val)))+geom_point(aes(color=sig,size=abs(pct.1-pct.2)))+theme_bw()+theme(panel.grid=element_blank())+xlab(paste("logFC",group1,"-",group2))+ylab("-log10(P)")+scale_color_manual(values=c("grey30","red"))
	p1 <- p1 + geom_hline(yintercept=0) + geom_vline(xintercept=0)
	png(volfile2,width=1800,height=1500,res=300)
	print(p1)
	dev.off()
	p1 <- ggplot(DA,aes(x=pct.1,y=pct.2))+geom_point(aes(color=sig,size=-log10(p_val)))+theme_bw()+theme(panel.grid=element_blank())+xlab(paste("Proportion ",group1))+ylab(paste("Proportion ",group2))+scale_color_manual(values=c("grey30","red"))
	p1 <- p1 + geom_abline(slope=1,intercept=0,linetype="dashed")
	png(pctfile,width=1800,height=1500,res=300)
	print(p1)
	dev.off()
	
}


## Correlate RNA and ATAC DE for each population (promoter) (https://satijalab.org/signac/articles/pbmc_vignette.html)
load("step8.integrated.RData")
celltype.list = unique(integrated$celltype)
for(x in celltype.list)
{
	print(x)
	group1 = paste(x,"_P1_MI",sep="")
	group2 = paste(x,"_P1_CT",sep="")
	filename = paste("ATAC_DA_",group2,group1,".txt",sep="")
	rnafilename = paste("RNA_DA_",group2,group1,".txt",sep="")
	outputfile = paste("RNA_",filename,sep="")
	scatterfile = paste("logFC_scatter_",group1,group2,".png",sep="")
	DA = read.table(filename,header=TRUE)
	DE = read.table(rnafilename,header=TRUE) 
	region = rownames(DA)
	region.gene.control = ClosestFeature(control, regions = region)
	region.gene.LMNA1 = ClosestFeature(mi, regions = region)
	DA = cbind(DA,region.gene.control[match(rownames(DA),region.gene.control$query_region),])
	DA$rna.logFC = DE$avg_log2FC[match(DA$gene_name,rownames(DE))]
	DA$rna.p_val_adj = DE$p_val_adj[match(DA$gene_name,rownames(DE))]
	
	### focus on promoter
	filtered = DA[!is.na(DA$gene_name) & !is.na(DA$avg_log2FC),]
	print(dim(filtered))
	##filtered = filtered[ filtered$p_val_adj<0.05 & filtered$rna.p_val_adj<0.05, ]
	filtered = filtered[ filtered$distance<10000, ]
	p1 <- ggplot(filtered,aes(x=avg_log2FC,y=rna.logFC))+geom_point()+theme_bw()+theme(panel.grid=element_blank())
	p1 <- p1 + xlab("ATAC logFC") + ylab("RNA logFC") + geom_hline(yintercept=0)+geom_vline(xintercept=0)
	png(scatterfile,width=2000,height=2000,res=300)
	print(p1)
	dev.off()
	write.table(DA,file=outputfile,sep="\t",quote=FALSE)
}


## plot heatmap
## extract the CM only
require(reshape)
require(scales)
require(ggplot2)

for(x in celltype.list)
{
	print(x)
	CM = subset(integrated, idents = x)
	group1 = paste(x,"_P1_MI",sep="")
	group2 = paste(x,"_P1_CT",sep="")
	filename = paste("ATAC_DA_",group2,group1,".txt",sep="")
	plotfile = paste("ATAC_DA_",group1,group2,"_heatmap.png",sep="")
	DA = read.table(filename,header=TRUE)
	DA = DA[ DA$p_val_adj<0.05, ]

	## export the count into matrix
	count = CM@assays$peaks@data[rownames(CM@assays$peaks@data) %in% rownames(DA), ]
	dim(count)

	## record the grouping of each cell
	plotdata = melt(as.matrix(count))
	names(plotdata) = c("peak","cell","count")
	plotdata$condition = CM@meta.data$dataset[match(plotdata$cell,rownames(CM@meta.data))]
	plotdata$logFC = DA$avg_log2FC[match(plotdata$peak,rownames(DA))]

	plotdata$count[ plotdata$count>4]=4
	p1 <- ggplot(plotdata,aes(x=cell,y=reorder(peak,logFC)))+geom_tile(aes(fill=count))+theme_bw()+theme(panel.grid=element_blank())
	p1 <- p1 + facet_grid(.~condition,space="free",scale="free") + theme(axis.text=element_blank(),axis.ticks=element_blank())
	p1 <- p1 + scale_fill_gradientn(values=rescale(c(0,1,2,3,4)),colors=c("LMNA1dnightblue","blue","turquoise","white","yellow"))
	png(plotfile,width=3000,height=3000,res=300)
	print(p1)
	dev.off()
}

## perform motif analysis
load("step9_final_labelled.object.RData")
DefaultAssay(integrated) = "peaks"

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(x = JASPAR2020,opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
save(pfm,file="PFM.JASPAR2020.RData")

# add motif information
Idents(integrated) = integrated$celltype

integrated <- AddMotifs(
  object = integrated,
  genome = BSgenome.Rnorvegicus.UCSC.rn6,
  pfm = pfm
)
save(integrated,file="step10_motif_final_labelled.object.RData")

## Read in DA analysis
pfm.data = data.frame(id=c(),symbol=c())
for(i in 1:length(pfm@listData))
{
	pfm.data = rbind(pfm.data,data.frame(id=c(names(pfm@listData)[i]),symbol=c(pfm@listData[[i]]@name)))
}
pfm.data


celltype.list = unique(integrated$celltype)
for(x in celltype.list)
{
	print(x)
	group1 = paste(x,"_P1_MI",sep="")
	group2 = paste(x,"_P1_CT",sep="")
	filename = paste("ATAC_DA_",group2,group1,".txt",sep="")
	DA = read.table(filename,header=TRUE)
	DA = DA[ DA$p_val_adj<0.05, ]
	## perform motif analysis without background
	DA.motifs <- FindMotifs(
	  object = integrated,
	  features = rownames(DA)
	)
	write.table(DA.motifs,file=paste(x,"_motif_nobg.txt"),sep="\t",quote=FALSE,row.names=F)
	
	p1 <- MotifPlot(object = integrated,motifs = head(rownames(DA.motifs)))
	png(paste(x,"topmotifs_nobackground.png",sep="_"),width=2500,height=1500,res=300)
	print(p1)
	dev.off()
	
	## Identify background
	# find peaks open in Pvalb or Sst cells
	open.peaks <- AccessiblePeaks(integrated, idents = c(x))
	
	# match the overall GC content in the peak set
	meta.feature <- GetAssayData(integrated, assay = "peaks", slot = "meta.features")
	peaks.matched <- MatchRegionStats(
	  meta.feature = meta.feature[open.peaks, ],
	  query.feature = meta.feature[rownames(DA), ],
	  n = 50000
	)
	DA.motifs <- FindMotifs(
	  object = integrated,
	  features = rownames(DA),
	  background=peaks.matched
	)
	p1 <- MotifPlot(object = integrated,motifs = head(rownames(DA.motifs)))
	png(paste(x,"topmotifs_background.png",sep="_"),width=2500,height=1500,res=300)
	print(p1)
	dev.off()
	write.table(DA.motifs,file=paste(x,"_motif_bg.txt"),sep="\t",quote=FALSE,row.names=F)
}


## Run ChromVar
pfm.data = data.frame(id=c(),symbol=c())
for(i in 1:length(pfm@listData))
{
        pfm.data = rbind(pfm.data,data.frame(id=c(names(pfm@listData)[i]),symbol=c(pfm@listData[[i]]@name)))
}
pfm.data

integrated <- RunChromVAR(
  object = integrated,
  genome = BSgenome.Rnorvegicus.UCSC.rn6
)
DefaultAssay(integrated) <- 'chromvar'
save(integrated,file="step11_chromvar_final_labelled.object.RData")

## run differential activity analysis
### change the identity of cells into MI/CT
Idents(integrated) = integrated$celltype
integrated$celltype.dataset <- paste(Idents(integrated), integrated$dataset, sep = "_")
Idents(integrated) <- "celltype.dataset"

celltype.list = unique(integrated$celltype)
for(x in celltype.list)
{
        print(x)
        filename = paste(x,"_DA.txt",sep="")
        group1 = paste(x,"_P1_MI",sep="")
        group2 = paste(x,"_P1_CT",sep="")
        filename = paste("MOTIF_",group2,"_",group1,".txt",sep="")
	motiffile = paste("MOTIF_",group2,"_",group1,"_motif.png",sep="")
	DA <- FindMarkers(
  		object = integrated,
		ident.1 = group1,
		ident.2 = group2,
		only.pos = TRUE,
		mean.fxn = rowMeans,
		fc.name = "avg_diff"
	)
	DA$symbol = pfm.data$symbol[match(rownames(DA),pfm.data$id)]
	print(head(DA))
	write.table(DA,file=filename,sep="\t",quote=FALSE)

	p1 <- MotifPlot(object = integrated, motifs = head(rownames(DA)),assay = 'peaks')
	png(motiffile,width=2500,height=2500,res=300)
	print(p1)
	dev.off()
}



## optionals if seurat or homer (https://satijalab.org/signac/articles/motif_vignette.html)

## run cicero ("https://satijalab.org/signac/articles/cicero.html")

## Perform ChromVAR analysis

## cross talk analysis

#### This script is actually to look at doublet rates

library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
require(data.table)

## change this accoarding to organism
library(EnsDb.Rnorvegicus.v79)
library(BSgenome.Rnorvegicus.UCSC.rn5)


## For cluster
options(bitmapType="cairo")

set.seed(1234)
# load the RNA and ATAC data

counts <- Read10X_h5("../outs/filtered_feature_bc_matrix.h5")
fragpath <- "../outs/atac_fragments.tsv.gz"


gtf = rtracklayer::import("/labs/joewu/wlwtan/annotation/hg38/arc/gencode.v38.filtered.annotation.biotype.gtf")
gene.coords <- gtf[gtf$type == 'gene']
seqlevelsStyle(gene.coords) <- 'UCSC'
annotation <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')


heart <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
heart[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(heart) = "RNA"
heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "MT-")

qc = heart@meta.data
write.table(qc,file="QC.txt",sep="\t",quote=FALSE)

save(heart,file="heart.RData")


## QC
## nCount_RNA  = total molecules (high means doublet)
## nFeature_RNA = total genes detected (low means dead cells)
DefaultAssay(heart) = "ATAC"
heart <- NucleosomeSignal(heart)
heart <- TSSEnrichment(heart)


png("QC_violin_seurat.png",width=3000,height=2000,res=300)
VlnPlot(object=heart,features=c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC","TSS.enrichment", "nucleosome_signal"),ncol=3,pt.size=0)
dev.off()

save(heart,file="heart_qc.RData")
quit()

heart <- subset(
  x = heart,
  subset = nCount_RNA > 500 &
    nCount_RNA < 100000 &
    nCount_ATAC > 400 &
    nCount_ATAC < 90000 &
    nFeature_RNA > 400 &
    nFeature_RNA < 10000 &
    nFeature_ATAC > 400 &
    nFeature_ATAC < 90000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1.5
)

save(heart,file="heart_3.RData")
quit()

#### peak-calling
# call peaks using MACS2
peaks <- CallPeaks(heart, macs2.path = "macs2")
save(peaks,file="raw_peaks.RData")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
save(peaks,file="filtered_peaks.RData")


# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(heart),
  features = peaks,
  cells = colnames(heart)
)
save(macs2_counts,file="heart_masc2.RData")

# create a new assay using the MACS2 peak set and add it to the Seurat object
heart[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

save(heart,file="final_heart_ATAC.RData")



### Process RNA
heart.rna = heart
DefaultAssay(heart.rna) <- "RNA"
heart.rna = NormalizeData(heart.rna)
heart.rna <- FindVariableFeatures(heart.rna)
heart.rna <- ScaleData(heart.rna)
heart.rna <- RunPCA(heart.rna)
heart.rna <- RunUMAP(heart.rna, dims = 1:30)
png("Elbow_RNA.png",width=2500,height=2500,res=300)
ElbowPlot(heart.rna)
dev.off()
heart.rna <- FindNeighbors(heart.rna, dims = 1:15)
heart.rna <- FindClusters(heart.rna, resolution = 1)
p1 <- DimPlot(heart.rna, group.by = "RNA_snn_res.1", label = TRUE) + NoLegend() + ggtitle("RNA")

png("UMAP_RNA.png",width=1500,height=1500,res=300)
print(p1)
dev.off()

save(heart.rna, file="RNA.UMAP.RData")


### Process ATAC
heart.atac = heart
DefaultAssay(heart.atac) <- "peaks"
heart.atac <- FindTopFeatures(heart.atac, min.cutoff = 5)
heart.atac <- RunTFIDF(heart.atac)
heart.atac <- RunSVD(heart.atac)
heart.atac <- RunUMAP(heart.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
heart.atac <- FindNeighbors(object = heart.atac, reduction = 'lsi', dims = 2:30)
heart.atac <- FindClusters(object = heart.atac, verbose = FALSE, algorithm = 3)
save(heart.atac,file="clustered_heart.atac.RData")
p2 <- DimPlot(heart.atac, group.by = "peaks_snn_res.0.8", label = FALSE) + NoLegend() + ggtitle("ATAC")

png("UMAP_ATAC.png",width=1500,height=1500,res=300)
print(p2)
dev.off()

png("UMAP_RNA_ATAC.png",width=2500,height=1800,res=300)
p1 + p2
dev.off()

save(heart.rna,heart.atac,heart,file="named_clustered_heart.RData")


### find anchors
# quantify gene activity
## make sure this is in "peaks" mode
gene.activities <- GeneActivity(heart.atac, features = VariableFeatures(heart.rna))

# add gene activities as a new assay
heart.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(heart.atac) <- "ACTIVITY"
heart.atac <- NormalizeData(heart.atac)
heart.atac <- ScaleData(heart.atac, features = rownames(heart.atac))

# find anchor
# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = heart.rna, query = heart.atac, features = VariableFeatures(object = heart.rna), reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")


## add annotation to ATAC metadata
heart.rna$annotations = Idents(heart.rna)
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = heart.rna$annotations,weight.reduction = heart.atac[["lsi"]], dims = 2:30)
heart.atac <- AddMetaData(heart.atac, metadata = celltype.predictions)




## replotting the UMAP
p2 <- DimPlot(heart.atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("ATAC")
png("labelled_UMAP.png",width=2500,height=1800,res=300)
p1 + p2
dev.off()

save(heart,heart.atac,heart.rna,file="heart.use.for.downstream.RData")

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(heart.rna)
refdata <- GetAssayData(heart.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = heart.atac[["lsi"]],dims = 2:30)
heart.atac[["RNA"]] <- imputation

coembed <- merge(x = heart.rna, y = heart.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

DimPlot(coembed, group.by = c("orig.ident", "seurat_annotations"))

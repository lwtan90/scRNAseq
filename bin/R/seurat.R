#### seurat.R #######
# This scripts takes raw output from cellranger and run standard seurat pipeline.
# It can also be adapted to run output from cellranger-arc.
# Section 1: Required packages.
# Section 2: Reading output from cellranger
# Section 3: QC
# Section 3.9: Thresholding and clean up Seurat object (Key step)
# Section 4: Peak calling for ATAC (if multiome)
# Section 5: Scale and Normalize / sctransform
# Section 6: UMAP / FindNeighbours / Cluster
# Section 7: Other potential functions
# The sections can be modified to suit your purpose.
# We rarely rely on the output from this analysis to come out with anything meaningful.
# You can use the RDS from Section 6 for integration, which is more meaningful.
#### This script is actually to look at doublet rates

#### Section 1: Package loading
## Core packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)
library(sctransform)
library(glmGamPoi)

## This is for the utility of some genomeDB functions for atac processing
## Note: Ideally you want to load corresponding packages for the organism you are working on, but it doesnt matter
library(EnsDb.Rnorvegicus.v79)
library(BSgenome.Rnorvegicus.UCSC.rn5)

## if multiome processing
library(Signac) # This is required for multiome processing

## Becuase th SCG cluster lacks PNG driver, this options flag is only for cluster. Ignore if using RStudio
options(bitmapType="cairo")

## This is to ensure that the random number generation is stable. 
set.seed(1234)

## Loading the gene annotation for scATAC-seq
## I prefer to load from custom gtf file to keep it consistent with the cellranger-arc database.
## You can try with the data from EnsDb R package, but I discourage it.
gtf = rtracklayer::import("/labs/joewu/wlwtan/annotation/hg38/arc/gencode.v38.filtered.annotation.biotype.gtf")
gene.coords <- gtf[gtf$type == 'gene']
seqlevelsStyle(gene.coords) <- 'UCSC'
annotation <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')



#### Section 2: Loading output from cellranger/cellranger-arc
# load the RNA (and/or) ATAC data from cellranger
# HDF5 file
counts <- Read10X_h5("../outs/filtered_feature_bc_matrix.h5")

## This path declaration is only for the ATAC/Multiome. can be commented out if running scRNA-seq
fragpath <- "../outs/atac_fragments.tsv.gz"

## Creation of Seurat object (Multiome).
## Note: I have yet to test this command on scRNA-seq alone.
# This will create RNA slot in seurat object
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



#### Section 3: Quality control (initial scan)
# This section is mainly to assess the general performance of each library.
# The actual threshold setting should be performed with the other replicates.
## First, we will assess the quality of the RNA component.
## nCount_RNA  = total molecules (high means doublet)
## nFeature_RNA = total genes detected (low means dead cells)
## percent.mt = proportion of reads from mitochondrial transcription (sign of nuclei damage if high)
# Some gene annotation uses other combination of naming for mitochondrial genes, eg "Mt-", "mt-".
# Please be advise to edit the pattern="MT-" accordingly.
DefaultAssay(heart) = "RNA"
heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "MT-")

### QC for ATAC (if multiome)
DefaultAssay(heart) = "ATAC"
heart <- NucleosomeSignal(heart)
heart <- TSSEnrichment(heart)

## Output for combined QC assessment with other replicates.
## I suggest saving the output in RDS, and continue with the Section 3.9 with the RDS object.
qc = heart@meta.data
write.table(qc,file="QC.txt",sep="\t",quote=FALSE)

png("QC_violin_seurat.png",width=3000,height=2000,res=300)
VlnPlot(object=heart,features=c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC","TSS.enrichment", "nucleosome_signal"),ncol=3,pt.size=0)
dev.off()

saveRDS(heart,file="heart_qc.RDS")

#### Section 3.9 Subset the seurat objects based on combined thresholding
# Use the figures generated by qc.R to decide on various parameters.
# The followings are tailored for one experiment, not universal.
# The DefaultAssay doesn't matter. 
heart = readRDS("heart_qc.RDS") ## unless you continue the qc.R analysis in a parallel server, just load the RDS object generated from step 3.
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
## This RDS will be used for Step 4 onwards.
saveRDS(heart,file="seurat_cleaned.RDS")


#### Section 4: Peak-calling for ATAC (if multiome)
# call peaks using MACS2
heart = readRDS("seurat_cleaned.RDS") ## unless you continue the step3.9 analysis in a parallel server, just load the RDS object generated from step 3.9.
DefaultAssay(heart) = "ATAC"
peaks <- CallPeaks(heart, macs2.path = "macs2")
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
saveRDS(peaks,file="filtered_peaks.RDS")

# quantify counts in each peak and save in "peaks" slot
macs2_counts <- FeatureMatrix(
  fragments = Fragments(heart),
  features = peaks,
  cells = colnames(heart)
)
saveRDS(macs2_counts,file="heart_masc2.RDS")
# create a new assay using the MACS2 peak set and add it to the Seurat object
heart[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
saveRDS(heart,file="final_heart_ATAC.RDS")



## Section 5: Scale and Normalize / sctransform
## Preferred Workflow (SCTransform)
## This is compatibale with all the downstream scripts.
DefaultAssay(heart) <- "RNA"
heart <- SCTransform(heart, method="glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(heart, file="SCTransformed.RDS")
## Alternative, to activate, just change if(0) to if(1)
## This is incompatible with the downstream scripts
if(0){
	heart = NormalizeData(heart)
	heart <- FindVariableFeatures(heart)
	heart <- ScaleData(heart)
	saveRDS(heart, file="scaled_normalized.RDS")
}

### Section 6: UMAP / FindNeighbours / Cluster
heart <- RunPCA(heart)
heart <- RunUMAP(heart, dims = 1:30)

## Elbow plot to determine number of PCs to detain
## Usually around 10-15
png("Elbow_RNA.png",width=2500,height=2500,res=300)
ElbowPlot(heart)
dev.off()

heart <- FindNeighbors(heart, dims = 1:15)
## Finding clusters. I find that resolution between 0.2 to 0.5 gives better clusters than 1 (too many).
heart <- FindClusters(heart, resolution = 0.5)
Idents(heart) = heart$seurat_clusters

## Visualize the UMAP with defined clusters
## By default, the Identity of cells will be displayed. To set, use Ident(heart) = heart$...
p1 <- DimPlot(heart, label = TRUE) + NoLegend() + ggtitle("RNA")
png("UMAP_RNA.png",width=1500,height=1500,res=300)
print(p1)
dev.off()

## If you are running scRNA-seq, this RDS can be used for integration analysis.
saveRDS(heart, file="SCTransformed.RDS")


#### Section 7: ATAC-analysis (Multiome)
## We will use LSI for the ATAC-seq integration, so it is important to run this step
DefaultAssay(heart) <- "peaks"
heart <- FindTopFeatures(heart, min.cutoff = 10)
heart <- RunTFIDF(heart)
heart <- RunSVD(heart)
heart <- RunUMAP(object = heart, reduction = 'lsi', dims = 2:30)
heart <- FindNeighbors(object = heart, reduction = 'lsi', dims = 2:30)
heart <- FindClusters(object = heart, verbose = FALSE, algorithm = 3)

## This can be used for Multiome integration analysis/
## Dont worry about the grouping because the final clustering from integration matters more.
saveRDS(heart,file="LSI_sctransformed_heart_RNA_ATAC.RDS")
p2 <- DimPlot(heart, label = TRUE) + NoLegend() + ggtitle("ATAC")
png("UMAP_RNA_ATAC.png",width=3500,height=1500,res=300)
p1 + p2
dev.off()

#########END OF SCRIPT###########
## The rest can be used for future purposes



### find anchors
# quantify gene activity
## make sure this is in "peaks" mode
gene.activities <- GeneActivity(heart, features = VariableFeatures(heart))

# add gene activities as a new assay
heart[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(heart) <- "ACTIVITY"
heart <- NormalizeData(heart)
heart <- ScaleData(heart, features = rownames(heart))

# find anchor
# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = heart, query = heart, features = VariableFeatures(object = heart), reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")


## add annotation to ATAC metadata
heart$annotations = Idents(heart)
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = heart$annotations,weight.reduction = heart[["lsi"]], dims = 2:30)
heart <- AddMetaData(heart, metadata = celltype.predictions)




## replotting the UMAP
p2 <- DimPlot(heart, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("ATAC")
png("labelled_UMAP.png",width=2500,height=1800,res=300)
p1 + p2
dev.off()

save(heart,heart,heart,file="heart.use.for.downstream.RData")

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(heart)
refdata <- GetAssayData(heart, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = heart[["lsi"]],dims = 2:30)
heart[["RNA"]] <- imputation

coembed <- merge(x = heart, y = heart)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

DimPlot(coembed, group.by = c("orig.ident", "seurat_annotations"))

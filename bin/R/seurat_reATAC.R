### This script will take the merged region from sctransform_ATAC.R script, and recount matrix for integeration purpose
## Places for change
# P1: path to the merged file
# P2: the genome annotation
# P3: the RDS file containing ATAC data

library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
require(data.table)

## change this accoarding to organism
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

options(bitmapType="cairo")

set.seed(1234)

## Key step: read in the merged features
#P1
merged.features = readRDS("../../../integration_atac/merged.feature.RDS")

## Check if the path to fragment file exists
fragpath <- "../outs/atac_fragments.tsv.gz"


## for annotation purpose
# P2
gtf = rtracklayer::import("/labs/joewu/wlwtan/annotation/hg38/arc/gencode.v38.filtered.annotation.biotype.gtf")
gene.coords <- gtf[gtf$type == 'gene']
seqlevelsStyle(gene.coords) <- 'UCSC'
annotation <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
annotation


### The RData here should be coming from the SCtrasnformed RNA data
# P3
heart = readRDS("sctransformed_reclustered.RDS")

### must change to ATAC
DefaultAssay(heart) <- "ATAC"

# start here
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(heart),
  # edit here to allow common range
  features = merged.features,
  cells = colnames(heart)
)
saveRDS(macs2_counts,file="merged_heart_masc2.RDS")

# create a new assay using the MACS2 peak set and add it to the Seurat object
heart[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
saveRDS(heart,file="merged_final_heart_ATAC.RDS")

#heart = readRDS("final_heart_ATAC.RDS")
DefaultAssay(heart) = "peaks"
heart <- FindTopFeatures(heart, min.cutoff = 10)
heart <- RunTFIDF(heart)
heart <- RunSVD(heart)
heart <- RunUMAP(object = heart, reduction = 'lsi', dims = 2:30)
heart <- FindNeighbors(object = heart, reduction = 'lsi', dims = 2:30)
heart <- FindClusters(object = heart, verbose = FALSE, algorithm = 3)
DimPlot(object = heart, label = TRUE) + NoLegend()


### This output should be used for all downstream analysis for ATAC-seq
saveRDS(heart,file="merged_LSI_final_heart_ATAC.RDS")


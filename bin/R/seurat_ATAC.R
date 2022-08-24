#### This script is actually to look at doublet rates

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

fragpath <- "../outs/atac_fragments.tsv.gz"

gtf = rtracklayer::import("/labs/joewu/wlwtan/annotation/hg38/arc/gencode.v38.filtered.annotation.biotype.gtf")
gene.coords <- gtf[gtf$type == 'gene']
seqlevelsStyle(gene.coords) <- 'UCSC'
annotation <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
annotation

### The RData here should be coming from the SCtrasnformed RNA data
heart = readRDS("sctransformed_reclustered.RDS")

### must change to ATAC
DefaultAssay(heart) <- "ATAC"

peaks <- CallPeaks(heart, macs2.path = "macs2")
saveRDS(peaks,file="raw_peaks.RDS")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
saveRDS(peaks,file="filtered_peaks.RDS")


# quantify counts in each peak
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

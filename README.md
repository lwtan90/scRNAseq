# single-cell Multiome / RNA analysis  
Author: Wilson Tan  
Date: Aug 2022  
Purpose: Document analysis workflow for scRNA-seq/Multiome  

  
## Alignment of FASTQ files  
Tool: cellranger (scRNA-seq) and cellrange-arc (multiome)  
Command:  
```
# for full command, see script bin/align_10x_cellrangerarc.sh

cellranger-arc count --id=$ID --reference=$REF --libraries=$LIB --localcores=8 --localmem=120

```  

  
## Quality control of the 10x sequencing data  
Tools: Seurat / cellrange output  
```
# for full command, see R script in bin/R/seurat.R

Rscript bin/R/seurat.R  

```  

### Description for this script  
This Rscript will read the output from 10x cellranger, and generate SCTransformed (normalized) data for analysis or integration.  

Usually, I prefer to set the threshold that is consistent across different replicates to ensure faire comparison, but this should be reviewed on a case-by-case basis.  

  
![Example QC of scRNA data](/images/QC.png)  
  
```  
Rscript bin/R/qc.R
```  
  
  
## Integrate Multiple Seurat Objects / Run for comparison (Downstream)  
This Rscript will read individual objects from different prep into an integrated object. After Integration, the data will be rescaled, normalized and dimension-reduced. This will be the starting point for all downstream analysis, and hypothesis generation.  

```  

Rscript bin/R/integrate.R

```  
  
All output should be saved in RDS / R data for other applications. Together in the R folder, includes various R script:  
- bin/R/integrate.R (integrate seurat objects, and normalize)  
- bin/R/findMarkers.R (identify cell type markers)  
- bin/R/DE.R (if pairwise comparison is needed, for differential analysis)  
- bin/R/dotplot.R (visualize the top 10 markers per cell type)  
- bin/R/aggregateCount.R (generate pseudobulk for DESeq2)  
- bin/R/monocle3.R (run pseudotime analysis)  
- bin/R/heatmap.R (plot heatmap)  
  
![Example QC of scRNA data](/images/UMAP.png)  



## Subset Cell Type of Interest for Further Analysis  
Sometimes you might want to look at specific cell type to decide if new clusters might be formed using a new sets of variable features.  
This script to use is in the bin/R/subset.R  
```  
Rscript subset.R
```    


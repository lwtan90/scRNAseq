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

 

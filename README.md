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

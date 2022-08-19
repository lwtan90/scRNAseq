#!/bin/bash

#### This pipeline will take FASTQ file only 
module load cellranger-arc


#### Global

ID=control
REF=/labs/joewu/wlwtan/annotation/hg38/arc/hg38
LIB=samplesheet.csv


cellranger-arc count --id=$ID --reference=$REF --libraries=$LIB --localcores=8 --localmem=120


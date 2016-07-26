#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N SORT_BEFORE_COMBINE
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=36gb
#PBS -q old

module load java
java -Xmx15g -jar /scratch/jandrews/bin/picard-tools-2.2.1/picard.jar SortVcf \
	I= /scratch/jandrews/Data/Variant_Calling/Coding/FINAL/merge_samtools.vcf  \
	O= /scratch/jandrews/Data/Variant_Calling/Coding/FINAL/merge_samtools.sorted.vcf \
	SEQUENCE_DICTIONARY= /scratch/jandrews/Ref/ucsc.hg19.dict

java -Xmx15g -jar /scratch/jandrews/bin/picard-tools-2.2.1/picard.jar SortVcf \
	I= /scratch/jandrews/Data/Variant_Calling/Coding/FINAL/merge_vs.vcf \
	O= /scratch/jandrews/Data/Variant_Calling/Coding/FINAL/merge_vs.sorted.vcf \
	SEQUENCE_DICTIONARY= /scratch/jandrews/Ref/ucsc.hg19.dict 
module remove java


#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N SORT_VCF_FOR_GATK
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb
#PBS -q old

module load java

for f in /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/MERGED/*.vcf; do
	java -Xmx8g -jar /scratch/jandrews/bin/picard-tools-2.2.1/picard.jar SortVcf \
	I="$f" \
	O="$f".sorted \
	SEQUENCE_DICTIONARY=/scratch/jandrews/Ref/hg19.dict 
	rename .vcf.sorted .sorted.vcf "$f".sorted
	rename .sorted.sorted.vcf .sorted.vcf *.sorted.sorted*
done

module remove java
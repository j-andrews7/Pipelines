#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N Var_Call_K_22
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=36gb

module load samtools-1.2
module load bcftools-1.2

for file in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/Batch22/*sorted.bam; do
	base=${file##*/}
	samtools mpileup -u -t DP -f /scratch/jandrews/Ref/hg19.fa "$file" | bcftools call -cv -O v | /scratch/jandrews/bin/vcfutils.pl varFilter -D100 > /scratch/jandrews/Data/Variant_Calling/Non_Coding/NEW_VCFs/${base%.*}.vcf 
	
done
wait
module remove samtools-1.2
module remove bcftools-1.2


#!/bin/sh
 

# give the job a name to help keep track of running jobs (optional)
#PBS -N Var_Call_BCF_NEW

#PBS -m e

#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=16gb

module load samtools-1.2
module load bcftools-1.2

for file in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch9/*.bam; do

	samtools mpileup -u -t DP -f /scratch/jandrews/Ref/hg19.fa $file | bcftools call -mv -O v - | /scratch/jandrews/bin/vcfutils.pl varFilter -D100 > /scratch/jandrews/Data/Variant_Calling/Coding/BCFTools/New/VCFs/$file.vcf &
	
done
wait
module remove samtools-1.2
module remove bcftools-1.2


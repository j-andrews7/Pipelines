#!/bin/sh
 

# give the job a name to help keep track of running jobs (optional)
#PBS -N Var_Call_VarScan

#PBS -m e

#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=24gb

module load samtools-1.2
module load bcftools-1.2

for file in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch29/*.bam; do
	base=${file##*/}
	samtools mpileup -t DP -f /scratch/jandrews/Ref/hg19.fa $file | java -Xmx15g -jar /scratch/jandrews/bin/VarScan.v2.3.9.jar mpileup2cns --min-coverage 5 --min-avg-qual 15 --variants 1 --output-vcf 1 > /scratch/jandrews/Data/Variant_Calling/Coding/VarScan/VCFs/"${base%.*}"_varscan.vcf  &
	
done
wait
module remove samtools-1.2
module remove bcftools-1.2


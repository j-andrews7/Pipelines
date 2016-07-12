#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N Check_mpileup

#PBS -m e

#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb

module load samtools-1.2
module load htslib
module load bcftools-1.2

samtools mpileup -v -f /scratch/jandrews/Ref/hg19.fa /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch5/DL191_accepted_hits.sorted.bam > mpileup_DL191.vcf.gz

module remove htslib
module remove samtools-1.2
module remove bcftools-1.2

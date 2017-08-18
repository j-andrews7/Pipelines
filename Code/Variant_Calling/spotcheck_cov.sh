#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N Spot_Check_Cov

#PBS -m e

#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb

module load samtools

samtools depth -q 15 -b /scratch/jandrews/Data/Variant_Calling/Coding/spot_check.bed /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch5/DL191_accepted_hits.sorted.bam > spotter_DL191.txt

module remove samtools


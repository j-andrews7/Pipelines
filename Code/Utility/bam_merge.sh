#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MERGE_BAMs_CC_CB
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=24gb

module load samtools

samtools merge /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/all_CC_CB_merged.bam -b /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/CC_CB_bam_list.txt

module remove samtools


#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PINDEL_200
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=168:00:00,vmem=64gb

module load pindel
module load samtools

pindel -T 8 -f /scratch/jandrews/Ref/hg19.fa -i /scratch/jandrews/Data/Indel_SV_Calling/Configs/pindel_config_rna_200bp_ins.txt -c ALL -o /scratch/jandrews/Data/Indel_SV_Calling/Results/RNA_200bp

module remove samtools
module remove pindel

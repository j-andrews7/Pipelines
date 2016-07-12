#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N Get_Fasta

#PBS -m e

#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb

module load bedtools2

for file in /scratch/jandrews/Data/Motif_Analyses/Genomic_Ranges_Bed/*_sigonly.bed; do

	base=${file##*/}
	bedtools getfasta -fi /scratch/jandrews/Ref/hg19.fa -bed "$file" -fo stdout | fold -w 50 > /scratch/jandrews/Data/Motif_Analyses/Results/${base%.*}.fa ;
	
done
wait

module remove bedtools2

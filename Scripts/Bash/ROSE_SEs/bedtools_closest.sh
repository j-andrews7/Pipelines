#!/bin/sh
 

# give the job a name to help keep track of running jobs (optional)
#PBS -N bedtools_closest

#PBS -m e

#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=8gb

module load bedtools2

for fold in /scratch/jandrews/Data/ChIP_Seq/MACS/ROSE_SEs_From_Ind_Sample_Peaks/RESULTS/*/; do

	cd "$fold" 
	bedtools closest -a ROSE_SuperEnhancers_sorted.bed -b /scratch/jandrews/Ref/gencode.v19.annotation_sorted_genes_only.bed -d > "$fold"${PWD##*/}_ROSE_SuperEnhancers_Sorted_Gencode_Annotated.bed &
	
done
wait
module remove bedtools2


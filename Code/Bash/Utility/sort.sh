#!/bin/sh
 

# give the job a name to help keep track of running jobs (optional)
#PBS -N sort

#PBS -m e

#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

for fold in /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks_Redo/RESULTS/*/; do

	sort -k 1,1 -k2,2n  "$fold"*SuperEnhancers.bed > "$fold"ROSE_SuperEnhancers_sorted.bed & 
	
done
wait


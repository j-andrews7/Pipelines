#!/bin/sh 

# give the job a name to help keep track of running jobs (optional)
#PBS -N move

#PBS -m e

#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

for fold in /scratch/jandrews/Data/ChIP_Seq/MACS/ROSE_SEs/*/; do

	cd "$fold"
	mv *_merged_sorted_annotations.bed /scratch/jandrews/Data/ChIP_Seq/MACS/ROSE_SEs/SE_Intersects/ & 
	
done
wait


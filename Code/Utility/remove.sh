#!/bin/sh 

# give the job a name to help keep track of running jobs (optional)
#PBS -N remove

#PBS -m e

#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

for fold in /scratch/jandrews/Data/ChIP_Seq/MACS/ROSE_SEs_From_Ind_Sample_Peaks/SE_INTERSECTS/*/; do

	cd "$fold"
	rm *SEs.bed & 
	
done
wait


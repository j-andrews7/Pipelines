#!/bin/sh
 

# give the job a name to help keep track of running jobs (optional)
#PBS -N bedops_unique

#PBS -m e

#PBS -l nodes=1:ppn=4,walltime=8:00:00,vmem=8gb

for fold in /scratch/jandrews/Data/ChIP_Seq/MACS/ROSE_SEs_From_Ind_Sample_Peaks/UNIQUE_*/; do

	cd "$fold"
	bedops --everything --ec *.bed \
    | bedmap --echo --count --fraction-both 0.25 - \
    | awk -F"|" '($2 == 1)'  > ./Results/All_Unique_SEs_FLDL_Subgroups.bed 
	
done
wait


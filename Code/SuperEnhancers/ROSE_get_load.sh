#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N GET_SE_LOAD_ALL
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=12gb

module load samtools
module load R

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

cd /scratch/jandrews/bin/rose/

for fold in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/B*/; do
	
	for file in "$fold"*.bam; do
		echo "$file"
		base=${file##*/}
		python ROSE_bamToGFF.py -b "$file" -m 1 -r -i /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks/UNIQUE_SES/BEDOPS_Overlap_Method/Final/Recurrent_Only/ALL_RECURRENT_SES_K27AC_LOAD/ALL_RECURRENT_SES_GFF/All_Recurrent_SEs_BEDOPS.gff -o /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks/UNIQUE_SES/BEDOPS_Overlap_Method/Final/Recurrent_Only/ALL_RECURRENT_SES_K27AC_LOAD/ALL_RECURRENT_SES_K27AC/${base%.*}_LOAD.gff ;
	done
done
wait
module remove samtools
module remove R

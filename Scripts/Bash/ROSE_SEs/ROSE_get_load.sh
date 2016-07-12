#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N GET_SE_LOAD_DL
#PBS -m e
#PBS -l nodes=1:ppn=4,walltime=8:00:00,vmem=12gb

module load samtools-1.2
module load R

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

cd /scratch/jandrews/bin/rose/

for fold in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/*/; do
	
	for file in "$fold"*.bam; do
		base=${file##*/}
		python ROSE_bamToGFF.py -b "$file" -m 1 -r -i /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks_Redo/UNIQUE_SES/Multiinter_Method/wOUT_VGA_VGR/UNIQUE_SES_GFF/DL_unique_SEs.gff -o /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks_Redo/UNIQUE_SES/Multiinter_Method/wOUT_VGA_VGR/UNIQUE_SES_K27AC_LOAD/UNIQUE_DL_SES_LOAD/${base%.*}_DL_SES_LOAD.gff ;
	done
done
wait
module remove samtools-1.2
module remove R

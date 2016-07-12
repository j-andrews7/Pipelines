#!/bin/sh
 

# give the job a name to help keep track of running jobs (optional)
#PBS -N ROSE_ind

#PBS -m e

#PBS -l nodes=1:ppn=4,walltime=8:00:00,vmem=12gb

module load samtools-1.2
module load R

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

cd /scratch/jandrews/bin/rose/

for file in /scratch/jandrews/Data/ChIP_Seq/K27AC/Batch22/*.bam; do
	
	echo "$file"
	base=${file##*/}
	python ROSE_main.py -g HG19 -t 2500 -r "$file" -i /scratch/jandrews/Data/ChIP_Seq/MACS/ROSE_SEs_From_Ind_Sample_Peaks/PEAKS_GFF/${base%.*}_peaks.gff -o /scratch/jandrews/Data/ChIP_Seq/MACS/ROSE_SEs_From_Ind_Sample_Peaks/RESULTS/${base%.*} &
	
done
wait
module remove samtools-1.2
module remove R


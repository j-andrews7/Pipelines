#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N IND_ROSE_B1
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=24gb

module load samtools-1.2
module load R
export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

cd /scratch/jandrews/bin/rose/

for file in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/NO_CTRL/Batch1/*sorted.bam; do
	
	echo "$file"
	base=${file##*/}
	python ROSE_main.py -g HG19 -t 2500 -r "$file" -i /scratch/jandrews/Data/ChIP_Seq/ROSE/STANDARD/PEAKS_GFF/${base/.sorted.bam/""}_peaks.gff -o /scratch/jandrews/Data/ChIP_Seq/ROSE/STANDARD/RESULTS/${base/.sorted.bam/""} 
	
done
wait
module remove samtools-1.2
module remove R


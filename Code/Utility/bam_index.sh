#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N INDEX_BAM

#PBS -m e

#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=12gb

module load samtools-1.2

for fold in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/Batch*/; do
	cd "$fold"
	for f in *BL_removed.bam; do
		samtools index "$f" 
	done
	wait
done 
	
module remove samtools-1.2


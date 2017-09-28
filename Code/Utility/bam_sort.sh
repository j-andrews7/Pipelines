#!/bin/sh
 

# give the job a name to help keep track of running jobs (optional)
#PBS -N BAM_SORT

#PBS -m e

#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb

module load samtools-1.2

for f in /scratch/jandrews/Data/ChIP_Seq/BAMs/SELECT/*.bam; do
		
	base=${f##*/}
	samtools sort -T $f.sorted -o /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/${base%.*}.sorted.bam "$f" 
	
done
wait
module remove samtools-1.2


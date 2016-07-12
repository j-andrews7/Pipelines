#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N RM_BLACKLIST_REGIONS

#PBS -m e

#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=32gb

module load samtools-1.2

for fold in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/Batch*/; do

	cd "$fold"
	for f in *.bam; do
		echo "$f"
		base=${f##/*}
		samtools view -b -t /home/jandrews/Ref/hg19.fa -L /home/jandrews/Ref/hg19_blacklist_regions_removed.bed -o "${base%.*}".BL_removed.bam "$f" &
	done
	wait
	
done
wait
module remove samtools-1.2


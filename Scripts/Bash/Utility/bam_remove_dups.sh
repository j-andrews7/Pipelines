#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N REMOVE_DUPS
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=32gb

module load java
module load R

for fold in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch*/; do

		cd "$fold"
		for f in *.sorted.bam; do
			echo "$f"
			base=${f%%.*}
			java -jar /export/picard-tools-2.0.1/picard.jar	MarkDuplicates INPUT="$f" OUTPUT="$base".dups_rmvd.bam REMOVE_DUPLICATES=true METRICS_FILE="$base".dup_metrics.txt ;
		wait
		done	
wait
done

module remove R
module remove java
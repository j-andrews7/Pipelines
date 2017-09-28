#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N bin_wigs
#PBS -m e
#PBS -q old
#PBS -l nodes=1:ppn=12,walltime=24:00:00,vmem=64gb

for fold in /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/WIGS/K27AC/TREAT/Batch*/; do

	cd "$fold"

	for f in "$fold"/*.wig; do
		perl /scratch/jandrews/bin/bin_whole_genome_wig.pl "$f" &
	done
	wait

	mv .bin /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/BIN_PROCESSING/K27AC/

done
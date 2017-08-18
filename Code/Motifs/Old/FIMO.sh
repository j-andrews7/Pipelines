#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N Motif_Hunt

#PBS -m e

#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

for file in /scratch/jandrews/Data/Motif_Analyses/Sequences/Batch2/*.fa; do

	base=${file##*/}
	fimo --o /scratch/jandrews/Data/Motif_Analyses/Results/${base%.*} /scratch/jandrews/bin/meme/db/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme "$file" &
	
done
wait


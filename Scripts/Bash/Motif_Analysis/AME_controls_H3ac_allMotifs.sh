#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N ame_motifs_w_controls_H3ac_allMotifs

#PBS -m e

#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

for file in /scratch/jandrews/Data/Motif_Analyses/Sequences/Sig_Only/H3ac*.fa; do

	base=${file##*/}
	ame --verbose 3 --pvalue-report-threshold 1 --o /scratch/jandrews/Data/Motif_Analyses/Results/${base%.*}_ALL --control /scratch/jandrews/Data/Motif_Analyses/Sequences/Ranked/Batch1/H3ac_pvalues_FLvDL_NonTSS_FAIRE_peaks_ranked.fa "$file" /scratch/jandrews/bin/meme/db/motif_databases/EUKARYOTE/jolma2013.meme &
	
done
wait


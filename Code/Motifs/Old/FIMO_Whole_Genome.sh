#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N Motif_Hunt_GW
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=168:00:00,vmem=64gb

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

fimo --o /scratch/jandrews/Data/Motif_Analyses/Results/HG19_FLvDL_Enriched /scratch/jandrews/Data/Motif_Analyses/Motif_Databases/enriched_NonTSS_FAIRE_sigdiff_H3AC_K27AC_FLvDL_2016.meme /scratch/jandrews/Ref/hg19.fa
	



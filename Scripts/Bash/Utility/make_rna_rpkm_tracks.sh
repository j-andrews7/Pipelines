#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MAKE_RNA_RPKM_TRACKS_1
#PBS -m e
#PBS -q old
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=36gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

batch=Batch1/
mark=_RNA
treat_suffix=.sorted.bam
treat=/scratch/jandrews/Data/RNA_Seq/BAMs/

for f in /scratch/jandrews/Data/RNA_Seq/BAMs/"$batch"/*"$treat_suffix"; do
	base=${f##*/}
	samp=${base%%_*}
	bamCoverage -p max -of bigwig --normalizeUsingRPKM  -bl /scratch/jandrews/Ref/ENCODE_Blacklist_hg19.bed -b "$f" -o "$treat""$batch""$samp""$mark".rpkm.bw ;
done
wait



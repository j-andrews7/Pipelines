#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MAKE_CHIP_RPKM_TRACKS_1
#PBS -m e
#PBS -q old
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=64gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

batch=Batch1/

# For K4ME3, set -e to 200, for FAIRE use 100, for other marks, use 300. This is just double the -shiftsize used for macs
# and is supposed to be the fragment length.

cd /scratch/jandrews/Data/ChIP_Seq/T_Cell/BAMs/"$batch"

for f in *.BL_removed.bam; do
	bamCoverage  -e 300 -p 8 -of bigwig -bs 10 --normalizeUsingRPKM  -bl /scratch/jandrews/Ref/ENCODE_Blacklist_hg19.bed -b "$f" -o "$f".bw ;
done
wait

rename .bam.bw .RPKM.bw *.bw 

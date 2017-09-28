#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MAKE_CHIP_RPKM_TRACKS_INPUT_SUBT_1
#PBS -m e
#PBS -q old
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=64gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

batch=Batch1/
mark=_K27AC
treat_suffix=.sorted.BL_removed.bam
treat=/scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/
control=/scratch/jandrews/Data/ChIP_Seq/BAMs/INPUT/

# For K4ME3, set -e to 200, for FAIRE use 100, for other marks, use 300. This is just double the -shiftsize used for macs
# and is supposed to be the fragment length.

for f in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/"$batch"/*"$treat_suffix"; do
    base=${f##*/}
    samp=${base%%_*}
    bamCompare -ratio subtract -e 300 -p max -of bigwig --normalizeUsingRPKM  -bl /scratch/jandrews/Ref/ENCODE_Blacklist_hg19.bed \
    -b1 "$f" -b2 "$control""$batch""$samp"_INPUT.sorted.bam -o "$treat""$batch""$samp""$mark".BL_removed.input_subt.rpkm.bw ;
done
wait



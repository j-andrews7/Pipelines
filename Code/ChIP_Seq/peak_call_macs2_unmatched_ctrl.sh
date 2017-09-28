#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PEAK_CALL_BL_RMVD_14
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate py2

# Idea here is to throw all the samples and the single input file for each of them in a given folder.
# So that this can be submitted in batches.

treat=/scratch/jandrews/Data/ChIP_Seq/T_Cell/BAMs/Batch10/

for f in "$treat"*C.sorted.BL_removed.bam; do
    base=${f##*/}
    macs2 callpeak -t "$f" -c "$treat"*INPUT.sorted.BL_removed.bam -n "$base" --outdir "$treat" --tempdir /scratch/jandrews/Scratch/MACS2 -q 0.01 -B --SPMR -m 10 50 &
done
wait

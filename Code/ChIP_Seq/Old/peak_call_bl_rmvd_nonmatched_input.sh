#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PEAK_CALL_BL_RMVD_UNMATCHED_1
#PBS -m e
#PBS -q old
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=36gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate py2

batch=Batch1/
mark=_K27AC
treat_suffix=.sorted.BL_removed.bam
treat=/scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/

for f in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/NO_CTRL/"$batch"/*"$treat_suffix"; do
    base=${f##*/}
    base=${base%.*}
    samp=${base%_*}
    macs14 -t "$f" -c /scratch/jandrews/Data/ChIP_Seq/BAMs/INPUT/CLL378_INPUT.sorted.bam -f BAM -g hs -n /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/"$samp""$mark"_Unmatched_Ctrl -w -S --nomodel --shiftsize=150 &
done
wait

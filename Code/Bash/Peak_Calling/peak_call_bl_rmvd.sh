#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PEAK_CALL_BL_RMVD_22
#PBS -m e
#PBS -q old
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=36gb

batch=Batch22/
mark=_K27AC
treat_suffix=.sorted.BL_removed.bam
control_suffix=_INPUT.sorted.bam
treat=/scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/
control=/scratch/jandrews/Data/ChIP_Seq/BAMs/INPUT/

for f in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/"$batch"/*"$treat_suffix"; do
	base=${f##*/}
	samp=${base%%_*}
	macs14 -t "$f" -c "$control""$batch""$samp""$control_suffix" -f BAM -g hs -n /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/"$samp""$mark" -w -S --nomodel --shiftsize=150 &
done
wait



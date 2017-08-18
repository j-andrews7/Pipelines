#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PEAK_CALL_NOCTRL_STANDARD
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb

batch=Batch4/
mark=_K27AC
treat_suffix=.sorted.bam
control_suffix=_INPUT.sorted.bam
treat=/scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/
control=/scratch/jandrews/Data/ChIP_Seq/BAMs/INPUT/

for f in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/NO_CTRL/"$batch"/*"$treat_suffix"; do
	base=${f##*/}
	base=${base%.*}
	samp=${base%_*}
	macs14 -t "$f" -f BAM -g hs -n /scratch/jandrews/Data/ChIP_Seq/MACS/STANDARD/"$samp""$mark" -w -S --nomodel --shiftsize=150 &
done
wait



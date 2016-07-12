#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PEAK_CALL_NO_CTRL
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

base=TS102214_NAIVE
base_add=_K27AC
treat_batch=Batch10/
control_batch=Batch14/
treat_suffix=.sorted.bam
control_suffix=_INPUT.sorted.bam
treat=/scratch/jandrews/Data/ChIP_Seq/BAMs/SELECT/
control=/scratch/jandrews/Data/ChIP_Seq/BAMs/INPUT/

macs14 -t "$treat""$base""$base_add""$treat_suffix" -f BAM -g hs -n /scratch/jandrews/Data/ChIP_Seq/MACS/SELECT/"$base""$base_add" -w -S --nomodel --shiftsize=150




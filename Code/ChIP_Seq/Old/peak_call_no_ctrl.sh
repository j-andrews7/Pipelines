#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PEAK_CALL_NO_CTRL
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

base=VGR11
base_add=_X
treat_batch=Batch18/
control_batch=Batch14/
treat_suffix=.sorted.bam
control_suffix=_INPUT.sorted.bam
treat=/scratch/jandrews/Data/ChIP_Seq/BAMs/FAIRE/
control=/scratch/jandrews/Data/ChIP_Seq/BAMs/INPUT/

macs14 -t "$treat""$treat_batch""$base""$base_add""$treat_suffix" -f BAM -g hs -n /scratch/jandrews/Data/ChIP_Seq/MACS/STANDARD/"$base""$base_add" -w -S --nomodel --shiftsize=50




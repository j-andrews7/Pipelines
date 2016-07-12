#!/bin/sh
 

# give the job a name to help keep track of running jobs (optional)
#PBS -N tophat_align_TS072111A

#PBS -m e

#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=12gb

module load tophat-1.4.1.1
export PATH=/scratch/jandrews/bin/samtools-0.1.19:$PATH

cd /scratch/jandrews/Data/RNA_Seq/Alignment_Workspace/FastQs

tophat -p 8 -o /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/TS072111A -G ../gencode.v19.annotation_sorted.txt hg19 TS072111A_1.fq.gz TS072111A_2.fq.gz

module remove tophat-1.4.1.1


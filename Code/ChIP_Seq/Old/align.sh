#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N align_chip
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=64gb

module load bowtie2
module load samtools-1.2

for f in /scratch/jandrews/Data/ChIP_Seq/T_Cell/*.fq.gz; do
  bowtie2 -p 8 -x /scratch/jandrews/Ref/hg19 -U "$f" | samtools view -bS - | samtools sort - > ${f%.*}.sorted.bam &
  samtools index ${f%.*}.sorted.bam
done

module remove bowtie2
module remove samtools-1.2

#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N align_chip
#PBS -m be
#PBS -l nodes=1:ppn=8,walltime=72:00:00,vmem=128gb

module load bowtie2
module load samtools

# Directory containing reference genome fasta, bowtie indices, blacklist, etc.
genome=/scratch/jandrews/Ref/hg19
# Base directory containing reads and where additional directories will be created.
base_dir=/scratch/jandrews/Data/ChIP_Seq/T_Cell/READS

cd "$base_dir"
mkdir BAMS

# Alignment.
for f in *.fq.gz; do
  echo Aligning "$f"
  bowtie2 -p 6 -x "$genome" -U "$f" | samtools view -@ 8 -b -S -u - | samtools sort -@ 8 -m 5G -o ./BAMS/"${f%.*}".sorted.bam -;
  samtools index /scratch/jandrews/Data/ChIP_Seq/T_Cell/READS/BAMS/"${f%.*}".sorted.bam
done

# Removal of blacklist reads.
cd BAMS
for f in *.bam; do
    echo Removing blacklisted reads from "$f"
    base=${f##/*}
    samtools view -b -t /home/jandrews/Ref/hg19.fa -L /home/jandrews/Ref/hg19_blacklist_regions_removed.bed -o "${base%.*}".BL_removed.bam "$f";
done

module remove bowtie2
module remove samtools

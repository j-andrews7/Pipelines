#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N VCF_TO_BED
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

for f in /scratch/jandrews/Data/Variant_Calling/Non_Coding/*.vcf; do
    base=${f##*/}
    vcf2bed < "$f" > ${base%.*}.bed ;
done
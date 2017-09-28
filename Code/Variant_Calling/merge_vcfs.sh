#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MERGE_VCFs

#PBS -m e

#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=36gb


module load bcftools-1.2

bcftools merge -O v -m none -i DP:sum,QV:avg,DP4:sum /scratch/jandrews/Data/Variant_Calling/Non_Coding/VCFs_With_Quals/*.gz > /scratch/jandrews/Data/Variant_Calling/Non_Coding/merged_noncoding_variants.vcf

module remove bcftools-1.2

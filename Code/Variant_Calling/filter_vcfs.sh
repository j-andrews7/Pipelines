#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N FILTER_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=36gb


module load bcftools-1.2

for f in /scratch/jandrews/Data/Variant_Calling/Coding/BCFTools/New/VCFs/*.vcf; do
    base=${f##*/}
    bcftools filter -i 'DP>=5 & QUAL>=15' "$f" > /scratch/jandrews/Data/Variant_Calling/Coding/BCFTools/New/VCFs/${base%.*}.filtered.vcf ;
done

module remove bcftools-1.2

#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MERGE_ALL_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=36gb


module load bcftools

bcftools merge -O v -m none --force-samples -i set:join /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/MERGED/*.gz > /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/MERGED/ALL_VARIANTS.MERGED.RNA_DP10.RNA_NODUPS.CHIP_MULTIMARK.vcf

module remove bcftools

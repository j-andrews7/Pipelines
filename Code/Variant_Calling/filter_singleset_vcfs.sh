#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N FILTER_SINGLESET_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb

# This script filters variants found by only one of the two callers used. A variant has to be found in two datasets (already done for chip_seq),
# so the set would have to be "bcftools-varscan", "chip_seq", "varscan-chip_seq", or "bcftools-chip_seq" to be included in the output.

module load bcftools

for f in /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/*Combined.vcf; do
    bcftools filter -e 'INFO/set="bcftools" || INFO/set="varscan"' "$f" > "$f".filtered ;
    rename .vcf.filtered .multiset.vcf "$f".filtered
done

module remove bcftools
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PROCESS_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

module load bcftools-1.2

for fold in /scratch/jandrews/Data/Variant_Calling/Non_Coding/VCFs_With_Quals/Ind_Samp_Merged/*/; do
	cd "$fold"
	bgzip -c ${PWD##*/}_variants.vcf > ${PWD##*/}_variants.vcf.gz
	tabix -p vcf ${PWD##*/}_variants.vcf.gz
    bcftools filter -i 'DP>=10' ${PWD##*/}_variants.vcf.gz > ${PWD##*/}_variants_filtered.vcf
done

module remove bcftools-1.2
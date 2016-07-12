#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MERGE_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb

module load bcftools-1.2
for fold in /scratch/jandrews/Data/Variant_Calling/Non_Coding/VCFs_With_Quals/Ind_Samp_Merged/*/; do
	cd "$fold"
	base=${PWD##*/}
	cd ${PWD##*/}_Ind_VCFs/
	bcftools merge -O v -m none -i DP:sum,QV:avg,DP4:sum *.gz > ../"$base"_variants.vcf
done

module remove bcftools-1.2
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N ISEC_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb

module load bcftools-1.2
module load bedtools2
for fold in /scratch/jandrews/Data/Variant_Calling/Non_Coding/VCFs_With_Quals/Ind_Samp_Merged/*/; do
	cd "$fold"
	bcftools isec -p ${PWD}/Full_isec -c none ${PWD##*/}_variants.vcf.gz /scratch/jandrews/Data/Variant_Calling/Non_Coding/merged_noncoding_multitype_variants.vcf.gz
	bedtools intersect -wa -wb -a ${PWD##*/}_variants.vcf -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/merged_noncoding_multitype_funseq_MAF0,01_sorted_positions.bed > ${PWD##*/}_funseq_isec.txt
done
module remove bcftools-1.2
module remove bedtools2

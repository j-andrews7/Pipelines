#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N ISEC_EVERYTHING
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb

module load bedtools2
for fold in /scratch/jandrews/Data/Variant_Calling/Non_Coding/VCFs_With_Quals/Ind_Samp_Merged/*/; do
	cd "$fold"
	rm ${PWD##*/}_isecs_summary.txt
	touch ${PWD##*/}_isecs_summary.txt
	echo -en ${PWD##*/}' filtered, multitype variants: ' >> ${PWD##*/}_isecs_summary.txt
	total_fil_vars=$(wc -l < ${PWD}/Full_isec/0002.vcf)
	echo -e $total_fil_vars >> ${PWD##*/}_isecs_summary.txt

	bedtools intersect -wa -a ${PWD}/Full_isec/0002.vcf -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/MMPIDS/MMPID_NonTSS_FAIRE_POSITIVE_outsideSEs.bed | uniq - > ${PWD##*/}_variants_inMMPIDS_outsideSEs.vcf

	# Too lazy to change the variable names. Deal with it.
	echo -en ${PWD##*/}' filtered, multitype variants in regular enhancers: ' >> ${PWD##*/}_isecs_summary.txt
	total_fil_vars=$(wc -l < ${PWD##*/}_variants_inMMPIDS_outsideSEs.vcf)
	echo -e $total_fil_vars >> ${PWD##*/}_isecs_summary.txt

	bedtools intersect -wa -a ${PWD}/Full_isec/0002.vcf -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/SES/ALL_SES_POSITIONS_SORTED.bed | uniq - > ${PWD##*/}_variants_inSEs.vcf

	echo -en ${PWD##*/}' filtered, multitype variants in super enhancers: ' >> ${PWD##*/}_isecs_summary.txt
	total_fil_vars=$(wc -l < ${PWD##*/}_variants_inSEs.vcf)
	echo -e $total_fil_vars >> ${PWD##*/}_isecs_summary.txt

	bedtools intersect -wa -a ${PWD}/Full_isec/0002.vcf -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/TSS/2kbTSStranscripts_gencode19_protein_coding_positions_uniq.bed | uniq - > ${PWD##*/}_variants_in2kbTSS.vcf

	echo -en ${PWD##*/}' filtered, multitype variants in 2kb TSSs: ' >> ${PWD##*/}_isecs_summary.txt
	total_fil_vars=$(wc -l < ${PWD##*/}_variants_in2kbTSS.vcf)
	echo -e $total_fil_vars >> ${PWD##*/}_isecs_summary.txt

	bedtools intersect -wa -a ${PWD}/Full_isec/0002.vcf -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/TSS/2kbTSStranscripts_gencode19_protein_coding_inSEs.bed | uniq - > ${PWD##*/}_variants_in2kbTSS_inSEs.vcf

	echo -en ${PWD##*/}' filtered, multitype variants in 2kb TSSs in SEs: ' >> ${PWD##*/}_isecs_summary.txt
	total_fil_vars=$(wc -l < ${PWD##*/}_variants_in2kbTSS_inSEs.vcf)
	echo -e $total_fil_vars >> ${PWD##*/}_isecs_summary.txt

	echo -e "" >> ${PWD##*/}_isecs_summary.txt

	rev ${PWD##*/}_funseq_isec.txt | cut -f1-3 | rev > ${PWD##*/}_funseq_isec_positions.bed
	sort -k1,1 -k2,2n ${PWD##*/}_funseq_isec_positions.bed | uniq - > ${PWD##*/}_funseq_isec_positions_sorted.bed

	echo -en ${PWD##*/}' funseq annotated variants: ' >> ${PWD##*/}_isecs_summary.txt
	total_fil_vars=$(wc -l < ${PWD##*/}_funseq_isec_positions_sorted.bed)
	echo -e $total_fil_vars >> ${PWD##*/}_isecs_summary.txt

	bedtools intersect -wa -a ${PWD##*/}_funseq_isec_positions_sorted.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/MMPIDS/MMPID_NonTSS_FAIRE_POSITIVE_outsideSEs.bed | uniq - > ${PWD##*/}_funseq_variants_inMMPIDS_outsideSEs.bed

	echo -en ${PWD##*/}' funseq annotated variants in regular enhancers: ' >> ${PWD##*/}_isecs_summary.txt
	total_fil_vars=$(wc -l < ${PWD##*/}_funseq_variants_inMMPIDS_outsideSEs.bed)
	echo -e $total_fil_vars >> ${PWD##*/}_isecs_summary.txt

	bedtools intersect -wa -a ${PWD##*/}_funseq_isec_positions_sorted.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/SES/ALL_SES_POSITIONS_SORTED.bed | uniq - > ${PWD##*/}_funseq_variants_inSEs.bed

	echo -en ${PWD##*/}' funseq annotated variants in super enhancers: ' >> ${PWD##*/}_isecs_summary.txt
	total_fil_vars=$(wc -l < ${PWD##*/}_funseq_variants_inSEs.bed)
	echo -e $total_fil_vars >> ${PWD##*/}_isecs_summary.txt

	bedtools intersect -wa -a ${PWD##*/}_funseq_isec_positions_sorted.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/TSS/2kbTSStranscripts_gencode19_protein_coding_positions_uniq.bed | uniq - > ${PWD##*/}_funseq_variants_in2kbTSS.bed

	echo -en ${PWD##*/}' funseq annotated variants in 2kb TSSs: ' >> ${PWD##*/}_isecs_summary.txt
	total_fil_vars=$(wc -l < ${PWD##*/}_funseq_variants_in2kbTSS.bed)
	echo -e $total_fil_vars >> ${PWD##*/}_isecs_summary.txt

	bedtools intersect -wa -a ${PWD##*/}_funseq_isec_positions_sorted.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/TSS/2kbTSStranscripts_gencode19_protein_coding_inSEs.bed | uniq - > ${PWD##*/}_funseq_variants_in2kbTSS_inSEs.bed

	echo -en ${PWD##*/}' funseq annotated variants in 2kb TSSs in SEs: ' >> ${PWD##*/}_isecs_summary.txt
	total_fil_vars=$(wc -l < ${PWD##*/}_funseq_variants_in2kbTSS_inSEs.bed)
	echo -e $total_fil_vars >> ${PWD##*/}_isecs_summary.txt

done
module remove bedtools2

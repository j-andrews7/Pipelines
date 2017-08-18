#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N ISEC_EVERYTHING
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb

module load bedtools2

rm /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
touch /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

for f in /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/*positions.bed; do

	base=${f##*/}

	sort -k1,1 -k2,2n "$f" | uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted.bed
	echo -en ${base%.*}' funseq annotated variants: ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < $f)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

	bedtools intersect -v -a /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted.bed -b /scratch/jandrews/Ref/Ig_Loci.bed | uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed

	echo -en ${base%.*}' funseq annotated variants (Ig filtered): ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

	bedtools intersect -wa -a /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/MMPIDS/MMPID_NonTSS_FAIRE_POSITIVE_outsideSEs.bed | uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_inMMPIDS_outsideSEs_noIgLoci.bed

	echo -en ${base%.*}' funseq annotated variants in regular enhancers (Ig filtered): ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_inMMPIDS_outsideSEs_noIgLoci.bed)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

	bedtools intersect -wa -a /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/SES/ALL_SES_POSITIONS_SORTED.bed | uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_inSEs_noIgLoci.bed

	echo -en ${base%.*}' funseq annotated variants in super enhancers (Ig filtered): ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_inSEs_noIgLoci.bed)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

	bedtools intersect -wa -a /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/TSS/2kbTSStranscripts_gencode19_protein_coding_positions_uniq.bed | uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_in2kbTSS_noIgLoci.bed

	echo -en ${base%.*}' funseq annotated variants in 2kb TSSs (Ig filtered): ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_in2kbTSS_noIgLoci.bed)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

	bedtools intersect -wa -a /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/TSS/2kbTSStranscripts_gencode19_protein_coding_inSEs.bed | uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_in2kbTSS_inSEs_noIgLoci.bed

	echo -en ${base%.*}' funseq annotated variants in 2kb TSSs in SEs (Ig filtered): ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_in2kbTSS_inSEs_noIgLoci.bed)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

done
module remove bedtools2
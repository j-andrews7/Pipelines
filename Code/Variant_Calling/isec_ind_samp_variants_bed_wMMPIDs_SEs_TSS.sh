#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N ISEC_EVERYTHING
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb

module load bedtools2

rm /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt
touch /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt

for f in /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/FUNSEQ_ANNOT/FLDL_ONLY/*noncodingOnly.positions.uniq.bed; do

	base=${f##*/}

	sort -k1,1 -k2,2n "$f" > /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.sorted.bed

	echo -en ${base%%.*}'_SNVs: ' >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt
	var_count=$(sed '/^\s*#/d' $f | wc -l)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt

	bedtools intersect -header -v -a /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.sorted.bed \
	-b /scratch/jandrews/Ref/Ig_Loci.bed > /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.sorted.noIgLoci.bed

	echo -en ${base%%.*}'_SNVs_Ig_Rmvd: ' >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt
	var_count=$(sed '/^\s*#/d' /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.sorted.noIgLoci.bed | wc -l)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt

	bedtools intersect -header -wa -a /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.sorted.noIgLoci.bed \
	-b /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/MMPIDS/NonTSS_FAIREpositive_MMPIDs.uniq.outsideSEs.bed \
	| uniq - > /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.inMMPIDS.outsideSEs.sorted.noIgLoci.bed

	echo -en ${base%%.*}'_SNVs_inCEs_Ig_Rmvd: ' >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt
	var_count=$(sed '/^\s*#/d' /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.inMMPIDS.outsideSEs.sorted.noIgLoci.bed | wc -l)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt

	bedtools intersect -header -wa -a /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.sorted.noIgLoci.bed \
	-b /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/SES/All_SEs.bed \
	| uniq - > /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.inSEs.sorted.noIgLoci.bed

	echo -en ${base%%.*}'_SNVs_inSEs_Ig_Rmvd: ' >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt
	var_count=$(sed '/^\s*#/d' /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.inSEs.sorted.noIgLoci.bed | wc -l)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt

	bedtools intersect -header -wa -a /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.sorted.noIgLoci.bed \
	-b /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/TSS/gencode.v19.annotation_sorted.linc_proteincoding.2kbTSS.uniq.outsideSEs.bed \
	| uniq - > /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.in2kbTSS.sorted.noIgLoci.bed

	echo -en ${base%%.*}'_SNVs_in2kbTSS_Ig_Rmvd: ' >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt
	var_count=$(sed '/^\s*#/d' /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.in2kbTSS.sorted.noIgLoci.bed | wc -l)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt

	bedtools intersect -header -wa -a /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.sorted.noIgLoci.bed \
	-b /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/TSS/gencode.v19.annotation_sorted.linc_proteincoding.2kbTSS.uniq.inSEs.bed \
	| uniq - > /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.in2kbTSS.inSEs.sorted.noIgLoci.bed

	echo -en ${base%%.*}'_SNVs_in2kbTSS_inSEs_Ig_Rmvd: ' >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt
	var_count=$(sed '/^\s*#/d' /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FUNSEQ_ANNO_IND_SAMPS/FLDL_SEP_ANNO/${base%.*}.in2kbTSS.inSEs.sorted.noIgLoci.bed | wc -l)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/INTERSECTS/FLDL_funseq_isecs_summary.txt

done
module remove bedtools2
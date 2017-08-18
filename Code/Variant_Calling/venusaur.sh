#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N VENUSAUR
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=36gb
#PBS -q old

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

base_dir=/scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/VENUSAUR/

for file in /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/VENUSAUR/COMMON_SNPS_RMVD/Batch1/*.vcf; do
	base=${file##*/}
	python /scratch/jandrews/bin/tfsites_checker.py -i "$file" -o "$base_dir"RESULTS/COMMON_SNPS_RMVD/${base%.*}.motifs.vcf -r /scratch/jandrews/Ref/hg19.fa -m "$base_dir"J2016CORE_fpr_01.txt -bp "$base_dir"bp_all.txt -ci "$base_dir"GM12878_TF151_names_final.bed -fm -fc &
	
done
wait
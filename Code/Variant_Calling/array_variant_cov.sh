#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N SNP_ARRAY_COV
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=8:00:00,vmem=16gb

module load bedtools2

for fold in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch*/; do

	cd "$fold"
	for f in *.bam; do
		echo "$f"
		base=${f%%_*}
		bedtools multicov -bams "$f" -bed /scratch/jandrews/Data/Variant_Calling/Coding/Array_Comparisons/SNP_Array_Ind_Samples/Het_Only/"${base}"_array.vcf  > /scratch/jandrews/Data/Variant_Calling/Coding/Array_Comparisons/Cov_At_Each_Array_Variant_In_RNASeq/Het_Only/"${base}"_array_coverage.vcf &
	done

done
wait
module remove bedtools2


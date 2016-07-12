#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N merge_SE_annotations

#PBS -m e

#PBS -l nodes=1:ppn=4,walltime=4:00:00,vmem=4gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for fold in /scratch/jandrews/Data/ChIP_Seq/MACS/ROSE_SEs_From_Ind_Sample_Peaks/RESULTS/*/; do
	
	cd "$fold"
	mv *_SuperEnhancers_ENHANCER_TO_GENE.txt ${PWD##*/}_ROSE_SuperEnhancers_Sorted_RefSeq_Annotated.txt 
	python3 /scratch/jandrews/bin/merge_SE_annotations.py *Gencode* *RefSeq* ${PWD##*/}_ROSE_SEs_merged_annotations.bed ;
	sort -k1,1 -k2,2n ${PWD##*/}_ROSE_SEs_merged_annotations.bed > ${PWD##*/}_ROSE_SEs_merged_sorted_annotations.bed ;
	rm ${PWD##*/}_ROSE_SEs_merged_annotations.bed &
	
done
wait


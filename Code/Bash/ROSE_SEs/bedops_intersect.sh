#!/bin/sh
 

# give the job a name to help keep track of running jobs (optional)
#PBS -N bedops_intersect

#PBS -m e

#PBS -l nodes=1:ppn=4,walltime=8:00:00,vmem=8gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for fold in /scratch/jandrews/Data/ChIP_Seq/MACS/ROSE_SEs_From_Ind_Sample_Peaks/SE_INTERSECTS/*/; do

	cd "$fold"
	bedops --everything *.bed \
    | bedmap --echo-map --ec --sweep-all --fraction-either 0.25 - \
    | sed 's/\:/\n/g' - \
    | sort-bed - \
    | uniq - \
    | python /scratch/jandrews/bin/parse_BEDOPS.py \
    | sort -k1,1 -k2,2n > ${PWD##*/}_SEs.bed 
	
done
wait


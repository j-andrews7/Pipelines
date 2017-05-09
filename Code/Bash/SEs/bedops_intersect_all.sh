#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N BOPS_ISEC_ALL
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for fold in /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks_Redo/SE_INTERSECTS/BEDOPS_Overlap_Method/FLDL_TS_CCCB_For_All/; do

	cd "$fold"
	bedops --everything *.bed \
    | bedmap --echo-map --ec --sweep-all --fraction-either 0.25 - \
    | sed 's/\:/\n/g' - \
    | sort-bed - \
    | uniq - \
    | python /scratch/jandrews/bin/parse_BEDOPS.py  > ${PWD##*/}_SEs.bed 
	
done
wait


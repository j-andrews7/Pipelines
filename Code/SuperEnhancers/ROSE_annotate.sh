#!/bin/sh
 

# give the job a name to help keep track of running jobs (optional)
#PBS -N ROSE_annotate

#PBS -m e

#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb

module load samtools-1.2
module load R

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

cd /scratch/jandrews/bin/rose/

for fold in /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Merged_BAMs/RESULTS/*/; do
	
	python ROSE_geneMapper.py -f -g HG19 -i "$fold"/*SuperEnhancers.table.txt -o "$fold" ;
	
done
wait
module remove samtools-1.2
module remove R


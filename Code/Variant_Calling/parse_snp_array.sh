#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PARSE_SNP_ARRAYS
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python3 /scratch/jandrews/bin/parse_snp_array.py /scratch/jandrews/Data/Variant_Calling/SNP_Arrays/GenomeWideSNP_6.na35.annot.csv /scratch/jandrews/Data/Variant_Calling/SNP_Arrays/LYMPHOMASNP77_GTYPE_2014.txt
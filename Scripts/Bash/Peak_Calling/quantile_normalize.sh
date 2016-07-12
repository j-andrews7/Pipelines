#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N quantile_normalize
#PBS -m e
#PBS -l nodes=1:ppn=1:haswell,walltime=24:00:00,vmem=96gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

cd /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/BIN_PROCESSING/K27AC

python /scratch/jandrews/bin/quantile_normalize_new.py master_table_RPKM.bed 5
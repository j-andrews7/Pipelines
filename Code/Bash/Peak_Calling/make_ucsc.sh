#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N MAKE_UCSC
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=64gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

cd /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/BIN_PROCESSING/K27AC

python /scratch/jandrews/bin/MakeUCSC.py -i QN_master_table_RPKM.bed -o QN_master_table_RPKM_final.bed
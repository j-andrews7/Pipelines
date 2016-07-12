#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N CALC_RPM
#PBS -m e
#PBS -q old
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb

cd /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/BIN_PROCESSING/K27AC

perl /scratch/jandrews/bin/Calc_RPM_ChIP_Seq.pl chrcoords_master_table.bed master_table_RPKM.bed
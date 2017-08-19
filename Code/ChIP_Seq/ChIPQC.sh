#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N ChIPQC
#PBS -m be
#PBS -l nodes=1:ppn=8,walltime=96:00:00,vmem=96gb

module load R

Rscript /scratch/jandrews/bin/ChIPQC_Cell_Lines.R

module remove R
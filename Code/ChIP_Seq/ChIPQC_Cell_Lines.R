library(ChIPQC)
setwd("/scratch/jandrews/Data/ChIP_Seq/T_Cell/ChIPQC")
samples = read.csv("CellLine_Experiment_Sheet.csv")

experiment = ChIPQC(samples)
ChIPQCreport(experiment)
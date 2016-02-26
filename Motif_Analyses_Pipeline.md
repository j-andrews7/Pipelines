Up to date as of 01/22/2015
An attempt to find motifs enriched in ChIP-Seq peaks (FAIRE in this case) located within REs that differ significantly for various ChIP-Seq between different cell types.
In this case, it was looking at FAIRE peaks located within NonTSS REs with significantly different K27AC or H3AC load between FL and DL samples.

All necessary scripts should be in Jared's code folder: N:\Bioinformatics\Jareds_Code

1.) Get correct columns from MMPID p-value comparisons.
awk -v OFS='\t' '{print $7, $8, $9, $6}' K27ac_pvalues_DLvCC_chrom_NonTSS_FAIREpositive_FAIREdata.txt > K27ac_pvalues_DLvCC_NonTSS_FAIRE_peaks.bed


2A.) Sort by pvalue, should be in fourth column (low to high).
sort -k4 -g K27ac_pvalues_FLvDL_NonTSS_FAIRE_peaks.bed > K27ac_pvalues_FLvDL_NonTSS_FAIRE_peaks_ranked.bed 


2B.) Shuffle randomly as well. Used as comparison to determine how ranking actually affects the program.
shuf K27ac_pvalues_FLvDL_NonTSS_FAIRE_peaks.bed > K27ac_pvalues_FLvDL_NonTSS_FAIRE_peaks_shuffled.bed


3.) Create file with only those that are significant from the ranked files (H3AC and K27AC - FLvDL).
cat K27ac_pvalues_FLvDL_NonTSS_FAIRE_peaks_ranked.bed | awk '$4 <= 0.05' > K27ac_pvalues_FLvDL_NonTSS_FAIRE_peaks_sigonly.bed
cat H3ac_pvalues_FLvDL_NonTSS_FAIRE_peaks_ranked.bed | awk '$4 <= 0.05' > H3ac_pvalues_FLvDL_NonTSS_FAIRE_peaks_sigonly.bed


4.) Get fasta sequences for each FAIRE peak.
Bash script (getfasta.sh):
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N Get_Fasta
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb

module load bedtools2

for file in /scratch/jandrews/Data/Motif_Analyses/Genomic_Ranges_Bed/*.bed; do

        base=${file##*/}
        bedtools getfasta -fi /scratch/jandrews/Ref/hg19.fa -bed "$file" -fo stdout | fold -w 50 > /scratch/jandrews/Data/Motif_Analyses/Results/${base%.*}.fa &

done
wait

module remove bedtools2


5.) Run AME (Analysis of Motif Enrichment - part of MEME Suite) using the fasta files as the sequence databases. Run on the cluster if you can, takes ~48 hours a file locally.
For tool docs: http://meme-suite.org/doc/ame.html?man_type=web

For ranked vs shuffled (NOT RECOMMENDED):
Bash script (AME.sh):
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N ame_motifs
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=168:00:00,vmem=48gb

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

for file in /scratch/jandrews/Data/Motif_Analyses/Sequences/Ranked/Batch1/*.fa; do

        base=${file##*/}
        ame --o /scratch/jandrews/Data/Motif_Analyses/Results/${base%.*} "$file" /scratch/jandrews/bin/meme/db/motif_databases/EUKARYOTE/jolma2013.meme &

done
wait

For significant vs rest (RECOMMENDED), as it's much quicker and gives more meaningful results. This will report p-vals for all motifs, not just those that are significant.
Can remove '--pvalue-report-threshold 1' to only report those that are significant.
Bash script (AME_controls_K27ac_allMotifs.sh):
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N ame_motifs_w_controls_K27ac_allMotifs
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

for file in /scratch/jandrews/Data/Motif_Analyses/Sequences/Sig_Only/K27ac*.fa; do

	base=${file##*/}
	ame --verbose 3 --pvalue-report-threshold 1 --o /scratch/jandrews/Data/Motif_Analyses/Results/${base%.*} --control /scratch/jandrews/Data/Motif_Analyses/Sequences/Ranked/Batch4/K27ac_pvalues_FLvDL_NonTSS_FAIRE_peaks_ranked.fa "$file" /scratch/jandrews/bin/meme/db/motif_databases/EUKARYOTE/jolma2013.meme &
	
done
wait


Bash script (AME_controls_H3ac_allMotifs.sh):
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N ame_motifs_w_controls_H3ac_allMotifs
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

for file in /scratch/jandrews/Data/Motif_Analyses/Sequences/Sig_Only/H3ac*.fa; do

	base=${file##*/}
	ame --verbose 3 --pvalue-report-threshold 1 --o /scratch/jandrews/Data/Motif_Analyses/Results/${base%.*} --control /scratch/jandrews/Data/Motif_Analyses/Sequences/Ranked/Batch1/H3ac_pvalues_FLvDL_NonTSS_FAIRE_peaks_ranked.fa "$file" /scratch/jandrews/bin/meme/db/motif_databases/EUKARYOTE/jolma2013.meme &
	
done
wait


NOTE: Initial runs for the shuffled sequences gave no results. May have to reduce number of sequences used (~120k currently).
NOTE: Ranked sequences gave results. Trying only sig sequences vs all (shuffled) to see if it gives better results.
NOTE: Running with control sequences (sig vs all) drastically improves runtime. Highly suggested.


6.) Get enriched motifs from output files, go into meme motif database used for AME, and copy all the enriched motifs into a new .meme file. Annoying and time-consuming. 
Could be coded probably.


7.) Search entire genome for enriched motifs with FIMO to get their locations.
Bash script (FIMO_Whole_Genome.sh):
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N Motif_Hunt_GW
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=168:00:00,vmem=64gb

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

fimo --o /scratch/jandrews/Data/Motif_Analyses/Results/FIMO/HG19_FLvDL_Enriched /scratch/jandrews/Data/Motif_Analyses/Motif_Databases/enriched_NonTSS_FAIRE_sigdiff_H3AC_K27AC_FLvDL_2016.meme /scratch/jandrews/Ref/hg19.fa


8.) Convert FIMO output to BED format for intersections. Must have BEDOPS installed (or write a script, find another program to do so, whatever).
gff2bed < fimo.gff > fimo.bed


9.) Remove garbage chromosomes, haplotypes, etc.
(sed '/_g/d' fimo.bed | sed '/chrM/d' | sed '/chrY/d') > fimo_clean.bed

File should now be ready for intersection with variant files.


###-To prepare the GM ChIP-seq file Jackie had made for downstream intersections with variants, etc.

1). Initial parsing to remove "Rep1/2/3" from TFs, remove 4th column with TF count at position, and remove replicate duplicates.
Python script (fix_gm_chip.py): python fix_gm_chip.py GM12878_TF151_merge1.bed GM12878_TF151_fixed.bed


2.) Go through output file and manually change the necessary TF names to corresponding genes that can be found in Sarah's FPKM table.


3.) Last runthrough to remove duplicate TFs for a given region (due to different antibodies, etc). Warns user of TFs not found in the FPKM table. Skips fully-numeric TF entries.
Python script (fix_gm_chip_round2.py): python fix_gm_chip_round2.py GM12878_TF151_fixed.bed All_GENEFPKMtracking_DLFLVGACC_genetype.txt GM12878_TF151_final.bed

GM file should now be ready to intersect with variant files.

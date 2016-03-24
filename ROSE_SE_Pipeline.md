Up to date as of 01/25/2015
Jared's SE Pipeline using ROSE - Young et al, 2013. Highly advise doing this on a computing cluster, as it will take you ages otherwise. Some steps were broken up into multiple
for the sake of clarity. Several can easy be combined/piped together.

All necessary scripts should be in Jared's code folder: N:\Bioinformatics\Jareds_Code
Some scripts (especially bash scripts) may be slightly different than listed here, but they should function in the same way. 
All python scripts use python 3, but ROSE requires python2.7. Be aware of which you need and use a virtualenv (anaconda on cluster) as needed.

Highly recommend using individual peak files to call SEs for each sample. Merging BAMs for a given cell type and calling SEs on that gives muddy results unless all your samples
sequenced very well (in which case you probably don't want to go that route anyway!).

1.) Download ROSE (https://bitbucket.org/young_computation/rose) and stick somewhere. Personal bin folder is always good.


2.) Place all BAMs (K27AC or TF (MED1, BRD4, etc) ChIP-Seq) into a directory.


3A.) Sort BAMs first. Can just be done from command line if you only have a few:
Bash script (bam_sort.sh)
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N bam_sort
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=4:00:00,vmem=8gb

module load samtools-1.2

for file in /scratch/jandrews/Data/ChIP-Seq/K27ME3/Batch2/*.bam; do

	base=${file##*/}
	samtools sort -T $file.sorted -o /scratch/jandrews/Data/ChIP-Seq/K27AC/Batch2/${base%.*}.sorted.bam $file &
	
done
wait
module remove samtools-1.2


3B.) Then index:
Bash script (bam_index.sh)
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N bam_index
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=4:00:00,vmem=4gb

module load samtools-1.2

for file in /scratch/jandrews/Data/ChIP-Seq/K27ME3/Batch2/*.bam; do

	samtools index $file &
	
done
wait
module remove samtools-1.2


4.) Place 5 sorted BAMs with indexes (K27AC or TF ChIP-Seq) into each batch directory.

TO-DO UPDATE
5.) Call peaks with MACS.
Bash script (peak_call.sh)- Could likely be parallelized if a bit of thought went into it, just annoying when some files have input controls and others don't: 
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N peak_call
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

macs14 -t /scratch/jandrews/Data/ChIP-Seq/K27AC/Batch8/FL202_K27AC.bam -c /scratch/jandrews/Data/ChIP-Seq/INPUT/Batch4/FL202_INPUT.bam -f BAM -g hs -n /scratch/jandrews/Data/ChIP-Seq/MACS/FL202_K27AC -w -S --nomodel --shiftsize=150


5.) Place all peaks.bed files from MACS into a single directory. 


6.) Remove garbage chromosomes and unnecessary columns. Run below command from within folder containing the peaks.bed files for each sample.

for F in *.bed; do
	base=${F##*/}
	(sed '/_g/d' "$F" | sed '/chrM/d' | sed '/chrY/d' | cut -f 1-4) > ${base%.*}.clean.bed ;
	cut -f 1-4 ${base%.*}.clean.bed > "$F"
	rm ${base%.*}.clean.bed
done


7.) Convert to gff format.
These files will be used as the "enhancers" that are used by ROSE.
If on cluster, set appropriate version of python as default, not necessary if done locally:
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python3  /scratch/jandrews/bin/ROSE_bed2gff.py <input.bed>

###For individual samples:
for F in *.bed; do
	python3 /scratch/jandrews/bin/ROSE_bed2gff.py "$F"
done

8.) Run ROSE. 
-t specifies areas around TSS to exclude peaks for stitching. Can be omitted if wanted.
NOTE: ROSE is stupid and won't run properly if the output folder isn't a new folder. Be sure to delete old results or specify a new output folder before running again if you want to play with the settings.

Bash script (ROSE_ind.sh)
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N ROSE_ind
#PBS -m e
#PBS -l nodes=1:ppn=4,walltime=8:00:00,vmem=12gb

module load samtools-1.2
module load R

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

cd /scratch/jandrews/bin/rose/

for file in /scratch/jandrews/Data/ChIP_Seq/K27AC/Batch20/*.bam; do
	
	echo "$file"
	base=${file##*/}
	python ROSE_main.py -g HG19 -t 2500 -r "$file" -i /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks/PEAKS_GFF/${base%.*}_peaks.gff -o /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks/RESULTS/${base%.*} &
	
done
wait
module remove samtools-1.2
module remove R


9A.) Annotate results (uses ref_seq annotations provided with ROSE).
Bash script (ROSE_annotate.sh): 
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N ROSE
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

module load samtools-1.2
module load R

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

cd /scratch/jandrews/bin/rose/

for fold in /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_All_Merged_Peaks/*/; do

	python ROSE_geneMapper.py -f -g HG19 -i "$fold"/*SuperEnhancers.table.txt -o "$fold" ;
	
done
wait
module remove samtools-1.2
module remove R


9B.) Annotate with Gencode (v19, genes only, from our master annotation files) to retain info on lincs, etc.
i.) First sort (bash script - sort.sh):
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N sort
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

for fold in /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_All_Merged_Peaks/*/; do

	sort -k 1,1 -k2,2n  "$fold"*SuperEnhancers.bed > "$fold"ROSE_SuperEnhancers_sorted.bed & 
	
done
wait

ii.) Then annotate (batch script - bedtools_closest.sh):
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N bedtools_closest
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=8gb

module load bedtools2

for fold in /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_All_Merged_Peaks/*/; do
	
	cd "$fold"
	bedtools closest -a ROSE_SuperEnhancers_sorted.bed -b /scratch/jandrews/Ref/gencode.v19.annotation_sorted_genes_only.bed -d -t all > "$fold"${PWD##*/}_ROSE_SuperEnhancers_Sorted_Gencode_Annotated.bed &
	
done
wait
module remove bedtools2


10.) Merge the two annotations into a single file.
Bash script (merge_SE_annotations.sh):
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N merge_SE_annotations
#PBS -m e
#PBS -l nodes=1:ppn=4,walltime=4:00:00,vmem=4gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for fold in /scratch/jandrews/Data/ChIP_Seq/MACS/ROSE_SEs_From_All_Merged_Peaks/*/; do
	
	cd "$fold"
	mv *_SuperEnhancers_ENHANCER_TO_GENE.txt ${PWD##*/}_ROSE_SuperEnhancers_Sorted_RefSeq_Annotated.txt
	python3 /scratch/jandrews/bin/merge_SE_annotations.py *Gencode* *RefSeq* ${PWD##*/}_ROSE_SEs_merged_annotations.bed ;
	sort -k1,1 -k2,2n ${PWD##*/}_ROSE_SEs_merged_annotations.bed > ${PWD##*/}_ROSE_SEs_merged_sorted_annotations.bed ;
	rm ${PWD##*/}_ROSE_SEs_merged_annotations.bed & 
	
done
wait


11.) Move files into a new folder. Create subfolders for each cell type and group appropriately (all CC in CC folder, DL in DL, etc.). Create one folder with all
files in it. Used for the first method below.


###-Two different methods to determine "unique" SEs below-###
The first takes into account overlap between SEs in different cell types, only calling those that overlap by less than 25% as unique.
The second method just takes all SEs for a given cell type, concatenates and merges them, and then does a multi-intersect with clustering.

#The first method (mine):


1.) This intersects all the SE files, merging those that overlap by >25% of either element. Has to be done twice to remove redundant overlaps. May have to cycle through a few more times
if you really crank the overlap up to like 0.99. Smaller overlap catches most of them the first time through. Idea is that A and B overlap, so it grabs the range for them. A and C also overlap, 
so it pulls the range for them as well. B and C don't overlap, so the range isn't pulled from them, but A+B and A+C overlap, so the second run through shores up that redundancy. Pulls the 
ranges of the overlaps first, then the sample column of the overlaps and pastes them together. Then parses output file to get SEs found in each cell type, unique to each cell type, and unique 
and recurrent to each cell type.

Bash script (bedops_full.sh):
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N BEDOPS_FULL
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for fold in /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks/SE_INTERSECTS/BEDOPS_Overlap_Method/All/; do

	cd "$fold"
    bedops --everything *.bed \
    | sort-bed - > all_files.bed ;
    bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 all_files.bed \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | uniq - > ranges.bed
    bedmap --echo-map-id-uniq --ec --sweep-all --fraction-either 0.25 ranges.bed all_files.bed > samples.bed ;
    paste ranges.bed samples.bed > ./Results/All_SEs.bed

    grep DL ./Results/All_SEs.bed > ./Results/DL_SEs.bed
    sed '/FL\|CC\|CB\|TS/d' ./Results/DL_SEs.bed > ./Results/Unique_DL_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_DL_SEs.bed ./Results/Recurrent_Unique_DL_SEs.bed ";" 
    sed '/TS\|CC\|CB/d' ./Results/DL_SEs.bed | sed '/FL/!d' - > ./Results/Unique_DL_FL_SEs.bed
    sed '/TS\|FL/d' ./Results/DL_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_DL_CC_CB_SEs.bed
    sed '/FL\|CC\|CB/d' ./Results/DL_SEs.bed | sed '/TS/!d' - > ./Results/Unique_DL_TS_SEs.bed
    grep FL ./Results/All_SEs.bed > ./Results/FL_SEs.bed
    sed '/DL\|TS\|CC\|CB/d' ./Results/FL_SEs.bed > ./Results/Unique_FL_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_FL_SEs.bed ./Results/Recurrent_Unique_FL_SEs.bed ";"
    sed '/DL\|CC\|CB/d' ./Results/FL_SEs.bed | sed '/TS/!d' - > ./Results/Unique_FL_TS_SEs.bed
    sed '/DL\|TS/d' ./Results/FL_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_FL_CC_CB_SEs.bed
    grep TS ./Results/All_SEs.bed > ./Results/TS_SEs.bed
    sed '/FL\|DL\|CC\|CB/d' ./Results/TS_SEs.bed > ./Results/Unique_TS_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_TS_SEs.bed ./Results/Recurrent_Unique_TS_SEs.bed ";"
    sed '/DL\|FL/d' ./Results/TS_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_TS_CC_CB_SEs.bed
    grep CC ./Results/All_SEs.bed > ./Results/CCCB_SEs.bed
    grep CB ./Results/All_SEs.bed >> ./Results/CCCB_SEs.bed
    sort-bed ./Results/CCCB_SEs.bed | uniq - > ./Results/CC_CB_SEs.bed
    sed '/FL\|TS\|DL/d' ./Results/CC_CB_SEs.bed > ./Results/Unique_CC_CB_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_CC_CB_SEs.bed ./Results/Recurrent_Unique_CC_CB_SEs.bed ";"

    rm ranges.bed
    rm samples.bed
    rm all_files.bed

done
wait



#The second method (Patrick's):

1.) For each cell type, concatenate the SEs for all samples into one file and sort.
cat *.bed | sort -k1,1 -k2,2n > CC_SEs_cat_sorted.bed


2.) Merge each SE list for each cell type.
module load bedtools2
mergeBed -i CC_CB_cat_sorted.bed -c 4 -o distinct > CC_CB_merged.bed


### (OPTIONAL, but RECOMMENDED) Isolate only those SEs that are recurrent within a given cell type and BETWEEN INDIVIDUALS 
(e.g. an SE in CC and CB from same individual does not count as recurrent). Actually only looks at CC/CB right now, edit script to extend functionality if needed.
##NOTE: THIS IS A MAJOR DRAWBACK OF THIS METHOD, AS YOU'RE REMOVING SEs THAT AREN'T RECURRENT BEFORE FILTERING FOR THOSE THAT ARE UNIQUE. THIS MEANS A NON-RECURRENT SE IN
A GIVEN CELL TYPE CAN'T BE USED TO FILTER NON-UNIQUE SEs IN OTHER CELL TYPES (IN WHICH IT MAY BE RECURRENT).##
	2B.) Run script on merged_SE files for each cell type to remove those only occurring in a single sample. 
	export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
	source activate anaconda
	for f in *.bed; do
		python /scratch/jandrews/bin/get_recurrent_SEs.py "$f" "$f".recurrent <delim of sample name column in quotes - ";">
	done
	rename .bed.recurrent _recurrent.bed *.bed


3.) Move all merged files into a single folder, then multiintersect. Clean up header if needed.
module load bedtools2
bedtools multiinter -cluster -header -i *.bed > All_SEs_multiinter.bed


4.) Parse each unique record to its specific sample file. Will produce a file for each data column containing the positions of the unique SEs for the column (each cell type in this case.)
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda
python parse_multiinter_output.py All_SEs_multiinter.bed


###-GET K27AC Signal For a set of SEs from BAMS-###
We try to use this to determine how "unique" our unique SEs really are. For instance, do DL-specific SEs as defined by our pipeline also have high K27AC signal in CC samples?
If so, are they really that unique? Outliers are easily identified by doing this as well, e.g., samples that are really driving the majority of the "unique" SE calls.

1.) Convert the resulting BED files for 'unique' SEs of a given cell type to GFF. 
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda
for f in *.bed; do
	python /scratch/jandrews/bin/ROSE_bed2gff.py "$f"
done


2.) Move the GFF files into a new folder and run ROSE_bamToGFF.py to get RPM'd K27AC load for each sample at each unique SE position. So will run on every sample multiple times, once for each cell type (DL, FL, TS, CC_CB in this case). '-r' options gives RPM results, '-m 1' is required for the script to run properly (sets the bin size to 1, they coded it lazily as hell).
This can be finicky and at times will randomly not run on the cluster if run as a job (literally no idea why). Just run from within the rose bin folder if so - copy and paste commands from here.
Batch script (ROSE_get_load.sh):
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N GET_SE_LOAD_DL
#PBS -m e
#PBS -l nodes=1:ppn=4,walltime=8:00:00,vmem=12gb

module load samtools-1.2
module load R

export PATH=~/Enthought/Canopy_64bit/User/bin:${PATH}
source ~/.bash_profile

cd /scratch/jandrews/bin/rose/

for fold in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/*/; do
	
	for file in "$fold"*.bam; do
		base=${file##*/}
		python ROSE_bamToGFF.py -b "$file" -m 1 -r -i /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks_Redo/UNIQUE_SES/Multiinter_Method/wOUT_VGA_VGR/UNIQUE_SES_GFF/DL_unique_SEs.gff -o /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks_Redo/UNIQUE_SES/Multiinter_Method/wOUT_VGA_VGR/UNIQUE_SES_K27AC_LOAD/UNIQUE_DL_SES_LOAD/${base%.*}_DL_SES_LOAD.gff ;
	done
done
wait
module remove samtools-1.2
module remove R


3.) Within the folder for K27AC load of unique SEs for each cell type (4 in this case), run script to calculate signal from the RPM values in each GFF file. (signal = RPM (density) * length of SE)
This is how ROSE calculates it. This script will take each file in the folder and essentially intersect them, yielding a single file in BED-type format with the K27AC signal for each sample
at each SE so that they may be easily plotted. I create a folder for 'Included' samples and 'Excluded' ones, as the previous step gets the load for all the BAMs, many of which I don't care about.

Python script (calc_SE_signal.py):
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda
python /scratch/jandrews/bin/calc_SE_signal.py <output.bed> <gff files>


4.) Create line graph with each SE as a data series, will plot a line for each SE showing signal in each sample that'll let you determine how "unique" they really are.
(Hint: not very). I like to take the signal from the Unique SEs of each cell type for each sample, copy them into excel, and calculate the log2(FC) of K27AC signal for each SE
over the median enhancer for that sample. Figure out the median enhancer from the all_enhancer output table from ROSE - includes signal. Can then smash these FC ratios for the SEs
unique to each cell type together and create a heatmap in R. This gives a better representation of how "unique" the calls really are, as it's comparing within the samples themselves.
I included 4 tonsil samples (bulk CD19+) in the unique analysis, but omitted them from the heatmap, as they're such a heterogeneous population that it causes a lot of noise.


5.) Create said heatmap in R. Keep Rowv = False to keep rows in order. I usually save several variations of each heatmap.  
R code:

data <- read.table("signal_table_from_excel.txt", header = TRUE)
data$START <- NULL
data$END <- NULL
data$CHROM <- NULL

data_dm <- data.matrix(data)

# Can drop rows at this point if needed: data_set <- data_dm[1:100,]

# Used to change the color scaling on the key. Can toy with the breaks to get the color scaling appropriate.
shades <- c(seq(-6,2.25,length=1),seq(2.26,6,length=500))

heatmap.2(no_21314_dm, density.info = "none", col=colorRampPalette(c("white","red4"))(500), margins = c(6,5), keysize = 1, cexRow = 0.6, cexCol = 0.7, trace="none", breaks=shades, main = "Unique Recurrent SEs vs Median Enhancer", key.xlab = "Log2(FC) of SE K27AC Signal over Median Enhancer", Rowv = FALSE, Colv=FALSE)



###-Determine amplifications/deletions that overlap SEs.
Need to run GISTIC to determine focal amplifications/deletions first. See CN_Analysis_Pipeline.txt for more details. 




###-Check K27AC of SEs using RPM'd, QN'd MMPID K27AC values for each sample that Liv generated.
1.) Get bed file of SE positions using cut or awk or whatever. Reorder the MMPID K27AC file if necessary.


2.) Intersect the SE positions with the MMPID K27AC values file. Probably have to take MMPID header off first, paste it into notepad or something to add back later.
cat MMPID_K27AC_RPM_QN_FINAL.bed ALL_SES_POSITIONS_SORTED.bed > SE_K27AC_VALS_CAT.bed
module load bedtools2
bedtools intersect -wa -wb -a ALL_SES_POSITIONS_SORTED.bed -b MMPID_K27AC_RPM_QN_FINAL_NOHEADER_CHR.bed > SE_MMPID_ISEC.bed


3.) Cut the columns with the MMPID positions, don't care about them. 
cut -f1-3,7- SE_MMPID_ISEC.bed > SE_MMPID_VALUES.bed


4.) Merge and sum all the values for the sample columns.
module load bedtools2
bedtools merge -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88 -o sum -i SE_MMPID_VALUES.bed > SE_SUMMED_MMPID_VALUES.bed


5.) Add back header and do whatever.

From here on out, can do whatever analyses you need. Actual documentation, what a wonder.

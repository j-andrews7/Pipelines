# Super Enhancer Calling

##### Up to date as of 03/25/2016

A super enhancer pipeline using ROSE - Young et al, 2013. Highly advise doing this on a computing cluster if possible. Some steps were broken up into multiple for the sake of clarity. Several can easy be combined/piped together.

This was done on the CHPC cluster, so all of the `export`, `source`, and `module load/remove` statements are to load the various software necessary to run the command(s) that follow. If you're running this locally and the various tools needed are located on your `PATH`, you can ignore these.

> Bash scripts are submitted on the cluster with the `qsub` command. Check the [CHPC wiki](http://mgt2.chpc.wustl.edu/wiki119/index.php/Main_Page) for more info on cluster commands and what software is available. All scripts listed here should be accessible to anyone in the Payton Lab, i.e., **you should be able to access everything in my scratch folder and run the scripts from there if you so choose.**

All necessary scripts should be here: **N:\Bioinformatics\Jareds_Code**  
They are also in `/scratch/jandrews/bin/` or `/scratch/jandrews/Bash_Scripts/` on the cluster as well as stored on my local PC and external hard drive.  

An _actual_ workflow (Luigi, Snakemake, etc) could easily be made for this with a bit of time, maybe I'll get around to it at some point.

**Software Requirements:**
- [BEDOPS](http://bedops.readthedocs.org/en/latest/index.html)
- [ROSE](https://bitbucket.org/young_computation/rose)
- [Samtools, BCFtools, and HTSlib](http://www.htslib.org/)  
  - These should be available on the CHPC cluster.
- [Python3](https://www.python.org/downloads/)
  - Use an [anaconda environment](http://mgt2.chpc.wustl.edu/wiki119/index.php/Python#Anaconda_Python) if on the CHPC cluster (also useful for running various versions of python locally).  
- [bedtools](http://bedtools.readthedocs.org/en/latest/)
  - Also available on the CHPC cluster.
  
---


##### 1A.) Sort BAMs first. 
**Bash script (bam_sort.sh)**
```Bash
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
```


##### 1B.) Then index:
**Bash script (bam_index.sh)**
```Bash
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
```


##### 2.) Place BAMs into batches.
Do the same with any INPUT controls, but in a separate directory. Be sure the INPUT BAM and it's corresponding
factor BAM are in a Batch directory with the same number, i.e. `./INPUT/BATCH1/sample1_input.sorted.bam` & 
`./K27AC/BATCH1/sample1_k27ac.sorted.bam`. Put the factor bams with no control in a separate directory, i.e. 
`./K27AC/NO_CTRL/BATCH1/sample_wNo_ctrl_k27ac.sorted.bam`.

This will make things much easier for peak calling later on. 


##### 3.) Remove ENCODE blacklist regions.
[ENCODE](https://sites.google.com/site/anshulkundaje/projects/blacklists) found that ChIP-seq experiments are inherently enriched for huge numbers of reads in specific regions of the genome (usually regions with tons of repeats like near centromeres, telomeres, etc). They suggest removing these reads before calling peaks, but this is even more important for super enhancer calling for which the signal of a given peak is taken into account. Regions with artificially high signal can lead to false positives that can affect the SE curves in a big way, introducing additional noise into data that is already often difficult to decipher.

**Bash script (bam_remove_blacklist.sh)**
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N RM_BLACKLIST_REGIONS

#PBS -m e

#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=32gb

module load samtools-1.2

for fold in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/Batch*/; do

	cd "$fold"
	for f in *.bam; do
		echo "$f"
		base=${f##/*}
		samtools view -b -t /home/jandrews/Ref/hg19.fa -L /home/jandrews/Ref/hg19_blacklist_regions_removed.bed -o "${base%.*}".BL_removed.bam "$f" &
	done
	wait
	
done
wait
module remove samtools-1.2
```


##### 4.) Index the blacklist removed BAMs.
Use the same script as above, just edit it slightly to grab the right files. 


##### 5.) Call peaks with MACS.
**Bash script (peak_call_bl_rmvd.sh and variations)**
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PEAK_CALL_BL_RMVD_22
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=36gb

batch=Batch22/
mark=_K27AC
treat_suffix=.sorted.BL_removed.bam
control_suffix=_INPUT.sorted.bam
treat=/scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/
control=/scratch/jandrews/Data/ChIP_Seq/BAMs/INPUT/

for f in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/"$batch"/*"$treat_suffix"; do
	base=${f##*/}
	samp=${base%%_*}
	macs14 -t "$f" -c "$control""$batch""$samp""$control_suffix" -f BAM -g hs -n /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/"$samp""$mark" -w -S --nomodel --shiftsize=150 &
done
wait
```
  
##### 6.) Consolidate peaks.bed files from MACS into a single directory. 
  
  
##### 7.) Scrub 'em.
Remove the garbage chromosomes and unnecessary columns. Run the below command from within folder containing the peaks.bed files for each sample.

```Bash
for F in *.bed; do
	base=${F##*/}
	(sed '/_g/d' "$F" | sed '/chrM/d' | sed '/chrY/d' | cut -f 1-4) > ${base%.*}.clean.bed ;
	cut -f 1-4 ${base%.*}.clean.bed > "$F"
	rm ${base%.*}.clean.bed
done
```
  
##### 8.) Convert to gff format.
These files will be used as the "enhancers" that are used by ROSE. If on the cluster, set appropriate version of python as default, not necessary if done locally.  
**Python script (ROSE_bed2gff.py)**
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for F in *.bed; do
	python3 /scratch/jandrews/bin/ROSE_bed2gff.py "$F"
done
```
  
  
##### 9.) Run ROSE. 
`-t` specifies areas around the TSS to exclude peaks for stitching. Can be omitted if wanted. ROSE is stupid and won't run properly if the output folder isn't a new folder. Be sure to delete old results or specify a new output folder before running again if you want to play with the settings.

**Bash script (ROSE_ind.sh)**
```Bash
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
```
  
  
##### 10A.) Annotate results.
Uses ref_seq annotations provided with ROSE.  
**Bash script (ROSE_annotate.sh)**
```Bash
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
```
  
  
##### 10B.) Annotate with Gencode. 
v19, genes only, from our master annotation files to retain info on lincs, etc.  
**i.) First sort.**  
**Bash script - (sort.sh)**
```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N sort
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

for fold in /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_All_Merged_Peaks/*/; do

	sort -k 1,1 -k2,2n  "$fold"*SuperEnhancers.bed > "$fold"ROSE_SuperEnhancers_sorted.bed & 
	
done
wait
```
  
**ii.) Then annotate.**  
**Bash script - (bedtools_closest.sh)**
```Bash
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
```
  
  
##### 11.) Merge the two annotations.  
This removes redundancies between the two annotations and creates a single file.  
**Bash script (merge_SE_annotations.sh)**
```Bash
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
```
  
  
##### 12.) Conglomerate files.  
Copy the merged annotation file for each sample into a single directory. This is used for the first method below.

---

### Two different methods to determine "unique" SEs  
The first takes into account overlap between SEs in different cell types, only calling those that overlap by **less than 25%** as unique. The second method just takes all SEs for a given cell type, concatenates and merges them, and then does a multi-intersect with clustering. If the SEs overlap at all between samples, they will be merged.

#### The first method (mine)  

##### 1.) Recurrently intersect, merge, and filter.  
This intersects all the SE files, **merging those that overlap by >25% of either element**. Has to be done several times to remove redundant overlaps. The idea is that SE_1 and SE_2 overlap, so it grabs the range for them. SE_1 and SE_3 also overlap, so it pulls the range for them as well. SE_2 and SE_3 don't overlap, so the range isn't pulled from them, but SE_1+SE_2 and SE_2+SE_3 overlap, so the second run through shores up that redundancy. It pulls the ranges of the overlaps first, then the sample column of the overlaps and pastes them together. The output file is then parsed to get SEs found in each cell type, unique to each cell type, and unique and recurrent to each cell type and pretty much every other comparison we could want here.

**Bash script (bedops_full.sh)**
```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N BEDOPS_FULL
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for fold in /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks/SE_INTERSECTS/BEDOPS_Overlap_Method/SELECT_SAMPLES/; do

    cd "$fold"
    bedops --everything *.bed \
    | sort-bed - > all_files.bed ;
    bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 all_files.bed \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | uniq - > ranges.bed
    bedmap --echo-map-id-uniq --ec --sweep-all --fraction-either 0.25 ranges.bed all_files.bed > samples.bed ;
    paste ranges.bed samples.bed > ./Results/All_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/All_SEs.bed ./Results/Recurrent_All_SEs.bed ";" 

    sed '/CC\|CB\|VG\|TS/!d' ./Results/All_SEs.bed > ./Results/NORMAL_SEs.bed
    sed '/FL\|DL\|CLL/d' ./Results/All_SEs.bed > ./Results/NORMAL_ONLY_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/NORMAL_ONLY_SEs.bed ./Results/Recurrent_NORMAL_ONLY_SEs.bed ";"
    sed '/FL\|DL\|CLL/!d' ./Results/All_SEs.bed > ./Results/TUMOR_SEs.bed
    sed '/CC\|CB\|VG\|TS/d' ./Results/All_SEs.bed > ./Results/TUMOR_ONLY_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/TUMOR_ONLY_SEs.bed ./Results/Recurrent_TUMOR_ONLY_SEs.bed ";"

    grep DL ./Results/All_SEs.bed > ./Results/DL_SEs.bed
    sed '/FL\|CC\|CB\|TS\|CLL\|VG/d' ./Results/DL_SEs.bed > ./Results/Unique_DL_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_DL_SEs.bed ./Results/Recurrent_Unique_DL_SEs.bed ";" 
    sed '/TS\|CC\|CB\|CLL\|VG/d' ./Results/DL_SEs.bed | sed '/FL/!d' - > ./Results/Unique_DL_FL_SEs.bed
    sed '/TS\|FL\|CLL\|VG/d' ./Results/DL_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_DL_CC_CB_SEs.bed
    sed '/FL\|CC\|CB\|CLL\|MEM\|VG/d' ./Results/DL_SEs.bed | sed '/NAIVE/!d' - > ./Results/Unique_DL_NAIVE_SEs.bed
    sed '/FL\|CC\|CB\|CLL\|NAIVE\|VG/d' ./Results/DL_SEs.bed | sed '/MEM/!d' - > ./Results/Unique_DL_MEM_SEs.bed
    sed '/FL\|CC\|CB\|CLL\|TS\|VGR/d' ./Results/DL_SEs.bed | sed '/VGA/!d' - > ./Results/Unique_DL_VGA_SEs.bed
    sed '/FL\|CC\|CB\|CLL\|TS\|VGA/d' ./Results/DL_SEs.bed | sed '/VGR/!d' - > ./Results/Unique_DL_VGR_SEs.bed
    sed '/FL\|CC\|CB\|VG\|TS/d' ./Results/DL_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_DL_CLL_SEs.bed


    grep FL ./Results/All_SEs.bed > ./Results/FL_SEs.bed
    sed '/DL\|CC\|CB\|TS\|CLL\|VG/d' ./Results/FL_SEs.bed > ./Results/Unique_FL_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_FL_SEs.bed ./Results/Recurrent_Unique_FL_SEs.bed ";"
    sed '/TS\|DL\|CLL\|VG/d' ./Results/FL_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_FL_CC_CB_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|MEM\|VG/d' ./Results/FL_SEs.bed | sed '/NAIVE/!d' - > ./Results/Unique_FL_NAIVE_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|NAIVE\|VG/d' ./Results/FL_SEs.bed | sed '/MEM/!d' - > ./Results/Unique_FL_MEM_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|TS\|VGR/d' ./Results/FL_SEs.bed | sed '/VGA/!d' - > ./Results/Unique_FL_VGA_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|TS\|VGA/d' ./Results/FL_SEs.bed | sed '/VGR/!d' - > ./Results/Unique_FL_VGR_SEs.bed
    sed '/DL\|CC\|CB\|VG\|TS/d' ./Results/FL_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_FL_CLL_SEs.bed


    grep NAIVE ./Results/All_SEs.bed > ./Results/NAIVE_SEs.bed
    sed '/DL\|CC\|CB\|MEM\|CLL\|VG\|FL/d' ./Results/NAIVE_SEs.bed > ./Results/Unique_NAIVE_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_NAIVE_SEs.bed ./Results/Recurrent_Unique_NAIVE_SEs.bed ";"
    sed '/FL\|DL\|CLL\|VG\|MEM/d' ./Results/NAIVE_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_NAIVE_CC_CB_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|FL\|VG/d' ./Results/NAIVE_SEs.bed | sed '/MEM/!d' - > ./Results/Unique_NAIVE_MEM_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|MEM\|VGR\|FL/d' ./Results/NAIVE_SEs.bed | sed '/VGA/!d' - > ./Results/Unique_NAIVE_VGA_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|MEM\|VGA\|FL/d' ./Results/NAIVE_SEs.bed | sed '/VGR/!d' - > ./Results/Unique_NAIVE_VGR_SEs.bed
    sed '/DL\|CC\|CB\|VG\|MEM\|FL/d' ./Results/NAIVE_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_NAIVE_CLL_SEs.bed


    grep MEM ./Results/All_SEs.bed > ./Results/MEM_SEs.bed
    sed '/DL\|CC\|CB\|NAIVE\|CLL\|VG\|FL/d' ./Results/MEM_SEs.bed > ./Results/Unique_MEM_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_MEM_SEs.bed ./Results/Recurrent_Unique_MEM_SEs.bed ";"
    sed '/FL\|DL\|CLL\|VG\|NAIVE/d' ./Results/MEM_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_MEM_CC_CB_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|NAIVE\|VGR\|FL/d' ./Results/MEM_SEs.bed | sed '/VGA/!d' - > ./Results/Unique_MEM_VGA_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|NAIVE\|VGA\|FL/d' ./Results/MEM_SEs.bed | sed '/VGR/!d' - > ./Results/Unique_MEM_VGR_SEs.bed
    sed '/DL\|CC\|CB\|VG\|NAIVE\|FL/d' ./Results/MEM_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_MEM_CLL_SEs.bed


    grep VGA ./Results/All_SEs.bed > ./Results/VGA_SEs.bed
    sed '/DL\|CC\|CB\|TS\|CLL\|VGR\|FL/d' ./Results/VGA_SEs.bed > ./Results/Unique_VGA_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_VGA_SEs.bed ./Results/Recurrent_Unique_VGA_SEs.bed ";"
    sed '/FL\|DL\|CLL\|VGR\|TS/d' ./Results/VGA_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_VGA_CC_CB_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|TS\|FL/d' ./Results/VGA_SEs.bed | sed '/VGR/!d' - > ./Results/Unique_VGA_VGR_SEs.bed
    sed '/DL\|CC\|CB\|VGR\|TS\|FL/d' ./Results/VGA_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_VGA_CLL_SEs.bed


    grep VGR ./Results/All_SEs.bed > ./Results/VGR_SEs.bed
    sed '/DL\|CC\|CB\|TS\|CLL\|VGA\|FL/d' ./Results/VGR_SEs.bed > ./Results/Unique_VGR_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_VGR_SEs.bed ./Results/Recurrent_Unique_VGR_SEs.bed ";"
    sed '/FL\|DL\|CLL\|VGA\|TS/d' ./Results/VGR_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_VGR_CC_CB_SEs.bed
    sed '/DL\|CC\|CB\|VGA\|TS\|FL/d' ./Results/VGR_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_VGR_CLL_SEs.bed


    grep CLL ./Results/All_SEs.bed > ./Results/CLL_SEs.bed
    sed '/DL\|CC\|CB\|TS\|VG\|FL/d' ./Results/CLL_SEs.bed > ./Results/Unique_CLL_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_CLL_SEs.bed ./Results/Recurrent_Unique_CLL_SEs.bed ";"
    sed '/FL\|DL\|VG\|TS/d' ./Results/CLL_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_CLL_CC_CB_SEs.bed


    grep CC ./Results/All_SEs.bed > ./Results/CCCB_SEs.bed
    grep CB ./Results/All_SEs.bed >> ./Results/CCCB_SEs.bed
    sort-bed ./Results/CCCB_SEs.bed | uniq - > ./Results/CC_CB_SEs.bed
    sed '/FL\|TS\|DL\|VG\|CLL/d' ./Results/CC_CB_SEs.bed > ./Results/Unique_CC_CB_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_CC_CB_SEs.bed ./Results/Recurrent_Unique_CC_CB_SEs.bed ";"

    rm ranges.bed
    rm samples.bed
    rm all_files.bed

done
wait
```


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

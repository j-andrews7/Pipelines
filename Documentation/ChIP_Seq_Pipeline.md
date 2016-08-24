# ChIP-SEQ Pipeline
**Last updated 08/04/2016**  
Author: jared.andrews07@gmail.com  

This document describes the bioinformatics pipeline used to analyze the Payton Lab's histone ChIP-seq data. This pipeline is pretty linear, but the `.wig` and `peaks.bed` files **can be handled separately** until the last portion of the pipeline in which they are normalized and merged. Additional file manipulations may be necessary (removal of headers, switching columns around, etc), though considerable effort has been made to minimize this as much as possible. **This is not the end-all, be-all, but it should be a good place to start.** The scripts were created/maintained **by 4 different people over several years**, though efforts have been made to streamline things recently.

This was done on the CHPC cluster, so all of the `export`, `source`, and `module load/remove` statements are to load the various software necessary to run the command(s) that follow. If you're running this locally and the various tools needed are located on your `PATH`, you can ignore these. They're more so I can just copy and paste when running through this myself.

> Bash scripts are submitted on the cluster with the `qsub` command. Check the [CHPC wiki](http://mgt2.chpc.wustl.edu/wiki119/index.php/Main_Page) for more info on cluster commands and what software is available. All scripts listed here should be accessible to anyone in the Payton Lab.

All necessary scripts should be here: **N:\Bioinformatics\Jareds_Code**  
They are also in `/scratch/jandrews/bin/` or `/scratch/jandrews/Bash_Scripts/` on the cluster as well as stored on my local PC and external hard drive.  

An _actual_ workflow (Luigi, Snakemake, etc) could easily be made for this with a bit of time, maybe I'll get around to it at some point.

**Software Requirements:**  
- [Samtools](http://www.htslib.org/)  
  - This should be available on the CHPC cluster.
- [Python3 & 2.7](https://www.python.org/downloads/)
  - Use an [anaconda environment](http://mgt2.chpc.wustl.edu/wiki119/index.php/Python#Anaconda_Python) if on the CHPC cluster (also useful for running various versions of python locally).  
  - MACS requires Python 2.7.
- [bedtools](http://bedtools.readthedocs.org/en/latest/)
  - Also available on the CHPC cluster.
- [Perl](https://www.perl.org/)
  - For old legacy scripts. Disgusting, I know. Will replace with python in time.
- [MACS 1.4](http://liulab.dfci.harvard.edu/MACS/)
  - MACS is our peak caller of choice. There are tons of them out there, but it's a pretty popular one.
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - This isn't necessary if you already have aligned BAMs that you're working from.
- [R](https://www.r-project.org/)
  - One or two of the python scripts invoke some R code.
  - May also need the preprocessCore package installed. Just use R studio to install.

  
#### Sections  
- [Preprocessing and Peak Calling](#peak-calling)
- [Peak Processing](#peak-processing)
- [Binning and Normalization](#binning-and-normalization)
- [Normalize Peak Regions Only](#normalize-peak-regions-only)
- [Making RPM Tracks Directly From Bams](#making-tracks)
- [Other Useful Scripts](#other-useful-scripts)


---

## Peak Calling  
This section walks through aligning the fastQ files to the genome, sorting and indexing the bam files, and actually calling the peaks.

#### 1.) Create bowtie indexes if necessary.  
These are necessary for bowtie to work properly. Download your `genome.fa` for your organism/reference of choice and run the build command on it. The second argument is just the prefix to attach to the index files.  
```Bash
bowtie2-build hg19.fa hg19
```

#### 2.) Align fastQs then sort and index BAMs.  
I like to use an interactive job on the CHPC cluster for this.  

> Histone mark files have relatively few duplicate reads in my experience (<5%), but TF ChIP is almost certainly a different story. Might pay to test a few files with/without duplicates removed to see how they compare in peak calling. I don't think it's a huge issue, but something to keep in mind at least. MACS is supposed to handle duplicates pretty well.

```Bash
qsub -I -l nodes=1:ppn=8,walltime=15:00:00,vmem=16gb
module load bowtie2
module load samtools-1.2

for f in *.fastq; do
  bowtie2 -p 8 -x /scratch/jandrews/Ref/hg19 -U "$f" | samtools view -bS - | samtools sort - > ${f%.*}.sorted.bam ;
  samtools index ${f%.*}.sorted.bam
done
```

#### 3.) Place BAMs into batches.  
Do the same with any INPUT controls, but in a separate directory. Be sure the INPUT BAM and it's corresponding
factor BAM are in a Batch directory with the same number, i.e. `./INPUT/BATCH1/sample1_input.sorted.bam` & 
`./K27AC/BATCH1/sample1_k27ac.sorted.bam`. Put the factor bams with no control in a separate directory, i.e. 
`./K27AC/NO_CTRL/BATCH1/sample_wNo_ctrl_k27ac.sorted.bam`.

This will make things much easier for peak calling later on. 


#### 4.) Remove ENCODE blacklist regions.  
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

#### 5.) (Optional) Remove PCR duplicates from BAMs.
MACs should really handle duplicates fine on its own, so I wouldn't worry too much about this step. This is how I'd do it though. Most of the files only have 3-5% duplicate reads, so I doubt it makes much of a difference.

**Bash script (bam_remove_dups.sh):**  
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N REMOVE_DUPS
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb

module load java
module load R

for fold in /scratch/jandrews/Data/ChIP_Seq/K27AC/Batch*/; do

		cd "$fold"
		for f in *.BL_removed.bam; do
			echo "$f"
			base=${f%%.*}
			java -jar /export/picard-tools-2.0.1/picard.jar	MarkDuplicates INPUT="$f" OUTPUT="$base".dups_rmvd.bam REMOVE_DUPLICATES=true METRICS_FILE="$base".dup_metrics.txt &
		wait
		rename .dups_rmvd.bam .BL_rmvd.dups_rmvd.sorted.bam *.dups_rmvd.bam
		done	
wait
done

module remove R
module remove java
```

#### 6.) Index the blacklist removed BAMs.  
Same command used previous.  


#### 7.) Call peaks with MACS.  
Ensure the `--shiftsize` option is appropriate for your data. For FAIRE I use `--shiftsize=50`, for K4ME3 `--shiftsize=100`, and for other histone marks I use `--shiftsize=150`. For TFs, I usually see things in the `--shiftsize=50 to 100` range, but it's probably worth trying to read more about it. Or trying to let MACS figure out the shiftsize itself (remove the `--nomodel` & `--shiftsize` options).

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

I **highly recommend** at least taking a *peek* at the `peaks.xls` files for each sample. A low number of peaks called is indicative of poor sequencing quality. You can also run the files through a QC package like [ChIPQC](http://bioconductor.org/packages/release/bioc/html/ChIPQC.html), which give additional metrics and are pretty easy to use. 
  
---
  
## Peak Processing  
This section explains how to handle the `peaks.bed` files that are output from MACs. It's really just a matter of cleaning the data up and merging the peaks for all the samples.

#### 1.) Consolidate peaks.bed files from MACS into a single directory for a given mark/TF. 
  
  
#### 2.) Scrub 'em.  
Remove the garbage chromosomes and unnecessary columns. Run the below command from within folder containing the peaks.bed files for each sample. The **python script (rename_columns.py)** replaces the 4th column with the actual sample name. 

```Bash
for F in *.bed; do
	base=${F##*/}
	sed '/_g\|chrM\|chrY\|chrX\|random\|chr23/d' "$F" | cut -f 1-4 | sort -k1,1 -k2,2n - > ${base%.*}.clean.bed ;
	cut -f 1-4 ${base%.*}.clean.bed > "$F"
	rm ${base%.*}.clean.bed
done

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python /scratch/jandrews/bin/rename_columns.py *.bed
```

#### 3.) Concatenate and merge the peak files.  
This merges overlapping peaks between samples so that we end up with just one set of genomic positions for the peaks.

```Bash
cat *.fixed.bed | sort -k1,1 -k2,2n - > catsort_peaks.bed

module load bedtools2
bedtools merge -c 4 -o distinct,count_distinct -i catsort_peaks.bed > merged_peaks.bed
```

#### 4.) Add the mark ID to the 4th column.  
This script adds a peak ID to the 4th column.  
```Bash
perl /scratch/jandrews/bin/Add_XID_tocol.pl merged_peaks.bed merged_peaks_ID.bed
args: <Mark_ID> 4 1
```

#### 4.b) (Optional) Add MMPID to the 4th column.  
If you want to merge these peaks with peaks from another mark or factor, you can assign an MMPID to them now. This script is quite old and probably needs to be edited due to changes made earlier in the pipeline, as it adds counts as well. This step may be extraneous at this point, keeping it here for historical purposes.

NOTE: Make sure delimiter for fourth column list is a "," and that format of each entry is "samplename_mark".

**Python script (add_MMPID_count_marks.py):**  
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python3 /scratch/jandrews/bin/add_MMPID_count_marks.py input.bed output.bed
```

#### 5.) Bin the genome.  
We need to create bins for the resolution we want to calculate RPM, do normalizations with, etc. Previously, we've done 200 bp bins, but I don't care how big the files are, so I'm doing 100 bp here. The `hg19.bed` file just contains the size of each chromosome - `chr1	249250621`, etc. 

**Perl script (binning_genome_binnumrestart_JAedit.pl):**  
```Bash
perl binning_genome_binnumrestart_JAedit.pl hg19.chrom.sizes hg19.100bp_bins.bed
Input: 50
```

#### 6.) Intersect the merged peak file with the binned genome.  
Also going to go ahead and cut a bunch of columns, grabbing only the ones for the position, bin #, and mark ID. Then sort lexicographically.

```Bash
module load bedtools2

bedtools intersect -wa -wb -a /scratch/jandrews/Ref/hg19.100bp_bins.bed -b merged_peaks_ID.bed | cut -f1-4,8 - | sort -k1,1 -k2,2n - > merged_peaks_bins_intersect.bed
```

Now we can set this file aside while we normalize the actual signal for each sample.

---

## Binning and Normalization  
This section shows how to process the wig files output from MACS to get the load for each sample across the genome.

#### 1.) Bin the wig files.  
This script breaks the wig files up into 100 bp bins with the signal for each sample. You can edit the `$bin_size` variable at the beginning of the script to adjust the bin size to match the bins you used for your `peaks.bed` file if it's not 100. We don't need to do this for any control/input samples we have. This will output a `.bin` file for each sample. Also be sure edit the path within the Perl script to your `chrom.sizes` file. 

This Perl script **sucks** and is incredibly slow. Run it overnight if you're doing it on a large number of files. I break it up into batches so that it can run on files in parallel.

**Bash script (bin_wigs.sh)**

```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N bin_wigs
#PBS -m e
#PBS -q old
#PBS -l nodes=1:ppn=12,walltime=24:00:00,vmem=64gb

for fold in /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/WIGS/K27AC/TREAT/Batch*/; do

	cd "$fold"

	for f in "$fold"/*.wig; do
		perl /scratch/jandrews/bin/bin_whole_genome_wig.pl "$f" &
	done
	wait

	# Move these guys to a differet folder for processing.
	mv *.bin /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/BIN_PROCESSING/K27AC/
done
```

#### 2.) Combine the `.bin` files.  
This script will combine the .bin files for a given mark into a table. It will also remove chromosomes and such. I run it in an **interactive session** because it eats memory. Sort it afterwards.

**Python script (combine_bins.py):**  
```Bash
qsub -I -l nodes=1:ppn=1,walltime=4:00:00,vmem=32gb
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python3 /scratch/jandrews/bin/combine_bins.py bin_directory master_table.bin

(head -n1 master_table.bin && tail -n+2 master_table.bin | sort -k1,1 -k2,2n) > master_table_sorted.bin
mv master_table_sorted.bin master_table.bin
``` 


#### 3.) Convert .bin table to .bed format.  
Now we want to add the chromosome positions of each bin in between the chromosome and bin columns for our master table. This uses the binned genome we created previously. It'll prompt you for some input, use the values below.

**Perl script (bin_to_chromcoords_JAedit.pl):**   
```Bash
perl /scratch/jandrews/bin/bin_to_chromcoords_JAedit.pl /scratch/jandrews/Ref/hg19.100bp_bins.bed master_table.bin
Input prompt args: 3 2 1
```

#### 4.) Get number of aligned and total reads for each sample.  
These numbers are necessary for the next step.

```Bash
module load samtools

for fold in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/Batch*; do
	cd "$fold";
	for f in *BL_removed.bam; do
		mapped="$(samtools view -c -F 4 "$f")"
		echo "${f%%.*}" '=>' "${mapped}",;
	done
done

# For those without control bams as well.
for fold in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/NO_CTRL/Batch*; do
	cd "$fold";
	for f in *BL_removed.bam; do
		mapped="$(samtools view -c -F 4 "$f")"
		echo "${f%%.*}" '=>' "${mapped}",;
	done
done
```

#### 5.) Edit the RPM normalization script.  
Now we're going to edit the next script with the total mapped reads from the previous step. For each sample you have, just add the number to the hash at the beginning of the script (**Calc_RPM_ChIP_Seq.pl**). *Yes, this is obnoxious, I know, but you should just be able to copy and paste the values in*. 


#### 6.) Normalize by reads for each sample.  
This script will use the header from the input bed file and compare it to the hash. It will throw errors if it can't find something it's looking for, so make sure everything you expect to be found is actually found. This will calculate RPKMs for each bin for every sample.

**Bash script (calc_rpm.sh)**

```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N CALC_RPM
#PBS -m e
#PBS -q old
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=64gb

cd /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/BIN_PROCESSING/K27AC

perl /scratch/jandrews/bin/Calc_RPM_ChIP_Seq.pl chrcoords_master_table.bed master_table_RPKM.bed
```

#### 7.) Quantile normalize samples.  
This script will quantile normalize the samples for each bin, ironing out differences in the statistical properties between distributions of the samples. Makes it easier to compare across samples. It's usually used in microarray analysis. You have to denote the column in which the *data to actually QN starts* (5th column here).

> This script is tricky to run. It's tough to run on the CHPC cluster due to the way R is installed, so running it locally is easiest, provided you've got the memory to do it (hint, you don't). The file is likely quite large at this point. If you have to run it on the cluster, you'll have to install R in your home directory. I recommend R-3.2.1 since it includes a bunch of libraries you need that are left out in later versions. R is not memory efficient and your file may be **10gb+**, so this sucker takes a while to run.

**Bash scripts (quantile_normalize.sh):**

```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N quantile_normalize
#PBS -m e
#PBS -l nodes=1:ppn=1:haswell,walltime=24:00:00,vmem=96gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

cd /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/BIN_PROCESSING/K27AC

python /scratch/jandrews/bin/quantile_normalize_new.py master_table_RPKM.bed 5
```

#### 8.) Make UCSC tracks.  
Now our data is normalized and ready to be uploaded to UCSC. This script will fix the last bin in each chromosome so that UCSC doesn't throw errors and will create a gzipped bedgraph file from each data column that can be directly uploaded to UCSC. 

**Bash script (make_ucsc.sh):**

```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N MAKE_UCSC
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=64gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

cd /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/BIN_PROCESSING/K27AC

python /scratch/jandrews/bin/MakeUCSC.py -i QN_master_table_RPKM.bed -o QN_master_table_RPKM_final.bed
```

---

## Normalize Peak Regions Only
This section will take our big RPKM table we just generated and the peaks for our marks and condense the bins for each peak. This is really only necessary if you want to compare the previously created tracks for the QN'd bins to how quantile normalizing only the peak regions across samples to each other looks. *I don't think there's much of a difference, honestly.*

#### 1.) Get the files.  
Remember that `.bed` file we made for the peaks way back when? Yeah, we're going to need that again, so copy `merged_peaks_bins_intersect.bed` and the `master_table_RPKM.bed` file into the same folder.

#### 2.) Grab the bins for each peak.  
Yeah, let's do that. 

**Perl script (condense_by_bin_JAedit.pl):**

```Bash
perl /scratch/jandrews/bin/condense_by_bin_JAedit.pl merged_peaks_bins_intersect.bed master_table_RPKM.bed master_peaks_RPKM_condensed.bed
```

#### 3.) Quantile normalize samples.  
This script will quantile normalize the samples for each bin, ironing out differences in the statistical properties between distributions of the samples. Makes it easier to compare across samples. It's usually used in microarray analysis. You have to denote the column in which the *data to actually QN starts* (6th column here).

> This script is tricky to run. It's tough to run on the CHPC cluster due to the way R is installed, so running it locally is easiest, provided you've got the memory to do it (hint, you don't). The file is likely quite large at this point. R is not memory efficient and your file may be **10gb+**, so this sucker takes a while to run.

**Bash scripts (quantile_normalize.sh):**

```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N quantile_normalize
#PBS -m e
#PBS -l nodes=1:ppn=1:haswell,walltime=24:00:00,vmem=96gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

cd /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/BIN_PROCESSING/K27AC

python /scratch/jandrews/bin/quantile_normalize_new.py master_table_RPKM.bed 6
```

#### 4.) Condense by peak IDs and sum RPM values.  
This will get a single value for each peak for each sample. Useful for direct comparisons, etc

**Python script (sum_RPMs_merge_peakIDs.py):**  
```Bash
python3 sum_RPMs_merge_peakIDs.py QN_master_peaks_RPKM_condensed.bed QN_master_peaks_RPKM_condensed_final.bed
```

#### 5.) Make UCSC tracks.  
Now our data is normalized and ready to be uploaded to UCSC. This script will fix the last bin in each chromosome so that UCSC doesn't throw errors and will create a gzipped bedgraph file from each data column that can be directly uploaded to UCSC. You might have to *edit the script* to change the data column start.

**Bash script (make_ucsc.sh):**

```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N MAKE_UCSC
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=64gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

cd /scratch/jandrews/Data/ChIP_Seq/MACS/BL_REMOVED/BIN_PROCESSING/K27AC

python /scratch/jandrews/bin/MakeUCSC.py -i QN_master_table_RPKM.bed -o QN_master_table_RPKM_final.bed
```

#### 6.) (Optional) Calculate averages and fold change by cell type for each QN file from R. 
This script will calculate avgs and fold changes by cell type for each peak and print them to the output file. It will also create UCSC tracks for the data columns (the same as the `Make_UCSC.py` script) if you want it to. I haven't used this guy in a long while, but it's here if that data is something you might want.

**Python script (avg_fc_vals_by_celltype.py):**

```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python3 /scratch/jandrews/bin/avg_fc_vals_by_celltype.py -i QN_master_table_RPKM_condensed_final.bed -o QN_master_table_RPKM_condensed_final_FC.bed -ucsc False
```

---

## Making Tracks
This script will allow you to make **RPM** (reads per million mapped reads) bigwig tracks directly from `.bam` files that you can then observe in UCSC. This script requires **pybedtools** - `pip install pybedtools` and bedtools in your path. I iterate through folders each containing a few files, but you can also just chuck them all in a folder. I just didn't feel like waiting that long. Once done, just throw the bigwig files into a folder that can be seen from outside your network and link to them from UCSC. Saves you the hassle of uploading large files, though you do need to store them.

**Python script (bam_to_RPM_bigwig.py):**  
```Bash
qsub -I -l nodes=1:ppn=1,walltime=24:00:00,vmem=32gb 

module load bedtools2

for fold in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/Batch8*; do 
	cd "$fold"; 
	for f in *.BL_removed.bam; do 
		echo "$f"; 
		python /scratch/jandrews/bin/bam_to_RPM_bigwig.py "$f" "${f%.*}".bw; 
	done; 
done
```

---

## Other Useful Scripts   
These scripts might be useful for other analyses, but they're pretty complex, so I'm not going to go into them much here. Just read their usage statements to figure out what they do/how to use them.

Get genes in range of a loci, given a GTF file containing genes. Can also filter out those that overlap a TSS, etc. Useful for identifying REs that aren't promoters. Again, may want to read the script itself, as it has a number of options.  
`get_genes_in_range.py`

To do a variety of analyses between samples, use below script. Reading it is your best bet, it has many options that grant it power and flexibility. Good for comparisons between groups of samples, as it can calculate fold changes, p-vals, etc.  
`Analyze_Loci.py`

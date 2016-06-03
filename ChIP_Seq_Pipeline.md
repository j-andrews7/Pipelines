# ChIP-SEQ Pipeline
**Last updated 05/31/2016**  
Author: jared.andrews07@gmail.com  

This document describes the bioinformatics pipeline used to analyze the Payton Lab's histone ChIP-seq data. This pipeline is pretty linear, but the `.wig` and `peaks.bed` files **can be handled separately** until the last portion of the pipeline in which they are normalized and merged. Additional file manipulations may be necessary (removal of headers, switching columns around, etc), though considerable effort has been made to minimize this as much as possible. **This is not the end-all, be-all, but it should be a good place to start and may warn you of possible caveats ahead of time.**

This was done on the CHPC cluster, so all of the `export`, `source`, and `module load/remove` statements are to load the various software necessary to run the command(s) that follow. If you're running this locally and the various tools needed are located on your `PATH`, you can ignore these. They're more so I can just copy and paste when running through this myself.

> Bash scripts are submitted on the cluster with the `qsub` command. Check the [CHPC wiki](http://mgt2.chpc.wustl.edu/wiki119/index.php/Main_Page) for more info on cluster commands and what software is available. All scripts listed here should be accessible to anyone in the Payton Lab, i.e., **you should be able to access everything in my scratch folder and run the scripts from there if you so choose.**

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
  - May also need the preprocessCore package installed. Just use R studio.
  
#### Sections  
- [Preprocessing and Peak Calling](#peak-calling)
- [Peak Processing](#peak-processing)
- [Binning and Normalization](#binning-and-normalization)
- [Intersect with broad K4ME3 Peaks](#intersect-with-broad-k4me3-peaks)

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


#### 5.) Index the blacklist removed BAMs.
Same command used previous.  


#### 6.) Call peaks with MACS.
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

I **highly recommend** at least taking a *peek* at the `peaks.xls` files for each sample. A low number of peaks called is indicative of crappy sequencing quality. You can also run the files through a QC package like [ChIPQC](http://bioconductor.org/packages/release/bioc/html/ChIPQC.html), which give additional metrics and are pretty easy to use. 
  
## Peaks Processing
This section explains how to handle the `peaks.bed` files that are output from MACs. It's really just a matter of cleaning the data up and merging the peaks for all the samples.

#### 1.) Consolidate peaks.bed files from MACS into a single directory for a given mark/TF. 
  
  
#### 2.) Scrub 'em.
Remove the garbage chromosomes and unnecessary columns. Run the below command from within folder containing the peaks.bed files for each sample. The python script replaces the 4th column with the actual sample name. 

```Bash
for F in *.bed; do
	base=${F##*/}
	sed '/_g\|chrM\|chrY\|chrX\|random\|chr23/d' "$F" | cut -f 1-4 | sort -k1,1 -k2,2n - > ${base%.*}.clean.bed ;
	cut -f 1-4 ${base%.*}.clean.bed > "$F"
	rm ${base%.*}.clean.bed
done

python rename_columns.py *.clean.bed
```

#### 3.) Concatenate and merge the peak files.
This merges overlapping peaks between samples so that we end up with just one set of genomic positions for the peaks.



## Binning and Normalization
This section process the wig files to get the load for each sample across the genome in bins.

####-Wig to Bigwig & Reverse for Visualization of Peaks on genome browser-####
#Done to spot check data quality

--Create genome.chrom.sizes
fetchChromSizes <db (hg19, etc)> > <db>.chrom.sizes

--Convert wig to bigwig
wigToBigWig <inputwig.gz> <genome.chrom.sizes> <myBigWig.bw>

--Convert bigwig to wig
bigWigToWig in.bigWig out.wig

--Trackline for typical bigwig file loading to UCSC
track name=<name> type=bigWig bigDataUrl=http://example.path.edu/filename.bw

--To create bedgraphs from normalized files for UCSC
#browser position is where track will start in genome by default. below is CD69
#color is R,G,B
#can't get priority to work for ordering tracks
#autoscale off so that each track is scaled to same y axis limits
#viewLimits y axis min and max

browser position chr12:9896666-9921913 
track type=bedGraph name=T174_K27ACQN description=T174_K27ACQN visibility=2 color=200,0,0 altColor=0,100,200 priority=10 autoScale=off alwaysZero=off gridDefault=on graphType=bar viewLimits=0:3
#Note: Alternatively, the scripts detailed below now have the capability of manually creating bedgraph files for you.

#Copy and paste data from file with no header here.



####-MACS Bin Files Processing-####
--bin wig peak file into 200bp bins using perl script
bin_whole_genome_jpedit_spedit_jaedit.pl <peaks.wig files> 
note: must have hg19.chrom_sizes_trimmed file in same folder

--combine .bin files into 1 txt file for a given mark
python3 combine_bins.py <input directory> <output.txt>
Note: This script requires no post-processing for sorting, removing junk chromosomes, etc. Must have natsort package installed to use though. 

--create binned genome
perl binning_genome_binnumrestart_JAedit.pl chrom_length.txt output.bed

--Turn bin files to .bed format, adding chromosome positions of each bin in between the chromosome and bin columns
perl bin_to_chromcoords_JAedit.pl <binned_chromosomes file> <any number of binned files>
input args: 3 2 1

--To get a count of number of mapped reads in a BAM file, which must be done for each sample to use for RPM:
samtools view -c -F 4 filename.bam

--To get a count of number of ummapped reads in a BAM file:
samtools view -c -f 4 filename.bam

--To get total number of reads for a BAM file
samtools view -c HG00173.chrom11.ILLUMINA.bwa.FIN.low_coverage.20111114.bam

--Use read counts from previous step to normalize by RPM for each sample name in the combined_bins file for each mark: use chrcoords_mark_BINS_sort_chr.bin files
perl Calculating_RPM_forchipseq_calc_JAedit.pl <input.bin with data starting in 4th column> <output file>
NOTE: Must edit the hash in the beginning of this script to add your sample names and the number of aligned reads in each. May also have to adjust the column in which data begins for the file within the script. 




####-MACS Bed Files Processing-####


##Plan to make a script to do the enclosed manual steps in one easy step eventually
--Rename 4th column in MACS14_peak.bed files to contain the sample name and mark. This is for those with additional info in the filename (TS081414_NAIVE_K27AC, for example.)
python3 rename_columns.py <MACS14_peak.bed files>
NOTE: filename format MUST be samplename_celltype,etc_histonemark... .bed
May just want to edit these as necessary

--Concatenate all bed files in a directory
cat *.bed > output.bed

--delete lines containing _g random genomic, M chromosome, Y chromosome loci
sed '/_g\|chrY\|chrX\|chr23/d' file.bed  > output.bed

--Keep only first 4 columns of bed file
cut -f 1-4 file.bed>file.bed

--Sort a bed file by chromosome and then by start
sort -k1.4,1.5 -k2 -V input.bed > output.bed

--Merge features in a bed file
mergeBed -c 4 -o collapse,count -i K27AC_PEAKS_catsort.bed > K27AC_PEAKS_merged.bed

--For merging of all merged peak files, add MMPID for each line to the fourth column. Place count of unique marks for the region in 6th column
python3 add_MMPID_count_marks.py <input.bed> <output.bed>
NOTE:Make sure delimiter for fourth column list is a "," and that format of each entry is "samplename_mark".

##If using ROSE to call SEs, convert the file to .gff format. This script adds IDs.
python3  convert_to_gff_ROSE.py <input.bed>

--For individual merged mark files, add ID to 4th column
perl Add_XID_tocol.pl <mergefile.bed> <outfile.bed>
args: <Mark_ID> 4 1

--create binned genome or copy one already made.
perl binning_genome_binnumrestart_JAedit.pl chrom_length.txt output.bed

--intersect each individual merged mark file with the binned genome file
bedtools intersect -a Binned_chromosomes.bed -b H3AC_PEAKS_final.bed -wa -wb -sorted > H3AC_PEAKS_intersect.bed

--cut all columns excepts chr#, start, stop, bin#, and mark_ID
cut -f1-4,8 H3AC_PEAKS_intersect.bed  > H3AC_PEAKS_intersect_cut.bed




####-BIN/BED CONVERGENCE-####
--move the <mark>_PEAKS_intersect_cut.bed files and <mark>_RPM.bed files into the same folder

--Condense the bins by comparing to the <mark>_PEAKS_intersect_cut.bed file for each mark. End up with an outfile file for each mark. Also includes peak IDs.
perl condense_by_bin_JAedit.pl <mark_PEAKS_intersect_cut.bed> <mark_RPM.bed> <output_filename>

--Condense by peak IDs and sum RPM for each sample for each peak
python3 sum_RPMs_merge_peakIDs.py <mark_RPM_condensed.bed> <output.bed (mark_RPM_condensed_summed.bed suggested)>




####-Quantile Normalization-####

#This actually uses R to perform the QN, but without the hassle of actually using R and its wonky output. Doesn't normalize to a known distribution (yet)
quantile_normalize.py <RPM table for multiple samples>

--Calculate averages and fold change by cell type for each QN file from R.
python3 avg_fc_vals_by_celltype.py -i <input.txt> -o <output.txt> -ucsc <True or False, False by default>
#Can upload gzipped bedgraph files directly to UCSC to visualize peaks.




####-Analysis-####

--Get genes in range of a loci, given a GTF file containing genes. Can also filter out those that overlap a TSS, etc. Useful for identifying REs that aren't promoters. Again, 
may want to read the script itself, as it has a number of options.
get_genes_in_range.py

--To do a variety of analyses between samples, use below script. Reading it is your best bet, it has many options that grant it power and flexibility.
Analyze_Loci.py



#####-UCSC VISUALIZATION-#####
#Can create bedgraph files from the original fold change values for viewing in UCSC. 

RECOMMENDED:
Just use the avg_fc_vals_by_celltype.py and MakeUCSC.py scripts on the condensed bins and uncondensed bins respectively.
They will create gzipped bedgraph files for each data column for you and save a lot of time.

#Create normalized tracks from QN on non-condensed bins for the samples/marks of interest. 
#These serve as useful comparators when viewing FC between cell types for each samples. Get the first three columns and sample column using awk.
awk -F'\t' '{print $1,$2,$3, $23}' OFS="\t" <input.txt> > <output.txt>

#Remove header from file in EmEditor, Textwrangler, using awk or sed, whatever. Chromosome end must also be corrected.
#Use the below script, as the binned end of chromosome isn't exact. UCSC requires exact ends. 
python3 chrom_bin_fix.py <input.bg> <output.bg>

#Then upload data into UCSC.
#For uploading, the track line NEEDS TO BE IN THE FILE. Otherwise it will just load as a BED file and won't display properly. 

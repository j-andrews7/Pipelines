# ChIP-SEQ Pipeline
**Last updated 08/18/2017**  
Author: jared.andrews07@gmail.com  

This document describes the bioinformatics pipeline used to analyze the Payton Lab's histone ChIP-seq data. This pipeline is pretty linear, but additional file manipulations may be necessary (removal of headers, switching columns around, etc), though considerable effort has been made to minimize this as much as possible. **This is not the end-all, be-all, but it should be a good place to start.**  This pipeline was originally created/maintained **by 4 different people over several years**, but recent advances in the field and development of new tools have allowed many of the homebrewed scripts to be removed. It's mostly composed of well-touted, commonly used tools and packages now.

This was done on the CHPC cluster, so all of the `export`, `source`, and `module load/remove` statements are to load the various software necessary to run the command(s) that follow. If you're running this locally and the various tools needed are located on your `PATH`, you can ignore these. They're more so I can just copy and paste when running through this myself.

> Bash scripts are submitted on the cluster with the `qsub` command. Check the [CHPC wiki](http://mgt2.chpc.wustl.edu/wiki119/index.php/Main_Page) for more info on cluster commands and what software is available. 

All necessary scripts should be [here](https://github.com/j-andrews7/Pipelines/tree/master/Code). They are also in `/scratch/jandrews/bin/` or `/scratch/jandrews/Bash_Scripts/` on the cluster as well as stored on my local PC and external hard drive.  

An _actual_ workflow (Luigi, Snakemake, etc) could easily be made for this with a bit of time, maybe I'll get around to it at some point.

**Software Requirements:**  
- [Samtools](http://www.htslib.org/)  
  - This should be available on the CHPC cluster.
- [Python3 & 2.7](https://www.python.org/downloads/)
  - Use an [anaconda environment](http://mgt2.chpc.wustl.edu/wiki119/index.php/Python#Anaconda_Python) if on the CHPC cluster (also useful for running various versions of python locally).  
  - MACS requires Python 2.7.
- [bedtools](http://bedtools.readthedocs.org/en/latest/)
  - Also available on the CHPC cluster.
- [MACS 2](https://github.com/taoliu/MACS)
  - MACS2 is our peak caller of choice. There are tons of them out there, but it's a pretty popular one.
  - Install with `pip install macs2`
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - This isn't necessary if you already have aligned BAMs that you're working from.
- [kentUtils](https://github.com/ENCODE-DCC/kentUtils)
  - Also on CHPC cluster.
- [R](https://www.r-project.org/)
  - Also need the DiffBind and ChIPQC packages installed.
```
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("ChIPQC", "DiffBind"))
```

  
#### Sections 
- [Spot Checking Files](#quick-spot-check)
- [Alignment](#alignment)
- [Peak Calling](#peak-calling)
- [QC with ChIPQC](#quality-control)
- [Normalize Peak Regions Only](#normalize-peak-regions-only)
- [Making Genome Browser Tracks](#making-tracks)
- [Other Useful Scripts](#other-useful-scripts)


---

## Quick Spot Check
First things first, check your actual sequence files to be sure your data isn't hot garbage before you go through all of this. I generally recommend [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), as it gives a good overview of how well your sequencing went. It's also dead easy to use and will indicate whether you might need to trim adaptors off your reads. It will tell you how many duplicate reads you have, which for histone ChIP-seq shouldn't be more than 20% or so. 

---

## Alignment

#### 1.) Create bowtie indexes if necessary.  
These are necessary for bowtie to work properly. Download your `genome.fa` for your organism/reference of choice and run the build command on it. The second argument is just the prefix to attach to the index files.  
```Bash
bowtie2-build hg19.fa hg19
```
#### 2.) Align and remove blacklisted reads.
This aligns the fastQ files to the genomes, sorts and indexes them, and removes reads overlapping [blacklisted regions](https://sites.google.com/site/anshulkundaje/projects/blacklists), which can screw with peak calling and comparisons.

**Bash script (align_remove_bl.sh)**
```
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N align_chip
#PBS -m be
#PBS -l nodes=1:ppn=8,walltime=72:00:00,vmem=128gb

module load bowtie2
module load samtools

# Directory containing reference genome fasta, bowtie indices, blacklist, etc.
genome=/scratch/jandrews/Ref/hg19
# Base directory containing reads and where additional directories will be created.
base_dir=/scratch/jandrews/Data/ChIP_Seq/T_Cell/READS

cd "$base_dir"
mkdir BAMS

# Alignment.
for f in *.fq.gz; do
  echo Aligning "$f"
  bowtie2 -p 6 -x "$genome" -U "$f" | samtools view -@ 8 -b -S -u - | samtools sort -@ 8 -m 5G -o ./BAMS/"${f%.*}".sorted.bam -;
  samtools index ./BAMS/"${f%.*}".sorted.bam
done

# Removal of blacklist reads.
cd BAMS
for f in *.bam; do
    echo Removing blacklisted reads from "$f"
    base=${f##/*}
    samtools view -@ 8 -b -t /scratch/jandrews/Ref/hg19.fa -L /scratch/jandrews/Ref/hg19_blacklist_regions_removed.bed -o "${base%.*}".BL_removed.bam "$f";
    samtools index "${base%.*}".BL_removed.bam;
done

module remove bowtie2
module remove samtools
```

---

## Peak Calling  
Now we're ready to call peaks with `MACS2`. First, move the `bam` files into batches - all the actual samples along with a single input into each folder that will be used for all of the samples in that folder.

#### 1.) Call peaks with MACS2.
`MACS2` has a lot of options and things you can tweak. Chief among them are the `--qvalue` and `--mfold` options. `qvalue` is the threshold for which to consider a peak significant and real - it's recommended to set it to `0.01` for ChIP-seq with relatively sharp peaks expected. Think transcription factors and sharp histone marks like H3K4me3, H3K27ac. It can be relaxed to `0.05` or `0.1` and the `--broad` option used if you expect broad peaks (H3K27me3, H3K9me3, etc). The `--mfold` option is used to select the regions within MFOLD range of high-confidence enrichment ratio against background to build the extension model. The regions must be lower than the upper limit, and higher than the lower limit of fold enrichment. DEFAULT:5,50 means using all regions not too low (>5) and not too high (<50) to build paired-peaks model. The `MACS2` author recommends playing with this value, and I tweak it here to be slightly more stringent with the regions used to build the paired-peaks model.

The `-B` option will output signal tracks in `bedGraph` format (which we'll convert to `bigWig` format later). The `--SPMR` option will scale this signal to signal per million reads, normalizing for read depth differences between samples. This will make our tracks look better and allow them to be much more appropriate for figures.

> The original `MACS` sometimes had trouble building an appropriate model that would extend reads to better represent the size of the actual DNA fragment and therefore, binding site for your protein. If this extension isn't done, you tend to get a pileup of reads on both sides of the actual peak from the forward and reverse strands, ending up with a dip in the middle. `MACS2` is *supposed* to be better at this, but you can still disable it with `--nomodel` and set the `--extsize` on your own if you want. For the old `MACS`, which used an option called `--shiftsize` instead of `--extsize`, for FAIRE I used `--shiftsize=50`, for H3K4me3 `--shiftsize=100`, and for other histone marks I used `--shiftsize=150`. For TFs, you should try to set it to the size of the typical binding site for that TF. **I don't mess with any of that here, I let `MACS2` build the model for me.**

**Be sure to adjust the effective genome size with `-g` if needed. Change to `mm` for mouse, as human is the default.**

**Bash script (peak_call_macs2_unmatched_ctrl.sh)**  
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PEAK_CALL_BL_RMVD_1
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate py2

# Idea here is to throw all the samples and the single input file for each of them in a given folder.
# So that this can be submitted in batches.

treat=/scratch/jandrews/Data/ChIP_Seq/T_Cell/BAMs/Batch1/

for f in "$treat"*C.sorted.BL_removed.bam; do
    base=${f##*/}
    macs2 callpeak -t "$f" -c "$treat"*INPUT.sorted.BL_removed.bam -n "$base" --outdir "$treat" -q 0.01 -B --SPMR -m 10 50 &
done
wait
```

This will spit out a bunch of files. I recommend reading the `MACS2` documentation to fully understand them all. For now, all we really care about are the `narrowPeak` and `bedGraph` files.

---

## Quality Control
You'll likely want to know how well your ChIP actually worked, right? Looking at the percentage of reads in peaks or in blacklisted regions can tell you a lot about how well the enrichment worked. And knowing how close to each other your replicates are doesn't hurt either. Fortunately, [ChIPQC](http://bioconductor.org/packages/release/bioc/vignettes/ChIPQC/inst/doc/ChIPQC.pdf) is a nifty `R` package that can do all that for you with little work on your part. The linked manual goes into way more detail than I will here, so check there if you get stuck.

#### 1.) Set up your experiment sample sheet.
This sample sheet will tell ChIPQC where to look for files, the different conditions or treatment between samples, which samples are recplicates, etc. An example is below.

| SampleID | Tissue | Factor  | Condition | Treatment  | Replicate | bamReads                                                       | ControlID | bamControl                                                   | Peaks                                                                                               | PeakCaller |
|----------|--------|---------|-----------|------------|-----------|----------------------------------------------------------------|-----------|--------------------------------------------------------------|-----------------------------------------------------------------------------------------------------|------------|
| HHu      | HH     | H3AC    | None      | UNTREATED  | 1         | ../BAMs/Batch11/HH.UNTREATED.H3AC.sorted.BL_removed.bam        | HHc       | ../BAMs/Batch11/HH.UNTREATED.INPUT.sorted.BL_removed.bam     | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HH.UNTREATED.H3AC.sorted.BL_REMOVED.peaks.narrowPeak           | bed        |
| HHr      | HH     | H3AC    | None      | ROMIDEPSIN | 1         | ../BAMs/Batch11/HH.ROMI.H3AC.sorted.BL_removed.bam             | HHc       | ../BAMs/Batch11/HH.UNTREATED.INPUT.sorted.BL_removed.bam     | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HH.ROMI.H3AC.sorted.BL_REMOVED.peaks.narrowPeak                | bed        |
| HUT78u   | HUT78  | H3AC    | None      | UNTREATED  | 1         | ../BAMs/Batch12/HUT78.UNTREATED.H3AC.sorted.BL_removed.bam     | HUT78c    | ../BAMs/Batch12/HUT78.UNTREATED.INPUT.sorted.BL_removed.bam  | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HUT78.UNTREATED.H3AC.sorted.BL_REMOVED.peaks.narrowPeak        | bed        |
| HUT78r   | HUT78  | H3AC    | None      | ROMIDEPSIN | 1         | ../BAMs/Batch12/HUT78.ROMI.H3AC.sorted.BL_removed.bam          | HUT78c    | ../BAMs/Batch12/HUT78.UNTREATED.INPUT.sorted.BL_removed.bam  | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HUT78.ROMI.H3AC.sorted.BL_REMOVED.peaks.narrowPeak             | bed        |
| OCILY7u  | OCILY7 | H3AC    | None      | UNTREATED  | 1         | ../BAMs/Batch14/OCILY7.UNTREATED.H3AC.sorted.BL_removed.bam    | OCILY7c   | ../BAMs/Batch14/OCILY7.UNTREATED.INPUT.sorted.BL_removed.bam | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/OCILY7.UNTREATED.H3AC.sorted.BL_REMOVED.peaks.narrowPeak       | bed        |
| HHu      | HH     | H3K27AC | None      | UNTREATED  | 1         | ../BAMs/Batch11/HH.UNTREATED.H3K27AC.sorted.BL_removed.bam     | HHc       | ../BAMs/Batch11/HH.UNTREATED.INPUT.sorted.BL_removed.bam     | ../MACS2/BL_REMOVED/NARROW_PEAK/H3K27AC/HH.UNTREATED.H3K27AC.sorted.BL_REMOVED.peaks.narrowPeak     | bed        |
| HHr      | HH     | H3K27AC | None      | ROMIDEPSIN | 1         | ../BAMs/Batch11/HH.ROMI.H3K27AC.sorted.BL_removed.bam          | HHc       | ../BAMs/Batch11/HH.UNTREATED.INPUT.sorted.BL_removed.bam     | ../MACS2/BL_REMOVED/NARROW_PEAK/H3K27AC/HH.ROMI.H3K27AC.sorted.BL_REMOVED.peaks.narrowPeak          | bed        |
| HUT78u   | HUT78  | H3K27AC | None      | UNTREATED  | 1         | ../BAMs/Batch12/HUT78.UNTREATED.H3K27AC.sorted.BL_removed.bam  | HUT78c    | ../BAMs/Batch12/HUT78.UNTREATED.INPUT.sorted.BL_removed.bam  | ../MACS2/BL_REMOVED/NARROW_PEAK/H3K27AC/HUT78.UNTREATED.H3K27AC.sorted.BL_REMOVED.peaks.narrowPeak  | bed        |
| HUT78r   | HUT78  | H3K27AC | None      | ROMIDEPSIN | 1         | ../BAMs/Batch12/HUT78.ROMI.H3K27AC.sorted.BL_removed.bam       | HUT78c    | ../BAMs/Batch12/HUT78.UNTREATED.INPUT.sorted.BL_removed.bam  | ../MACS2/BL_REMOVED/NARROW_PEAK/H3K27AC/HUT78.ROMI.H3K27AC.sorted.BL_REMOVED.peaks.narrowPeak       | bed        |
| OCILY7u  | OCILY7 | H3K27AC | None      | UNTREATED  | 1         | ../BAMs/Batch14/OCILY7.UNTREATED.H3K27AC.sorted.BL_removed.bam | OCILY7c   | ../BAMs/Batch14/OCILY7.UNTREATED.INPUT.sorted.BL_removed.bam | ../MACS2/BL_REMOVED/NARROW_PEAK/H3K27AC/OCILY7.UNTREATED.H3K27AC.sorted.BL_REMOVED.peaks.narrowPeak | bed        |

#### 2.) 

#### 9.) Make UCSC tracks from the peaks.bed files.
Uploading individual BED files to UCSC is annoying when you have dozens of samples. Use a [track hub](https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html) to keep track of your samples. It's a bit of a pain to set up, but it'll make your life much easier.

Anyway, you'll want to make bigBed files from the `peaks.bed` file for each sample. To do so, you'll need to round the 5th column of each file to an integer and make sure it's within the range of 0-1000 or the `bedToBigBed` utility will throw a fit. I made a python script to do just that. If the score is > 1000, it will just cap it at 1000. Additional commands to get rid of chromosomes we don't care about and convert to `bigBed`.

**Python script (div_peaks_bed_scores.py):**
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda
module load kentUtils
fetchChromSizes hg19 > hg19.chrom.sizes

for f in *.bed; do
	python /scratch/jandrews/bin/div_peaks_bed_scores.py "$f" "$f".round
	(sed '/_g\|chrM\|chrY\|chrX\|chr23/d' "$f".round) > "$f".final 
	bedToBigBed "$f".final hg19.chrom.sizes "$f".bb
	rm "$f"
done

rename .bed.bb .bb *.bed.bb
```

Now you're ready to hook it up in your track hub.

---
  

---

## Making Tracks
This script will allow you to make **RPKM** (reads per kilobase per million mapped reads) bigwig tracks directly from `.bam` files that you can then observe in UCSC. This script requires the Python package **deeptools** - `pip install deeptools`. I iterate through folders each containing a few files, but you can also just chuck them all in a folder. I just didn't feel like waiting that long. Once done, just throw the bigwig files into a folder that can be seen from outside your network and link to them from UCSC. Saves you the hassle of uploading large files, though you do need to store them.

*You can also **subtract input reads** from these if you want, see the other 'input_subtracted' version of this script for that.*

**If doing it for RNA-seq, remove the `-e` option.**

**Bash script (make_chip_rpkm_tracks.sh & variants):**  
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MAKE_CHIP_RPKM_TRACKS_1
#PBS -m e
#PBS -q old
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=64gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

batch=Batch1/
mark=_K27AC
treat_suffix=.sorted.BL_removed.bam
treat=/scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/

# For K4ME3, set -e to 200, for FAIRE use 100, for other marks, use 300. This is just double the -shiftsize used for macs
# and is supposed to be the fragment length.

for f in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/"$batch"/*"$treat_suffix"; do
	base=${f##*/}
	samp=${base%%_*}
	bamCoverage  -e 300 -p 3 -of bigwig --normalizeUsingRPKM  -bl /scratch/jandrews/Ref/ENCODE_Blacklist_hg19.bed -b "$f" -o "$treat""$batch""$samp""$mark".BL_removed.rpkm.bw ;
done
wait
```

---

## Other Useful Scripts   
These scripts might be useful for other analyses, but they're pretty complex, so I'm not going to go into them much here. Just read their usage statements to figure out what they do/how to use them.

Get genes in range of a loci, given a GTF file containing genes. Can also filter out those that overlap a TSS, etc. Useful for identifying REs that aren't promoters. Again, may want to read the script itself, as it has a number of options.  
`get_genes_in_range.py`

To do a variety of analyses between samples, use below script. Reading it is your best bet, it has many options that grant it power and flexibility. Good for comparisons between groups of samples, as it can calculate fold changes, p-vals, etc.  
`Analyze_Loci.py`

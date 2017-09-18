# ChIP-SEQ Pipeline
**Last updated 09/18/2017**  
Author: jared.andrews07@gmail.com  

This document describes the bioinformatics pipeline used to analyze the Payton Lab's histone ChIP-seq data. This pipeline is pretty linear, but additional file manipulations may be necessary (removal of headers, switching columns around, etc), though considerable effort has been made to minimize this as much as possible. **This is not the end-all, be-all, but it should be a good place to start.**  This pipeline was originally created/maintained **by 4 different people over several years**, but recent advances in the field and development of new tools have allowed many of the homebrewed scripts to be removed. It's mostly composed of well-touted, commonly used tools and packages now.

This was done on the CHPC cluster, so all of the `export`, `source`, and `module load/remove` statements are to load the various software necessary to run the command(s) that follow. If you're running this locally and the various tools needed are located on your `PATH`, you can ignore these. They're more so I can just copy and paste when running through this myself.

> Bash scripts are submitted on the cluster with the `qsub` command. Check the [CHPC wiki](http://mgt2.chpc.wustl.edu/wiki119/index.php/Main_Page) for more info on cluster commands and what software is available. 

All necessary scripts should be [here](https://github.com/j-andrews7/Pipelines/tree/master/Code/ChIP_Seq). They are also in `/scratch/jandrews/bin/` or `/scratch/jandrews/Bash_Scripts/` on the cluster as well as stored on my local PC and external hard drive.  

An _actual_ workflow (Luigi, Snakemake, etc) could easily be made for this with a bit of time, maybe I'll get around to it at some point.

**Software Requirements:**  
- [Samtools](http://www.htslib.org/)  
  - This should be available on the CHPC cluster.
- [Python2.7](https://www.python.org/downloads/)
  - Use an [anaconda environment](http://mgt2.chpc.wustl.edu/wiki119/index.php/Python#Anaconda_Python) if on the CHPC cluster (also useful for running various versions of python locally).  
  - MACS2 requires Python 2.7.
- [bedtools](http://bedtools.readthedocs.org/en/latest/)
  - Also available on the CHPC cluster.
- [MACS2](https://github.com/taoliu/MACS)
  - MACS2 is our peak caller of choice. There are tons of them out there, but it's a pretty popular one.
  - Install with `pip install macs2`
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - This isn't necessary if you already have aligned BAMs that you're working from.
- [kentUtils](https://github.com/ENCODE-DCC/kentUtils)
  - Also on CHPC cluster.
- [R](https://www.r-project.org/)
  - Version 3.3.3 is probably what you'll want - I've had issues with some of the below packages with R 3.4+
  - Also need the DiffBind, BiocParallel and ChIPQC packages installed. I had trouble installing these on the CHPC cluster.
```
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("ChIPQC", "DiffBind", "BiocParallel"))
```

  
#### Sections 
- [Spot Checking Files](#quick-spot-check)
- [Alignment](#alignment)
- [Peak Calling](#peak-calling)
- [QC with ChIPQC](#quality-control)
- [Making Genome Browser Tracks](#making-tracks)
- [Differential Binding Analyses](#differential-binding-analyses)


---

## Quick Spot Check
First things first, check your actual sequence files to be sure your data isn't hot garbage before you go through all of this. I generally recommend [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), as it gives a good overview of how well your sequencing went. It's also dead easy to use and will indicate whether you might need to trim adaptors off your reads. It will tell you how many duplicate reads you have, which for histone ChIP-seq shouldn't be more than 20% or so. For TF ChIP, this might get higher, as the peaks are usually sharper. Any peakcaller worth it's weight will handle these appropriately, so it's usually not necessary to remove them, though it likely wouldn't hurt anything if you did. I've heard arguments for both choices.

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

The `-B` option will output signal tracks in `bedGraph` format (which we'll convert to `bigWig` format later). The `--SPMR` option will scale this signal to signal per million reads, normalizing for read depth differences between samples. However, these `bedGraph` files look like crap in browsers, so we'll make our own `bigWig` files from the `bam` files a bit later.

> The original `MACS` sometimes had trouble building an appropriate model that would extend reads to better represent the size of the actual DNA fragment and therefore, binding site for your protein. If this extension isn't done, you tend to get a pileup of reads on both sides of the actual peak from the forward and reverse strands, ending up with a dip in the middle. `MACS2` is *supposed* to be better at this, but you can still disable it with `--nomodel` and set the `--extsize` on your own if you want. For the old `MACS`, which used an option called `--shiftsize` instead of `--extsize`, for FAIRE I used `--shiftsize=50`, for H3K4me3 `--shiftsize=100`, and for other histone marks I used `--shiftsize=150`. For TFs, you should try to set it to the size of the typical binding site for that TF. **I don't mess with any of that here, I let `MACS2` build the model for me at least the first time around.**

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

This will spit out a bunch of files. I recommend reading the `MACS2` documentation to fully understand them all. For now, all we really care about are the `narrowPeak`.

---

## Quality Control
You'll likely want to know how well your ChIP actually worked, right? Looking at the percentage of reads in peaks or in blacklisted regions can tell you a lot about how well the enrichment worked. And knowing how close to each other your replicates are doesn't hurt either. Fortunately, [ChIPQC](http://bioconductor.org/packages/release/bioc/vignettes/ChIPQC/inst/doc/ChIPQC.pdf) is a nifty `R` package that can do all that for you with little work on your part. The linked manual goes into way more detail than I will here, so check there if you get stuck.

#### 1.) Set up your experiment sample sheet.
This sample sheet will tell ChIPQC where to look for files, the different conditions or treatment between samples, which samples are recplicates, etc. An example is below.

| SampleID        | Tissue | Factor | Condition | Treatment  | Replicate | bamReads                                                         | ControlID       | bamControl                                                       | Peaks                                                                                         | PeakCaller |
|-----------------|--------|--------|-----------|------------|-----------|------------------------------------------------------------------|-----------------|------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|------------|
| HH_unt_3        | HH     | H3AC   | None      | UNTREATED  | 1         | ../BAMs/Batch11/HH.UNTREATED.H3AC.sorted.BL_removed.bam          | HHc_unt_inp     | ../BAMs/Batch11/HH.UNTREATED.INPUT.sorted.BL_removed.bam         | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HH.UNTREATED.H3AC.sorted.BL_removed.peaks.narrowPeak     | narrow     |
| HH_romi_3       | HH     | H3AC   | None      | ROMIDEPSIN | 1         | ../BAMs/Batch11/HH.ROMI.H3AC.sorted.BL_removed.bam               | HHc_romi_inp    | ../BAMs/Batch11/HH.UNTREATED.INPUT.sorted.BL_removed.bam         | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HH.ROMI.H3AC.sorted.BL_removed.peaks.narrowPeak          | narrow     |
| HUT78_unt_3     | HUT78  | H3AC   | None      | UNTREATED  | 1         | ../BAMs/Batch12/HUT78.UNTREATED.H3AC.sorted.BL_removed.bam       | HUT78c_unt_inp  | ../BAMs/Batch12/HUT78.UNTREATED.INPUT.sorted.BL_removed.bam      | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HUT78.UNTREATED.H3AC.sorted.BL_removed.peaks.narrowPeak  | narrow     |
| HUT78_romi_3    | HUT78  | H3AC   | None      | ROMIDEPSIN | 1         | ../BAMs/Batch12/HUT78.ROMI.H3AC.sorted.BL_removed.bam            | HUT78c_romi_inp | ../BAMs/Batch12/HUT78.UNTREATED.INPUT.sorted.BL_removed.bam      | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HUT78.ROMI.H3AC.sorted.BL_removed.peaks.narrowPeak       | narrow     |
| OCILY7_unt_3    | OCILY7 | H3AC   | None      | UNTREATED  | 1         | ../BAMs/Batch14/OCILY7.UNTREATED.H3AC.sorted.BL_removed.bam      | OCILY7c_unt_inp | ../BAMs/Batch14/OCILY7.UNTREATED.INPUT.sorted.BL_removed.bam     | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/OCILY7.UNTREATED.H3AC.sorted.BL_removed.peaks.narrowPeak | narrow     |
| HHc_romi_inp    | HH     | Input  | None      | ROMIDEPSIN | c1        | ../BAMs/Batch11/HH.UNTREATED.INPUT.sorted.BL_removed.bam         | HHc_romi_inp    | ../BAMs/Batch11/HH.UNTREATED.INPUT.sorted.BL_removed.bam         | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HH.ROMI.H3AC.sorted.BL_removed.peaks.narrowPeak          | narrow     |
| HUT78c_romi_inp | HUT78  | Input  | None      | ROMIDEPSIN | c2        | ../BAMs/Batch12/HUT78.UNTREATED.INPUT.sorted.BL_removed.bam      | HUT78c_romi_inp | ../BAMs/Batch12/HUT78.UNTREATED.INPUT.sorted.BL_removed.bam      | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HUT78.ROMI.H3AC.sorted.BL_removed.peaks.narrowPeak       | narrow     |
| OCILY7c_unt_inp | OCILY7 | Input  | None      | UNTREATED  | c3        | ../BAMs/Batch14/OCILY7.UNTREATED.INPUT.sorted.BL_removed.bam     | OCILY7c_unt_inp | ../BAMs/Batch14/OCILY7.UNTREATED.INPUT.sorted.BL_removed.bam     | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/OCILY7.UNTREATED.H3AC.sorted.BL_removed.peaks.narrowPeak | narrow     |
| HHc_unt_inp     | HH     | Input  | None      | UNTREATED  | c4        | ../BAMs/Batch11/HH.UNTREATED.INPUT.sorted.BL_removed.Copy.bam    | HHc_unt_inp     | ../BAMs/Batch11/HH.UNTREATED.INPUT.sorted.BL_removed.Copy.bam    | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HH.UNTREATED.H3AC.sorted.BL_removed.peaks.narrowPeak     | narrow     |
| HUT78c_unt_inp  | HUT78  | Input  | None      | UNTREATED  | c5        | ../BAMs/Batch12/HUT78.UNTREATED.INPUT.sorted.BL_removed.Copy.bam | HUT78c_unt_inp  | ../BAMs/Batch12/HUT78.UNTREATED.INPUT.sorted.BL_removed.Copy.bam | ../MACS2/BL_REMOVED/NARROW_PEAK/H3AC/HUT78.UNTREATED.H3AC.sorted.BL_removed.peaks.narrowPeak  | narrow     |

It won't display your inputs on the same plot as the corresponding sample (contrary to what they show in their vignette) if you don't specifically list them as samples as well. Kind of annoying, as a given `bam` file can only correspond to one peakset, so if two different peaksets have the same input (as is the case here), you *have* to make a copy of the input `bam` file. If you don't care if the inputs are on the same plot as the other samples, feel free to leave them off. If you get an error like `'names' attribute must be the same length as the vector` make all of your input IDs unique and it should clear things up.

#### 2.) Load your sample sheet into `ChIPQC`.
Pretty simple. Don't type the `>`, it's just the code prompt. If you get the error stated above despite having all unique names, make sure `biocParallel` is loaded and enter `register(SerialParam())`. That will usually fix the error.

```R
> library(ChIPQC)
> samples = read.csv("example_QCexperiment.csv")

```

#### 3.) Create the ChIPQC Experiment.
Also really simple. This may take several hours to run. It will display some summary statistics about your experiment and create an interactive `html` summary report for your viewing. It includes a whole host of images and tables with nice captions that explain the various statistics and metrics for each sample. You can change the `facetBy` parameter if you want to group your samples on some other metadata like `Condition`.

```R
> experiment = ChIPQC(samples)
> experiment
> ChIPQCreport(experiment, facetBy=c("Tissue", "Treatment"))
```

#### 4.) Interpret your results.
Read the report, learn the metrics, and determine if your samples are of good quality. I pay particular attention to `RiP%` (Reads in Peaks) metric, which is a good measure of enrichment. If any need to be resubmitted for sequencing, now's the time. Save this so you can include some QC figures in the supplement of your fascinating future paper.

---

## Making Tracks
Differential peaks in a table are great and all, but it'd be better to *see* them, right? So let's whip up some tracks and put them in a place viewable to UCSC. Continous ChIP-seq data files are big, even when compressed, so uploading them to a genome browser isn't really feasible. Fortunately, there's an easy way around this - using a [track hub](https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html). 
These take some time to set up, but are easy to add tracks to and allow for easier viewing of your data. I won't go over how to create them here, but there are a few things to note about them. Most importantly, **your data has to be located in an area accessible outside your network**. This can be tricky if you don't have control over your firewall, and you might have to get IT involved. The genome browser needs to access these files to view the data for *the area that you currently want to view*. So the whole track isn't loaded at once, which really reduces memory and performance needs.
Secondly, the data files have to be in compressed formats, **so the `narrowPeak` and `bedGraph` formats won't work for this**. We need to convert them to `bigBed` and `bigWig` formats respectively, which is pretty easy.

#### 1.) Convert the `narrowPeak` files to `bigBed` format.
Easy enough with the handy [`narrowpeak2bb.sh` script](https://github.com/j-andrews7/Pipelines/blob/Master/Code/ChIP_Seq/narrowPeak2bb.sh). You'll need to edit this script to point to the location of wherever you put the `narrowPeak.as` file as well. Then just navigate to your directory containing the `narrowPeak` files, get your genome chromosome sizes, and run the script on the files.

```Bash
fetchChromSizes hg19 > hg19.sizes
for f in *.narrowPeak; do
	narrowpeak2bb.sh "$f" hg19.sizes
done
```

#### 2.) Convert the `bam` files to `bigWig` format.
Again, pretty easy. I use [deepTools](https://deeptools.github.io/) to create these files. The `-e` option should be set to your average fragment size (which can be guessed or determined from the output of `ChIPQC`). The `-bs` options sets the bin size, so it can be decreased for increased resolution or increased for smaller files. The `-bl` option allows you to specify a blacklist, though we've already removed these reads in this case. The `-b` option is for your input file, and `-o` is the output file. The `-p` options allows you to set the number of processing cores to be utilized. There are *many* other options also available, but you can figure those out on your own. Can normalize to 1x genome coverage, subtract input reads from your sample tracks, etc.

**Bash script (make_chip_rpkm_tracks.sh):
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MAKE_CHIP_RPKM_TRACKS_1
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=64gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

batch=Batch1/

cd /scratch/jandrews/Data/ChIP_Seq/T_Cell/BAMs/"$batch"

for f in *.BL_removed.bam; do
	bamCoverage  -e 200 -p 8 -of bigwig -bs 10 --normalizeUsingRPKM  -bl /scratch/jandrews/Ref/ENCODE_Blacklist_hg19.bed -b "$f" -o "$f".bw ;
done
wait

rename .bam.bw .RPKM.bw *.bw 
```

Now you just stick both these sets of files somewhere UCSC or another genome browser can access them and create your track hub.

---

## Differential Binding Analyses
Calling peaks is great, but presence/absence of peaks isn't the best metric to compare across samples. You'd likely rather compare normalized signal across samples for a consensus peak set. Depending on your experimental design, there are multiple ways/tools to do this. 

### With MAnorm
The most straight forward is by comparing two samples, say a treated and untreated sample, and comparing their signal for a peak set that merges all the peaks of those samples. You can then determine which peaks have significantly different signal between the two samples. The [MAnorm R/bash code](http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/MAnorm.htm) can do this for you easily enough. 

>I actually had to modify this code to get it running, as it hasn't been updated in some time. The `bedTools` arguments were out of date for a few commands (mostly just input file arguments were swapped for `coverageBed`). It also had a really high memory footprint for large files, which was taken care of by adding a few extra sorting steps. This code can be found [here](https://github.com/j-andrews7/MAnorm_updated).

#### 1.) Run MAnorm.
```Bash
./MAnorm.sh  sample1_peaks.bed  sample2_peaks.bed  sample1_read.bed  sample2_read.bed  150  150
```

The first 4 parameters should be input files in `bed` format with no header lines. The first 2 files have **ONLY 3 columns**: chromosome, start, end. The next 2 files should have 4 columns: chromosome, start, end, strand (+/-). Use `bamtobed` from `bedTools` to convert the `BAM` files to `bed` files. Use `cut` to generate the correct file format from the `narrowPeak` and `bed` files for peaks and reads, respectively. 
The last 2 parameters are the number of bp to be shifted for each read. These two parameters are found from the MACS `xls` peak files after "# d =".
MAnorm.r is called from MAnorm.sh, and there is no need to run it separately.  Check the file `Rcommand.out` for the output file from running the R script for error tracking.

This will create two tables - one with only the common merged peaks, and one with the common and unique peaks. I use the table with all of them and filter by the `-LOG10(pvalue)` column (>5, equivalent to p-value of 0.00001) and the `M-value` column. The M-value essentially refers to the magnitude of the change between the two samples for a given peak. **I use M-values of 1 and -1 in conjunction with the p-value filter to identify the peaks enriched in my two different conditions**. You can just do this in Excel, or use the `classfy_by_peaks.sh` script from MAnorm. 

#### 2.) Plot results and compare signal at differential peaks.
Now you've got a bunch of lists, but you probably want to compare those lists for different conditions and visualize those differences.

You've really got options at this point. You can use [plot.ly](https://plot.ly/) with R/Python/Javascript/MatLab/a web editor to create some really pretty, *interactive* graphs and plots. Another new program I'm very fond of is [EaSeq](easeq.net), which makes generating genome-wide occupancy maps and manipulating your data in fairly complex ways a breeze. It can also help you annotate your peaks by nearest gene, etc. It is, unfortunately, a **Windows only program** though.

#### 3.) Pathway and motif enrichment analyses.
You may also be interested in if the genes your differential peaks are near happen to be similar in any way. This is where pathway enrichment analyses come in handy. [GREAT](http://bejerano.stanford.edu/great/public/html/index.php) is my personal favorite tool for this, as all it requires is a list of genomic regions and will perform all the gene ontology (GO) analyses you'll ever need. It hits a lot of databases.

[AME](http://meme-suite.org/tools/ame) is what I typically use for motif enrichment analyses, which are useful if you want some idea of what transcription factors may be enacting the changes in binding you've found. It's part of the enormous MEME suite, which also has a ton of other tools if you're doing TF ChIP or trying to determine the binding motif for a transcription factor/protein that doesn't have a known one.

### With DiffBind
`DiffBind` is a nifty R package written by the same group that did `ChIPQC` - it even uses the same sample sheet. Remove any poor quality samples before this step and be sure to add any metadata you may need to the sample sheet (Treatment, Conditions, etc).



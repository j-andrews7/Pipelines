# CN Analysis 

**Up to date as of 05/09/2016.**  
jared.andrews07@gmail.com

> This version differs from the "old" version because it **treats the CNVs on a sample-by-sample basis for the MMPID/Circuit analysis** in order to keep them as small as possible (i.e. it doesn't merge them together for samples of a given cell type). It also tries to find the **"golden ticket"** CNVs by searching for **minimal common regions**. Also gets into creating some figures to actually try to show the effects of the CNVs on the activity of SEs, lncRNAs, and MMPIDs.

The aim of this pipeline is to get all copy number changes for all samples for which we have SNP arrays. These are then intersected with CNAs identified in other publications or with our SE, lincRNA, and MMPID data. As with most things, this started off relatively simple and grew to be more complicated as results were analyzed and additional approaches tried.

This was done on the CHPC cluster, so all of the `export`, `source`, and `module load/remove` statements are to load the various software necessary to run the command(s) that follow. If you're running this locally and the various tools needed are located on your `PATH`, you can ignore these.

> Bash scripts are submitted on the cluster with the `qsub` command. Check the [CHPC wiki](http://mgt2.chpc.wustl.edu/wiki119/index.php/Main_Page) for more info on cluster commands and what software is available. All scripts listed here should be accessible to anyone in the Payton Lab, i.e., **you should be able to access everything in my scratch folder and run the scripts from there if you so choose.**

All necessary scripts should be here: **N:\Bioinformatics\Jareds_Code**  
They are also in `/scratch/jandrews/bin/` or `/scratch/jandrews/Bash_Scripts/` on the cluster as well as stored on my local PC and external hard drive.  

An _actual_ workflow (Luigi, Snakemake, etc) could easily be made for this with a bit of time, maybe I'll get around to it at some point.

#### Software Requirements
- [BEDOPS](http://bedops.readthedocs.org/en/latest/index.html)
- [Samtools](http://www.htslib.org/)  
  - This should be available on the CHPC cluster.
- [Python3](https://www.python.org/downloads/)
  - Use an [anaconda environment](http://mgt2.chpc.wustl.edu/wiki119/index.php/Python#Anaconda_Python) if on the CHPC cluster (also useful for running various versions of python locally).  
    -Some of the scripts use NumPy, matplotlib, seaborn, pybedtools, and pandas. I'm too lazy to post links, but just activating your virtualenv and using `pip install numpy`, etc should get these all installed very easily.
- [bedtools](http://bedtools.readthedocs.org/en/latest/)
  - Also available on the CHPC cluster.
- [Affymetrix Genotyping Console](http://www.affymetrix.com/estore/browse/level_seven_software_products_only.jsp?productId=131535#1_1)

#### Sections  
- [Calling CNVs from SNP6 Arrays](#cnv-calling)
- [Finding Minimal Common Regions Across CNVs](#finding-minimal-common-regions-across-cnvs)
- [Integrating SEs with CNVs by Cell Type](#integrate-se-data)
- [Integrating Circuit Table Data to Filter MMPIDs in CNVs](#integrating-circuit-table-data-for-mmpids)
- [Creating Boxplots of SE Signal on a Sample-by-Sample CNV Basis](#observe-se-signals-inside-and-outside-cnvs)
- [Comparing lncRNA Expression in CNVs](#comparing-lincrna-expression-in-cnvs)

---

## CNV Calling
This section focuses on generating the CNVs for each sample that was run on an Affymetrix SNP6 array.


#### 1.) Set up environment.
Download the [Affymetrix Genotyping Console](http://www.affymetrix.com/estore/browse/level_seven_software_products_only.jsp?productId=131535#1_1) program and the na32 annotation db. You'll have to set a library folder when opening the program for the first time - this is where you should stick annotation/databse files. Download the na32 library from within the program and stick it in your library folder, along with the na32 CN annotation file from Affy's site. Again, this is all for the **SNP 6 arrays**, and I used the older annotations (na32 instead of na35) for continuity with other analyses. 


#### 2.) Conglomerate files.
Take all cel & arr files you want to use and stick them in a folder. Load the data into the program. Run the Copy Number/LOH analysis, it'll generate a CNCHP file for each sample. Go to Workspace->Copy Number/LOH Analysis and export them, making sure to include the Log2 ratio for each marker (along with it's chromosome), marker ID, and position.


#### 3.) Run the Segment Reporting Tool. 
Be sure to click the option to generate a summary file, which will have the segments for each sample. 


#### 4.) Scrub the summary file.
The first column of the summary file will have the array ID, which should be replaced with the sample name. In addition, the `#` header lines can be removed and the only columns that need to be kept are those for `sample, Chr, Start_Linear_Position, End_Linear_Position, #Markers, Copy_Number_State`. `awk/sed/cut/paste` make swapping the columns around pretty easy.

#### 5.) Convert to bed-like format.
Just need to order columns - `Chr, Start, End, Sample, #Markers, Copy_Number_State` and add `chr` to first column. This also removes segments with fewer than 5 markers in them if they haven't been removed already.

```Bash
awk -v OFS='\t' '{print $2, $3, $4, $1, $5, $6}' Lymphoma_SNP77_CNCHP_011615_segment_summary_cleaned.txt \
| awk '{print "chr" $0}' \
| awk '{
if ($5 >=5) 
print $0
}' > Lymphoma_SNP77_CNCHP_011615_segment_summary_cleaned.bed
```

#### 6.) Grab segments specific to each cell-type.
Go ahead and delete those on the Y chromosome as well.
```Bash
grep 'DL' Lymphoma_SNP77_CNCHP_011615_segment_summary_cleaned.bed | sed '/chrY/d' > DL_CNVs.bed
grep 'FL' Lymphoma_SNP77_CNCHP_011615_segment_summary_cleaned.bed | sed '/chrY/d' > FL_CNVs.bed
grep 'CLL' Lymphoma_SNP77_CNCHP_011615_segment_summary_cleaned.bed | sed '/chrY/d' > CLL_CNVs.bed
```

#### 7.) Break each CNV file into amps/dels
Before we merge the CNVs for each cell-type, need to break the amps and dels into separate, sorted files.

**Amplifications**
```Bash
awk '{
if ($6 >= 2)
print $0
}' DL_CNVs.bed | sort -k1,1 -k2,2n > DL_AMPS.bed
```

**Deletions**
```Bash
awk '{
if ($6 <= 1)
print $0
}' DL_CNVs.bed | sort -k1,1 -k2,2n > DL_DELS.bed
```

#### 8.) Merge the amps/dels for each cell-type.
This also creates a list of the samples in which each CNV is found as well as the actual copy number of the regions merged, which can *sometimes* be helpful for determining if there are any high-level samples (multiple copy gain/loss) for that CNV, though you can't determine which sample it's necessarily in without going back and looking in each. It will pull the max marker numbers for the merged region as well.

```Bash
module load bedtools2
mergeBed -c 4,5,6 -o distinct,max,collapse -i CLL_AMPS.bed > CLL_AMPS_MERGED.bed
```

#### 9.) Annotate known CNVs.
While CNAs unique to the samples are likely to be most interesting, a CNV found in other cells and cancers may still have a functional impact. As such, we **don't remove common CNVs**, rather we simply make a note that the CNV is found elsewhere. There are many lists of common CNVs (HapMap, DGV, and a [2015 Nature Genetics Review](http://www.nature.com/nrg/journal/v16/n3/full/nrg3871.html) among them). Here, I use the [DGV gold standard list](http://dgv.tcag.ca/dgv/app/downloads?ref=) after some funky parsing to remove extraneous info. 


```Bash
sed -n '2~3p' DGV.GoldStandard.July2015.hg19.gff3 \
| python /scratch/jandrews/bin/parse_DGV_gff.py > DGV.GoldStandard.July2015.hg19.parsed.gff3

module load bedtools2

# Variable to make switching between cell types easy.
type=CLL
bedtools intersect -loj -a "$type"_AMPS_MERGED.bed -b ../DGV.GoldStandard.July2015.hg19.parsed.gff3 \
| cut -f1-6,15 \
| sort -k1,1 -k2,2n \
| uniq - \
| bedtools merge -c 4,7 -o distinct -i - > "$type"_AMPS_MERGED_ANNOT.bed

bedtools intersect -loj -a "$type"_DELS_MERGED.bed -b ../DGV.GoldStandard.July2015.hg19.parsed.gff3 \
| cut -f1-6,15 \
| sort -k1,1 -k2,2n \
| uniq - \
| bedtools merge -c 4,7 -o distinct -i - > "$type"_DELS_MERGED_ANNOT.bed

# Do for non-merged amps and such as well.
bedtools intersect -loj -a "$type"_AMPS.bed -b ../DGV.GoldStandard.July2015.hg19.parsed.gff3 \
| cut -f1-6,15 \
| sort -k1,1 -k2,2n \
| uniq - > "$type"_AMPS_ANNOT.bed

bedtools intersect -loj -a "$type"_DELS.bed -b ../DGV.GoldStandard.July2015.hg19.parsed.gff3 \
| cut -f1-6,15 \
| sort -k1,1 -k2,2n \
| uniq - > "$type"_DELS_ANNOT.bed
```

#### 10.) Filter for recurrence.
This can easily be done by grepping for the `,` delimiter in our sample column for the merged files.

```Bash
# Variable to make switching between cell types easy.
type=CLL
cat "$type"_AMPS_MERGED_ANNOT.bed | python /scratch/jandrews/bin/get_recurrent_CNAs.py > RECURRENT_"$type"_AMPS_MERGED_ANNOT.bed
cat "$type"_DELS_MERGED_ANNOT.bed | python /scratch/jandrews/bin/get_recurrent_CNAs.py > RECURRENT_"$type"_DELS_MERGED_ANNOT.bed
```

#### 11.) Filter for true CNAs.
Useful to do this for both those that are recurrent and those that are not. The period must be escaped as linux will otherwise consider '.' as any character. I like to move the different cell types into their own folders at this point as well.

```Bash
type=CLL
grep '\.' "$type"_AMPS_MERGED_ANNOT.bed > "$type"_AMPS_MERGED_ANNOT_CNA.bed
grep '\.' RECURRENT_"$type"_AMPS_MERGED_ANNOT.bed > RECURRENT_"$type"_AMPS_MERGED_ANNOT_CNA.bed
grep '\.' "$type"_DELS_MERGED_ANNOT.bed > "$type"_DELS_MERGED_ANNOT_CNA.bed
grep '\.' RECURRENT_"$type"_DELS_MERGED_ANNOT.bed > RECURRENT_"$type"_DELS_MERGED_ANNOT_CNA.bed
```

---

## Finding Minimal Common Regions Across CNVs
This utilizes the files generated in the [CNV calling section](#cnv-calling) to try to identify the minimal common regions (MCRs) between the CNVs in each sample. This is *kind of* what [GISTIC](http://www.broadinstitute.org/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=216&p=t) tries to do.

#### 1.) Bin the Genome.
Chopping the genome into bins will allow for greater resolution of which segments of the CNVs actually overlap between samples. I used 5k bins here, you can use whatever you'd like. But for **plotting**, you should use **at least 25kb bin** if you're trying to plot the entire genome. Smaller bins can be used if you make multiple plots (one for each chromosome). The genome file should just have the start (0) and end of each chromosome.

```Bash
bedops --chop 5000 hg19.bed > hg19.5kb_bins.bed
```

#### 2.) Intersect CNVs with the Binned Genome.
This will create a new column for each sample with a count for the number of elements that overlap each bin. After iterating through all the files, we'll end up with a matrix that we can use to plot after a bit of tweaking. This can be done with CNVs for **only one cell type** or **combining them all together**. I show doing it for both DL and FL samples here. 

```Bash
module load bedtools2

cp hg19.5kb_bins.bed FLDL_AMPS_MATRIX_5KB.bed
cp hg19.5kb_bins.bed FLDL_DELS_MATRIX_5KB.bed

# First the amps.
for f in *AMPS_ANNOT*; do
	echo "$f"
	bedtools intersect -loj -c -a FLDL_AMPS_MATRIX_5KB.bed -b "$f" > FLDL_AMPS_MATRIX_5KB.bed.temp
	mv FLDL_AMPS_MATRIX_5KB.bed.temp FLDL_AMPS_MATRIX_5KB.bed
done

# Then the dels.
for f in *DELS_ANNOT*; do
	echo "$f"
	bedtools intersect -loj -c -a FLDL_DELS_MATRIX_5KB.bed -b "$f" > FLDL_DELS_MATRIX_5KB.bed.temp
	mv FLDL_DELS_MATRIX_5KB.bed.temp FLDL_DELS_MATRIX_5KB.bed
done
```

#### 3.) Scrub and condense amp/del matrices.
We need to do a few things before we can plot the data. First, the bins with multiple amps/dels overlapping need these values >1 to be reduced to 1. Second, the dels need their values converted to negatives. Lastly, the the amp/dels need to be merged for each bin and the file needs to be sorted numerically rather than lexicographically. This script does all of that. **Pay attention to file order, it *must* go amps then dels.** This uses a fair amount of memory, so I ran it in an interactive job on the cluster.

**Python script (condense_cn_matrices.py):**
```Bash
qsub -I -l nodes=1:ppn=4,walltime=15:00:00,vmem=16gb

python /scratch/jandrews/bin/condense_cn_matrices.py \ /scratch/jandrews/Data/CN_Analyses/CNVs_By_Cell_Type/CNV_MCRs/FLDL_AMPS_MATRIX_5KB.bed \ /scratch/jandrews/Data/CN_Analyses/CNVs_By_Cell_Type/CNV_MCRs/FLDL_DELS_MATRIX_5KB.bed \ /scratch/jandrews/Data/CN_Analyses/CNVs_By_Cell_Type/CNV_MCRs/FLDL_CN_MATRIX_5KB.bed \

```

#### 4.) Break up into chromosomes.
The entire genome doesn't look so great on a single figure. So let's break it up into chromosomes.

```Bash
mkdir 5KB_SPLIT_RESULTS
for chr in `cut -f 1 FLDL_CN_MATRIX_5KB.bed | uniq`; do
        grep -w $chr FLDL_CN_MATRIX_5KB.bed > 5KB_SPLIT_RESULTS/$chr.FLDL_CN_MATRIX_5KB.bed
done
```

#### 5.) Add a header.
These will be used as columns for plotting. It's also nice to know what's in a file.
```Bash
for f in *MATRIX*; do
	{ printf 'CHR\tSTART\tEND\tDL135\tDL166\tDL188\tDL191\tDL237\tDL252\tDL273\tDL3A193\tFL120\tFL125\tFL139\tFL153\tFL174\tFL202\tFL238\tFL255\tFL301\tFL313\tFL3A145\n'; cat "$f"; } > "$f".temp
done
mv "$f".temp "$f"

for f in ./5KB_SPLIT_RESULTS/*.bed; do
	{ printf 'CHR\tSTART\tEND\tDL135\tDL166\tDL188\tDL191\tDL237\tDL252\tDL273\tDL3A193\tFL120\tFL125\tFL139\tFL153\tFL174\tFL202\tFL238\tFL255\tFL301\tFL313\tFL3A145\n'; cat "$f"; } > "$f".temp
	mv "$f".temp "$f"
done
```

The next few steps are going to be a bit wonky in that they don't __*necessarily*__ have to be done in order. We're going to do **plotting** next, but the files we just created *will be used to find the MCRs* in a later step.


#### 6.) Plot the data.  
The first three columns will be used as the row labels. The rest of the data will be used as the column labels. This script can be edited to change the aesthetics and size. Note that if you try to make enormous figs, you may run into a segfault that results in your column labels not being printed, though the rest of the figure will display correctly. Adding lines between the columns/chromosomes, etc, in photoshop or powerpoint is helpful, and you'll probably want to relabel the columns.

*Note: I had trouble getting the Seaborn package to run on the cluster, so I ran this locally.*

**Note 2.0: The arrays didn't have probes for portions of certain chromosomes (chr13-15, 21, 22), so don't freak if these chunks are empty. Just crop them out of the plot since they're usually at the top.**

**Python script (plot_cn_bins.py):**
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python /scratch/jandrews/bin/plot_cn_bins.py DL_AMPS_MATRIX_5KB.bed 

cd 5KB_SPLIT_RESULTS
for f in *.bed; do
	python ../plot_cn_bins.py "$f"
done
```

#### 7.) Go back and find the MCRs.
I use bin sizes of 1kb for this, though you could go even smaller if you wanted even greater resolution. This is done for the amp and del matrices seperately. You can set the cutoff for the percentage of samples that must have the CNV for a given bin for it to be considered part of an MCR.

**Python script (find_cn_mcrs.py):**
```Bash
"""
Given a matrix of CN counts for a binned genome, identify bins that have the cnv for a given percentage of samples.
Merge them to identify minimal common regions between all the samples. Assumes the input file has a header.

Usage: python3 find_cn_mcrs.py -i <input.bed> -o <output.bed> -s <percentage as decimal>

Args:
    -i input.bed (required) = A matrix containing the bins for the genome and a column for each sample defining whether
    	or not the sample has a cnv for the bin.
    -o output.bed (required) = Name of output file.
    -p (optional) = Percentage of samples that must contain the cnv for a bin to be reported as a MCR. Default = 0.25.
"""
```

**Actual use:**
```Bash
# Script requires bedtools, so load it up.
module load bedtools2

python /scratch/jandrews/bin/find_cn_mcrs.py -i FLDL_AMPS_MATRIX_1KB.bed -o FLDL_AMPS_1KB.bed
python /scratch/jandrews/bin/find_cn_mcrs.py -i FLDL_DELS_MATRIX_1KB.bed -o FLDL_DELS_1KB.bed
```

Can take these output and intersect with lncRNAs, SEs, MMPIDs, etc.

---

## Integrate SE Data
This utilizes the files generated in the [CNV calling section](#cnv-calling) to look at the SEs located in CNVs on a cell-type basis.

#### 1.) Intersect CNVs with SEs and Enhancers.
Do every intersect you could ever really want. 

```Bash
# Use a variable to make it easy to switch cell types.
type=CLL

# Get the MMPIDs that lie outside the SEs.
module load bedtools2
bedtools intersect -v -a MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b All_SEs.bed > MMPID_NonTSS_FAIREPOS_OUTSIDE_ALL_SEs.bed
bedtools intersect -v -a MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b "$type"_SEs.bed > MMPID_NonTSS_FAIREPOS_OUTSIDE_"$type"_SEs.bed


# Intersect with the amps & dels.
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/All_SEs.bed -b "$type"_AMPS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/ALL_SEs_IN_"$type"_AMPS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/"$type"_SEs.bed -b "$type"_AMPS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/"$type"_SEs_IN_"$type"_AMPS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/All_SEs.bed -b RECURRENT_"$type"_AMPS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/ALL_SEs_IN_RECURRENT_"$type"_AMPS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/"$type"_SEs.bed -b RECURRENT_"$type"_AMPS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/"$type"_SEs_IN_RECURRENT_"$type"_AMPS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/All_SEs.bed -b "$type"_DELS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/ALL_SEs_IN_"$type"_DELS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/"$type"_SEs.bed -b "$type"_DELS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/"$type"_SEs_IN_"$type"_DELS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/All_SEs.bed -b RECURRENT_"$type"_DELS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/ALL_SEs_IN_RECURRENT_"$type"_DELS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/"$type"_SEs.bed -b RECURRENT_"$type"_DELS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/"$type"_SEs_IN_RECURRENT_"$type"_DELS.bed

# Then the CNAs.
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/All_SEs.bed -b "$type"_AMPS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/ALL_SEs_IN_"$type"_AMPS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/"$type"_SEs.bed -b "$type"_AMPS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/"$type"_SEs_IN_"$type"_AMPS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/All_SEs.bed -b RECURRENT_"$type"_AMPS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/ALL_SEs_IN_RECURRENT_"$type"_AMPS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/"$type"_SEs.bed -b RECURRENT_"$type"_AMPS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/"$type"_SEs_IN_RECURRENT_"$type"_AMPS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/All_SEs.bed -b "$type"_DELS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/ALL_SEs_IN_"$type"_DELS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/"$type"_SEs.bed -b "$type"_DELS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/"$type"_SEs_IN_"$type"_DELS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/All_SEs.bed -b RECURRENT_"$type"_DELS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/ALL_SEs_IN_RECURRENT_"$type"_DELS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/"$type"_SEs.bed -b RECURRENT_"$type"_DELS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/"$type"_SEs_IN_RECURRENT_"$type"_DELS_CNA.bed

# Then MMPIDs outside the SEs and in the amps/dels.
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_ALL_SEs.bed -b "$type"_AMPS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_ALL_SEs_IN_"$type"_AMPS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_"$type"_SEs.bed -b "$type"_AMPS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_"$type"_SEs_IN_"$type"_AMPS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_ALL_SEs.bed -b RECURRENT_"$type"_AMPS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_ALL_SEs_IN_RECURRENT_"$type"_AMPS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_"$type"_SEs.bed -b RECURRENT_"$type"_AMPS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_"$type"_SEs_IN_RECURRENT_"$type"_AMPS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_ALL_SEs.bed -b "$type"_DELS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_ALL_SEs_IN_"$type"_DELS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_"$type"_SEs.bed -b "$type"_DELS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_"$type"_SEs_IN_"$type"_DELS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_ALL_SEs.bed -b RECURRENT_"$type"_DELS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_ALL_SEs_IN_RECURRENT_"$type"_DELS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_"$type"_SEs.bed -b RECURRENT_"$type"_DELS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_"$type"_SEs_IN_RECURRENT_"$type"_DELS.bed

# Then MMPIDs outside the SEs and in the CNAs.
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_ALL_SEs.bed -b "$type"_AMPS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_ALL_SEs_IN_"$type"_AMPS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_"$type"_SEs.bed -b "$type"_AMPS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_"$type"_SEs_IN_"$type"_AMPS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_ALL_SEs.bed -b RECURRENT_"$type"_AMPS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_ALL_SEs_IN_RECURRENT_"$type"_AMPS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_"$type"_SEs.bed -b RECURRENT_"$type"_AMPS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_"$type"_SEs_IN_RECURRENT_"$type"_AMPS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_ALL_SEs.bed -b "$type"_DELS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_ALL_SEs_IN_"$type"_DELS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_"$type"_SEs.bed -b "$type"_DELS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_"$type"_SEs_IN_"$type"_DELS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_ALL_SEs.bed -b RECURRENT_"$type"_DELS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_ALL_SEs_IN_RECURRENT_"$type"_DELS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIREPOS_OUTSIDE_"$type"_SEs.bed -b RECURRENT_"$type"_DELS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_OUTSIDE_"$type"_SEs_IN_RECURRENT_"$type"_DELS_CNA.bed

# Then MMPIDs in the amps/dels.
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b "$type"_AMPS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_IN_"$type"_AMPS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b RECURRENT_"$type"_AMPS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_IN_RECURRENT_"$type"_AMPS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b "$type"_DELS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_IN_"$type"_DELS.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b RECURRENT_"$type"_DELS_MERGED_ANNOT.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_IN_RECURRENT_"$type"_DELS.bed

# Then MMPIDs and in the CNAs.
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b "$type"_AMPS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_IN_"$type"_AMPS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b RECURRENT_"$type"_AMPS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_IN_RECURRENT_"$type"_AMPS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b "$type"_DELS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_IN_"$type"_DELS_CNA.bed
bedtools intersect -wa -wb -a ./INTERSECTS/FLDL_CCCB_ONLY/SEs_MMPIDs/MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b RECURRENT_"$type"_DELS_MERGED_ANNOT_CNA.bed > ./INTERSECTS/FLDL_CCCB_ONLY/MMPIDs_IN_RECURRENT_"$type"_DELS_CNA.bed
```

#### 2.) Get overlapping genes in CNVs/CNAs.
Now we grab the genes in each CNV/CNA in each of these files. You can use whatever annotations you want for this, I used only the genes in Gencode v19 (includes LINCs). The `-size` option allows you to add wings to your input file to get genes within a certain range of a feature, but I only want those that directly overlap. Could also go back and run this script on the original CNV/CNA files for each cell type if wanted. Be sure to set the chromosome, start, and end column positions properly. 

```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

# MMPID intersects folder.
for file in *.bed; do
  python /scratch/jandrews/bin/get_genes_in_range.py -size 0 -i "$file" -C 4 -S 5 -E 6 -g /scratch/jandrews/Ref/gencode.v19.annotation_sorted_genes_only.bed -o ${file%.*}_GENES.bed
done

# SE intersects folders.
for file in *.bed; do
  python /scratch/jandrews/bin/get_genes_in_range.py -size 0 -i "$file" -C 5 -S 6 -E 7 -g /scratch/jandrews/Ref/gencode.v19.annotation_sorted_genes_only.bed -o ${file%.*}_GENES.bed
done

# Original CNV/CNA files.
for file in *.bed; do
  python /scratch/jandrews/bin/get_genes_in_range.py -size 0 -i "$file" -g /scratch/jandrews/Ref/gencode.v19.annotation_sorted_genes_only.bed -o ${file%.*}_GENES.bed
done
```

#### 3.) Get SE target gene(s).
Now we can also get the target genes for the SEs using the same script. This time I add 100 kb wings to the start/stop of the SE though, as a SE's targets genes may be quite a ways away. Will then move the SE_Genes column to immediately after the SE positions within the file.

```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for file in *GENES.bed; do
  python /scratch/jandrews/bin/get_genes_in_range.py -size 100000 -i "$file" -C 1 -S 2 -E 3 -g /scratch/jandrews/Ref/gencode.v19.annotation_sorted_genes_only.bed -o ${file%.*}_SE_GENES.bed
  awk -F"\t" -v OFS="\t" '{print $1, $2, $3, $4, $11, $5, $6, $7, $8, $9, $10}' ${file%.*}_SE_GENES.bed > ${file%.*}.100KB_SE_GENES.bed
  rm ${file%.*}_SE_GENES.bed
done
```

#### 4.) Rank CNVs containing SEs by recurrence.
This will rank the amps/dels containing SEs by recurrence of the CNV. 

```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

type=CLL

python3 /scratch/jandrews/bin/rank_CNV_by_recurrence.py ALL_SEs_IN_"$type"_AMPS_GENES.100KB_SE_GENES.bed RANKED_BY_CNV_RECUR_ALL_SEs_IN_"$type"_AMPS_GENES.100KB_SE_GENES.bed 
python3 /scratch/jandrews/bin/rank_CNV_by_recurrence.py ALL_SEs_IN_"$type"_DELS_GENES.100KB_SE_GENES.bed RANKED_BY_CNV_RECUR_ALL_SEs_IN_"$type"_DELS_GENES.100KB_SE_GENES.bed 

python3 /scratch/jandrews/bin/rank_CNV_by_recurrence.py "$type"_SEs_IN_"$type"_AMPS_GENES.100KB_SE_GENES.bed RANKED_BY_CNV_RECUR_"$type"_SEs_IN_"$type"_AMPS_GENES.100KB_SE_GENES.bed 
python3 /scratch/jandrews/bin/rank_CNV_by_recurrence.py "$type"_SEs_IN_"$type"_DELS_GENES.100KB_SE_GENES.bed RANKED_BY_CNV_RECUR_"$type"_SEs_IN_"$type"_DELS_GENES.100KB_SE_GENES.bed 
```

At this point, you should have more or less whatever you need to do whatever you want. I like to compare the SE genes to a list of important B cell genes to see if anything really pops out. 

#### 5.) Compare SE genes to gene list.  
The gene list is just a text file with a gene symbol on each line.  

**Python script (pull_interesting_genes.py)**
```Bash
"""Given a gene list and an input file, prints lines from the input file that have genes found in the gene list in the 
specified gene column. Gene list should have one gene symbol per line.

Usage: python3 pull_interesting_genes.py -i <input.bed> -o <output.bed>

Args:
    -i input.bed = Name of input file to process.
    -o output.bed = Name of output file.
    -g gene_list.txt = Name of gene list file.
    -c Gene column = Column in which gene(s) reside in input file. 
    -d Delimiter = Delimiter of the gene column in the input file if they are lists. Should be quoted (e.g. ";", not just ;).
"""
```

**Actual use:**
```Bash
python /scratch/jandrews/bin/pull_interesting_genes.py \
-i RANKED_BY_CNV_RECUR_ALL_SEs_IN_DL_DELS_GENES.100KB_SE_GENES.bed \
-o RANKED_BY_CNV_RECUR_ALL_SEs_IN_DL_DELS_GENES.100KB_SE_GENES.GC_B_GENES.bed \
-g /scratch/jandrews/Ref/GC_B_CELL_GENES.txt \
-c 5 \
-d ";"
```

---

## Integrating Circuit Table Data for MMPIDs
As the circuit tables for FL and DL contain fold-change values for FAIRE, H3AC, and K27AC for each sample at each MMPID relative to the average of the CCs, they can be useful for us to determine which MMPIDs are actually affected by the CNVs. **This section doesn't use the merged CNVs.**

#### 1.) Process the circuit tables.
The circuit tables need a bit of work before they can be used directly (though I left these lying around somewhere after parsing them, so searching around might save you some trouble). The **DL table** erroneously has a sample called DL140, which is actually a CLL sample and needs to be removed. In addition, the MMPIDs are in the first column, but we want a psuedo-bed format for intersecting so we move them to the 4th column.

The **FL table** has some extraneous columns we don't need as well as having the MMPID in the first column.

```Bash
# DL table.
sed '/DL140/d' DL_circuits_May13_2014.txt | awk -v OFS='\t' '{print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10, $11}' > DL_circuits_May13_2014.processed.txt

# FL table.
cut -f1-4,10-16 FL_circuits_May2014.txt | awk -v OFS='\t' '{print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10, $11}' > FL_circuits_May2014.final.txt
```

There was also a mixup with the expression of one sample - DL3A538, as it's expression was listed for DL3B538 (not an actual sample) on a different line. As such, this had to be remedied with a one-off script that grabbed the expression value for each gene for DL3B538 and stuck it in the line for DL3A538.  
**Python script (fix_dl_circuit_expression.py):**

```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python /scratch/jandrews/bin/fix_dl_circuit_expression.py DL_circuits_May13_2014.processed.txt DL_circuits_May13_2014.final.txt
```

#### 2.) Condense the AMPS/DELS.
This will fix lines being repeated if they overlapped multiple common CNV annotations (line for gain, another for loss, etc, will all be condensed into one).  
**Python script (condense_cnv_annot_by_samp.py):**

```Bash
python /scratch/jandrews/bin/condense_cnv_annot_by_samp.py -i DL_AMPS_ANNOT_GENES.bed -o DL_AMPS_ANNOT_GENES_CONDENSED.bed
```

#### 3.) Intersect with AMPS/DELS tables for the appropriate cell type.
Copy and paste the header from the circuit table somewhere and then remove it or bedtools will fight you. We'll have to add headers after this anyway. 

```Bash
# To remove header after copy and pasting it into a text editor.
sed -i -e "1d" DL_circuits_May13_2014.final.txt

module load bedtools2
bedtools intersect -wa -wb -a DL_circuits_May13_2014.final.bed -b ../DL_CNVs/DL_AMPS_ANNOT_GENES.bed > DL_CIRCUIT_MMPIDs_IN_DL_AMPS_MERGED_ANNOT_GENES.bed
bedtools intersect -wa -wb -a DL_circuits_May13_2014.final.bed -b ../DL_CNVs/DL_DELS_ANNOT_GENES.bed > DL_CIRCUIT_MMPIDs_IN_DL_DELS_MERGED_ANNOT_GENES.bed

bedtools intersect -wa -wb -a FL_circuits_May2014.final.bed -b ../FL_CNVs/FL_AMPS_ANNOT_GENES.bed > FL_CIRCUIT_MMPIDs_IN_FL_AMPS_MERGED_ANNOT_GENES.bed
bedtools intersect -wa -wb -a FL_circuits_May2014.final.bed -b ../FL_CNVs/FL_DELS_ANNOT_GENES.bed > FL_CIRCUIT_MMPIDs_IN_FL_DELS_MERGED_ANNOT_GENES.bed
module remove bedtools2
```

#### 4.) Add headers.
Files are quite complicated at this point, so now we add headers since we should be done intersecting.

```Bash
{ printf 'MMPID_CHR\tMMPID_ST\tMMPID_END\tMMPID\tSAMPLE\tFC_X\tFC_H3AC\tFC_K27AC\tABS_K4\tGENE\tEXPRESSION\tCNV_CHR\tCNV_ST\tCNV_END\tCNV_SAMPLES\tANNOT\tCNV_GENES\n'; cat DL_CIRCUIT_MMPIDs_IN_DL_AMPS_ANNOT_GENES_CONDENSED.bed; } > DL_CIRCUIT_MMPIDs_IN_DL_AMPS_ANNOT_GENES_CONDENSED_HEADER.bed

# File names are long enough as is. Try to keep them short.
mv DL_CIRCUIT_MMPIDs_IN_DL_AMPS_ANNOT_GENES_CONDENSED_HEADER.bed DL_CIRCUIT_MMPIDs_IN_DL_AMPS_ANNOT_GENES_CONDENSED.bed
mv DL_CIRCUIT_MMPIDs_IN_DL_DELS_ANNOT_GENES_CONDENSED_HEADER.bed DL_CIRCUIT_MMPIDs_IN_DL_DELS_ANNOT_GENES_CONDENSED.bed
```

#### 5.) Parse for potentially interesting MMPIDs.
This script matches samples between the MMPID and CNVs and then checks if the sample meets the FC cutoffs for either H3AC or K27AC at the MMPID that we'd expect with a deletion or amplification as specified. With the `-r` option, it also only prints MMPIDs that meet the cutoffs recurrently, i.e., in at least two samples. It removes lines with no expression values for the gene.

```
Tries to match MMPIDs that hit the FC cutoff specified in the circuit table for a given sample to a CNV occurring in that sample as well.

Usage: python match_FC_to_CNVs.py -t <amp or del> -i <input.bed> -o <output.bed> [OPTIONS]

Args:
	(required) -t <amp or del> = Type of CNVs that are in the file. 
 	(required) -i <input.bed> = Name of locus list file to process.
 	(required) -o <output.bed> = Name of output file to be created.
 	(optional) -cut <value> = Linear FC magnitude cutoff for H3AC/K27AC used to filter out uninteresting MMPIDs (default=2). Results must meet cutoff to be included.
 	(optional) -r = If included, only includes MMPIDs that are recurrent after applying the cutoff filters and such.
```

**Actual use:**
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python /scratch/jandrews/bin/match_FC_to_CNVs.py -t amp \
-i DL_CIRCUIT_MMPIDs_IN_DL_AMPS_ANNOT_GENES_CONDENSED.bed \
-o DL_CIRCUIT_MMPIDs_IN_DL_AMPS_ANNOT_GENES_CONDENSED_2FC.bed \
-cut 2 \
-r

python /scratch/jandrews/bin/match_FC_to_CNVs.py -t amp \
-i DL_CIRCUIT_MMPIDs_IN_DL_DELS_ANNOT_GENES_CONDENSED.bed \
-o DL_CIRCUIT_MMPIDs_IN_DL_DELS_ANNOT_GENES_CONDENSED_2FC.bed \
-cut 2 \
-r
```

#### 6.) Remove non-FAIRE positive MMPIDs.
Pretty much any enhancer has a FAIRE peak in it, so we can remove those that don't. Have to remove the header and add it back after.

```Bash
module load bedtools2

type=DL
# Remove header.
sed -i -e "1d" RECUR_"$type"_CIRCUIT_MMPIDs_IN_"$type"_AMPS_ANNOT_GENES_CONDENSED_2FC.bed
sed -i -e "1d" RECUR_"$type"_CIRCUIT_MMPIDs_IN_"$type"_DELS_ANNOT_GENES_CONDENSED_2FC.bed

# Intersect
bedtools intersect -wa -a RECUR_"$type"_CIRCUIT_MMPIDs_IN_"$type"_AMPS_ANNOT_GENES_2FC.bed -b MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed > RECUR_"$type"_CIRCUIT_FAIRE_POS_MMPIDs_IN_"$type"_AMPS_ANNOT_GENES_CONDENSED_2FC.bed
bedtools intersect -wa -a RECUR_"$type"_CIRCUIT_MMPIDs_IN_"$type"_DELS_ANNOT_GENES_2FC.bed -b MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed > RECUR_"$type"_CIRCUIT_FAIRE_POS_MMPIDs_IN_"$type"_DELS_ANNOT_GENES_CONDENSED_2FC.bed

# Add back header.
{ printf 'MMPID_CHR\tMMPID_ST\tMMPID_END\tMMPID\tSAMPLE\tFC_X\tFC_H3AC\tFC_K27AC\tABS_K4\tGENE\tEXPRESSION\tCNV_CHR\tCNV_ST\tCNV_END\tCNV_SAMPLES\tANNOT\tCNV_GENES\n'; cat RECUR_"$type"_CIRCUIT_FAIRE_POS_MMPIDs_IN_"$type"_AMPS_ANNOT_GENES_CONDENSED_2FC.bed; } > RECUR_"$type"_CIRCUIT_FAIRE_POS_MMPIDs_IN_"$type"_AMPS_ANNOT_GENES_CONDENSED_2FC_HEADER.bed
{ printf 'MMPID_CHR\tMMPID_ST\tMMPID_END\tMMPID\tSAMPLE\tFC_X\tFC_H3AC\tFC_K27AC\tABS_K4\tGENE\tEXPRESSION\tCNV_CHR\tCNV_ST\tCNV_END\tCNV_SAMPLES\tANNOT\tCNV_GENES\n'; cat RECUR_"$type"_CIRCUIT_FAIRE_POS_MMPIDs_IN_"$type"_DELS_ANNOT_GENES_CONDENSED_2FC.bed; } > RECUR_"$type"_CIRCUIT_FAIRE_POS_MMPIDs_IN_"$type"_DELS_ANNOT_GENES_CONDENSED_2FC_HEADER.bed

# File names are long enough as is. Try to keep them short.
mv RECUR_"$type"_CIRCUIT_FAIRE_POS_MMPIDs_IN_"$type"_AMPS_ANNOT_GENES_CONDENSED_2FC_HEADER.bed RECUR_"$type"_CIRCUIT_FAIRE_POS_MMPIDs_IN_"$type"_AMPS_ANNOT_GENES_CONDENSED_2FC.bed
mv RECUR_"$type"_CIRCUIT_FAIRE_POS_MMPIDs_IN_"$type"_DELS_ANNOT_GENES_CONDENSED_2FC_HEADER.bed RECUR_"$type"_CIRCUIT_FAIRE_POS_MMPIDs_IN_"$type"_DELS_ANNOT_GENES_CONDENSED_2FC.bed
```

---

## Observe SE Signals Inside and Outside CNVs 
This section breaks the amp/del lists up *by sample* and then uses the signal from the B cell SEs to try to show an **increase** in SE signal in amps and a **decrease** in dels. It uses output from the **[CNV-calling section](#cnv-calling)** and also requires that the [SE Pipeline](https://github.com/j-andrews7/Pipelines/blob/master/ROSE_SE_Pipeline.md) has been completed.

#### 1.) Break the CNVs up by sample.
We can use the files already generated when looking at the CNVs on a cell-type basis. More specifically, we want the files containing the amps/dels for each sample **without** merging, but with annotations and genes *already added*. This script will create a file for each sample within the CNV list and stick the CNVs for that sample in the file. These files are also useful for looking at [lincRNA expression changes](name=#comparing-lincrna-expression-in-cnvs).

**Python script (get_cnvs_by_sample.py)**
```Bash
"""
Given a list of unmerged amps or dels with samples in the 4th column, prints the cnvs for every sample in the file to an 
individual file.

Usage: python3 get_cnvs_by_sample.py -i <input.bed> -o <output suffix>

Args:
    -i input.bed = Name of input file to process.
    -o output suffix = Suffix to add to append to the sample names to be used as file names.
"""
```

**Actual use:**
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python /scratch/jandrews/bin/get_cnvs_by_sample.py -i DL_AMPS_ANNOT_GENES_CONDENSED.bed -o _AMPS_ANNOT_GENES_CONDENSED.bed
python /scratch/jandrews/bin/get_cnvs_by_sample.py -i DL_DELS_ANNOT_GENES_CONDENSED.bed -o _DELS_ANNOT_GENES_CONDENSED.bed

python /scratch/jandrews/bin/get_cnvs_by_sample.py -i FL_AMPS_ANNOT_GENES_CONDENSED.bed -o _AMPS_ANNOT_GENES_CONDENSED.bed
python /scratch/jandrews/bin/get_cnvs_by_sample.py -i FL_DELS_ANNOT_GENES_CONDENSED.bed -o _DELS_ANNOT_GENES_CONDENSED.bed

python /scratch/jandrews/bin/get_cnvs_by_sample.py -i CLL_AMPS_ANNOT_GENES_CONDENSED.bed -o _AMPS_ANNOT_GENES_CONDENSED.bed
python /scratch/jandrews/bin/get_cnvs_by_sample.py -i CLL_DELS_ANNOT_GENES_CONDENSED.bed -o _DELS_ANNOT_GENES_CONDENSED.bed
```

#### 2.) Get signal for SEs in/outside CNVs.
This script intersects the SEs with the amps and dels for a given sample. The signal for the SEs found **in** the amp/del will be output to one file, while those found **outside** the amp/del will be output to another file. It also finds the SEs that are "unchanged" in a given sample (i.e., **not found in either the amps or dels**). It's kind of lazily coded, so it *assumes* the sample name will be the first thing in the input file name (e.g., <sample>_moreinfo). 

**Python script (get_sample_se_cnv_loads.py):**
```Bash
"""
Given lists of unmerged amps and dels for a sample, gets the signal for all SEs in that sample that lie within the amps/dels, 
outside them, and for those that are unchanged and spits this info out to multiple files.

Usage: python3 get_sample_se_cnv_loads.py <amps.bed> <dels.bed> <SE_signal.bed> <Overlap percentage as decimal>

Args:
    amps.bed = Name of amps file to process.
    dels.bed = Name of dels file to process.
    SE_signal.bed = Name of SE signal file.
    overlap_percentage = Percentage as decimal for determining if an SE should be considered "overlapping" a CNV or not.
"""
```

I cheat a bit here and bank on my files being named with the sample first and the directory being otherwise empty for anything that would match this pattern.

**Actual use:**
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for file in *AMPS_ANNOT*; do
	samp="$(echo "$file" | cut -d'_' -f1)"
	python /scratch/jandrews/bin/get_sample_se_cnv_loads.py "$samp"* FLDL_CCCB_ONLY_SES_SIGNAL.bed 0.25
done
```

#### 3.) Copy data into table.
Use excel (or write a script, I'm a guideline, not a cop), to get all of the signals into a format like so for all the comparisons you'd like to see:

| DL135      |             | DL188      |             |
|------------|-------------|------------|-------------|
| In amps    | Not in amps | In amps    | Not in amps |
| 3911.606   | 2203.344    | 11907.221  | 8764.877    |
| 16398.0495 | 2948.2706   | 8582.6538  | 3584.8002   |
| 12176.5078 | 2866.2588   | 7796.7622  | 881.8416    |

The columns will likely not be the same length. Can also make other tables like this, like "In amps" vs "In dels", etc.

#### 4.) Create box plots.
I used GraphPad Prism for this since it looks good, is pretty easy to use, and can do any stats you may want. I did unpaired, two-tailed, Welch's t-tests and got decent results. 

---

## Comparing lincRNA Expression in CNVs

This analysis is similar to the one for SEs above, but looking at lincRNA expression rather than SE signal. It does the analysis two ways, by comparing to the **median** expression across samples for each linc and also comparing to the **average** expression.

#### 1.) Cut unwanted columns and rows.
The original table has VGA/CC samples, but I only want to normalize to samples for which we also have CN data (so just the tumors). Also removes the sex chromosomes.

```Bash
cut -f4,7,8,26-30 --complement GENCODE_NOVEL_LINC_FPKMS_2SAMPS_OVER1.txt | sed '/chrX\|chrY\|chr23\|_g/d' - > GENCODE_NOVEL_LINC_FPKMS_2SAMPS_OVER1_CUT.txt
```

#### 2.) Calculate linc expression FCs for each sample in CNVs.
The mash everything together and get stats script. It calculates the log2 FC for each linc in each sample compared to the median and average of the linc FPKMs for all samples in the expression file. Spits out these values for all lincs inside/out the CNVs for each sample. You set the percentage that the linc must overlap the CNV to be considered "in" it. I usually try a few different percentages, but anything >= 25% will get rid on those fringe cases where the linc is barely overlapping the CNV.

**Python script (get_sample_linc_cnv_loads.py):**
```Bash
"""Given lists of unmerged amps and dels for a sample, gets the expression for all lincs in that sample that lie within the amps/dels, outside them, and for those that are unchanged and spits this info out to multiple files. The linc expression is log2 FC to the median/average of the other samples in the expression file. Also yields a few summary stats that may be helpful.

Linc RNA expression file should contain a header and the format should be:
CHR	START	STOP	GENE_ID	GENE_SHORT_NAME	<SAMPLE_FPKM_COLUMNS>

Usage: python3 get_sample_linc_cnv_loads.py <amps.bed> <dels.bed> <linc_expression.txt> <overlap percentage as decimal>

Args:
    amps.bed = Name of amps file to process.
    dels.bed = Name of dels file to process.
    linc_expression.txt = Name of linc_expression file.
    overlap_percentage = Percent of linc that must overlap the CNV for it to be considered 'overlapping'.
"""
```

I cheat a bit here and bank on my files being named with the sample first and the directory being otherwise empty for anything that would match this pattern.

**Actual use:**
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for file in *AMPS_ANNOT*; do
	samp="$(echo "$file" | cut -d'_' -f1)"
	python /scratch/jandrews/bin/get_sample_linc_cnv_loads.py "$samp"* GENCODE_NOVEL_LINC_FPKMS_2SAMPS_OVER1_CUT.txt 0.25
done
```

#### 3.) Copy data into table.
Use excel (or write a script, I'm a guideline, not a cop), to get all of the signals into a format like so for all the comparisons you'd like to see:

| DL135   |             |             |             |
|---------|-------------|-------------|-------------|
| In amps |  In amps    | Not in amps | Not in amps | 
| Median  | Average     | Median      |  Average    |
| 1.0785  |  0.8923     | 0.2133      | -0.2176     |
| 0.1346  | -0.2342     | -0.5646     | -0.3435     |

The columns will likely not be the same length. Can also make other tables like this, like "In amps" vs "In dels", etc.

#### 4.) Visualize.
Again, box plots in prism are likely the best way to look at this data, and it allows for easy t-tests, etc.

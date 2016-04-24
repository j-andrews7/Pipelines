# CN Analysis 

**Up to date as of 04/22/2016.**  
jared.andrews07@gmail.com

> This version differs from the "old" version because it **treats the CNVs on a sample-by-sample basis for the MMPID/Circuit analysis** in order to keep them as small as possible (i.e. it doesn't merge them together for samples of a given cell type). It also tries to find the **"golden ticket"** CNVs by searching for minimal recurrent regions. Also gets into creating some figures to actually try to show the effects of the CNVs on the activity of SEs and MMPIDs.

The aim of this pipeline is to get all copy number changes for all samples for which we have SNP arrays. These are then intersected with CNAs identified in other publications or with our SE and MMPID data. As with most things, this started off relatively simple and grew to become more complicated as results were viewed and additional approaches tried.

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
- [bedtools](http://bedtools.readthedocs.org/en/latest/)
  - Also available on the CHPC cluster.
- [Affymetrix Genotyping Console](http://www.affymetrix.com/estore/browse/level_seven_software_products_only.jsp?productId=131535#1_1)

#### Sections  
- [Calling CNVs from SNP6 Arrays](name=#segment)
- [Integrating SEs with CNVs by Cell Type](name=#SEsCellType)
- [Integrating Circuit Table Data to Filter MMPIDs in CNVs](name=#MMPIDS)
- [Creating Boxplots of SE Signal on a Sample-by-Sample CNV Basis](name=#SEBoxPlots)

---

## <a name="segment"></a>CNV Calling
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

## <a name="SEsCellType"></a>Integrate SE Data for CNVs by Cell Type
This utilizes the files generated in the above section to look at the SEs located in CNVs on a cell-type basis.

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

## <a name="MMPIDS"></a> Integrating Circuit Table Data for MMPIDs in CNVs
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

## <a name="SEBoxPlots"></a> Observe SE Signals Inside/Outside CNVs on a Sample-by-Sample Basis
This section breaks the amp/del lists up by sample and then uses the signal from the B cell SEs to try to show an **increase** in SE signal in amps and a **decrease** in dels. It uses output from the **first section** and also requires that the [SE Pipeline](https://github.com/j-andrews7/Pipelines/blob/master/ROSE_SE_Pipeline.md) has been completed.

#### 1.) Break the CNVs up by sample.
We can use the files already generated when looking at the CNVs on a cell-type basis. More specifically, we want the files containing the amps/dels for each sample **without** merging, but with annotations and genes *already added*. This script will create a file for each sample within the CNV list and stick the CNVs for that sample in the file. These files are also useful for looking at [lincRNA expression changes](name=#) and 

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

**Python script (compare_SE_signal_cnvs.py):**
```Bash
"""
Given lists of unmerged amps and dels for a sample, gets the signal for all SEs in that sample that lie within the amps/dels, 
outside them, and for those that are unchanged and spits this info out to multiple files.

Usage: python3 compare_SE_signal_cnvs.py <amps.bed> <dels.bed> <SE_signal.bed>

Args:
    amps.bed = Name of amps file to process.
    dels.bed = Name of dels file to process.
    SE_signal.bed = Name of SE signal file.
"""
```

I cheat a bit here and bank on my files being named with the sample first and the directory being otherwise empty for anything that would match this pattern.

**Actual use:**
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for file in *AMPS*; do
	samp="$(echo "$file" | cut -d'_' -f1)"
	python /scratch/jandrews/bin/compare_SE_signal_cnvs.py "$samp"* FLDL_CCCB_ONLY_SES_SIGNAL.bed
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

## <a name="lincs"></a> Comparing lincRNA Expression in CNVs

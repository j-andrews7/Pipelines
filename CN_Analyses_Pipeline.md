# CN Analysis 
The aim of this pipeline is to get all copy number changes for all samples for which we have SNP arrays. These can then be intersected with CNAs identified in other publications or with our SE and MMPID data.

**Up to date as of 03/29/2016.**  
jared.andrews07@gmail.com

---

##### 1.) Download the [Affymetrix Genotyping Console](http://www.affymetrix.com/estore/browse/level_seven_software_products_only.jsp?productId=131535#1_1) program and the na32 annotation db. 
You'll have to set a library folder when opening the program for the first time - this is where you should stick annotation/databse files. Download the na32 library from within the program and stick it in your library folder, along with the na32 CN annotation file from Affy's site. Again, this is all for the **SNP 6 arrays**, and I used the older annotations (na32 instead of na35) for continuity with other analyses. 


##### 2.) Take all cel & arr files you want to use and stick them in a folder. 
Load the data into the program. Run the Copy Number/LOH analysis, it'll generate a CNCHP file for each sample. Go to Workspace->Copy Number/LOH Analysis and export them, making sure to include the Log2 ratio for each marker (along with it's chromosome), marker ID, and position.


##### 3.) Run the Segment Reporting Tool. 
Be sure to click the option to generate a summary file, which will have the segments for each sample. 


##### 4.) Scrub the summary file.
The first column of the summary file will have the array ID, which should be replaced with the sample name. In addition, the `#` header lines can be removed and the only columns that need to be kept are those for `sample, Chr, Start_Linear_Position, End_Linear_Position, #Markers, Copy_Number_State`. `awk/sed/cut/paste` make swapping the columns around pretty easy.

##### 5.) Convert to bed-like format.
Just need to order columns - `Chr, Start, End, Sample, #Markers, Copy_Number_State` and add `chr` to first column. This also removes segments with fewer than 5 markers in them if they haven't been removed already.

```Bash
awk -v OFS='\t' '{print $2, $3, $4, $1, $5, $6}' Lymphoma_SNP77_CNCHP_011615_segment_summary_cleaned.txt \
| awk '{print "chr" $0}' \
| awk '{
if ($5 >=5) 
print $0
}' > Lymphoma_SNP77_CNCHP_011615_segment_summary_cleaned.bed
```

##### 6.) Grab segments specific to each cell-type.
Go ahead and delete those on the Y chromosome as well.
```Bash
grep DL Lymphoma_SNP77_CNCHP_011615_segment_summary_cleaned.bed | sed '/chrY/d' > DL_CNVs.bed
grep FL Lymphoma_SNP77_CNCHP_011615_segment_summary_cleaned.bed | sed '/chrY/d' > FL_CNVs.bed
grep CLL Lymphoma_SNP77_CNCHP_011615_segment_summary_cleaned.bed | sed '/chrY/d' > CLL_CNVs.bed
```

##### 7.) Break each CNV file into amps/dels
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

##### 8.) Merge the amps/dels for each cell-type.
This also creates a list of the samples in which each CNV is found as well as the actual copy number of the regions merged, which can *sometimes* be helpful for determining if there are any high-level samples (multiple copy gain/loss) for that CNV, though you can't determine which sample it's necessarily in without going back and looking in each. It will pull the max marker numbers for the merged region as well.

```Bash
module load bedtools2
mergeBed -c 4,5,6 -o distinct,max,collapse -i CLL_AMPS.bed > CLL_AMPS_MERGED.bed
```

##### 9.) Annotate known CNVs.
While CNAs unique to the samples are likely to be most interesting, a CNV found in other cells and cancers may still have a functional impact. As such, we **don't remove common CNVs**, rather we simply make a note that the CNV is found elsewhere. There are many lists of common CNVs (HapMap, DGV, and a [2015 Nature Genetics Review](http://www.nature.com/nrg/journal/v16/n3/full/nrg3871.html) among them). Here, I use the [DGV gold standard list](http://dgv.tcag.ca/dgv/app/downloads?ref=).


```Bash
module load bedtools2
bedtools intersect -wa -wb -a CLL_AMPS_MERGED.bed -b DGV.GoldStandard.July2015.hg19.gff3 > CLL_AMPS_MERGED_ANNOT.bed
```

TO-DO - ORDERING COLUMNS, REMOVING UNWANTED SAMPLES, BREAKING INTO AMP/DEL LISTS AND MERGING


6.) Filter segmentation file, only keeping segments with at least 5 markers present in them. Can be adjusted if you want to be more stringent.
awk -F'\t' -v OFS='\t' {(if $5 >= 5) print $1, $2, $3, $4, $5, $6}' Segs_For_GISTIC.txt > Segs_For_Gistic_Min5.txt


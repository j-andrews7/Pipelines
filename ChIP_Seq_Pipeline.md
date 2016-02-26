ChIP-SEQ Pipeline

This document describes the bioinformatics pipeline used to analyze the Payton Lab's normal B cell ChIP-Seq samples for a variety of chromatin marks (mostly K27AC and k4ME3). The .bin and peaks.bed files can be handled separately until the last portion of the pipeline in which they are normalized and merged. Some of these commands are not strictly necessary (viewing bigwig files on UCSC, fastq-dump, etc), but may be useful, particularly for handling/verifying external data. Additional file manipulations may be necessary (removal of headers, switching columns around, etc), though considerable effort has been made to minimize this as much as possible. Learn sed & awk to do these. This is not the end-all, be-all, but it should be a good place to start and may warn you of possible caveats ahead of time.


####-File Conversions-####
--Download and convert SRA file to fastq file
fastq-dump <SRA accession>

--Convert fastq file to bam file
bowtie2 -p 3 -x <path to genome index files and prefix> <fastq file.fastq> > <output file.sam>

--Convert sam to bam
samtools view -bS <file.sam> > <file.bam>

--Sort bam files
samtools sort <file.bam> <file_sorted>

--Index bam files
samtools index <test.bam> 

--Pipeline from SRA to BAM combined commands
fastq-dump <SRA accession> && bowtie2 -p3 -x <path to genome index files and prefix> <fastqfile.fastq> > <output.sam> && samtools view -bS <file.sam> > <file.bam> && samtools sort <file.bam> <file> && samtools index <file.bam>


***IMPORTANT***
If the data is from the SOLiD platform, it is garbage and you will not be able to get squat from it. Don't bother even trying.



####-Peak Calling-####

###MACS
--For other marks - calls peaks
macs14 -t <treated file> -c <control file> -f BAM -g <genome prefix - hs, mm> -n NAME -w -S --nomodel --shiftsize=150

--For K4ME3 - calls peaks
macs14 -t <treated file> -c <control file> -f BAM -g <genome prefix - hs, mm> -n NAME -w -S --nomodel --shiftsize=100

--For FAIRE - calls peaks
macs14 -t <treated file> -c <control file> -f BAM -g hs -n NAME -w -S --nomodel --shiftsize=50


###HOMER
#Make tag directories from all alignment files. Must have samtools installed if using .bam format.
batchMakeTagDirectory.pl <keyFile> -genome hg19 [options]

# make UCSC files for all directories
batchParallel.pl makeUCSCfile none -o auto -f *TagDir/

# find peaks for all directories
batchParallel.pl findPeaks none -o auto -style histone -size 1000 -minDist 2500 -fdr 0.005 -f *TagDir/
-Can alter minDist (size required between peaks) and size (of peaks) in order to get data in a form that you want
-These parameters give pretty good peak overlap with the MACS settings above.

#Normalization is a problem with HOMER. Can use it for peak annotations, etc, but not recommended for actual peak calling.



##-ROSE for SE Calling (CHPC)-##
See ROSE_SE_pipeline.txt.



####-Wig to Bigwig & Reverse for Visualization of Peaks on genome browser-####
#Done to spot check data quality

--Create genome.chrom.sizes
fetchChromSizes <db (hg19, etc)> > <db>.chrom.sizes

--Convert wig to bigwig
wigToBigWig <inputwig.gz> <genome.chrom.sizes> <myBigWig.bw>

--Convert bigwig to wig
bigWigToWig in.bigWig out.wig

--Trackline for typical bigwig file loading to UCSC
track name=<name> type=bigWig bigDataUrl=http://pathresearch.wustl.edu/paytonlab/H3K4ME3/TsB/filename.bw

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
NOTE: Must edit the hash in the beginning of this script to add your sample names and the number of aligned reads in each. May also have to adjust 
the column in which data begins for the file within the script. 




####-MACS Bed Files Processing-####


##Plan to make a script to do the enclosed manual steps in one easy step eventually
--Rename 4th column in MACS14_peak.bed files to contain the sample name and mark. This is for those with additional info in the filename (TS081414_NAIVE_K27AC, for example.)
python3 rename_columns.py <MACS14_peak.bed files>
NOTE: filename format MUST be samplename_celltype,etc_histonemark... .bed
May just want to edit these as necessary

--Rename 4th column in MACS14_peak.bed files to contain the sample name and mark. 
python3 rename_columns_second.py <MACS14_peak.bed files>
NOTE: filename format MUST be samplename_histonemark... .bed
May just want to edit these as necessary

--Concatenate all bed files in a directory
cat *.bed > output.bed

--delete lines containing _g random genomic, M chromosome, Y chromosome loci
(sed '/_g/d' file.bed | sed '/chrM/d' | sed '/chrY/d') > output.bed

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

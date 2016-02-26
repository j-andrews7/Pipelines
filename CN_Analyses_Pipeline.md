CN Analysis - Aim is to detect focal amplifications (or deletions) with GISTIC.
Up to date as of 02/08/2106
jared.andrews07@gmail.com

###-To determine amplifications/deletions associated with the tumors, use GISTIC.
Need a segmentation file and markers file to run. Segmentation file needs 6 columns: sample ID, chrom, start, end, # of markers, log2(CN)-1. Markers file needs 3: marker_ID, chrom, position. Be sure to use the CN probes from the SNP6 arrays from Affy. This was actually a huge hassle to figure out, so here's a basic breakdown. If you're running it on new samples at least, otherwise find the markers and segmentation files I made for this previously. 


1.) Download affy Genotype Console program and the na32 annotation db. Will have to stick it in the library folder you specify when opening program for first time. Will also have to download
the na32 library and stick it in the same place. And the na32 CN annotation file. 


2.) Take all cel & arr files you want to use and stick them in a folder. Point program to it to load the data. Run the Copy Number/LOH analysis, it'll generate a CNCHP file for each sample. 
Go to Workspace->Copy Number/LOH Analysis and export them, making sure to include the Log2 ratio for each marker (along with it's chromosome), marker ID, and position.


3.) Run the Segment Reporting Tool. Be sure to click the option to generate a summary file, which will have the segments for each sample. This takes a bit, so can do next step in meantime.


4.) Copy all of the CNCHP files to an empty directory. 


5.) Create a file containing the array number located in the CNCHP file along with the actual sample name. One sample per line, tab delimited. 
For example, for a file named 18117.CN.CNCHP.txt:
18117	RAJI


5.) Run this script to generate a segmentation file containing all segments found in each sample along with a markers file for all the markers used.
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python3  /scratch/jandrews/bin/prep_for_gistic_new.py <sample_names.txt> <output.txt> <input files - *.CNCHP.txt>


6.) Filter segmentation file, only keeping segments with at least 5 markers present in them. Can be adjusted if you want to be more stringent.
awk -F'\t' -v OFS='\t' {(if $5 >= 5) print $1, $2, $3, $4, $5, $6}' Segs_For_GISTIC.txt > Segs_For_Gistic_Min5.txt


7.) Take output files and run GISTIC from within the GISTIC directory. Can adjust '-maxseg' if you don't want to include files with tons of segments. I increased it to 6000 from default 2500
to include a tumor sample. 

./gp_gistic2_from_seg -seg ~/Analyses/Transfer/GISTIC/Segs_For_Gistic_Min5.txt -mk ~/Analyses/Transfer/GISTIC/Segs_For_GISTIC.txt.markers -refgene ./refgenefiles/hg19.mat -maxseg 6000 -b ~/Analyses/Transfer/GISTIC/RESULTS

# CN Analysis 
The aim of this pipeline is to get all copy number changes for all samples for which we have SNP arrays. These can then be intersected with CNAs identified in other publications or with our SE and MMPID data.

**Up to date as of 03/24/2016.**  
jared.andrews07@gmail.com

---

##### 1.) Download the [Affymetrix Genotyping Console](http://www.affymetrix.com/estore/browse/level_seven_software_products_only.jsp?productId=131535#1_1) program and the na32 annotation db. 
You'll have to set a library folder when opening the program for the first time - this is where you should stick annotation/databse files. Download the na32 library from within the program and stick it in your library folder, along with the na32 CN annotation file from Affy's site. Again, this is all for the **SNP 6 arrays**, and I used the older annotations (na32 instead of na35) for continuity with other analyses. 


##### 2.) Take all cel & arr files you want to use and stick them in a folder. 
Load the data into the program. Run the Copy Number/LOH analysis, it'll generate a CNCHP file for each sample. Go to Workspace->Copy Number/LOH Analysis and export them, making sure to include the Log2 ratio for each marker (along with it's chromosome), marker ID, and position.


##### 3.) Run the Segment Reporting Tool. 
Be sure to click the option to generate a summary file, which will have the segments for each sample. 


##### 4.) Scrub the summary file.
The first column of the summary file will have the array ID, which should be replaced with the sample name.
In addition, the `#` header lines can be removed and the only columns that need to be kept are those for
`sample, Chr, Start_Linear_Position, End_Linear_Position, #Markers, Copy_Number_State`. 'awk/sed/cut/paste' make
swapping the columns around pretty easy.

TO-DO - ORDERING COLUMNS, REMOVING UNWANTED SAMPLES, BREAKING INTO AMP/DEL LISTS AND MERGING


6.) Filter segmentation file, only keeping segments with at least 5 markers present in them. Can be adjusted if you want to be more stringent.
awk -F'\t' -v OFS='\t' {(if $5 >= 5) print $1, $2, $3, $4, $5, $6}' Segs_For_GISTIC.txt > Segs_For_Gistic_Min5.txt


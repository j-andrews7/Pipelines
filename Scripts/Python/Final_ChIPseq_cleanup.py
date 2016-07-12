#!/usr/bin/env python3
"""Cleans up files from wonky R output. Rounds numbers to 4 decimals, fixes chromosomes being converted to random numbers, 
makes sure ending bin of each chrom doesn't extend past actual end of chromosome, etc. Also has option to create bedgrapth files from each data column.

Usage: python3 Final_ChIPseq_cleanup.py <input.txt> <output.txt>

Args:
    input.txt = Name of input file to process 
    output.txt = Name of output file to be created
"""

import sys
from math import round

inputfile = sys.argv[1]
outputfile = sys.argv[2]

chrom_correct = {"1":"1", "12":"2", "16":"3", "17":"4", "18":"5","19":"6", "20":"7", "21":"8", "22":"9", "2":"10",
	"3":"11", "4":"12", "5":"13", "6":"14", "7":"15", "8":"16", "9":"17", "10":"18", "11":"19", "13":"20", "14":"21", "15":"22", "23":"X"}

#Dictionary of the offending bins and the correct chromosome stop position.
bin_correct = {("chr1", "249250800"):"249250621", ("chr2", "243199400"):"243199373", ("chr3","198022600"):"198022430", 
	("chr4", "191154400"):"191154276", ("chr5","180915400"):"180915260",("chr6","171115200"):"171115067", ("chr7", "159138800"):"159138663", 
	("chr8", "146364200"):"146364022", ("chr9", "141213600"):"141213431", ("chr10", "135534800"):"135534747",
	("chr11", "135006600"):"135006516", ("chr12", "133852000"):"133851895", ("chr13", "115170000"):"115169878", ("chr14", "107349600"):"107349540", 
	("chr15", "102531400"):"102531392", ("chr16", "90354800"):"90354753", ("chr17", "81195400"):"81195210", ("chr18", "78077400"):"78077248", 
	("chr19", "59129000"):"59128983", ("chr20", "63025600"):"63025520", 
	("chr21", "48130000"):"48129895", ("chr22", "51304600"):"51304566", ("chrX", "155270600"):"155270560"}

#Open input file
with open(inputfile) as f:

	#Open output file, "w" to make it writable
	output_file = open(outputfile, "w")

	#Get header and print it to the output file
	header = f.readline().strip()
	print(header, file=output_file)

	#Iterate through each line in input file
	for line in f:

		#The float function apparently can't handle the small e for scientific notation
		line = line.replace("e","E")

		#Strip newline character and split line by white space
		line = line.strip().split()

		#Add "chr" back to the chromosome element and correct for R changing the chromosome to a random number. 
		line[0] = "chr" + chrom_correct[line[0]]

		#Create tuple from chromosome and stop position for comparison to the bin_correct dictionary 
		tup = (line[0], str(line[2]))

		#Check if the tup is in the bin_correct dictionary
		if tup in bin_correct:
			line[2] = bin_correct[tup]

		#Convert start, stop, peak_ID of line to ints because R output is sometimes silly and will write chromosomal positions in scientific notation
		for i in line[1:4]:
			i_index = line.index(i)
			line[i_index] = int(float(i))

		for i in line[4:]:
			i_index = line.index(i)
			line[i_index] = round(i[, 4])

		#Print results
		print(*line[0:], sep="\t", file=output_file)

	output_file.close()
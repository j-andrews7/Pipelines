#!/usr/bin/env python3
"""Fixes end of actual chromosome not matching the last bin for bedgraph files of RPM'd, QN'd ChIP-seq samples.

Usage: python3 chrom_bin_fix.py <input.bg> <output.bg>

Args:
    input.bg = Name of input file to process 
    output.bg = Name of output file to be created
"""

import sys

inputfile = sys.argv[1]
outputfile = sys.argv[2]

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

	#Iterate through each line in input file
	for line in f:

		#Strip newline character and split line by white space
		line = line.strip().split()

		#Create tuple from chromosome and stop position for comparison to the bin_correct dictionary 
		tup = (line[0], str(line[2]))

		#Check if the tup is in the bin_correct dictionary
		if tup in bin_correct:
			line[2] = bin_correct[tup]

		#Print line
		result = "\t".join(line)
		print(result,file=output_file)

	#Close the output file
	output_file.close()
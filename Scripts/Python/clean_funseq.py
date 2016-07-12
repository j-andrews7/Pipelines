#!/usr/bin/env python3
"""Remove coding variants and those on X chrom after funseq annotation. Assumes bed output option from funseq.
Also prints just the positions to another output file.

Usage: python3 clean_funseq.py <input.bed> <output.bed>

Args:
    input.txt = Name of input file to process 
    output.txt = Name of output file to be created
"""

import sys

####-Variables-####

#Store file names
input_file = sys.argv[1]
output_name = sys.argv[2]
output_base = output_name.split(".")[0]

#open input file
with open(input_file) as f:

	#Grab header
	header = f.readline().strip()

	#Open output file
	outputfile = open(output_name, "w")

	positions_file = open(output_base+"_positions.bed", "w")

	#Print header to outputfile
	print(header, file=outputfile)

	#Iterate through each line of file
	for line in f:

		#Split and strip line and store as list
		line = line.strip()
		new_line = line.split("\t")

		chrom, start, stop = new_line[0], new_line[1], new_line[2]
		coding = new_line[6].split(";")[1]

		if coding == "No" and chrom != "chrX" and chrom != "chr23":
			print(chrom, start, stop, sep="\t", file=positions_file)
			print(line,file=outputfile)

	#Close output file
	outputfile.close()
	positions_file.close()



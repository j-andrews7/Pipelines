#!/usr/bin/env python3
"""For a .bed file after calling peaks, merging, cutting, and sorting, convert this file to gff format for use with ROSE.

Usage: python3 ROSE_bed2gff.py <input.bed> 

Args:
    input.bed = Name of input file to process 
"""

import sys

#Store file names
input_file = sys.argv[1]
base = input_file.split(".")[0]
output_file = base + ".gff"

with open(input_file) as f:

	#Open output file, "w" to make it writable
	output_f = open(output_file, "w")

	#Initiate variable to hold what next ID should be
	ID = 1

	#Iterate through each line in input file
	for line in f:

		#Strip newline character and split line by white space
		line = line.strip().split()

		#Print to outputfile
		print(line[0],ID,"RE",str(int(line[1]) + 1),line[2],".",".",".",ID,sep="\t", file=output_f)

		#Increase ID by one for next line
		ID += 1

	#Close output file
	output_f.close()

print("COMPLETE")

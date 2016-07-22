#!/usr/bin/env python3
"""Subtract 1 from the start of each line

Usage: python3 correct_start_stop.py <input.txt> <output.txt>

Args:
    input.txt = Name of input file to process 
    output.txt = Name of output file to be created
"""

import sys

####-Variables-####

#Store file names
input_file = sys.argv[1]
output_name = sys.argv[2]

#open input file
with open(input_file) as f:

	#Grab header
	header = f.readline().strip()

	#Open output file
	outputfile = open(output_name, "w")

	#Print header to outputfile
	print(header, file=outputfile)

	#Iterate through each line of file
	for line in f:

		#Split and strip line and store as list
		line = line.strip().split()

		#Subtract 1 from start and stop
		line[1] = int(line[1]) - 1

		print(*line[0:], sep="\t",file=outputfile, end="")
		print("",file=outputfile)

	#Close output file
	outputfile.close()



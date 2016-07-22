#!/usr/bin/env python3
"""For given .txt file(s), adds the chr prefix to each line for each column.

Usage: python3 add_chr_to_firstcol.py <input.txt file(s)> 

Args:
    input.bin = Name of .txt file to process 
    can use *.bin to run on all .txt files in the directory
"""

import sys

#Store file names
input_files = sys.argv[1:]

#Open input file
for item in input_files:
	with open(item) as f:

		output_name = (item.split(".")[0] + "_chr.txt")

		#Open output file, "w" to make it writable
		output_file = open(output_name, "w")

		#print out header
		header = f.readline().strip()
		print(header, file=output_file)

		#Iterate through each line in input file
		for line in f:

			#Strip newline character and get first column
			line = line.strip()
			line = line.split()
			chr_number = ("chr" + str(line[0]))

			#Print to outputfile
			print(chr_number, *line[1:], sep='\t', file=output_file)

		#Close output file
		output_file.close()

print("COMPLETE")
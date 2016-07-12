#!/usr/bin/env python3
"""Chops the long file paths for the sample column headers. Also removes the file extension for them.

Usage: python3 shorten_vcf_header.py <input.vcf> <output.vcf>

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

	#Open output file
	outputfile = open(output_name, "w")

	#Iterate through each line of file
	for line in f:

		if line.startswith("#C"):
			line = line.strip().split("\t")
			samples = line[8:]
			good_samples = []

			# Remove the excess crap for each sample name. Change to upper, remove unnecessary underscores, etc.
			for item in samples:
				item = item.split("/")[-1]
				item = item.split(".")[0]
				item = item.upper()
				if item.count("_") == 2:
					item = item.split("_")
					new = item[0] + item[1] + "_" + item[2]
					good_samples.append(new)
				else:
					good_samples.append(item)

			# Join samples, then rest of line.
			new_samples = "\t".join(good_samples)
			new_line = "\t".join(line[0:8])
			# Print results.
			print(new_line + "\t" + new_samples, file=outputfile)

		else:
			line = line.strip()
			print(line, file=outputfile)

	#Close output file
	outputfile.close()

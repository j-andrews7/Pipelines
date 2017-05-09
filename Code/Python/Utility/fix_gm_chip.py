#!/usr/bin/env python3
"""
For GM ChIP-seq data from ENCODE, cut last column with number of TFs binding that location. Removes Rep1/2/3 from end of TFs if present. 
Then removes repeats from those replicates if needed.

Usage: python3 fix_gm_chip.py <GM12878_TF151_merged1.txt> <output.txt>

01/15/2016 - jared.andrews07@gmail.com
"""

import sys

####-Main-####
# Grab args.
input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as f:

	output = open(output_file, "w")

	for line in f:

		line = line.strip().split("\t")
		TFs = line[3].split(";")

		fixed_TFs = []

		# Remove Reps, change to uppercase, add to list if not already present.
		for item in TFs:
			if item.endswith("Rep1") or item.endswith("Rep2") or item.endswith("Rep3"):
				item = item[:-4]
			item = item.strip().upper()
			if item not in fixed_TFs:
				fixed_TFs.append(item)

		# Add back delimiter.
		out_TFs = ";".join(fixed_TFs)
		fixed_TFs[:] = []

		print(line[0], line[1], line[2], out_TFs, sep="\t", file=output)

print("Complete")
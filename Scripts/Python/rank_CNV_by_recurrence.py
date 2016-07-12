#!/usr/bin/env python3
"""
Ranks SE-CNV intersects by recurrence. Assumes CNV samples in 8th column delimited by commas.

Usage: python3 rank_CNV_by_recurrence.py <input.bed> <output.bed>

Args:
    (required) <input.bed> = Name of locus list file to process.
    (required) <output.bed> = Name of output file to be created.
"""

import sys
from operator import itemgetter

input_file = sys.argv[1]
output_file = sys.argv[2]

####-MAIN-####

# Open input file. Lazy.
with open(input_file) as f:

	# Get all lines and break each into a list.
	lines = [line.strip().split() for line in f]

	# Add counts for CNV samples in each line to end of line.
	for i in lines:
		idx = lines.index(i)
		samps = i[8].split(",")
		i.append(len(samps))
		lines[idx] = i

	# Sort.
	lines.sort(key=itemgetter(11), reverse=True)

	# Open output.
	output = open(output_file, "w")

	# Make header and print to output.
	header = "\t".join(["SE_CHR","SE_ST","SE_END","SE_SAMPS", "SE_GENES", "CNV_CHR","CNV_ST","CNV_END","CNV_SAMPS","CNV_ANNOT","CNV_GENES"])
	print(header, file=output)

	# Print
	for line in lines:
		wanted = line[0:-1]
		print("\t".join(wanted), file=output)
	

	# Close output file.
	output.close()
#!/usr/bin/env python3
"""
Grab expression for DL3A538 by grabbing value for DL3B538 (which isn't a sample) and moving it to DL3A538 lines.

Usage: python3 fix_dl_circuit_expression.py <input.txt> <output.txt> 
"""

import sys

####-Variables-####

input_name = sys.argv[1]
output_name = sys.argv[2]

# Use to hold expression value for DL3B538 for each gene.
exp_dict = {}

# Open input file.
with open(input_name) as f:

	# Skip header.
	f.readline()

	for line in f:

		line = line.strip().split("\t")
		gene = line[-2]
		val = line[-1]

		if "DL3B538" in line:
			exp_dict[gene] = val

# Open output.
output = open(output_name, "w")

# Open input file. Lazy.
with open(input_name) as f:

	# Get header and print to output.
	header = f.readline().strip()
	print(header, file=output)

	for line in f:
		line = line.strip().split("\t")
		gene = line[-2]

		if "DL3A538" in line:
			line[-1] = exp_dict[gene]
			out_line = "\t".join(line)
			print(out_line, file=output)

		# Don't print lines with fake sample.
		elif "DL3B538" in line:
			continue

		else:
			out_line = "\t".join(line)
			print(out_line, file=output)

# Close output file.
output.close()

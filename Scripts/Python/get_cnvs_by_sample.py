#!/usr/bin/env python3
"""Given a list of unmerged amps or dels with samples in the 4th column, prints the cnvs for every sample in the file to an 
individual file.

Usage: python3 get_cnvs_by_sample.py -i <input.bed> -o <output suffix>

Args:
    -i input.bed = Name of input file to process.
    -o output suffix = Suffix to add to append to the sample names to be used as file names.
"""

import sys
import argparse

parser = argparse.ArgumentParser(usage=__doc__)

# Options
parser.add_argument("-o", "--output", dest = "output_suf", required=True)
parser.add_argument("-i", "--input", dest = "input_file", required=True)


args = parser.parse_args()

input_file = args.input_file
output_suf = args.output_suf

# Lists to hold sample names and the lines for each sample.
sample_names = []
sample_data = []

# Open input.
with open(input_file) as f:

	for line in f:
		line = line.strip()
		data = line.split()

		sample = data[3]

		if sample not in sample_names:
			sample_names.append(sample)
			sample_data.append([])
		else:
			idx = sample_names.index(sample)
			sample_data[idx].append(line)

# Go through samples and print to output.
for entry in sample_names:
	ind = sample_names.index(entry)
	out = open(entry + output_suf, "w")

	# Get data corresponding to the sample name.
	lines = sample_data[ind]

	for l in lines:
		print(l, file=out)

	out.close()
#!/usr/bin/env python3
"""
For the output of a multiintersect file with clustering, output all the unique elements for each data column into a new file based on the header.
Must have a header.

Usage: python3 parse_multiinter_output.py <input.bed>

01/21/2016 - jared.andrews07@gmail.com
"""

import sys

####-Main-####
# Grab args.
input_file = sys.argv[1]

# Use to hold list for each sample.
master_list = []
counts = []
with open(input_file) as fi:

	header = fi.readline().strip()
	# Get sample names.
	samples = header.split("\t")[5:]

	# Add new list to master list for each sample column.
	for entry in samples:
		master_list.append([])
		counts.append(0)

	for line in fi:

		line = line.strip()
		line_elements = line.split("\t")
		pos = "\t".join(line_elements[0:3])
		s_list = line_elements[4].split(",")

		if len(s_list) == 1:
			# Get index of sample in master_list to add entry to.
			samp_idx = int(s_list[0]) - 1
			counts[0] = counts[0] + 1
			master_list[samp_idx].append(pos)
		else:
			count_idx = len(s_list) - 1
			counts[count_idx] = counts[count_idx] + 1

# Iterate through each column list in the master list.
for x in master_list:
	samp_name = samples[master_list.index(x)]
	output = open(samp_name + "_unique.bed", "w")
	for line in x:
		print(line,file=output)

for y in counts:
	idx = int(counts.index(y)) + 1
	print("Present in " + str(idx) + " cell type(s):" + str(y))

print("Complete")
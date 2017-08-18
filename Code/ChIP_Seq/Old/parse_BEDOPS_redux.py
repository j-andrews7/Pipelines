#!/usr/bin/env python3
"""Parse output from BEDOPS unions and intersects for appropriate output format.

Usage: python3 parse_BEDOPS.py > output 
"""

import sys

#Dictionary to hold info for each unique SE.
entries = {}

for line in sys.stdin:

	line = line.strip().split(";")

	#If only only sample in the line, just print it out.
	if len(line) == 1:
		out_line = "\t".join(line) + "\n"
		sys.stdout.write(out_line)

	else:
		first_samp = line[0].split("\t")
		chrom = first_samp[0]
		start = first_samp[1]
		end = first_samp[2]
		genes = first_samp[4]
		samples = [first_samp[3]]

		for next_sample in line[1:]:
			next_sample = next_sample.split("\t")
			#Get new sample name
			samples.append(next_sample[3])
			#Check start and stop positions and change if needed to get largest range.
			if int(next_sample[1]) < int(start):
				start = next_sample[1]
			if int(next_sample[2]) > int(end):
				end = next_sample[2]
			#Check genes and add them if necessary.
			curr_genes = next_sample[4]
			curr_genes = curr_genes.lstrip("[").rstrip("]").split(", ")
			#Lazy addition to gene string. Could make a list and keep track, but lazy.
			for item in curr_genes:
				if item not in genes:
					genes = genes.rstrip("]")
					genes = genes + ", " + item + "]"

		pos = (chrom,(start,end))

		#Check if position exists in dict, if not, add it. If yes, check number of samples and if new line has more, use it instead.
		if pos in entries:
			old_samps = entries[pos][0]
			#Check number of samples in each line and replace dictionary samples for the position if necessary.
			if len(samples) > len(old_samps):
				entries[pos] = (samples,genes)
		else:
			entries[pos] = (samples,genes)

for data in entries:
	line_chrom = data[0]
	line_pos = data[1]
	line_start = line_pos[0]
	line_end = line_pos[1]
	sample_list = ";".join(entries[data][0])
	genes = entries[data][1]
	o_line = line_chrom + "\t" + line_start + "\t" + line_end + "\t" + sample_list + "\t" + genes + "\n"
	sys.stdout.write(o_line)
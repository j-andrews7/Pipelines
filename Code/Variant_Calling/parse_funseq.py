#!/usr/bin/env python3
"""
Remove coding variants after funseq annotation and prints each sample to its own file. Assumes bed output option from funseq.
Also prints just the positions to another output file for each sample in the file, so end up with four files per sample in the file. Bases file names
on the ID in the 6th column of the input file. 

'-uniq' flag will only output the variants for each sample that aren't found in any other samples.

Usage: python3 parse_funseq.py -i <input.bed> 

Args:
    -i input.bed = Name of input file to process.
    -uniq (optional) = Only output variants for each sample that are not found in any other samples.
"""

import sys
import argparse
parser = argparse.ArgumentParser(usage=__doc__)

# Potential options
parser.add_argument("-uniq", "--unique", action="store_true")
parser.add_argument("-i", "--input", dest = "input_file", required=True)

args = parser.parse_args()

uniq = args.unique
input_file = args.input_file

with open(input_file) as f:

	header = f.readline().strip()

	# Use to determine if line is a new sample or not. Keep track of previous line's sample ID.
	new_samp = True
	prev_ID = ""

	for line in f:

		line = line.strip()
		if line.startswith("chrX") or line.startswith("chrY"):
			continue
		new_line = line.split("\t")

		chrom, start, stop, ID = new_line[0], new_line[1], new_line[2], new_line[5]
		coding_var = new_line[6].split(";")[1]

		samps = new_line[6].split(";")[14]

		if ID != prev_ID:
			new_samp = True

		if new_samp:
			new_samp = False

			if uniq:
				output = open(ID+".noncodingOnly.uniq.bed", "w")
				positions = open(ID+".noncodingOnly.positions.uniq.bed", "w")
				coding = open(ID+".codingOnly.uniq.bed", "w")
				coding_pos = open(ID+".codingOnly.positions.uniq.bed", "w")

				print(header, file=output)
				print(header, file=coding)

				if coding_var == "No" and samps == ".":
					print(chrom, start, stop, sep="\t", file=positions)
					line = line.replace('&',',')
					print(line,file=output)
					prev_ID = ID

				# Get those that are coding.
				elif samps == ".":
					print(chrom, start, stop, sep="\t", file=coding_pos)
					line = line.replace('&',',')
					print(line,file=coding)
					prev_ID = ID

			else:
				output = open(ID+".noncodingOnly.bed", "w")
				positions = open(ID+".noncodingOnly.positions.bed", "w")
				coding = open(ID+".codingOnly.bed", "w")
				coding_pos = open(ID+".codingOnly.positions.bed", "w")

				print(header, file=output)
				print(header, file=coding)

				if coding_var == "No":
					print(chrom, start, stop, sep="\t", file=positions)
					line = line.replace('&',',')
					print(line,file=output)
					prev_ID = ID

				# Get those that are coding.
				else:
					print(chrom, start, stop, sep="\t", file=coding_pos)
					line = line.replace('&',',')
					print(line,file=coding)
					prev_ID = ID

		elif coding_var == "No":
			if uniq:
				if samps == ".":
					print(chrom, start, stop, sep="\t", file=positions)
					line = line.replace('&',',')
					print(line,file=output)
					prev_ID = ID
			else:
				print(chrom, start, stop, sep="\t", file=positions)
				line = line.replace('&',',')
				print(line,file=output)
				prev_ID = ID

		elif coding_var == "Yes":
			if uniq:
				if samps == ".":
					print(chrom, start, stop, sep="\t", file=coding_pos)
					line = line.replace('&',',')
					print(line,file=coding)
					prev_ID = ID
			else:
				print(chrom, start, stop, sep="\t", file=coding_pos)
				line = line.replace('&',',')
				print(line,file=coding)
				prev_ID = ID

	output.close()
	positions.close()
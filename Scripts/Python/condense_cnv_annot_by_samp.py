#!/usr/bin/env python3
"""Takes an annotated amps or del bed file (with genes) and condenses the common CNV annotations so the same positions
for a sample aren't repeated multiple times if they overlapped several annotations.

Usage: python3 condense_cnv_annot_by_samp.py -i <input.bed> -o <output.bed>

Args:
    -i input.bed = Name of input file to process.
    -o output.bed = Name of output file.
"""

import sys
import argparse
from subprocess import call
parser = argparse.ArgumentParser(usage=__doc__)

# Potential options
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-i", "--input", dest = "input_file", required=True)

args = parser.parse_args()

input_file = args.input_file
output_file = args.output_file

# Dict to hold each line and combine the annotations.
results = {}

# Open input.
with open(input_file) as f:

	# These files shouldn't have headers so no need to skip lines.
	for line in f:
		line = line.strip().split("\t")
		annot = line[6]

		# Create list for data that stays the same across lines except for the common cnv annotation.
		static_data = line[0:6]
		static_data.append(line[7])
		static_data = "\t".join(static_data)

		if static_data not in results:
			results[static_data] = [annot]
		elif annot not in results[static_data]:
			annot_list = results[static_data]
			annot_list.append(annot)
			results[static_data] = annot_list

temp_out = open(output_file+".temp", "w")

# Create new line and print it out to a temp file.
for item in results:
	annots = results[item]
	annots = ";".join(annots)

	item = item.split()
	new_line = "\t".join([item[0],item[1],item[2],item[3],item[4],item[5],annots,item[6]])
	print(new_line, file=temp_out)

temp_out.close()

# Sort temp file and print to real output file.
cmd = 'sort -k1,1 -k2,2n ' + output_file + ".temp"
cmd = cmd.split()

with open(output_file,'w') as fout: 
    call(cmd,stdout=fout)

# Remove temp file.
cmd = 'rm ' + output_file + ".temp"
cmd = cmd.split()
call(cmd)

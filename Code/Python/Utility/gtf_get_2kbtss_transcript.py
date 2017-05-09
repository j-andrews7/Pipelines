#!/usr/bin/env python3
"""For a .gtf file, grab lines for transcripts and only retain columns for chr, TSS-2kb, TSS+2kb, and gene name.
End up with a .bed file containing these columns.

Usage: python3 gtf_get_2kbtss_transcript.py <input.gtf> 

Args:
    input.gtf = Name of input file to process.
"""

import sys

input_file = sys.argv[1]
base = input_file.split(".")[0]
output_file = base + "_2kbTSS.bed"

with open(input_file) as f:

	output_f = open(output_file, "w")
	print("Processing...")

	for line in f:

		# Skip header/comment lines.
		if line.startswith("#"):
			continue

		line = line.strip().split("\t")

		# Parse if line is a transcript.
		if line[2] == "transcript":
			gene_name = line[8].split(";")[4].split()[1].strip('"')

			if line[6] == "-":
				tss = line[4]
			else:
				tss = int(line[3]) - 1

			tss_start = int(tss) - 2000
			tss_end = int(tss) + 2000

			print(line[0],str(tss_start),str(tss_end),gene_name,sep="\t", file=output_f)


	#Close output file
	output_f.close()

print("COMPLETE")
#!/usr/bin/env python3
"""For a .gtf file, parse out lines for genes and only retain columns for chr, start, stop, and gene name.
End up with a .bed file containing these columns.

Usage: python3 parse_anno_gtf_get_genes.py <input.gtf> 

Args:
    input.gtf = Name of input file to process.
"""

import sys

#Store file names
input_file = sys.argv[1]
base = input_file.split(".")[0]
output_file = base + "_genes_only.bed"

with open(input_file) as f:

	#Open output file, "w" to make it writable
	output_f = open(output_file, "w")
	print("Processing...")

	#Iterate through each line in input file
	for line in f:

		#Skip lines that start with #
		if line.startswith("#"):
			print("Skipping line.")
			continue

		#Strip newline character and split line by tab.
		line = line.strip().split("\t")

		#Check if third column is "gene", skip if not.
		if line[2] == "gene":
			gene_name = line[8].split(";")[4].split()[1].strip('"')

			#Print to outputfile
			print(line[0],str((int(line[3])-1)),line[4],gene_name,sep="\t", file=output_f)


	#Close output file
	output_f.close()

print("COMPLETE")
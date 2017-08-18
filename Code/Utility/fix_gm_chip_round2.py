#!/usr/bin/env python3
"""
For partially processed GM ChIP-seq data from ENCODE, check if each TF/Gene name can be found in the FPKM table. 
If not found, will print that it's not found so that it can be changed/removed as needed. Also removes duplicate TFs
from each line that may have used different antibodies, etc, so they weren't removed in the first script.

Usage: python3 fix_gm_chip_round2.py <GM12878_TF151_fixed.bed> <FPKM_table.txt> <output.txt>

01/15/2016 - jared.andrews07@gmail.com
"""

import sys

####-Main-####
# Grab args.
input_file = sys.argv[1]
fpkm_file = sys.argv[2]
output_file = sys.argv[3]

fpkm_genes = []
missed_genes = []

# Get all genes in FPKM file.
with open(fpkm_file) as f:
	for line in f:
		line = line.strip().split("\t")
		gene = line[2]
		if gene not in fpkm_genes:
			fpkm_genes.append(gene)

# For each TF in the GM ChIP-seq file, check if it's name is located in the FPKM gene list.
# If not, scream at user.
with open(input_file) as fi:

	output = open(output_file, "w")

	for line in fi:

		line = line.strip().split("\t")
		TFs = line[3].split(";")

		prev_TFs = []

		# Remove Reps, change to uppercase, add to list if not already present.
		for item in TFs:
			item = item.strip().upper()
			# Skip those that are only numeric.
			if item.isdigit():
				continue
			if item not in fpkm_genes and item not in missed_genes:
				missed_genes.append(item)
			if item not in prev_TFs:
				prev_TFs.append(item)

		# Add back delimiter.
		if len(prev_TFs) > 0:
			out_TFs = ";".join(prev_TFs)
			prev_TFs[:] = []

			print(line[0], line[1], line[2], out_TFs, sep="\t", file=output)

	print("Missed these, dingus: " + str(missed_genes))

print("Complete")
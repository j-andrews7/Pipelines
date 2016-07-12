#!/usr/bin/env python3
"""Given a gene list and an input file, prints lines from the input file that have genes found in the gene list in the 
specified gene column. Gene list should have one gene symbol per line.

Usage: python /scratch/jandrews/bin/pull_interesting_genes.py -i <input.bed> -o <output.bed> -g <gene_list.txt> -c <1-based position of gene column> -d <"delimiter">

Args:
    -i input.bed = Name of input file to process.
    -o output.bed = Name of output file.
    -g gene_list.txt = Name of gene list file.
    -c Gene column = 1-based column in which gene(s) reside in input file. 
    -d Delimiter = Delimiter of the gene column in the input file if they are lists. Should be quoted (e.g. ";", not just ;).
"""

import sys
import argparse

parser = argparse.ArgumentParser(usage=__doc__)

# Options
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-i", "--input", dest = "input_file", required=True)
parser.add_argument("-g", "--gene", dest = "gene_file", required=True)
parser.add_argument("-c", "--col", dest = "gene_col", type=int, required=True)
parser.add_argument("-d", "--delim", dest = "delim", required=False, default="$")

args = parser.parse_args()

# Parser variables.
input_file = args.input_file
output_file = args.output_file
gene_file = args.gene_file
gene_col = args.gene_col - 1
delim = args.delim

# Open output.
out = open(output_file, "w")

# Get genes from gene list.
int_genes = [line.strip() for line in open(gene_file)]

# Open input.
with open(input_file) as f:

	# Get and print header.
	header = f.readline().strip()
	print(header, file=out)

	# Iterate through input file line by line.
	for line in f:
		line = line.strip()
		new_line = line.strip().split("\t")
		genes = new_line[gene_col].split(delim)

		# Check if any of the line's genes are in the gene list.
		for item in genes:
			if item in int_genes:
				print(line, file=out)
				break

# Close output
out.close()
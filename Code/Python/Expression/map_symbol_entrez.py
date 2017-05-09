#!/usr/bin/env python3
"""
Maps gene symbols to entrez gene IDs from Affy NetAffx batch queries containing both version 1 and 2 U133 array genes.

Usage: map_symbol_entrez.py <array_gene_list.txt> <input.txt> <output.txt>
"""

#Necessary packages
import sys

array_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

gene_dict = {}

with open(array_file) as f:

	f.readline()

	for line in f:
		line = line.strip().split("\t")
		genes = line[1].split("///")
		clean_genes = []
		entrez_ids = line[3].split("///")

		for gene in genes:
			gene = gene.strip()
			if gene not in clean_genes:
				clean_genes.append(gene)

		for gene in clean_genes:
			idx = clean_genes.index(gene)

			if gene not in gene_dict and len(clean_genes) == len(entrez_ids):
				e_id = entrez_ids[idx]
				e_id = e_id.strip()
				gene_dict[gene] = e_id
				print("Not skipped.")
			else:
				print("Skipped.")


output = open(output_file,"w")

with open(input_file) as f:

	out_list = []

	for line in f:

		id_list = line.strip().split("///")

		for entry in id_list:
			entry = entry.strip()
			if entry not in out_list and entry in gene_dict:
				out_list.append(gene_dict[entry])

	for item in out_list:
		print(entry, file=output)

output.close()
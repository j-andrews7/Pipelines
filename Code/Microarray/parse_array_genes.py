#!/usr/bin/env python3
"""
Grabs the entrez ID from output from the Affy NetAffx batch query tool.

Usage: parse_array_genes.py <input.txt> <output.txt>
"""

#Necessary packages
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

output = open(output_file,"w")

with open(input_file) as f:

	f.readline()

	out_list = []

	for line in f:

		ids = line.strip().split("\t")[3]

		id_list = ids.split("///")

		for entry in id_list:
			entry = entry.strip()
			if entry not in out_list and entry != "" and entry != "---":
				out_list.append(entry)

	for item in out_list:
		print(item, file=output)

output.close()
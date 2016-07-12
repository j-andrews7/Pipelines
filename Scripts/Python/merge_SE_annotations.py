#!/usr/bin/env python3
"""Merge Gencode and Refseq annotations for SEs after calling with ROSE. Adds sample name to 4th column and annotations to 5th column.

Usage: python3 merge_SE_annotations.py <gencode.bed> <refseq.txt> <output.bed> 

Args:
    gencode.bed = Name of SE file annotated with gencode from our master annotation files using bedtools closest.
    refseq.txt = Name of SE file annotated with hg19 refseq from ROSE geneMapper.py.
    output.bed = Name of output file.
"""

import sys
from collections import defaultdict

#Store file names
gencode_file = sys.argv[1]
refseq_file = sys.argv[2]
output_file = sys.argv[3]

sample = gencode_file.split("_")[0]

#Dictionary to hold SE_IDs and corresponding annotations for SE. Values are placed in a list, a new key is created if necessary, otherwise values are appended to the list.
SE_dict = defaultdict(list)

#Open gencode bed file.
with open(gencode_file) as f:

	for line in f: 

		#Strip newline character and split line by white space.
		line = line.strip().split("\t")

		#Hold genes for line.
		gene = line[9]

		#Grab SE position.
		chrom = line[0]
		start = line[1]
		end = line [2]

		SE_pos = (chrom,(start,end))

		SE_dict[SE_pos].append(gene)

#Now refseq
with open(refseq_file) as f:

	#Iterate through each line in input file
	for line in f:

		#Hold genes for line.
		genes = []

		#Strip newline character and split line by white space.
		line = line.strip().split("\t")

		#Grab SE position.
		chrom = line[1]
		start = line[2]
		end = line [3]

		SE_pos = (chrom,(start,end))

		for item in line[9:]:
			items = item.split(",")

			for thing in items:
				if thing not in genes and thing is not "":
					genes.append(thing)

		#Add each gene to the dictionary.
		for entry in genes:
			if entry not in SE_dict[SE_pos]:
				SE_dict[SE_pos].append(entry)

#Open output file
output = open(output_file, "w")

#Print results
for result in SE_dict:
	chrom = result[0]
	pos = result[1]
	genes = SE_dict[result]
	print(chrom,pos[0],pos[1],sample,genes, sep="\t", file=output)

output.close()

print("COMPLETE")

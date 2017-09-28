#!/usr/bin/env python3

"""
For a given vcf file, copy the QUAL value into an INFO field (QV), thereby allowing operations to be performed on it (avg, sum, etc)
when merging VCF files with bcftools, as it only reports the MAX QUAL value for each variant from the given samples. 

Usage: python3 add_qual_to_info.py <VCF_file> 
"""
import sys

####-Functions-####
def get_qual(line):
	""" Return QUAL value for a given line. """

	line = line.strip().split()
	qual = line[5]
	return qual
	

####-Main-####
# Grab input files
in_file = sys.argv[1]
out_file = in_file.split(".")[0]
print("Adding QUAL value to INFO field (QV).")

output = open(out_file+"_wQV.vcf","w")
info_line = '##INFO=<ID=QV,Number=1,Type=Float,Description="Quality value from QUAL field for sample.">'

# Open input file.
with open(in_file) as f:
	header = []
	head_print = False

	# Iterate through lines.
	for line in f:
		line = line.strip()
		# Store header lines in a list.
		if line.startswith('##'):
			header.append(line)
		elif line.startswith('#'):
			header.append(info_line)
			header.append(line)
		else:
			# Check if header lines have been printed, print them if not.
			if head_print == False:
				for item in header:
					print(item,file=output)
					head_print = True
			# Print data lines with qual val in info fields.
			else:
				qual_val = get_qual(line)
				line = line.split("\t")
				line[7] = line[7] + ";QV=" + qual_val
				new_line = "\t".join(line)
				print(new_line,file=output)
# Close output
output.close()
print("Complete.")
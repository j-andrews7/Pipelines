#!/usr/bin/env python3
"""
For a merged BCFTools vcf and a merged VarScan vcf, concatenate the two files. Duplicate records will print the record from whichever file has the variant
called in more samples. If they are the same, the record from the first file will be used. This is lazily coded, so it eats memory.
If each file is large (several Gb), may have trouble running locally.

Ensure column order in header of each file is the same before use.

Usage: python3 cat_merged_coding_vcfs.py <merged_bcf.vcf> <merged_varscan.vcf> <output.vcf>

01/14/2016 - jared.andrews07@gmail.com
"""

import sys
import copy

####-Functions-####
def check_num(s):
	""" Check if a number is in the string, returns True if so, otherwise False."""
	return any(i.isdigit() for i in s)


def parse_files(file1, file2):
	""" 
	Grabs the header lines from each file and creates two dicts, one for each file containing the position of each variant and 
	the number of samples called at the position for that file. 

	parameters:
		file1 = First file supplied.
		file2 = Second file supplied.

	returns:
		header = A list of the header lines from both files.
		file1_dict = Dict containing the positions and number of samples containing each variant for file1.
		file2_dict = Dict containing the positions and number of samples containing each variant for file2.
	"""

	print("Pulling records from each file.")

	header = []
	file1_dict = {}
	file2_dict = {}

	with open(file1) as f:

		for line in f:

			# Count number of samples with variant called.
			count = 0

			# Check for header lines and add them to header list.
			if line.startswith("##"):
				line = line.strip()
				header.append(line)
				continue

			line = line.strip().split("\t")
			# Have to use position and ref allele to take care of multiallelic sites. 
			var_pos = ((line[0], line[1]), line[3])
			data = line[9:]

			for item in data:
				if check_num(item):
					count += 1

			file1_dict[var_pos] = count

		header.append('##INFO=<ID=BOTH,Number=1,Type=String,Description="Present in BCFTools (BCF), VarScan (VS), or BOTH">')

	with open(file2) as f:

		for line in f:

			# Count number of samples with variant called.
			count = 0

			# Check for header lines and add them to header list.
			if line.startswith("##"):
				line = line.strip()
				header.append(line)
				continue

			if line.startswith("#"):
				line = line.strip()
				header.append(line)
				continue

			line = line.strip().split("\t")
			# Have to use position and ref allele to take care of multiallelic sites. 
			var_pos = ((line[0], line[1]), line[3])
			data = line[9:]

			for item in data:
				if check_num(item):
					count += 1

			file2_dict[var_pos] = count

	return header, file1_dict, file2_dict


def compare_files(file1, file2, output_file):
	""" 
	Compares the two files, printing  

	parameters:
		file1 = First file supplied.
		file2 = Second file supplied.
	"""

	# Get headers and dicts for each position in each file.
	header, f1_dict, f2_dict = parse_files(file1, file2)

	print("Comparing files and printing to output.")

	output = open(output_file, "w")

	for entry in header:
		print(entry, file=output)

	# Open first file.
	with open(file1) as f:

		for line in f:

			if line.startswith("#"):
				continue

			line = line.strip()
			curr_line = line.split("\t")

			# Get variant position and ref allele.
			data = ((curr_line[0], curr_line[1]), curr_line[3])

			# Check if it's present in both files.
			if data in f1_dict and data in f2_dict:
				f1_count = f1_dict[data]
				f2_count = f2_dict[data]

				# Mark in info field if variant found in both files
				curr_line[7] = curr_line[7] + ";BOTH=Yes"
				out_line = "\t".join(curr_line)

				# Check number of samples called in each if it's present in both.
				if f1_count >= f2_count:
					print(out_line, file=output)

			# If not, go ahead and print to output.
			else:
				curr_line[7] = curr_line[7] + ";BOTH=BCF"
				out_line = "\t".join(curr_line)
				print(out_line, file=output)

	# Then second file.
	with open(file2) as f:

		for line in f:

			if line.startswith("#"):
				continue

			line = line.strip()
			curr_line = line.split("\t")

			# Get variant position and ref allele.
			data = ((curr_line[0], curr_line[1]), curr_line[3])

			# Check if it's present in both files.
			if data in f1_dict and data in f2_dict:
				f1_count = f1_dict[data]
				f2_count = f2_dict[data]

				# Mark in info field if variant found in both files
				curr_line[7] = curr_line[7] + ";BOTH=Yes"
				out_line = "\t".join(curr_line)

				# Check number of samples called in each if it's present in both.
				if f2_count > f1_count:
					print(out_line, file=output)

			# If not, go ahead and print to output.
			else:
				curr_line[7] = curr_line[7] + ";BOTH=VS"
				out_line = "\t".join(curr_line)
				print(out_line, file=output)

	output.close()


####-Main-####
# Grab args.
input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
output_file = sys.argv[3]

compare_files(input_file1, input_file2, output_file)

print("Complete")





#!/usr/bin/env python3

"""
For a merged non_coding VCF file, filter calls that are only called for one mark type.
Format of samples in header must be 'sample_mark'.

Usage: python3 filter_nc_vcf.py <input.vcf> <output.vcf>
"""

import sys
import copy

####-Functions-####
def Get_Marks_From_Header(header_line):
	""" For the header with the samples, determine which columns contain what mark. Return a dict with the index for each sample column and 
	the mark it represents.

	parameters:
		header_line = Line containing the sample names and marks for each column.

	returns:
		marks_dict = A dict with the index for each sample column as the key and the mark it represents as the value.
	"""

	print("Determining mark for each sample column.")

	# Initialize dict to hold index of sample columns and associated mark.
	marks_dict = {}

	# Split header and store as a list, pulling the elements pertaining to sample columns.
	header_line = header_line.strip().split("\t")
	samples = header_line[9:]

	# Get mark and index of each sample in the sample list.
	for item in samples:
		idx = samples.index(item)
		mark = item.split("_")[-1]
		marks_dict[idx] = mark

	return marks_dict

def Check_Num(s):
	""" Check if a number is in the string, returns True if so, otherwise False."""
	return any(i.isdigit() for i in s)


####-Main-####
# Grab args.
input_file = sys.argv[1]
output_file = sys.argv[2]

# Open input VCF.
with open(input_file) as f:

	# Set variables to hold current and previous line position and sample marks.
	curr_pos = 0
	curr_marks = []
	prev_pos = 0
	prev_marks = []
	curr_line = None
	prev_line = None

	# Open output file.
	output = open(output_file, "w")

	# Iterate through lines of input.
	for line in f:

		# Use to keep track of element index, as those that are blank won't be unique so the index() function won't help.
		count = 0

		# Check for and print header lines.
		if line.startswith("##"):
			line=line.strip()
			print(line,file=output)
			continue

		# Get marks for each sample column in a dict based on index from true header.
		if line.startswith("#"):
			line=line.strip()
			print(line, file=output)
			mark_dict = Get_Marks_From_Header(line)
			print("Removing variants found in only one data type.")
			continue

		# Strip and split line.
		line_l = line.strip().split("\t")
		data = line_l[9:]  # Get data in a list.
		curr_pos = line_l[1]
		curr_line = line.strip()

		# Iterate through data
		for item in data:
			# Check if entry isn't blank and get mark for the column if not.
			if Check_Num(item):
				col_mark = mark_dict[count]

				# Check if mark isn't already in mark list and add it if not.
				if col_mark not in curr_marks:
					curr_marks.append(col_mark)
			count += 1

		# Handle multiallelic sites. If one allele found in one data type and another from the next line in another data type, still want those to be printed.
		if (prev_line is not None) and (curr_pos == prev_pos) and (len(prev_marks) > 0):
			for entry in prev_marks:
				if entry not in curr_marks:
					curr_marks.append(entry)

			if len(curr_marks) > 1:
				print(prev_line, file=output)

		if (curr_pos != prev_pos) and (len(prev_marks) > 1) and (prev_line is not None):
			print(prev_line, file=output)

		# Set current variables to previous variables.
		prev_line = copy.copy(curr_line)
		prev_pos = copy.copy(curr_pos)
		prev_marks = copy.copy(curr_marks)

		# Reset variables.
		curr_line = None
		curr_pos = 0
		curr_marks = []

	output.close()

print("Complete")





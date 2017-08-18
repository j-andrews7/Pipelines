#!/usr/bin/env python3
"""For CNCHP files from SNP arrays, create a single table containing the segments in each sample in appropriate format for GISTIC input.
Also creates a file for the markers used. 

Usage: python3 prep_for_gistic.py <sample_names.txt> <output.txt> CNCHP_files

Args:
    sample_names.txt = Name of file containing array number from CNCHP files and a corresponding sample name, one per line, tab-delimited.
    	Only need the number that begins the array file, e.g., 18117 for a file named 18117.CN5.CNCHP. Uses '.' as delim to identify this. 
    output.txt = Name of output segmentation file to create.
    CNCHP_files = CNCHP files created by Affymetrix Genotyping Console Copy Number/LOH Analysis on .cel files. Can list them or use wildcard to pass in.
"""
#02/05/2016
#jared.andrews07@gmail.com
#Check CN_Analysis_Pipeline.txt for more info. 

import sys

####-Functions-####
def Make_Markers_File(CNCHP_File, output_file_name):
	""" Creates markers file (Marker_ID, chrom, position)."""

	print("Creating markers file: " + output_file_name + ".markers")
	markers_file = open(output_file_name + ".markers", "w")

	with open(CNCHP_File) as f:
		for line in f:
			# Skip header lines
			if line.startswith("#"):
				continue
			line = line.strip().split("\t")
			print(str(line[0]),str(line[1]),str(line[2]), sep="\t", file=markers_file)


def Get_Samples(sample_file):
	""" Get the CNCHP numbers and corresponding sample names from the sample_file."""
	print("Connecting CNCHP IDs with sample names.")
	sample_names = {}

	with open(sample_file) as f:
		for line in f:
			line = line.strip().split("\t")
			if len(line) > 2:
				print("Sample names file has too many elements in a line. Fix it, dingus.")
				sys.exit()
			CN_ID, sample = line[0], line[1]
			sample_names[CN_ID] = sample

	return sample_names


def Process_CNCHP_File(CNCHP_file, samp_name, output_file):
	""" Process a CNCHP file to create segments and get CN log2ratios."""

	marker_count = 0
	value_total = 0
	curr_CN = None
	seg_CN = None
	seg_start = None
	seg_end = None
	curr_chrom = None
	seg_chrom = None
	prev_pos = None

	first = True

	with open(CNCHP_file) as f:
		for line in f:
			# Skip header.
			if line.startswith("#") or line.startswith("P"):
				continue

			line = line.strip().split("\t")
			curr_chrom, pos, curr_CN, val = line[1], int(line[2]), int(line[3]), float(line[4])

			# If first line in file, set things appropriately.
			if first:
				seg_start = pos
				prev_pos = pos
				seg_chrom = curr_chrom
				seg_CN = curr_CN
				marker_count += 1
				value_total += val
				first = False

			# Check if segment ends.	
			elif curr_chrom != seg_chrom or curr_CN != seg_CN:
				print(curr_chrom, seg_chrom, str(curr_CN), str(seg_CN), sep="\t")
				seg_end = prev_pos
				val_avg = round(value_total/marker_count, 4)
				print(samp_name, str(seg_chrom), str(seg_start), str(seg_end), str(marker_count), str(val_avg), sep="\t", file=output_file)
				marker_count = 1
				value_total = val
				seg_chrom = curr_chrom
				seg_start = pos
				seg_CN = curr_CN
				prev_pos = pos


			# If a running segment, take care of things as needed.
			else:
				marker_count += 1
				value_total += val
				prev_pos = pos
		
		# Take care of last seg.
		val_avg = round(value_total/marker_count, 4)
		print(samp_name, str(seg_chrom), str(seg_start), str(prev_pos), str(marker_count), str(val_avg), sep="\t", file=output_file)




####-Variables-####
# Store file names
samp_file = sys.argv[1] 
output_name = sys.argv[2]
CN_file_list = sys.argv[3:]

####-MAIN-####
# Run functions defined above
samples = Get_Samples(samp_file)
output = open(output_name, "w")

# Use to check if it's the first file so markers file can be made.
first_file_done = False

# Print header to output.
print("Sample", "Chrom", "Start", "End", "#_Markers", "Avg_Log2", sep="\t", file=output)

# Cycle through file list
for CNCHP in CN_file_list:

	CN_ID = CNCHP.split(".")[0]
	try:
		sample_name = samples[CN_ID]
	except Exception as e:
		print("No sample found for ", CN_ID ," in the sample names file. Check it.")
		sys.exit()

	# Make markers file if necessary.
	if not first_file_done:
		Make_Markers_File(CNCHP, output_name)
		first_file_done = True

	print("Processing: ", CNCHP)

	Process_CNCHP_File(CNCHP, sample_name, output)

output.close()


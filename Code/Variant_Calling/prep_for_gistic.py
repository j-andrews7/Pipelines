#!/usr/bin/env python3
"""For CNCHP files from SNP arrays, create a single table containing the segments in each sample in appropriate format for GISTIC input.
Also creates a file for the markers used. Utilizes segment summary table as well, so both Copy Number/LOG analysis and the segmenting tool
within the Affymetrix Genotyping Console must be used - be sure to check the box to create the summary table!

Usage: python3 prep_for_gistic.py <segment_summary.txt> <sample_names.txt> <output.txt> CNCHP_files

Args:
	segment_summary.txt = Segment summary table from Affymetrix Genotyping Console with samples named according to the second column in the sample_names file.
		Must have following column order: sample, chrom, start, end, #markers, copy number state (optional, really)
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
			print(str(line[0]),str(line[1]),str(line[2]), sep="\t", file=)


def Get_Samples(sample_file):
	""" Get the CNCHP numbers and corresponding sample names from the sample_file."""
	print("Connecting CNCHP IDs with sample names.")
	samples_names = {}

	with open(sample_file) as f:
		for line in f:
			line = line.strip().split("\t")
			if len(line) > 2:
				print("Sample names file has too many elements in a line. Fix it, dingus.")
				sys.exit()
			CN_ID, sample = line[0], line[1]
			samples_name[CN_ID] = sample

	return sample_names


def Process_Segment(seg_loc, markers_list):
	""" 
	For a gained/lost segment in the segment file, calculate the avg log2ratio for the markers located in the seg
	and check to be sure the # of markers actually located in the seg matches what's listed in the seg file. 

	Args:
		seg_pos = Position of the segment in a tuple - (chrom,(start,end)).
		markers_list = List of markers ((chrom,position),log2ratio value) for the sample.

	Returns:
		mark_count = Count of markers located in seg.
		value_avg = Average of log2ratio for markers.
	"""
	# Store segment chrom, start, and end.
	seg_chrom = seg_loc[0]
	seg_pos = seg_loc[1]
	seg_start = seg_pos[0]
	seg_end = seg_pos[1]

	# Count and running value total.
	mark_count = 0
	total_value = 0

	# Get only markers on correct chromosome to save time.
	proper_chrom_markers = [marker for marker in markers_list if (marker[0])[0] == seg_chrom]

	for item in proper_chrom_markers:
		item_pos = (item[0])[1]
		item_val = float(item[1])

		# Check if marker is in segment, add to count and log2 value if so.
		if int(item_pos) >= int(seg_start) and int(item_pos) <= int(seg_end):
			mark_count += 1 
			total_value += item_val

	# Calculate the average.
	value_avg = total_value/mark_count

	return mark_count, total_value


def Get_CNCHP_File(prefix, file_list):
	""" Given the ID prefix, return the appropriate file name."""

	for item in file_list:
		if item.contains(prefix):
			return item
		else:
			print("CNCHP prefix " + str(prefix) + " not matching any CNCHP files in the sample names list. Check sample names file.")
			print(file_list)
			sys.exit()


def Parse_CNCHP_File(CNCHP_file):
	""" Parse a CNCHP file to return a list of the marker positions and log2ratio value for each."""

	markers_list = []

	with open(CNCHP_file) as f:
		for line in f:
			# Skip headers.
			if line.startswith("#") or line.startswith("P"):
				continue

			line = line.strip().split("\t")
			chrom, pos, value = line[1], line[2], line[4]
			# Add to the dictionary.
			markers_list.append(((chrom,pos),value))

	return markers_list



####-Variables-####
# Store file names
seg_file = sys.argv[1]
samp_file = sys.argv[2] 
output_name = sys.argv[3]
CN_file_list = sys.argv[4:]

####-MAIN-####
# Run functions defined above
samples = Get_Samples(samp_file)

output = open(output_name, "w")

with open(seg_file) as f:
	# Get header and print to output
	header = f.readline().strip()
	print(header, file=output)

	# Determine if first sample has been processed or not.
	first_samp_done = False
	new_samp = True
	prev_samp = None  # Hold previous sample
	create_seg = True  # Whether to try to create a new seg or not between established segments.

	for line in f: 
		line = line.strip().split("\t")
		curr_sample = line[0]

		# Determine if new sample or not.
		if curr_sample == prev_samp:
			new_samp = False
		else:
			new_samp = True
			create_seg = True

		if new_samp:
			# Match sample name to CNCHP prefix in samples dict.
			CN_prefix = samples[curr_sample]
			CNCHP_file = Get_CNCHP_File(CN_prefix, CN_file_list)  # Grab appropriate file.

			print("Parsing data from " + CNCHP_file + " for sample " + curr_sample)

			# Get markers data for the new sample.
			curr_markers = Parse_CNCHP_File(CNCHP_file)

		# If first sample hasn't been processed, make markers file.
		if not first_samp_done:
			Make_Markers_File(CNCHP_file)
			first_samp_done = True

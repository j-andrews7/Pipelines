#!/usr/bin/env python3
"""For GFF files with K27AC RPM values for a given SE, calculate the K27AC signal for each file at the position and print results to a single file.
Each sample will be a new column.

Usage: python3 calc_SE_signal.py <output.bed> <input.gff files> 
"""

import sys
import subprocess

# Store file names
output_file = sys.argv[1]
input_files = sys.argv[2:]

temp_output = open(output_file+".temp", "w")

# Dict to hold positions of each SE and the data for each position (in a list).
data = {}
first_file = True
# List to hold sample names from each file.
sample_names = []

#Open input file
for item in input_files:

	name = item.split(".")[0]
	sample_names.append(name)

	with open(item) as f:

		# Skip header.
		f.readline()

		# Initialize each entry in the dict if first file.
		if first_file:
			for line in f:
				line = line.strip().split("\t")
				value = float(line[2])
				chrom = line[1].split("(")[0]
				positions = line[1].split(":")[1]
				start, end = int(positions.split("-")[0]) - 1, int(positions.split("-")[1])
				# Calc actual signal.
				signal = round(((end - start) * value),4)
				entry = (chrom, (start,end))
				data[entry] = [str(signal)]

			first_file = False

		else:
			for line in f:
				line = line.strip().split("\t")
				value = float(line[2])
				chrom = line[1].split("(")[0]
				positions = line[1].split(":")[1]
				start, end = int(positions.split("-")[0]) - 1, int(positions.split("-")[1])
				# Calc actual signal.
				signal = round(((end - start) * value),4)
				entry = (chrom, (start,end))
				data[entry].append(str(signal))

samples = "\t".join(sample_names)

# Save the header for later.
header = "CHROM" + "\t" + "START" + "\t" + "END" + "\t" + "SE_ID" + "\t" + samples


# Now print data dictionary.
for item in data:
	values = data[item]
	chrom, position = item[0], item[1]
	start, end = position[0], position[1]
	out_values = "\t".join(values)
	print(chrom, str(start), str(end), out_values, sep="\t", file=temp_output)

temp_output.close()

# Sort output.
cmd = "sort -k1,1 -k2,2n " + output_file + ".temp"
cmd = cmd.split()

with open(output_file+".sort",'w') as fout: 
    subprocess.call(cmd,stdout=fout)

# Sort output.
subprocess.call("rm " + output_file + ".temp", shell=True)

out = open(output_file, "w")

# Add SE_IDs. This is so lazy.
with open(output_file+".sort") as f:
	print(header, file=out)

	se_id = 1

	for line in f:
		line = line.strip().split()
		pos, data = line[0:3], line[3:]
		new_line = "\t".join(pos) + "\t" + str(se_id) + "\t" + "\t".join(data)
		print(new_line, file=out)
		se_id += 1

subprocess.call("rm " + output_file + ".sort", shell=True)

print("COMPLETE")
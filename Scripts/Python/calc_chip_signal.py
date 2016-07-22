#!/usr/bin/env python3
"""For GFF files with ChIP RPM values for given genomic regions, calculate the signal for each file at the position and print results to a single file.
Each sample will be a new column.

Usage: python3 calc_chip_signal.py <output.bed> <input.gff files> 
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
				ID = line[0]
				value = float(line[2])
				chrom = line[1].split("(")[0]
				positions = line[1].split(":")[1]
				start, end = int(positions.split("-")[0]) - 1, int(positions.split("-")[1])
				# Calc actual signal.
				signal = round(((end - start) * value),4)
				entry = ((chrom,ID), (start,end))
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
				entry = ((chrom,ID), (start,end))
				data[entry].append(str(signal))

samples = "\t".join(sample_names)

# Save the header for later.
header = "CHROM" + "\t" + "START" + "\t" + "END" + "\t" + "ID" + "\t" + samples


# Now print data dictionary.
for item in data:
	values = data[item]
	chrom_id, position = item[0], item[1]
	chrom, ID = chrom_id[0], chrom_id[1]
	start, end = position[0], position[1]
	out_values = "\t".join(values)
	print(chrom, str(start), str(end), ID, out_values, sep="\t", file=temp_output)

temp_output.close()

# Sort output.
cmd = "sort -k1,1 -k2,2n " + output_file + ".temp"
cmd = cmd.split()

with open(output_file,'w') as fout: 
    subprocess.call(cmd,stdout=fout)

# Sort output.
subprocess.call("rm " + output_file + ".temp", shell=True)

print("COMPLETE")
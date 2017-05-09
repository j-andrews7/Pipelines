#!/usr/bin/env python3
"""
08/01/2016
jared.andrews07@gmail.com
-------------------------

Given the amps and dels CN matrices for a binned genome, condense them into a single matrix. 
Bins with both amps and dels are changed to -1 so that they can be highlighted during plotting. 

Usage: python3 condense_cn_matrices.py <amps_matrix.bed> <dels_matrix.bed> <output.bed> 

Args:
    amps_matrix.bed = Name of amps matrix file to process.
    dels_matrix.bed = Name of dels matrix file to process.
"""

import sys
from subprocess import call

amps_file = sys.argv[1]
dels_file = sys.argv[2]
out_file = sys.argv[3]

# Initialize dict to hold each bin.
bin_dict = {}

# First the amps file.
with open(amps_file) as f:
	print("Collecting amps data.")

	for line in f:
		line = line.strip().split("\t")
		data = line[3:]  # Stick the data in a list.
		chrom, start, stop = line[0], line[1], line[2]  # Get position info.
		bin_dict[(chrom, (start, stop))] = data  # Store in the dict.

# Then the dels file.
with open(dels_file) as f:
	print("Collecting dels data and processing.")

	for line in f:
		line = line.strip().split("\t")
		data = line[3:]  # Stick the data in a list.
		chrom, start, stop = line[0], line[1], line[2]  # Get position info.

		amp_vals = bin_dict[(chrom, (start, stop))]  # Get the values for the position from the amps file.

		# Make sure the same number of data elements are in both the amp/del files. Kill script if not.
		if len(amp_vals) != len(data):
			print("Differing numbers of columns in the input files! Check 'em out.")
			sys.exit()

		final_vals = []  # List to hold eventual final vals for bin.
		i = 0  # Used at counter to keep index.

		for item in data:
			amp = amp_vals[i]  
			if int(item) < 2 and int(amp) > 2:  # If both an amp and del is in the bin, give it a value of 2 for later plotting.
				final_vals.append("-1")
				i += 1  # Increase counter.
				continue  # Move to next iteration.
			elif int(item) == 2 and int(amp) == 2:
				final_vals.append(str(item))
				i += 1
				continue
			elif int(amp) > 2:  # If the del value is not < 2, use the amp value. 
				final_vals.append(str(amp))
				i += 1
				continue
			elif int(item) < 2:  # If the amp value is not > 2, use the del value.
				final_vals.append(str(item))
				i += 1
				continue
		
		bin_dict[(chrom, (start, stop))] = final_vals  # Update the dict.

out = open(out_file + ".tmp", "w")  # Open temp output file.
tmp_out = out_file + ".tmp"

# Unpack dict and write to output.
print("Printing to temp file.")
for entry in bin_dict:
	chrom, pos = entry[0], entry[1]
	start, stop = str(pos[0]), str(pos[1])
	data = "\t".join(bin_dict[entry])
	line = "\t".join([chrom, start, stop, data])  # Create string for output.
	print(line, file=out)

out.close()

print("Sorting.")
# Sort and remove temp file.
with open(out_file, "w") as fout:
	cmd = ["sort", "-V", "-k1,1", "-k2,2", tmp_out]
	call(cmd, stdout=fout)

cmd = ["rm", tmp_out]
call(cmd)

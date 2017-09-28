#!/usr/bin/env python3
"""For a merged_SE file for a given cell type, only print SEs occurring in one or more samples. An SE found in multiple samples (CC & CB) from a given individual but not in any 
other samples is not considered recurrent.  

Usage: $​ ​python3 get_recurrent_SEs.py <merged_SEs.txt> <output.txt> <delim>
"""

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
delim = sys.argv[3]

output = open(output_file,"w")

with open(input_file) as f:
	for line in f:
		# List to hold numbers in sample names.
		samp_numbers = []
		line = line.strip()
		data = line.split("\t")[3]
		samples = data.split(delim)
		# Use to check if recurrence within an individual is the only recurrence. Or something. I'm half asleep.
		same_ind_recurrence = False

		# Check if CC and CB for same individual are causing recurrence.
		if len(samples) == 2:
			for item in samples:
				item_num = item.lstrip('CB').lstrip('CC')
				samp_numbers.append(item_num)
			if samp_numbers[0] == samp_numbers[1]:
				same_ind_recurrence = True
				print("Only recurrent within one individual.")

		if len(samples) > 1 and same_ind_recurrence == False:
			print(line, file=output)
output.close()
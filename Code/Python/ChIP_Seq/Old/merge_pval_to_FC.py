#!/usr/bin/env python3
"""Adds pvals from ANOVA as last column of the fold change file for a mark.

Usage: python3 merge_pval_to_FC.py <p_val_file.text> <FC_file.txt> <output.txt>

Args:
    p_val_file.txt = Name of file containing p_values to process 
    FC_file.txt = Name of file containing fold changes
    output.txt = Name of output file to be created

NOTE: p values must be in last column of the p_value file. Both files must be sorted, most easily done by peak ID. P vals are appended
to the last column of the FC file. 
"""

import sys

####-Variables-####

#Store file names
pval_input_file = sys.argv[1]
FC_input_file = sys.argv[2]
output_name = sys.argv[3]

#Store pvals in a list
p_vals = []

#open p_value input file
with open(pval_input_file) as f:

	#Skip header
	f.readline()

	#Iterate through each line of file
	for line in f:

		#Split and strip line and store as list
		line = line.strip().split()

		#Append last value in line to p_vals list
		p_vals.append(line[-1])

#Open FC input file
with open(FC_input_file) as fc:

	#Grab header and print to output file with addition of p_val column
	header = fc.readline().strip()

	#Open output file
	outputfile = open(output_name, "w")

	#Print header to output file
	print(header + "\tp_val", file=outputfile)

	#Used to keep track of pval to print
	count = 0

	#Iterate through each line of file
	for line in fc:

		#Add appropriate p_val to line
		line = line.strip() + "\t" + p_vals[count]
		count += 1

		#Print line to output file
		print(line, file=outputfile)

	#Close output file
	outputfile.close()

#Check if pval list was iterated all the way through
if len(p_vals) == count:
	print("Complete.")
else:
	print("Number of p values not consistent with number of lines in FC file, check output!")
	print(str(len(p_vals)) + "\t" + str(count))
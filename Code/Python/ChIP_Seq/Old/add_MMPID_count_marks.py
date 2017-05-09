#!/usr/bin/env python3
"""For a .bed file, add MMPID for each line to the fourth column. Place count of unique marks for the region in 6th column

Usage: python3 add_MMPID_count_marks.py <input.bed> <output.bed>

Args:
    input.bed = Name of input file to process 
    output.bed = Name of output file to be created

NOTE: Sort the input .bed file first. Command: sort -k1.4,1.5 -k2 -V input.bed > output.bed
"""

import sys

#Store file names
input_file = sys.argv[1]
output_name = sys.argv[2]

print("PERFORMING INCREDIBLY COMPLEX FUNCTIONS, PLEASE WAIT SIR/MADAM")

#Open input file
with open(input_file) as f:

	#Open output file, "w" to make it writable
	output_file = open(output_name, "w")

	#Initiate variable to hold what next MMPID should be
	MMPID = 1

	#Iterate through each line in input file
	for line in f:

		#Initiate list to hold marks
		marks_list = []

		#Strip newline character and split line by white space
		line = line.strip()
		line = line.split()

		#Get column with each sample_mark delimited by ","
		marks_col = line[3].split(",")

		#For each member in the list, split by "_" and if the mark isn't already in the marks list, add it
		for entry in marks_col:

			mark = entry.split("_")[1]

			if mark not in marks_list:
				marks_list.append(mark)

		#Get number of marks in the marks list
		marks_count = len(marks_list)

		#Print to outputfile
		print(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + str(MMPID) + "\t" + line[3] + "\t" + str(marks_count), file=output_file)

		#Increase MMPID by one for next line
		MMPID += 1

	#Close output file
	output_file.close()

print("COMPLETE")




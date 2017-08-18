#!/usr/bin/env python3
"""
Append last column of each file given as input to an output file. 

Usage: python3 lazy_outer_join.py <output.txt> <input files>

Data to be joined must be in last column. Grabs all data in the first file for chromosomal positions, etc. 
"""

import sys

####-Variables-####

#Store file names
output_name = sys.argv[1]

#List of input files
input_files = sys.argv[2:]

#Check if first file, if so, get all the columns
first_file = True

#List to hold each column to print to output
output_data = []

#List to hold eventual output header
output_header = ["CHR", "START", "END"]

print("Getting data from each file.")

#Iterate through each input file
for entry in input_files:

	#Split to get the sample name
	file_title = entry.strip().split("_")[-1]
	sample_name = file_title.split(".")[0]

	output_header.append(sample_name)

	#open input file
	with open(entry) as f:

		#Get header, doing nothing with it cause it's a trackline in this case
		header = f.readline()

		#Check if first file, grab entire line if so
		if first_file:
			first_file = False

			#Iterate through each line of file
			for line in f:

				#Split and strip line and store as list
				line = line.strip().split()

				#Append line to output_data list
				output_data.append(line)

		#Else only get last column
		else:

			#Keep track of line number
			count = 0

			#Iterate through each line of file
			for line in f:

				#Split and strip line and store as list
				line = line.strip().split()

				#Get appropriate list (line) to append data to
				output_line = output_data[count]

				#Append data in last column to the appropriate list within the output_data
				output_line.append(line[-1])

				count +=1

#Open output
output = open(output_name, "w")

#Print header
print(*output_header[0:], sep="\t", file=output)

#print data
for item in output_data:
	print(*item[0:], sep="\t", file=output)

#Close output file
output.close()

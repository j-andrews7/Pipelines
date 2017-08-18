#!/usr/bin/env python3
"""Subtract 1 from the start of each line

Usage: python3 correct_start_stop.py <input.txt> <output.txt>

Args:
    input.txt = Name of input file to process 
    output.txt = Name of output file to be created
"""

import sys

####-Variables-####

#Store file names
input_file = sys.argv[1]
output_name = sys.argv[2]

#Dictionary of the offending bins and the correct chromosome stop position.
bin_correct = {("chr1", "249250800"):"249250621", ("chr2", "243199400"):"243199373", ("chr3","198022600"):"198022430", 
	("chr4", "191154400"):"191154276", ("chr5","180915400"):"180915260",("chr6","171115200"):"171115067", ("chr7", "159138800"):"159138663", 
	("chr8", "146364200"):"146364022", ("chr9", "141213600"):"141213431", ("chr10", "135534800"):"135534747",
	("chr11", "135006600"):"135006516", ("chr12", "133852000"):"133851895", ("chr13", "115170000"):"115169878", ("chr14", "107349600"):"107349540", 
	("chr15", "102531400"):"102531392", ("chr16", "90354800"):"90354753", ("chr17", "81195400"):"81195210", ("chr18", "78077400"):"78077248", 
	("chr19", "59129000"):"59128983", ("chr20", "63025600"):"63025520", 
	("chr21", "48130000"):"48129895", ("chr22", "51304600"):"51304566", ("chrX", "155270600"):"155270560"}

#open input file
with open(input_file) as f:

	#Grab header
	header = f.readline().strip()

	#Open output file
	outputfile = open(output_name, "w")

	#Print header to outputfile
	print(header, file=outputfile)

	#Iterate through each line of file
	for line in f:

		#Split and strip line and store as list
		line = line.strip().split()

		#Create tuple from chromosome and stop position for comparison to the bin_correct dictionary 
		tup = (line[0], str(line[2]))

		#Check if the tup is in the bin_correct dictionary
		if tup in bin_correct:
			line[2] = bin_correct[tup]

		#Subtract 1 from start
		line[1] = int(line[1]) - 1

		print(*line[0:], sep="\t",file=outputfile, end="")
		print("",file=outputfile)

	#Close output file
	outputfile.close()



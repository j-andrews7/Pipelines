#!/usr/bin/env python3
"""CLeans up files from wonky R output. Rounds numbers to 2 decimals, fixes chromosomes being converted to random numbers, etc.

Usage: python3 R_cleanup.py <input.txt> <output.txt>

Args:
    input.txt = Name of input file to process 
    output.txt = Name of output file to be created
"""

import sys

inputfile = sys.argv[1]
outputfile = sys.argv[2]

chrom_correct = {"1":"1", "12":"2", "16":"3", "17":"4", "18":"5","19":"6", "20":"7", "21":"8", "22":"9", "2":"10",
	"3":"11", "4":"12", "5":"13", "6":"14", "7":"15", "8":"16", "9":"17", "10":"18", "11":"19", "13":"20", "14":"21", "15":"22", "23":"X"}

#Open input file
with open(inputfile) as f:

	#Open output file, "w" to make it writable
	output_file = open(outputfile, "w")

	#Get header and print it to the output file
	header = f.readline().strip()
	print(header, file=output_file)

	#Iterate through each line in input file
	for line in f:

		#The float function apparently can't handle the small e for scientific notation
		line = line.replace("e","E")

		#Strip newline character and split line by white space
		line = line.strip().split()

		#Add "chr" back to the chromosome element and correct for R changing the chromosome to a random number. 
		line[0] = "chr" + chrom_correct[line[0]]

		#Convert start, stop, peak_ID of line to ints because R output is sometimes silly and will write chromosomal positions in scientific notation
		for i in line[1:4]:
			i_index = line.index(i)
			line[i_index] = int(float(i))

		#Print chr, start, stop, & peak_ID values.
		print(line[0], line[1], line[2], line[3], sep="\t", end="\t",file=output_file)

		#Print the rest of the values
		for val in line[4:]:
			print("{0:.2f}".format(float(val)) + "\t", end="", file=output_file)

		#Print newline
		print("",file=output_file)

	output_file.close()
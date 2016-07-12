#!/usr/bin/env python3
"""
For MACS14 .bed peaks files, replaces the forth column containing generic MACS_PEAK
with the sample name and type of mark

Usage: python3 rename_columns.py <MACS14_peak.bed files>

Args:
    MACS14_peak.bed files = Path to MACS14_peak.bed file(s) 
"""
import sys
import os

#Store each filename in a list
filenames = sys.argv[1:]

#Iterate through each file
for item in filenames:

	#Open file, store as variable
	with open(item) as f:

		#Get all of the filename before the extension
		basename = os.path.splitext(os.path.basename(item))[0]

		#Set the file to write out to. Change "_fixed.bed" if wanted.
		outputname = basename + ".fixed.bed"

		#Split the filename by "_" and get the appropriate elements.
		split_basename = basename.split("_")[0:-1]
		replacement = "_".join(split_basename)

		#Open the file to output to
		outfile = open(outputname, "w")

		#Iterate through each line in input file, split by tab and replace fourth column element
		for line in f:
			line = line.split("\t")
			line[3] = replacement

			#Join elements in the line list
			complete_line = "\t".join(line)

			#Write line to the output file
			print(complete_line, file=outfile)

		outfile.close()
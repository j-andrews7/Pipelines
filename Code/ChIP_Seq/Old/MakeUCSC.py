#!/usr/bin/env python3
"""
Updated 06/07/2016
Author: jared.andrews07@gmail.com

Makes sure ending bin of each chrom doesn't extend past actual end of chromosome, etc. Creates a bedgraph file for each column that can be directly uploaded to UCSC. 

Usage: python3 MakeUCSC.py -i <input.txt> -o <output.txt> 

Args:
    (required) -i input.txt = Name of input file to process. Format must be chrom, start, end, bin or peak ID, then data columns.
    (required) -o output.txt = Name of output file to be created.
"""

import sys
#Used for command line options and such
import argparse
#Used to compress output bedgraph files
import gzip

parser = argparse.ArgumentParser(usage=__doc__)

#Dictionary of the offending bins and the correct chromosome stop position.
bin_correct = {("chr1", "249250800"):"249250621", ("chr2", "243199400"):"243199373", ("chr3","198022600"):"198022430", 
	("chr4", "191154400"):"191154276", ("chr5","180915400"):"180915260",("chr6","171115200"):"171115067", ("chr7", "159138800"):"159138663", 
	("chr8", "146364200"):"146364022", ("chr9", "141213600"):"141213431", ("chr10", "135534800"):"135534747",
	("chr11", "135006600"):"135006516", ("chr12", "133852000"):"133851895", ("chr13", "115170000"):"115169878", ("chr14", "107349600"):"107349540", 
	("chr15", "102531400"):"102531392", ("chr16", "90354800"):"90354753", ("chr17", "81195400"):"81195210", ("chr18", "78077400"):"78077248", 
	("chr19", "59129000"):"59128983", ("chr20", "63025600"):"63025520", 
	("chr21", "48130000"):"48129895", ("chr22", "51304600"):"51304566", ("chrX", "155270600"):"155270560", ("chrY", "59373600"):"59373566"}


####-Functions-####


def Make_UCSC(column_index):
	"""
	Creates a gzipped bedgraph file from the specified data column of the output file. Bases the trackline and file name on the column header.

	Parameters:
		column_index = The index of the column to be used as data for the 4th column of the resulting output file
	"""

	#Open the input (original output) file
	with open(args.output_file) as f:

		#Get the header and split into a list
		header=f.readline().strip().split()

		#Get the appropriate column
		basename=header[column_index]

		#Set the output name
		output_name=basename+".bedgraph.gz"

		#Tell user what's going on
		print("Writing to " + output_name)

		#Open the bedgraph output file
		output = gzip.open(output_name, "wt")

		gzip_header = ("track type=bedGraph name=QN_" + basename + " description=QN_" + basename + 
			" visibility=2 color=200,0,0 altColor=0,100,200 priority=10 autoScale=off alwaysZero=off gridDefault=on graphType=bar viewLimits=0:20\n")

		output.write(gzip_header)

		#Iterate through file and print to the output file
		for line in f:
			line = line.strip().split()

			#Get data
			keep = line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[column_index] + "\n"
			output.write(keep)

		output.close()


####-Variables-####

#-i <input.txt> -o <output.txt> -ucsc <True or False>
parser.add_argument("-i", "--input", dest = "input_file", required=True)
parser.add_argument("-o", "--output", dest = "output_file", required=True)

args = parser.parse_args()

data_cols = 0


####-Main-####
#Run the main loop
#Open input file
with open(args.input_file) as f:

	#Open output file, "w" to make it writable
	output = open(args.output_file, "w")

	#Get header and print it to the output file
	header = f.readline().strip()
	print(header, file=output)

	header = header.split("\t")
	data_cols=len(header[4:])

	#Iterate through each line in input file
	for line in f:

		#Strip newline character and split line by white space
		line = line.strip().split()

		#Create tuple from chromosome and stop position for comparison to the bin_correct dictionary 
		tup = (line[0], str(line[2]))

		#Check if the tup is in the bin_correct dictionary
		if tup in bin_correct:
			line[2] = bin_correct[tup]

		line = "\t".join(line)

		#Print results
		print(line, file=output)

	output.close()

	for sample in header[4:]:
		samp_index = header.index(sample)
		Make_UCSC(samp_index)



	
#!/usr/bin/env python3
"""
05/09/2016
jared.andrews07@gmail.com
-------------------------

Given a matrix of CN counts for a binned genome, identify bins that have the cnv for a given percentage of samples.
Merge them to identify minimal common regions between all the samples. Assumes the input file has a header.

Usage: python3 find_cn_mcrs.py -i <input.bed> -o <output.bed> -p <percentage as decimal>

Args:
    -i input.bed (required) = A matrix containing the bins for the genome and a column for each sample defining whether
    	or not the sample has a cnv for the bin.
    -o output.bed (required) = Name of output file.
    -p (optional) = Percentage of samples that must contain the cnv for a bin to be reported as a MCR. Default = 0.25.
"""

import sys
import argparse
import pybedtools as pbt

####-Variables-####
# Create argument parser
parser = argparse.ArgumentParser(usage=__doc__)

# Add accepted arguments
parser.add_argument("-i", "--input", dest = "input_file", required=True)
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-p", "--perc", dest ="perc", type=float, default=0.25)

# Grab arguments from the parser
args = parser.parse_args()


# Output so user can double check options
print( "Input file: {}\nOutput file: {}\nPercentage cutoff (as decimal): {} \n \n".format(
        args.input_file,
        args.output_file,
        args.perc
        ))

# Assign arguments to more straightforward variables
input_file = args.input_file
output_file = args.output_file
perc = args.perc

####-Functions-####

def GetMCRBins(input_file, output_file, perc):
	"""
	Given a matrix of genomic bins with cnv values in columns, find the bins with >= the specified percentage
	of samples with a cnv and print to the output file with the sum of the columns. Returns the output file name.
	"""

	# Create output name.
	output_file_full = output_file.split(".")[0] + ".MCR_sums.bed"

	out_file = open(output_file_full, "w")  # Open output file name.

	with open(input_file) as f:
		header = f.readline()  # Get header
		header = header.strip().split("\t")
		total = len(header[3:])  # Count total number of samples in file.
		for line in f:
			line = line.strip().split("\t")
			pos = "\t".join(line[0:3])  # Get bin position.
			data = [int(i) for i in line[3:]]  # Convert the list to ints.
			count = abs(sum(data))  # Get count in bin.
			cn_perc = count/total  # Calc percentage of samples with cnv.
			if cn_perc >= perc:
				print(pos + "\t" + str(count), file=out_file)  # Print to output file if percentage of samples with cnv is over the cutoff.
	out_file.close()

	return output_file_full


def MergeMCRBins(mcr_file, output_file):
	"""
	Given a file with bins that met the mcr cutoff, merge them and print to a final output file.
	"""

	# Create output name.
	output_file_full = output_file.split(".")[0] + ".MCRs.bed"

	a = pbt.BedTool(mcr_file)  # Get the bed file.
	b = a.merge().saveas(output_file_full)  # Merge and save it.

	return 


####-Main-####

def main(input_file, output_file, perc):
	"""
	For each bin in input_file, determine if enough samples have the cnv so as to meet the perc cutoff. If so, print 
	these bins to an output file with the sum of the columns. Then take this file and merge the regions to create condensed regions
	pertaining to the actual MCRs.
	"""

	mcr_file = GetMCRBins(input_file, output_file, perc)  # Find bins that meet the cutoff.
	MergeMCRBins(mcr_file, output_file)  # Merge the bins and print to output.


# Actually run.
if __name__ == '__main__':
	main(input_file, output_file, perc)

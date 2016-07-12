#!/usr/bin/env python3
"""
05/09/2016
jared.andrews07@gmail.com
-------------------------

Given a bed file, rank each feature by length and print the specified top percentage to an output file.

Usage: python3 rank_by_length.py -i <input.bed> -o <output.bed> -p <percentage as decimal>

Args:
    -i input.bed (required) = A bed file.
    -o output.bed (required) = Name of output file.
    -p (optional) = Percentage of top hits to print to another output file. Default = 0.05.
"""

import sys
import argparse
import operator

####-Variables-####
# Create argument parser
parser = argparse.ArgumentParser(usage=__doc__)

# Add accepted arguments
parser.add_argument("-i", "--input", dest = "input_file", required=True)
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-p", "--perc", dest ="perc", type=float, default=0.05)

# Grab arguments from the parser
args = parser.parse_args()


# Output so user can double check options
print( "Input file: {}\nOutput file: {}\nTop hits to pull (as decimal): {} \n \n".format(
        args.input_file,
        args.output_file,
        args.perc
        ))

# Assign arguments to more straightforward variables
input_file = args.input_file
output_file = args.output_file
perc = args.perc

####-Functions-####

def GetLengths(input_file, output_file):
	"""
	Given a bed file, create a dictionary for each row with the length of the feature as the key. Returns this dict
	and prints the sorted dict positions to the output file.
	"""

	range_dict = {}  # Make dictionary.

	output = open(output_file, "w")  # Open output file.

	with open(input_file) as f:
		for line in f:
			line = line.strip().split("\t")
			pos = "\t".join(line[0:3])  # Get bin position.
			length = int(line[2]) - int(line[1])
			range_dict[pos] = length
		sorted_lengths = sorted(range_dict.items(), key=operator.itemgetter(1), reverse=True)  # Sort by the values of the dict.
		# Returns a list of tuples that we iterate through and print out.
		for item in sorted_lengths:
			print(item[0] + "\t" + str(item[1]), file=output)

	output.close()

	return sorted_lengths


def PrintTopHits(sorted_lengths, output_file, perc):
	"""
	Given a sorted list of dict elements, print the specified top percentage.
	"""

	perc_out = str(int((perc * 100)))  # Get a string for the percentage.

	# Create output name.
	output_file_full = output_file.split(".")[0] + ".TOP" + perc_out + "PERC.bed"
	out = open(output_file_full, "w")
	num_feats = len(sorted_lengths)  # Get total number of elements in list.
	num_out = int(perc * num_feats)  # Get number of elements to output.

	count = 0  # To count.

	for x in sorted_lengths:
		if count <= num_out:  # If the number we've printed so far is less than the number to output, keep printing.
			print(str(x[0]) + "\t" + str(x[1]), file=out)
			count += 1
		else:
			break

	return 


####-Main-####

def main(input_file, output_file, perc):
	"""
	Get lengths of all elements in the input file and rank by them in descending order. Create another output file
	with only the top x percentage of hits.
	"""

	sorted_lengths = GetLengths(input_file, output_file)  # Calc lengths, rank, and print to output in order.
	PrintTopHits(sorted_lengths, output_file, perc)  # Get top x% hits and print to separate output file.


# Actually run.
if __name__ == '__main__':
	main(input_file, output_file, perc)

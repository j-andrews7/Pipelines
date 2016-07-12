#!/usr/bin/env python3
"""
For a given vcf file (1-based), get reference sequence in which it lies and determine if it matches motifs given in a motif file.

Usage: python3 motif_checker.py -i <input.vcf> -r <reference.fa> -m <motif.txt> -size <bases on each side of variant to grab> -o <output.txt> [OPTIONS]

Args:
    (required) -i <input.vcf> = Name of variant file to process.
    (required) -r <reference.fa> = Name of reference sequence file to get surrounding bases from variant. 3 by default. 
    (required) -o <output.txt> = Name of output file to be created.
    (optional) -m <motif.txt> = Tab-delimited key file containing motif name in first column and sequence in second column in rows.
    (optional) -size <wing size> = An integer to define the distance from each side of the variant position to get sequence.
"""

import sys
# Used for command line options and such
import argparse
from pyfaidx import Fasta
parser = argparse.ArgumentParser(usage=__doc__)


####-Functions-####
def get_surrounding_seq(chromo, var_pos, size):
	""" Return sequence containing variant base + specified number of bases on each side from reference sequence file.

	Args:
		chromo = Chromosome of variant.
		var_pos = Position of variant within chromosome.
		size = Number of bases on each side of variant to return as full sequence.

	Returns:
		ref_seq = Sequence containing the variant base + specified number of bases on each side from reference sequence file.
	"""

	# Generate range to grab sequence from.
	seq_range = (var_pos-size,var_pos+size)

	ref_seq = fa_ind[chromo][seq_range].seq

	return ref_seq

				 


def get_motifs(motif_f):
	"""Return motif sequences with names as a list of tuples from the motif file."""

	motif_list = []

	# Open provided motif file
	with open(motif_f) as f:

		# Iterate through file line by line.
		for line in f:
			line = line.strip().split("\t")
			name = line[0]
			seq = line[1]
			tup = (name,seq)
			motif_list.append(tup)

	return motif_list


def check_motif(motifs, ref_sequence):
	"""Check if reference sequence around variant contains the user-specified motif.

	Args: 	
		motifs = List of motif sequences for a given protein, TF, etc.
		ref_sequence = Sequence to check for motifs in.

	Return:
		result = Boolean for whether motif found or not.
	"""

	result = False

	# Iterate through motif list for given protein/TF and check if it's in the reference sequence containing the variant.
	for item in motifs:
		if item in ref_sequence:
			result = True

	return result


####-PARSER-####

# Create arguments and options
parser.add_argument("-i", "--input", dest = "input_file", required=True)
parser.add_argument("-r", "--ref", dest = "ref_file", required=True)
parser.add_argument("-m", "--motif", dest = "motif_file", required=False, default=None)
parser.add_argument("-size", "--wingsize", dest ="wing_size", type=int, default=3)
parser.add_argument("-o", "--output", dest = "output_file", required=True)

args = parser.parse_args()

# Output so user can double check options
print( "Input file {}\n Reference file {}\n Motif file {}\nWing size {}\n Output file {}\n\n".format(
        args.input_file,
        args.ref_file,
        args.motif_file,
        args.wing_size,
        args.output_file
        ))

# Easier to use argument variables
inp_file = args.input_file
ref_file = args.ref_file
motif_file = args.motif_file
wing_size = args.wing_size
out_file = args.output_file

####-Main-####
# Grab motif list from motif file.
if motif_file is not None:
	print("Creating motif list from " + motif_file)
	motif_l = get_motifs(motif_file)
	# Get all motif names and prep for header
	header_addition = "\tRef_Sequence"
	for item in motif_l:
		header_addition = header_addition + "\t" + item[0] + "_motif?"
else:
	header_addition = "\tRef_Sequence"

# Create index file from input fasta for quick searching
print("Creating index from reference sequence for efficient searching.")
fa_ind = Fasta(ref_file)

print("Analyzing variants. This may take a while.")

# turtle = """\                                      ___-------___
#                                    _-~~             ~~-_
#                                 _-~                    /~-_
#              /^\__/^\         /~  \                   /    \
#            /|  O|| O|        /      \_______________/        \
#           | |___||__|      /       /                \          \
#           |          \    /      /                    \          \
#           |   (_______) /______/     SCIENCE IS SLOW    \_________ \
#           |         / /         \                      /            \
#            \         \^\         \                  /               \     /
#              \         ||           \______________/      _-_       //\__//
#                \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
#                  ~-----||====/~     |==================|       |/~~~~~
#                   (_(__/  ./     /                    \_\      \.
#                          (_(___/                         \_____)_) 
#                          """

# Open output file.
output_f = open(out_file,"w")

# Open VCF file.
with open(input_file) as f:
	for i in range(115):  # Skip all the info lines.
		f.readline()

	# Create appropriate header.
	header = f.readline()
	header = header.strip() + header_addition

	print(header, file=output_f)

	# Iterate through each line of VCF.
	for line in f:

		# Skip header lines.
		if line.startswith("##"):
			next(f)

		line = line.strip()
		line_list = line.split("\t")

		# Skip lines that contain an INDEL.
		if len(line_list[3]) > 1 or len(line_list[4]) > 1:
			next(f)

		chrom = line_list[0]
		v_pos = line_list[1]

		# Get sequence surrounding the variant from the reference sequence.
		surr_seq = get_surrounding_seq(chrom, v_pos, wing_size)
		line = line + "\t" + surr_seq

		# Check if motifs are in the reference sequence.
		for motif in motif_l:
			motif_seqs = motif[1].split(";")
			if check_motif(motif_seqs,surr_seq):
				line = line + "\tYES"
			else:
				line = line + "\tNO"

		print(line, file=output_f)

# Close output file.
output_f.close()

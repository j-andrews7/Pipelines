#!/usr/bin/env python3
"""
Tries to match MMPIDs that hit the FC cutoff specified in the circuit table for a given sample to a CNV occurring in that sample as well.
These tables are a mess, this tries to really cut down the size and make them easier to manage.

Usage: python3 match_FC_to_CNVs.py -t <amp or del> -i <input.bed> -o <output.bed> [OPTIONS]

Args:
	(required) -t <amp or del> = Type of CNVs that are in the file. 
    (required) -i <input.bed> = Name of locus list file to process.
    (required) -o <output.bed> = Name of output file to be created.
    (optional) -cut <value> = Linear FC magnitude cutoff for H3AC/K27AC used to filter out uninteresting MMPIDs (default=2). Results must meet cutoff to be included.
    (optional) -r = If included, only includes MMPIDs that are recurrent after applying the cutoff filters and such.
"""

import sys
from math import log2
from subprocess import call

# Used for command line options and such
import argparse
parser = argparse.ArgumentParser(usage=__doc__)

####-Functions-####

def MatchSamples(line):
	"""Checks if the MMPID sample is found in the CNV samples list. Returns True if so, False if not."""

	line = line.strip().split("\t")
	mmpid_sample = line[4]
	cnv_samps = line[14].split(",")

	if mmpid_sample in cnv_samps:
		return True
	else:
		return False


def CheckFC(line, cutoff, cn_type):
	"""
	Checks if either K27AC or H3AC meet the FC cutoff for a given line. Returns true if so, False if not.
	Must be a loss if the type of CNV is a del, a gain if the CNV is an amp.
	"""
	
	# Log space makes this all easier.
	log2_cut = log2(cutoff)

	line = line.strip().split("\t")

	# Account for no data (0).
	if float(line[6]) != 0:
		h3_fc = log2(float(line[6]))
	else:
		h3_fc = 0

	if float(line[7]) != 0:
		k27_fc = log2(float(line[7]))
	else:
		k27_fc = 0

	# Check if expression values are available for the gene, if not, skip it.
	if str(line[10]) == "0":
		return False

	# Actually if either the H3AC or K27AC values meet the cutoff.
	if cn_type == "AMP":
		if h3_fc >= log2_cut or k27_fc >= log2_cut:
			return True
		else:
			return False

	elif cn_type == "DEL":
		if h3_fc <= -log2_cut or k27_fc <= -log2_cut:
			return True
		else:
			return False
	else:
		print("CNV type must be 'amp' or 'del', please set properly.")
		sys.exit()


def GetMMPIDs(inp_file):
	"""
	Goes through file and creates a dict for each MMPID with the samples as the value. 
	Creates a list for all MMPIDs that had multiple samples. This list is returned.
	"""

	print("Getting recurrent MMPIDs.")

	id_dict = {}

	with open(inp_file) as f:

		# Skip header.
		f.readline()

		for line in f:
			line = line.strip().split("\t")
			mmpid = str(line[3])
			sample = line[4]

			# Add MMPID to dict initially.
			if mmpid not in id_dict:
				id_dict[mmpid] = [sample]

			else:
				s_list = id_dict[mmpid]

				# If sample isn't already in the value list for the MMPID, add it.
				if sample not in s_list:
					s_list.append(sample)
					id_dict[mmpid] = s_list

	# List of MMPIDs to delete from the dict.
	del_list = []

	for item in id_dict:

		# Check for recurrence, delete MMPID if necessary.
		if len(id_dict[item]) == 1:
			del_list.append(item)

	for entry in del_list:	
		del id_dict[entry]

	results = list(id_dict)
	return results


def RecurFilter(inp_file, mmpid_list):
	"""
	Goes through file and removes lines that aren't recurrent for MMPID positions in multiple samples after applying cutoffs and other filters.
	"""

	print("Filtering file for recurrence.")

	with open(inp_file) as f:

		# Open output file and print header.
		output = open("RECUR_"+inp_file, "w")
		header = f.readline().strip()
		print(header, file=output)

		for line in f:
			line = line.strip()
			mmpid = str(line.split("\t")[3])

			if mmpid in mmpid_list:
				print(line, file=output)



####-Variables-####

# Create arguments and options
parser.add_argument("-i", "--input", dest = "input_file", required=True)
parser.add_argument("-t", "--type", dest = "type", required=True)
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-cut", "--cutoff", dest ="cutoff", type=float, default=2)
parser.add_argument("-r", "--recurrent", action="store_true")

args = parser.parse_args()


# Output so user can double check options
print( "Input file: {}\nOutput file: {}\nCNV type: {}\nFC cutoff: {}\nRecurrence filter: {}\n\n".format(
        args.input_file,
        args.output_file,
        args.type,
        args.cutoff,
        args.recurrent
        ))

# Easier to use argument variables.
inp_file = args.input_file
cn_type = args.type.upper()
out_file = args.output_file
fc_cut = args.cutoff
recur = args.recurrent

####-MAIN-####

# Open input file. Lazy.
with open(inp_file) as f:

	# Open output.
	output = open(out_file, "w")

	# Get header and print to output.
	header = f.readline().strip()
	print(header, file=output)

	for line in f:
		line = line.strip()
		
		# Check if the MMPID and CNV samples match.
		if MatchSamples(line):

			# Check if the FC changes are met and print if so.
			if CheckFC(line, fc_cut, cn_type):
				print(line, file=output)

	# Close output file.
	output.close()

# Filter for recurrence if necessary.
if recur:
	mmpids = GetMMPIDs(out_file)
	RecurFilter(out_file, mmpids)


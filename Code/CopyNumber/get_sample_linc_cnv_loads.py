#!/usr/bin/env python3
"""
Given lists of unmerged amps and dels for a sample, gets the expression for all lincs in that sample that lie within the amps/dels, 
outside them, and for those that are unchanged and spits this info out to multiple files. The linc expression is log2 FC to the median/average
of the other samples in the expression file. Also yields a few summary stats that may be helpful.

Linc RNA expression file should contain a header and the format should be:
CHR	START	STOP	GENE_ID	GENE_SHORT_NAME	<SAMPLE_FPKM_COLUMNS>

Usage: python3 get_sample_linc_cnv_loads.py <amps.bed> <dels.bed> <linc_expression.txt> <overlap_percentage>

Args:
    amps.bed = Name of amps file to process.
    dels.bed = Name of dels file to process.
    linc_expression.txt = Name of linc_expression file.
    overlap_percentage = Percent of linc that must overlap the CNV for it to be considered 'overlapping'.
"""

import sys
from math import log2
from statistics import median
from statistics import mean

amps_file = sys.argv[1]
dels_file = sys.argv[2]
linc_exp = sys.argv[3]
ovlp_perc = sys.argv[4]

####-FUNCTIONS-####

def GetSampleIdx(header, sample):
	""" 
	Get sample index in a given header. 
	"""
	
	header = header.strip().split("\t")

	for entry in header:
		if sample in entry:
			return header.index(entry)

	# Kill script if the sample isn't found.
	print("Sample not found in header. Check input files if you expect it.")
	sys.exit()


def GetSampleLincData(sample, linc_exp):
	""" 
	Get the data for each linc for that sample. Get log2 fold change for FPKM compared to average and median of all samples for the linc.

	Args:
		sample = Sample from the input file name.
		linc_exp = Name of file containing the expression data for each linc in every sample.

	Returns:
		linc_dict = Dict containing signal of every SE position for the sample {(chr, (start, stop)): ((linc_id, linc_name), signal)}
	"""

	# Dict to hold data.
	linc_dict = {}

	with open(linc_exp) as f:

		# Get the sample index.
		header = f.readline().strip()
		sample_idx = GetSampleIdx(header, sample)

		for line in f:
			line = line.strip().split("\t")
			data = [float(x) for x in line[5:]]  # Convert all data to floats.
			linc_med = log2(float(median(data)))  # Get log2 median of list.
			linc_avg = log2(float(mean(data)))  # Get log2 average of the list.
			linc_val = log2(float(line[sample_idx]))  # Get log2 of the linc FPKM for the sample.
			linc_med_FC = linc_val - linc_med
			linc_avg_FC = linc_val - linc_avg

			# Grab data and add to the dict.
			chrom, start, stop, linc_id, linc_name  = line[0], int(line[1]), int(line[2]), line[3], line[4]
			linc_dict[(chrom, (start, stop))] = ((linc_id, linc_name), (linc_med_FC, linc_avg_FC))

	return linc_dict


def GetSampleCNVData(input_file):
	""" 
	Get the positions of each CNV in the input file and return them as a list. 
	"""
	cnvs = []

	with open(input_file) as f:

		for line in f:
			line = line.strip().split("\t")
			chrom, start, stop, num = line[0], line[1], line[2], line[5]
			pos = ((chrom, (start, stop)), num)
			cnvs.append(pos)

	return cnvs


def Calc_CNVS_Length(cn_list):
	""" 
	Gets total length of all CNVs in the list. Returns this length. 
	"""

	length = 0  # Use to hold running total for length.

	# Iterate through each CN, grab the start/end, calculate the length and multiply it by the CN if necessary.
	for item in cn_list:
		cnv_pos = item[0]
		entry_pos = cnv_pos[1]
		cnv_chr = cnv_pos[0]
		cnv_start = int(entry_pos[0])
		cnv_end = int(entry_pos[1])
		cn_num = int(item[1])  # Whether one copy or more were lost/gained.

		# Calculate the length.
		mult = abs(2 - cn_num)
		cn_length = (cnv_end - cnv_start) * mult

		# Add it to the running length.
		length += cn_length

	return length


def Overlap(start1, end1, start2, end2, ovlp_perc):
    """ 
    Does the range (start1, end1) overlap with (start2, end2)? 
    If so, does it do so by the specified percentage?
    """
    
    linc_len = int(end1) - int(start1)  # Get length of linc.

    ovlp_len = max(0, min(int(end1), int(end2)) - max(int(start1), int(start2)))  # Get overlap between features.

    if ovlp_len != 0:
    	ovlp_ratio = ovlp_len/linc_len

    	# Check to see if ovlp percentage is hit.
    	if ovlp_ratio >= float(ovlp_perc):
    		return True

    return False


####-MAIN-####


def main(amps_file, dels_file, linc_exp, ovlp_perc):
	""" 
	Compares each SE to the CNV list and prints data for those that don't overlap to one file and those that do to another. 
	"""

	print("Amps file: " + amps_file)
	print("Dels file: " + dels_file)
	print("Expression data: " + linc_exp)
	ovlp_pr = int(float(ovlp_perc) * 100)
	print("Linc overlap percentage required (1 bp used if 0): " + str(ovlp_pr))
	# Get sample name.
	sample = amps_file.split("_")[0]
	# Get cnv positions.
	print("Grabbing CNV data.")
	amp_pos = GetSampleCNVData(amps_file)
	amp_len = Calc_CNVS_Length(amp_pos)/1000000
	print("MB amps: " + str(amp_len))
	del_pos = GetSampleCNVData(dels_file)
	del_len = Calc_CNVS_Length(del_pos)/1000000
	print("MB dels: " + str(del_len))
	# Get linc data.
	print("Grabbing linc data.")
	linc_data = GetSampleLincData(sample, linc_exp)
	print("Comparing lincs to CNVs.")

	# Open output files.
	amp_out_in = open("LOG2_FC_FPKM_LINCS_IN_AMPS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.bed", "w")
	amp_out_out = open("LOG2_FC_FPKM_LINCS_OUTSIDE_AMPS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.bed", "w")
	amp_sig_in = open("LOG2_FC_FPKM_ONLY_LINCS_IN_AMPS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.txt", "w")
	amp_sig_out = open("LOG2_FC_FPKM_ONLY_LINCS_OUTSIDE_AMPS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.txt", "w")

	del_out_in = open("LOG2_FC_FPKM_LINCS_IN_DELS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.bed", "w")
	del_out_out = open("LOG2_FC_FPKM_LINCS_OUTSIDE_DELS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.bed", "w")
	del_sig_in = open("LOG2_FC_FPKM_ONLY_LINCS_IN_DELS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.txt", "w")
	del_sig_out = open("LOG2_FC_FPKM_ONLY_LINCS_OUTSIDE_DELS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.txt", "w")

	unc = open("LOG2_FC_FPKM_LINCS_UNCHANGED_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.bed", "w")
	unc_sig = open("LOG2_FC_FPKM_ONLY_LINCS_UNCHANGED_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.txt", "w")

	# Print headers.
	header = "\t".join(["CHROM", "START", "END", "LINC_ID", "LINC_NAME", "LOG2_FC_MEDIAN", "LOG2_FC_AVERAGE"])
	print(header, file=amp_out_in)
	print(header, file=amp_out_out)
	print(header, file=del_out_in)
	print(header, file=del_out_out)

	print(header, file=unc)

	# Counts to determine how many lincss lie in amps and how many in dels.
	linc_amps_count = 0
	linc_dels_count = 0

	# Check if each linc overlaps with a CNV, print to appropriate files.
	for item in linc_data:

		# Bools to determine whether linc is found in cnvs or not.
		in_amps = False
		in_dels = False

		# Get all linc data. Tuple unpacking hell.
		linc_chr = str(item[0])
		linc_pos = item[1]
		linc_start = str(linc_pos[0])
		linc_end = str(linc_pos[1])
		id_info = linc_data[item][0]
		vals = linc_data[item][1]
		linc_id = str(id_info[0])
		linc_name = str(id_info[1])
		linc_med_FC = str(vals[0])
		linc_avg_FC = str(vals[1])

		for entry in amp_pos:
			pos = entry[0]
			cn_num = entry[1]
			entry_pos = pos[1]
			amp_chr = pos[0]
			amp_start = str(entry_pos[0])
			amp_end = str(entry_pos[1])

			# Check chroms.
			if linc_chr == amp_chr:
				if Overlap(linc_start, linc_end, amp_start, amp_end, ovlp_perc):
					linc_amps_count += 1
					line = "\t".join([linc_chr, linc_start, linc_end, linc_id, linc_name, linc_med_FC, linc_avg_FC])
					print(line, file=amp_out_in)
					print(linc_med_FC + "\t" + linc_avg_FC, file=amp_sig_in)
					in_amps = True  # Indicate we've found it in the amps.
					break  # Stop checking cnvs.

			else:
				continue

		for entry in del_pos:
			pos = entry[0]
			cn_num = entry[1]
			entry_pos = pos[1]
			del_chr = pos[0]
			del_start = str(entry_pos[0])
			del_end = str(entry_pos[1])

			# Check chroms.
			if linc_chr == del_chr:
				if Overlap(linc_start, linc_end, del_start, del_end, ovlp_perc):
					linc_dels_count += 1
					line = "\t".join([linc_chr, linc_start, linc_end, linc_id, linc_name, linc_med_FC, linc_avg_FC])
					print(line, file=del_out_in)
					print(linc_med_FC + "\t" + linc_avg_FC, file=del_sig_in)
					in_dels = True  # Indicate we've found it in the dels.
					break  # Stop checking cnvs.

			else:
				continue

		# Print to appropriate files for which SE is not in a CNV.
		if not in_amps:
			line = "\t".join([linc_chr, linc_start, linc_end, linc_id, linc_name, linc_med_FC, linc_avg_FC])
			print(line, file=amp_out_out)
			print(linc_med_FC + "\t" + linc_avg_FC, file=amp_sig_out)

		if not in_dels:
			line = "\t".join([linc_chr, linc_start, linc_end, linc_id, linc_name, linc_med_FC, linc_avg_FC])
			print(line, file=del_out_out)
			print(linc_med_FC + "\t" + linc_avg_FC, file=del_sig_out)
			
		if not in_amps and not in_dels:
			line = "\t".join([linc_chr, linc_start, linc_end, linc_id, linc_name, linc_med_FC, linc_avg_FC])
			print(line, file=unc)
			print(linc_med_FC + "\t" + linc_avg_FC, file=unc_sig)

	# Get the lincss altered per MB amp/del.
	linc_amps_MB = linc_amps_count/amp_len
	linc_dels_MB = linc_dels_count/del_len

	stats_out = open(str(ovlp_pr) + "PERC_OVLP" + "_SUMMARY_STATS_ " + sample + ".txt", "w")
	print("LINCS_IN_AMPS" + "\t" + "LINCS_IN_DELS" + "\t" + "MB_AMPS" + "\t" + "MB_DELS" + "\t" + "LINCS/MB_AMP" + "\t" + "LINCS/MB_DEL", file=stats_out)
	print(str(linc_amps_count) + "\t" + str(linc_dels_count) + "\t" + str(amp_len) + "\t" + str(del_len) + "\t" + str(linc_amps_MB) + "\t" + str(linc_dels_MB), file=stats_out)

# Actually run.
if __name__ == '__main__':
	main(amps_file, dels_file, linc_exp, ovlp_perc)

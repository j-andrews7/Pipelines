#!/usr/bin/env python3
"""
Given lists of unmerged amps and dels for a sample, gets the signal for all SEs in that sample that lie within the amps/dels, 
outside them, and for those that are unchanged and spits this info out to multiple files. Also yields a few summary stats that may be helpful.

Usage: python3 get_sample_se_cnv_loads.py <amps.bed> <dels.bed> <SE_signal.bed> <overlap_percentage>

Args:
    amps.bed = Name of amps file to process.
    dels.bed = Name of dels file to process.
    SE_signal.bed = Name of SE signal file.
    overlap_percentage = Percent of SE that must overlap the CNV for it to be considered 'overlapping'.
"""

import sys

amps_file = sys.argv[1]
dels_file = sys.argv[2]
SE_signal = sys.argv[3]
ovlp_perc = sys.argv[4]

####-FUNCTIONS-####

def GetSampleIdx(header, sample):
	""" Get sample index in SE_signal header. """
	
	header = header.strip().split("\t")

	for entry in header:
		if sample in entry:
			return header.index(entry)

	# Kill script if the sample isn't found.
	print("Sample not found in signal file. Check input files if you expect it.")
	sys.exit()


def GetSampleSEData(sample, SE_signal):
	""" 
	Get the data for each SE for that sample. 

	Args:
		sample = Sample from the input file name.
		SE_signal = Name of file containing the signal data for each SE in every sample.

	Returns:
		SE_dict = Dict containing signal of every SE position for the sample {(chr,(start, stop)): (se_id, signal)}
	"""

	# Dict to hold data.
	SE_dict = {}

	with open(SE_signal) as f:

		# Get the sample index.
		header = f.readline().strip()
		sample_idx = GetSampleIdx(header, sample)

		for line in f:
			line = line.strip().split("\t")

			# Grab data and add to the dict.
			chrom, start, stop, se_id, data = line[0], int(line[1]), int(line[2]), line[3], line[sample_idx]
			SE_dict[(chrom, (start, stop))] = (se_id, data)

	return SE_dict


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
    
    SE_len = int(end1) - int(start1)  # Get length of SE.

    ovlp_len = max(0, min(int(end1), int(end2)) - max(int(start1), int(start2)))  # Get overlap between features.

    if ovlp_len != 0:
    	ovlp_ratio = ovlp_len/SE_len

    	# Check to see if ovlp percentage is hit.
    	if ovlp_ratio >= float(ovlp_perc):
    		return True

    return False


####-MAIN-####


def main(amps_file, dels_file, SE_signal, ovlp_perc):
	""" 
	Compares each SE to the CNV list and prints data for those that don't overlap to one file and those that do to another. 
	"""

	print("Amps file: " + amps_file)
	print("Dels file: " + dels_file)
	ovlp_pr = int(float(ovlp_perc) * 100)
	print("SE overlap percentage required (1 bp used if 0): " + str(ovlp_pr))
	print("SE signal data: " + SE_signal)
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
	# Get SE data.
	print("Grabbing SE data.")
	SE_data = GetSampleSEData(sample, SE_signal)
	print("Comparing SEs to CNVs.")

	# Open output files.
	amp_out_in = open("SIGNAL_SES_IN_AMPS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.bed", "w")
	amp_out_out = open("SIGNAL_SES_OUTSIDE_AMPS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.bed", "w")
	amp_sig_in = open("SIGNAL_ONLY_SES_IN_AMPS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.txt", "w")
	amp_sig_out = open("SIGNAL_ONLY_SES_OUTSIDE_AMPS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.txt", "w")

	del_out_in = open("SIGNAL_SES_IN_DELS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.bed", "w")
	del_out_out = open("SIGNAL_SES_OUTSIDE_DELS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.bed", "w")
	del_sig_in = open("SIGNAL_ONLY_SES_IN_DELS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.txt", "w")
	del_sig_out = open("SIGNAL_ONLY_SES_OUTSIDE_DELS_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.txt", "w")

	unc = open("SIGNAL_SES_UNCHANGED_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.bed", "w")
	unc_sig = open("SIGNAL_ONLY_SES_UNCHANGED_" + sample + "_" + str(ovlp_pr) + "PERC_OVLP.txt", "w")

	# Print headers.
	header = "\t".join(["CHROM", "START", "END", "SE_ID", "SIGNAL"])
	print(header, file=amp_out_in)
	print(header, file=amp_out_out)
	print(header, file=del_out_in)
	print(header, file=del_out_out)

	print(header, file=unc)

	# Counts to determine how many SEs lie in amps and how many in dels.
	SE_amps_count = 0
	SE_dels_count = 0

	# Check if each SE overlaps with a CNV, print to appropriate files.
	for item in SE_data:

		# Bools to determine whether SE is found in cnvs or not.
		in_amps = False
		in_dels = False

		# Get all SE data
		se_chr = str(item[0])
		se_pos = item[1]
		se_start = str(se_pos[0])
		se_end = str(se_pos[1])
		se_id = str(SE_data[item][0])
		signal = str(SE_data[item][1])

		for entry in amp_pos:
			pos = entry[0]
			cn_num = entry[1]
			entry_pos = pos[1]
			amp_chr = pos[0]
			amp_start = str(entry_pos[0])
			amp_end = str(entry_pos[1])

			# Check chroms.
			if se_chr == amp_chr:
				if Overlap(se_start, se_end, amp_start, amp_end, ovlp_perc):
					SE_amps_count += 1
					line = "\t".join([se_chr, se_start, se_end, se_id, signal])
					print(line, file=amp_out_in)
					print(signal, file=amp_sig_in)
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
			if se_chr == del_chr:
				if Overlap(se_start, se_end, del_start, del_end, ovlp_perc):
					SE_dels_count += 1
					line = "\t".join([se_chr, se_start, se_end, se_id, signal])
					print(line, file=del_out_in)
					print(signal, file=del_sig_in)
					in_dels = True  # Indicate we've found it in the dels.
					break  # Stop checking cnvs.

			else:
				continue

		# Print to appropriate files for which SE is not in a CNV.
		if not in_amps:
			line = "\t".join([se_chr, se_start, se_end, se_id, signal])
			print(line, file=amp_out_out)
			print(signal, file=amp_sig_out)

		if not in_dels:
			line = "\t".join([se_chr, se_start, se_end, se_id, signal])
			print(line, file=del_out_out)
			print(signal, file=del_sig_out)
			
		if not in_amps and not in_dels:
			line = "\t".join([se_chr, se_start, se_end, se_id, signal])
			print(line, file=unc)
			print(signal, file=unc_sig)

	# Get the SEs altered per MB amp/del.
	SE_amps_MB = SE_amps_count/amp_len
	SE_dels_MB = SE_dels_count/del_len

	stats_out = open(str(ovlp_pr) + "PERC_OVLP" + "_SUMMARY_STATS_" + sample + ".txt", "w")
	print("SE_IN_AMPS" + "\t" + "SE_IN_DELS" + "\t" + "MB_AMPS" + "\t" + "MB_DELS" + "\t" + "SEs/MB_AMP" + "\t" + "SEs/MB_DEL", file=stats_out)
	print(str(SE_amps_count) + "\t" + str(SE_dels_count) + "\t" + str(amp_len) + "\t" + str(del_len) + "\t" + str(SE_amps_MB) + "\t" + str(SE_dels_MB), file=stats_out)

# Actually run.
if __name__ == '__main__':
	main(amps_file, dels_file, SE_signal, ovlp_perc)

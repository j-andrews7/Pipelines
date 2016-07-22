#!/usr/bin/env python3
"""
For a given directory of VCF files (1-based), count the number of SNVs for 96 trinucleotides for each VCF file and print the matrix to output.

Usage: python3 make_trinuc_matrix.py -i <input_dir> -r <reference.fa> -o <output.txt>

Args:
    (required) -i <input_dir> = Name of input directory containing VCFs.
    (required) -r <reference.fa> = Name of reference sequence file to get surrounding bases from variant. 
    (required) -o <output.txt> = Name of output file to be created.
"""

import sys
import os
# Used for command line options and such
import argparse
from pyfaidx import Fasta
parser = argparse.ArgumentParser(usage=__doc__)

tri_nucs = ["TGTT","TGGT","TGCT","TGAT","TGTG","TGGG","TGCG","TGAG","TGTC","TGGC"
	,"TGCC","TGAC","TGTA","TGGA","TGCA","TGAA","TCTT","TCGT","TCCT","TCAT","TCTG"
	,"TCGG","TCCG","TCAG","TCTC","TCGC","TCCC","TCAC","TCTA","TCGA","TCCA","TCAA"
	,"TATT","TAGT","TACT","TAAT","TATG","TAGG","TACG","TAAG","TATC","TAGC","TACC"
	,"TAAC","TATA","TAGA","TACA","TAAA","CAAA","CAAC","CAAG","CAAT","CACA","CACC"
	,"CACG","CACT","CAGA","CAGC","CAGG","CAGT","CATA","CATC","CATG","CATT","CGAA"
	,"CGAC","CGAG","CGAT","CGCA","CGCC","CGCG","CGCT","CGGA","CGGC","CGGG","CGGT"
	,"CGTA","CGTC","CGTG","CGTT","CTAA","CTAC","CTAG","CTAT","CTCA","CTCC","CTCG"
	,"CTCT","CTGA","CTGC","CTGG","CTGT","CTTA","CTTC","CTTG","CTTT"]

tri_nuc_counts = []

####-Functions-####
def get_surrounding_seq(chromo, var_pos, fa_ind):
	""" Return sequence containing variant base + specified number of bases on each side from reference sequence file.

	Args:
		chromo = Chromosome of variant.
		var_pos = Position of variant within chromosome.
		fa_ind = Index to reference genome.

	Returns:
		ref_seq = Sequence containing the variant base + specified number of bases on each side from reference sequence file.
	"""

	# Generate range to grab sequence from.
	seq_start = int(var_pos)-2
	seq_end = int(var_pos)+1
	ref_seq = fa_ind[chromo][seq_start:seq_end].seq

	return ref_seq

#Returns reverse complement of a DNA sequence
def Reverse_Complement(sequence):
    rcSeq = ""
 
    #define complement dictionary
    compDict = {"A": "T", "T":"A", "G":"C", "C":"G", "N":"N"}
 
    #Declare variable to hold reverse complement and assign to it the reversed sequence
    revSeq = sequence.upper()[::-1]
 
    #Create the complement of the reverse sequence by comparing each character of the sequence to the compDict
    for char in revSeq:
        rcSeq+=compDict[char]
 
    #return the reverse complement
    return rcSeq

####-PARSER-####

# Create arguments and options
parser.add_argument("-i", "--input", dest = "input_dir", required=True)
parser.add_argument("-r", "--ref", dest = "ref_file", required=True)
parser.add_argument("-o", "--output", dest = "output_file", required=True)

args = parser.parse_args()

# Output so user can double check options
print( "Input directory {}\n Reference file {}\n Output file {}\n\n".format(
        args.input_dir,
        args.ref_file,
        args.output_file
        ))

# Easier to use argument variables
inp_dir = args.input_dir
ref_file = args.ref_file
out_file = args.output_file

file_list = os.listdir(inp_dir)

####-Main-####

# Create index file from input fasta for quick searching
print("Creating index from reference sequence for efficient searching.")
fa_ind = Fasta(ref_file)
print("Analyzing variants. Creating count matrix.")

header = []

for item in tri_nucs:
	tri_nuc_counts.append([])

num_samps = len(file_list)

s_count = 0
for entry in tri_nuc_counts:
	i = 0
	while i < num_samps:
		entry.append(0)
		tri_nuc_counts[s_count] = entry
		i += 1
	s_count += 1

# Open VCF file.
#Iterate through each file
for item in file_list:

	#Get complete path of file
	fullpath = inp_dir + item

	#Get sample name and store in a list
	sample_name = item.split("_")[0]
	header.append(sample_name)
	samp_idx = header.index(sample_name)

	with open(fullpath) as f:

		# Iterate through each line of VCF.
		for line in f:

			# Skip header lines and those that are INDELS.
			if line.startswith("#") or "INDEL" in line:
				continue

			line_list = line.strip().split("\t")

			chrom = line_list[0]
			v_pos = line_list[1]
			ref_allele = line_list[3].upper()
			var_allele = line_list[4].split(",")

			# Check multiallelic sites
			if len(var_allele) > 1:
				if len(ref_allele) > 1 or len(var_allele[0]) > 1 or len(var_allele[1]) > 1:
					continue
				elif len(ref_allele) > 1 or len(var_allele[0]) > 1:
					continue

			# Get sequence surrounding the variant from the reference sequence.
			try:
				surr_seq = get_surrounding_seq(chrom, v_pos, fa_ind)
			except:  # Skip if it's too close to the end of a chromosome.
				continue

			for item in var_allele:
				if ref_allele.upper() == 'C' or ref_allele.upper() == 'T':
					tri = ref_allele.upper() + item.upper() + surr_seq[0].upper() + surr_seq[2].upper()
				else:
					surr_seq_rc = Reverse_Complement(surr_seq)
					tri = ref_allele.upper() + item.upper() + surr_seq[0].upper() + surr_seq[2].upper()
					tri = Reverse_Complement(ref_allele.upper()) + Reverse_Complement(item.upper()) + surr_seq_rc[0].upper() + surr_seq_rc[2].upper()
				try:
					idx = tri_nucs.index(tri)
					counts = tri_nuc_counts[idx]
					counts[samp_idx] += 1
					tri_nuc_counts[idx] = counts
				except:
					continue


# Open output file.
output_f = open(out_file,"w")

print("\t".join(header), file=output_f)
for x in tri_nucs:
	x_idx = tri_nucs.index(x)
	counts = tri_nuc_counts[x_idx]
	out_counts = "\t".join(str(x) for x in counts)
	print(x + "\t" + out_counts, file=output_f)


# Close output file.
output_f.close()

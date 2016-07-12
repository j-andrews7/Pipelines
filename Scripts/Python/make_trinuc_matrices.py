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
import numpy as np
import pandas as pd
parser = argparse.ArgumentParser(usage=__doc__)

# T>G mutation in TTT context represented as 'TG T.T' in dict below and output.
tri_nucs = ["TG T.T", "TG G.T", "TG C.T", "TG A.T", "TG T.G", "TG G.G", "TG C.G", "TG A.G", "TG T.C", "TG G.C", "TG C.C", "TG A.C", "TG T.A", "TG G.A", "TG C.A", "TG A.A", "TC T.T", "TC G.T", "TC C.T", "TC A.T", "TC T.G", "TC G.G", "TC C.G", "TC A.G", "TC T.C", "TC G.C", "TC C.C", "TC A.C", "TC T.A", "TC G.A", "TC C.A", "TC A.A", "TA T.T", "TA G.T", "TA C.T", "TA A.T", "TA T.G", "TA G.G", "TA C.G", "TA A.G", "TA T.C", "TA G.C", "TA C.C", "TA A.C", "TA T.A", "TA G.A", "TA C.A",
            "TA A.A", "CA A.A", "CA A.C", "CA A.G", "CA A.T", "CA C.A", "CA C.C", "CA C.G", "CA C.T", "CA G.A", "CA G.C", "CA G.G", "CA G.T", "CA T.A", "CA T.C", "CA T.G", "CA T.T", "CG A.A", "CG A.C", "CG A.G", "CG A.T", "CG C.A", "CG C.C", "CG C.G", "CG C.T", "CG G.A", "CG G.C", "CG G.G", "CG G.T", "CG T.A", "CG T.C", "CG T.G", "CG T.T", "CT A.A", "CT A.C", "CT A.G", "CT A.T", "CT C.A", "CT C.C", "CT C.G", "CT C.T", "CT G.A", "CT G.C", "CT G.G", "CT G.T", "CT T.A", "CT T.C", "CT T.G", "CT T.T"]

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


def Reverse_Complement(sequence):
    """
    Returns reverse complement of a nucleotide sequence given as a string.
    """

    rcSeq = ""

    # Define complement dictionary
    compDict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

    # Declare variable to hold reverse complement and assign to it the
    # reversed sequence
    revSeq = sequence.upper()[::-1]

    # Create the complement of the reverse sequence by comparing each
    # character of the sequence to the compDict
    for char in revSeq:
        rcSeq += compDict[char]

        # Return the reverse complement
    return rcSeq


def CreateFreqMatrix(matrix_file):
    """
    Creates a frequency matrix for each sample from a mutation counts matrix.

    Args:
            matrix_file = File containing mutation counts matrix (trinucleotide motifs x samples).
    """

    # Read in file as a dataframe, using the first row as the header and the
    # first column as the index.
    data = pd.read_table(matrix_file, sep="\t", header=0, index_col=0)
    out_file = ".".join(matrix_file.split(".")[:-1]) + ".freq.txt"
    # Sum each column and get frequency of each element.
    data.loc[:, :] = data.loc[:, :].div(data.sum(axis=0), axis=1)
    data.to_csv(out_file, sep="\t")  # Write to output.


####-PARSER-####

# Create arguments and options
parser.add_argument("-i", "--input", dest="input_dir", required=True)
parser.add_argument("-r", "--ref", dest="ref_file", required=True)
parser.add_argument("-o", "--output", dest="output_file", required=True)

args = parser.parse_args()

# Output so user can double check options
print("Input directory: {}\n Reference file: {}\n Output file: {}\n\n".format(
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
# Iterate through each file
for item in file_list:

    # Get complete path of file
    fullpath = inp_dir + item

    # Get sample name and store in a list
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
                    tri = ref_allele.upper() + item.upper() + " " + \
                        surr_seq[0].upper() + "." + surr_seq[2].upper()
                else:
                    surr_seq_rc = Reverse_Complement(surr_seq)
                    tri = Reverse_Complement(ref_allele.upper()) + Reverse_Complement(
                        item.upper()) + " " + surr_seq_rc[0].upper() + "." + surr_seq_rc[2].upper()
                try:
                    if tri.startswith("TC"):  # Debug
                        continue  # Debug
                    idx = tri_nucs.index(tri)
                    counts = tri_nuc_counts[idx]
                    counts[samp_idx] += 1
                    tri_nuc_counts[idx] = counts
                except:
                    continue


# Open output file.
output_f = open(out_file, "w")

print("\t".join(header), file=output_f)
for x in tri_nucs:
    x_idx = tri_nucs.index(x)
    counts = tri_nuc_counts[x_idx]
    out_counts = "\t".join(str(x) for x in counts)
    print(x + "\t" + out_counts, file=output_f)


# Close output file.
output_f.close()

# Create frequency matrix.
CreateFreqMatrix(out_file)

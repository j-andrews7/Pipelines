#!/usr/bin/env python3
"""
For specific Gene Symbols, get the chromosomal positions from a different file and stick in
columns 1-3 of original file.

Usage: python3 get_genesymb_coords.py <gene_data.txt> <gene_coords.txt> <output.txt>

Args:
    gene_data.txt = Name of input file with Gene_Symbols in first column followed by values from
    	samples.
    gene_coords.bed = Name of input file with chromosomal positions followed by Gene Symbols in 4th
    	column. Should not have a header.
    output.txt = Name of output file to be created.
"""

import sys

####-Functions-####


def get_gene_pos(gene_coords):
    """
    Grabs SE and chromosomal positions from a text file.

    Args:
            gene_coords = gene_coords file in txt format

    Returns:
            gene_dict = Dictionary in following format {gene_sym:(chrom, (start, end))}
    """

    gene_dict = {}

    with open(gene_coords) as f:

        for line in f:
            line = line.strip().split()

            gene_sym = line[3]
            chrom = line[0]
            start = line[1]
            stop = line[2]

            gene_dict[gene_sym] = (chrom, (start, stop))

    return gene_dict


def add_positions(gene_data_file, gene_pos_dict, output):
    """
    Adds chromosomal position to the gene_data_file for each gene_sym and writes output to the
    output file.

    Args:
        gene_data_file = data file with gene symbol in first column and data following.
        gene_pos_dict = Dictionary containing the chromosomal positions of each gene in
                the gene_coords file.
        output = Name of output file to write to.
    """

    with open(gene_data_file) as f:
        header = f.readline().strip().split()
        output_file = open(output, "w")
        print("CHR", "START", "STOP", "GENE", *header[1:], sep="\t", file=output_file)

        for line in f:
            line = line.strip().split()
            gene = line[0]
            data = line[1:]

            if gene in gene_pos_dict:
                positions = gene_pos_dict[gene]
                chrom = positions[0]
                start_stop = positions[1]
                start = start_stop[0]
                stop = start_stop[1]

                print(chrom, start, stop, gene, *data[0:], sep="\t", file=output_file)

            else:
                print(gene + " not found in gene position file!")

        output_file.close()


def main(gene_coords_file, gene_data_file, output_name):
    gene_dict = get_gene_pos(gene_coords_file)
    add_positions(gene_data_file, gene_dict, output_name)


if __name__ == '__main__':
    gene_data_file = sys.argv[1]
    gene_coords_file = sys.argv[2]
    output_name = sys.argv[3]
    main(gene_coords_file, gene_data_file, output_name)

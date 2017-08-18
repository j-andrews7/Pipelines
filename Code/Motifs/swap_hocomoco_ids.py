#!/usr/bin/env python3
"""
Replace TF Uniprot IDs in HOCOMOCO PWM file with the gene names from a tsv file from the site.

Usage: python3 swap_hocomoco_ids.py <input_motifs.txt> <metadata_file.tsv> <output_motifs.txt>

Args:
    input_motifs.txt (str): Name of HOCOMOCO PWM file to process - JASPAR format.
    metadata_file.tsv (str): Name of tsv file from HOCOMOCO with TF names in second column and model names in first.
    output_motifs.txt (str): Name of output file with the TF Uniprot IDs replaced with gene names.
"""

import sys


def main(motifs_file, metadata_file, output_file):
    """
    Replace TF Uniprot IDs in HOCOMOCO PWM file with the gene names from a tsv file from the site.

    Usage: python3 parse_hocomoco.py <input_motifs.txt> <metadata_file.tsv> <output_motifs.txt>

    Args:
        input_motifs.txt (str): Name of HOCOMOCO PWM file to process - JASPAR format.
        metadata_file.tsv (str): Name of tsv file from HOCOMOCO with TF names in second column and models in first col.
        output_motifs.txt (str): Name of output file with the TF Uniprot IDs replaced with gene names.
    """
    names_dict = {}
    with open(metadata_file) as f:

        f.readline()  # Skip header.
        for line in f:
            line = line.strip().split()
            model = line[0].split("_")[0]
            factor = line[1]
            names_dict[model] = factor

    output = open(output_file, "w")
    with open(motifs_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                line = line.split()
                tf_id = line[1]
                line[1] = names_dict[tf_id]
                print(*line, sep="\t", file=output)
            else:
                print(line, file=output)
    output.close()


if __name__ == '__main__':
    input_fi = sys.argv[1]
    meta_fi = sys.argv[2]
    output_fi = sys.argv[3]
    main(input_fi, meta_fi, output_fi)

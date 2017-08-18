#!/usr/bin/env python3
"""
For a given bed-like table of QN'd SE signals, utilize the corresponding copy number data to normalize
them appropriately (i.e. 1.5x if CN is 3, 0.5x if CN is 1, etc). User may set the overlap threshold for
the normalization to apply.

Usage: python3 filter_single_SEs.py -i <input.bed> -o <output.bed> [OPTIONS]

Args:
    (required) -se <se_signal.bed> = Name of bed-like table with QN'd SE signal data starting in the 6th
        column.
    (required) -cn = Path to folder containing CNV files for each sample to be normalized.
    (required) -o <output.bed> = Name of output file to be created.
    (optional) -t = Percentage of SE that must be overlapped by the CNV for normalization to apply.
        0.5 (50) by default.
"""

import os
import argparse

# Classes


class Position(object):
    """
    Use to represent and handle genomic ranges more easily.

    Args:
        chrom (str): Chromosome.
        start (int): Start position.
        end (int): End position.
    """

    def __init__(self, chrom, start_pos, end_pos):
        self.chrom = chrom
        self.start = start_pos
        self.end = end_pos

    def overlaps(self, pos_b):
        """
        Return whether self overlaps Position pos_b.

        Args:
            pos_b (Position): Another Position.

        Returns:
            overlap (bool): True if self overlaps with Position pos_b. False if not.
            percent_ovlp (float): Fraction of Position that overlaps Position pos_b.
        """

        overlap = False
        percent_ovlp = None

        if self.chrom == pos_b.chrom:
            start1, start2, end1, end2 = (self.start, pos_b.start, self.end, pos_b.end)

            if end1 >= start2 and end2 >= start1:
                overlap = True
                overlap_amount = get_overlap([self.start, self.end], [pos_b.start, pos_b.end])
                percent_ovlp = (self.end - self.start) / overlap_amount

        return overlap, percent_ovlp


class SuperEnhancer(object):
    """
    Use to hold/process information for each super enhancer and make output easy.
    """

    def __init__(self, line, all_sample_names):
        self.line_list = line.strip().split("\t")
        self.pos = Position(self.line_list[0], int(self.line_list[1]),
                            int(self.line_list[2]))
        self.iden = self.line_list[3]
        self.called = self.line_list[4]
        self.data = self.line_list[5:]
        # Create a dict from the sample names and their corresponding values for easy access.
        self.sample_vals = {sample: self.data[all_sample_names.index(sample)] for sample in all_sample_names}

    def get_output(self, header):
        """
        Return all SuperEnhancer info as a list for easy output.

        Args:
            header (str): Header of SE file used to order sample values in output.

        Returns:
            out_list (list): List of all SuperEnhancer info/data in order.
        """
        header_samps = [x.split('_')[0] for x in header.strip().split()[5:]]

        out_list = self.line_list[0:5]

        for item in header_samps:  # Add all data values to list.
            out_list.append(str(self.sample_vals[item]))

        return out_list

    def normalize_by_cn(self, sample, cn):
        """
        Normalize a SuperEnhancer value for a given sample based on a CNV in the sample.

        Args:
            sample (str): Sample name.
            cn (int): Binary copy number value used to scale data value for the SE in the appropriate sample.
        """

        se_val = float(self.sample_vals[sample])
        # print("OG SE VAL: " + str(se_val))
        if cn == 2:
            scale = 1
        elif cn == 1:
            scale = 2
        elif cn == 0:
            scale = 4
        elif cn == 3:
            scale = (2 / 3)
        elif cn == 4:
            scale = 0.5
        elif cn == 5:
            scale = (2 / 5)
        elif cn == 6:
            scale = (1 / 3)
        # print("CN: " + str(cn))
        # print("Scale: " + str(scale))
        se_val = se_val * scale  # Normalize the value.
        # print("New SE VAL: " + str(se_val))
        self.sample_vals[sample] = se_val


class CopyNumberLocus(object):
    """
    Use to hold/process information for each CNV.
    """

    def __init__(self, line):
        self.line_list = line.strip().split("\t")
        self.pos = Position(self.line_list[0], int(self.line_list[1]),
                            int(self.line_list[2]))
        self.sample = self.line_list[3]
        self.cn = int(self.line_list[5])


# Functions

def get_overlap(a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def parse_se_header(header):
    """
    Parse the se_file header to return sample names with their column indices as a list of tuples.

    Args:
        header (str): Header line of the se_file.

    Returns:
        samples (list of tuples): [(sample_name, column index)]
    """

    header = header.strip().split()
    samples = [x.split("_")[0] for x in header[5:]]

    return samples


def parse_se_file(se_file):
    """
    Parse the se_file into a dictionary, using the SE position and info as keys and dictionaries of the
    sample data as values.

    Args:
        se_file (str): Name of file containing SE data.

    Returns:
        se_list (list): List containing all Super_Enhancer objects created from file..
        header (str): Original header to be printed to output.
    """

    se_list = []

    with open(se_file) as f:
        header = f. readline().strip()
        samples = parse_se_header(header)

        for line in f:
            se_list.append(SuperEnhancer(line.strip(), samples))

    return se_list, header


def parse_all_cnv_files(cn_folder):
    """
    Create a list of CopyNumberLocus objects from all of the CNV files.

    Args:
        cn_folder (str): Path to directory containing all of the CNV files.

    Returns:
        cnv_list (list of CopyNumberLocus objects): List containing all CopyNumberLocus objects from
            all of the CNV files.
    """

    cnv_list = []

    for filename in os.listdir(cn_fold):
        with open(cn_fold + '/' + filename) as f:
            for line in f:
                new_cnv = CopyNumberLocus(line.strip())
                cnv_list.append(new_cnv)

    return cnv_list


def main(se_file, cn_folder, output_f, threshold):
    """
    Iterate through input file, determine those that meet the threshold, print to output file.
    """

    print("Parsing super enhancer file, creating SuperEnhancer objects.")
    se_list, header = parse_se_file(se_f)
    num_ses = len(se_list)
    count = 1
    print(str(num_ses) + " SuperEnhancer objects created.\n")

    print("Parsing all copy number files, creating CopyNumberLocus objects.")
    cnv_list = parse_all_cnv_files(cn_folder)
    print(str(len(cnv_list)) + " CopyNumberLocus objects created.\n")

    print("Iterating through SuperEnhancers to find overlaps and normalize.")

    warning_list = []  # Used to hold warnings about CNV samples not being present in SE file, so they're skipped.

    for se in se_list:
        if count % 50 == 0:
            print("SuperEnhancer " + str(count) + " of " + str(num_ses) + " being processed.")

        se_pos = se.pos
        count += 1

        for cnv in cnv_list:
            overlap, ovlp_perc = se_pos.overlaps(cnv.pos)

            if overlap and ovlp_perc >= threshold:
                try:
                    se.normalize_by_cn(cnv.sample, cnv.cn)
                    break  # Only let it be normalized by one CNV. Those right next to each other are usually the same.
                except:
                    issue = (cnv.sample + " has CNV data, but no SE data! It isn't being used for normalization, " +
                             "but if you're trying CN/SE comparisons, it will likely lead to issues down the line.")
                    if issue not in warning_list:
                        print(issue)
                        warning_list.append(issue)

    print("Printing SuperEnhancer records to output.")
    output = open(out_f, "w")
    print(header, file=output)
    for x in se_list:
        out_line = "\t".join(x.get_output(header))
        print(out_line, file=output)
    output.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-se", "--se_input", dest="se_file", required=True)
    parser.add_argument("-cn", "--cn_folder", dest="cn_folder", required=True)
    parser.add_argument("-o", "--output", dest="output_f", required=True)
    parser.add_argument(
        "-t", "--threshold", dest="threshold", required=False, default=0.5, type=float)

    args = parser.parse_args()

    # Easier to use argument variables
    se_f = args.se_file
    cn_fold = args.cn_folder
    out_f = args.output_f
    thresh = args.threshold

    main(se_f, cn_fold, out_f, thresh)

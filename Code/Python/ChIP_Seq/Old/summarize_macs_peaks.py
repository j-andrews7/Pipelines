#!/usr/bin/env python
"""
Summarizes number of peaks and average FDR rate (%) for a folder of MACS peaks.xls files.
Prints info to an easily copy and pasted file.

Usage: summarize_macs_peaks.py -i input_folder -o <output.txt> [OPTIONS]

Args:
    -i (str): Path to folder containing MACS peaks.xls files to process.
    -o (str): Path to output file.
    -f (float): FDR threshold to count number of peaks below (as percentage).
        (e.g. 10 will count number of peaks with <10% FDR.)
"""

import os
import argparse


def summarize_peaks_file(peaks_file, fdr_thresh):
    """
    Returns string that specifies whether sample had an input control run with it, the number of peaks,
    and the average FDR (%) of the peaks.

    Args:
        peaks_file (str): Path to peaks.xls file.
        fdr_thresh (float): FDR threshold to count number of peaks below (as percentage).

    Returns:
        summary (str): The summary info.
    """

    sample = os.path.basename(peaks_file).split('_')[0]
    peak_count = 0
    fdr_total = 0
    under_thresh = 0

    print("Processing " + peaks_file)

    with open(peaks_file) as f:

        # Skip comments and header.
        line = f.readline().strip()
        while line.startswith("#") or len(line.strip()) == 0:
            line = f.readline().strip()

        line = line.strip().split("\t")

        if len(line) == 8:
            fdr = False
        else:
            fdr = True

        for line in f:
            line = line.strip().split("\t")
            if fdr:
                fdr_total += float(line[8])
                peak_count += 1
                if fdr_thresh and float(line[8]) < fdr_thresh:
                    under_thresh += 1
            else:
                peak_count += 1

    if peak_count > 0:
        avg_fdr = fdr_total / peak_count
    else:
        peak_count = "Missing"
        avg_fdr = "Missing"

    if fdr:
        if fdr_thresh:
            summary = sample + "\t Yes\t" + str(peak_count) + "\t" + str(avg_fdr) + "\t" + str(under_thresh)
        else:
            summary = sample + "\t Yes\t" + str(peak_count) + "\t" + str(avg_fdr)
    else:
        summary = sample + "\t No\t" + str(peak_count) + "\tNA\tNA"

    return summary


def main(peaks_dir, output_file, fdr_thresh):
    """
    Args:
        peaks_dir (str): Path to directory containing the peaks.xls files.
        output_file (str): Path to output file to be created.
        fdr_thresh (float): FDR threshold to count number of peaks below (as percentage).
    """
    output = []

    # Iterate through files and summarize them.
    for filename in os.listdir(peaks_dir):
        if filename.endswith(".xls"):
            summary = summarize_peaks_file(peaks_dir + filename, fdr_thresh)
            output.append(summary)
            continue

    # Print to output.
    out_file = open(output_file, "w")
    print("SAMPLE", "INPUT_CTRL?", "PEAKS", "AVG_FDR%", "PEAKS<" + str(fdr_thresh) + "_FDR", sep="\t", file=out_file)
    print("\n".join(sorted(output)), file=out_file)
    out_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="input_dir", required=True)
    parser.add_argument("-o", "--output", dest="output_file", required=True)
    parser.add_argument("-f", "--fdr", dest="fdr_thresh", required=False, type=float, default=None)
    args = parser.parse_args()

    main(args.input_dir, args.output_file, args.fdr_thresh)

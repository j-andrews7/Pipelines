#!/usr/bin/env python3
"""
For a given bed-like table of QN'd SE signals, creates distribution plots for each SE. Outputs a new
table with additional columns specifying the number of SDs from the mean each sample is for each
SE. Makes a figure for each row.

Usage: python3 SE_dists.py -i <input.bed> -o <output.bed>

Args:
    (required) -i <input.bed> = Name of bed-like table with QN'd SE signal data starting in the 6th
        column. A unique ID should be present for each row in the 4th column.
    (required) -o <output.bed> = Name of output file to be created.
"""

import sys
import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from statistics import mean, median, pstdev
from statistics import pstdev

sns.set(style="white", palette="muted", color_codes=True)


####-PARSER-####

parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument("-i", "--input", dest="input_f", required=True)
parser.add_argument("-o", "--output", dest="output_f", required=True)

args = parser.parse_args()

# Easier to use argument variables
inp_f = args.input_f
out_f = args.output_f


####-FUNCTIONS-####

def calc_mean_stdev(data):
    """
    Given a list of data, return the mean and standard deviation for the population.
    """

    pop_stdev = pstdev(data)
    pop_mean = mean(data)

    return pop_mean, pop_stdev


def calc_median_mad(data):
    """
    Given a 1-d array (list) of data, return the median and median
    absolute deviation (MAD) for the population.
    """
    pop_median = median(data)
    mad = median([abs(x - pop_median) for x in data])

    return pop_median, mad


def get_sample_idx(sample, header):
    """
    For a given sample, return the column index of the sample in the header.
    """

    for item in header:
        if sample in item:
            return header.index(item)

    print(sample + " not found in header, check input files.")
    sys.exit()


def make_distplot(data, output_f, title, xlabel, prefix):
    """
    For a 1-d array of data, create a rugplot showing where the data lie.

    Args:
        data = 1-d array (list) of data.
        output_f = Name of output table. Figure will be named the same as the output table,
            but with `title` prepended to it and a `.pdf` extension.
        title = Title to use for plot.
        xlabel = Label to use for x-axis.
        prefix = Prefix to add to file name.
    """

    plt.figure(figsize=(8, 8), dpi=1200)
    displot = sns.distplot(data, hist=False, rug=True, color="b")
    out_name = prefix + "_" + title + "_" + output_f.split(".")[0] + ".pdf"
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Density')
    displot.figure.savefig(out_name)
    plt.close()


def main(inp_f, out_f):
    """
    Iterate through input file, determine those that meet the threshold, print to output file.
    """

    output = open(out_f, "w")

    with open(inp_f) as f:

        print("Processing and creating plots. This may take a while depending "
              "on the number of lines of input.")

        # Piece together new header.
        header = f.readline().strip().split("\t")
        new_header_samps = "\t".join([x.split("_")[0] + "_SDs_AROUND_MEAN" for x in header[
                                     5:]]) + "\t" + "\t". join([x.split("_")[0] + "_MADs_AROUND_MEDIAN" for x in header[5:]])
        new_header = "\t".join(header) + "\t" + "MEAN" + "\t" + "SD" + \
            "\t" + "MEDIAN" + "\t" + "MAD" + "\t" + new_header_samps
        print(new_header, file=output)

        for line in f:
            line = line.strip().split("\t")
            se_id = "SE " + str(line[3])
            # Assumes actual data starts in 6th column.
            data = [float(x) for x in line[5:]]
            se_samples = line[4].split(";")

            # Mean and std_dev of signal for that line.
            p_mean, p_stdev = calc_mean_stdev(data)
            p_median, p_mad = calc_median_mad(data)

            # Get MADs from mean for each sample.
            if p_mad == 0:
                mad_diffs = ["NA" for x in data]
                out_mad_diffs = "\t".join(mad_diffs)
            else:
                mad_diffs = [((x - p_median) / p_mad) for x in data]
                out_mad_diffs = "\t".join(
                    ["{0:.4f}".format(float(x)) for x in mad_diffs])

                # Make plot. 
                x_label = "MADs around Mean"
                prefix = "SD"
                make_distplot(mad_diffs, out_f, se_id, x_label, prefix)

            # Get SDs from mean for each sample.
            if p_stdev == 0:
                sd_diffs = ["NA" for x in data]
                out_sd_diffs = "\t".join(sd_diffs)
            else:
                sd_diffs = [((x - p_mean) / p_stdev) for x in data]
                out_sd_diffs = "\t".join(
                    ["{0:.4f}".format(float(x)) for x in sd_diffs])
                
                # Make plot. 
                x_label = "SDs around Mean"
                prefix = "SD"
                make_distplot(sd_diffs, out_f, se_id, x_label, prefix)

            new_line = ("\t".join(line) + "\t" + "{0:.4f}".format(
                float(p_mean)) + "\t" + "{0:.4f}".format(float(p_stdev)) +
                "\t" + "{0:.4f}".format(float(p_median)) + "\t" +
                "{0:.4f}".format(float(p_mad)) + "\t" + out_sd_diffs +
                "\t" + out_mad_diffs)

            print(new_line, file=output)
            
    output.close()


if __name__ == '__main__':
    main(inp_f, out_f)

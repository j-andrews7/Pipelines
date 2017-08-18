#!/usr/bin/env python3
"""
For a given bed-like table of QN'd SE signals, determine which of those called in single samples
are actually much different from the mean signal across all the samples in the file. Outputs a new
table with an additional column specifying the number of SDs from the mean for each SE called in
only one sample. Those that don't meet the user-sepcified cutoff can be excluded from output if
wanted. It will also create a distribution plot for the SEs in the output file.

Usage: python3 filter_single_SEs.py -i <input.bed> -o <output.bed> [OPTIONS]

Args:
    (required) -i <input.bed> = Name of bed-like table with QN'd SE signal data starting in the 6th
        column. Requires the samples in which the SE was called in a semi-colon delimited list in
        the 5th column.
    (required) -o <output.bed> = Name of output file to be created.
    (optional) -e = Option that when used, will not print lines that don't meet the specified
        cutoff. default: False.
    (optional) -t <threshold> = Numeric value for SDs above/below mean to use with exclude option.
        default: 2.0
"""

import sys
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from statistics import mean
from statistics import pstdev


def get_mean_stdev(data):
    """
    Given a list of data, return the mean and standard deviation for the population.
    """

    pop_stdev = pstdev(data)
    pop_mean = mean(data)

    return pop_mean, pop_stdev


def get_sample_idx(sample, header):
    """
    For a given sample, return the column index of the sample in the header.
    """

    for item in header:
        if sample in item:
            return header.index(item)

    print(sample + " not found in header, check input files.")
    sys.exit()


def make_rugplot(data, output_f):
    """
    For a 1-d array of data, create a rugplot showing where the data lie.
    """

    plt.figure(figsize=(8, 10), dpi=1200)
    displot = sns.distplot(data, hist=False, rug=True, color="r")
    out_name = output_f.split(".")[0] + ".png"
    plt.xlabel('Average SDs Above/Below Mean')
    plt.ylabel('Density')
    displot.figure.savefig(out_name)


def main():
    """
    Iterate through input file, determine those that meet the threshold, print to output file.
    """

    sns.set(style="white", palette="muted", color_codes=True)

    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="input_f", required=True)
    parser.add_argument("-o", "--output", dest="output_f", required=True)
    parser.add_argument(
        "-e", "--exclude", dest="exclude", action="store_true", required=False)
    parser.add_argument(
        "-t", "--threshold", dest="threshold", required=False, default=2.0, type=float)

    args = parser.parse_args()

    # Easier to use argument variables
    inp_f = args.input_f
    out_f = args.output_f
    exclude = args.exclude
    thresh = args.threshold

    # Catch threshold of 0.
    if thresh == 0:
        print(
            "Threshold cannot be set to 0. A value greater than or less than 0 must be used.")
        sys.exit()

    output = open(out_f, "w")

    with open(inp_f) as f:

        header = f.readline().strip().split("\t")
        new_header = "\t".join(header) + "\t" + "MEAN" + "\t" + "SD" + "\t" + "SDs_FROM_MEAN"
        print(new_header, file=output)

        avg_diff = []  # Used for plotting distribution of avg samp diff of each SE.

        for line in f:
            line = line.strip().split("\t")
            data = [float(x) for x in line[5:]]  # Assumes actual data starts in 6th column.
            samples = line[4].split(";")

            if len(samples) == 1:
                samp = samples[0]
                samp_idx = get_sample_idx(samp, header)
                samp_data = float(line[samp_idx])
                p_mean, p_stdev = get_mean_stdev(data)

                samp_diff = (samp_data - p_mean) / p_stdev

                # Exclude lines that don't meet the threshold if option is set.
                if exclude:
                    if thresh < 0:
                        if samp_diff < thresh:
                            avg_diff.append(samp_diff)
                            print(
                                "\t".join(line) + "\t" + "{0:.4f}".format(float(p_mean)) + "\t"
                                + "{0:.4f}".format(float(p_stdev)) + "\t"
                                + "{0:.4f}".format(float(samp_diff)), file=output)
                        else:
                            continue
                    if thresh > 0:
                        if samp_diff > thresh:
                            avg_diff.append(samp_diff)
                            print(
                                "\t".join(line) + "\t" + "{0:.4f}".format(float(p_mean)) + "\t"
                                + "{0:.4f}".format(float(p_stdev)) + "\t"
                                + "{0:.4f}".format(float(samp_diff)), file=output)
                        else:
                            continue
                else:
                    avg_diff.append(samp_diff)
                    print(
                        "\t".join(line) + "\t" + "{0:.4f}".format(float(p_mean)) + "\t"
                        + "{0:.4f}".format(float(p_stdev)) + "\t"
                        + "{0:.4f}".format(float(samp_diff)), file=output)

            else:
                avg_samp_diff = 0  # Used to hold running total of sample SD differences.
                p_mean, p_stdev = get_mean_stdev(data)

                for s in samples:
                    samp_idx = get_sample_idx(s, header)
                    samp_data = float(line[samp_idx])
                    samp_diff = (samp_data - p_mean) / p_stdev
                    avg_samp_diff += samp_diff

                avg_samp_diff = avg_samp_diff / len(samples)

                if exclude:
                    if thresh < 0:
                        if avg_samp_diff < thresh:
                            avg_diff.append(avg_samp_diff)
                            print(
                                "\t".join(line) + "\t" + "{0:.4f}".format(float(p_mean)) + "\t"
                                + "{0:.4f}".format(float(p_stdev)) + "\t"
                                + "{0:.4f}".format(float(avg_samp_diff)), file=output)
                        else:
                            continue
                    if thresh > 0:
                        if avg_samp_diff > thresh:
                            avg_diff.append(avg_samp_diff)
                            print(
                                "\t".join(line) + "\t" + "{0:.4f}".format(float(p_mean)) + "\t"
                                + "{0:.4f}".format(float(p_stdev)) + "\t"
                                + "{0:.4f}".format(float(avg_samp_diff)), file=output)
                        else:
                            continue
                else:
                    avg_diff.append(avg_samp_diff)
                    print(
                        "\t".join(line) + "\t" + "{0:.4f}".format(float(p_mean)) + "\t"
                        + "{0:.4f}".format(float(p_stdev)) + "\t"
                        + "{0:.4f}".format(float(avg_samp_diff)), file=output)

        make_rugplot(avg_diff, out_f)

    output.close()


if __name__ == '__main__':
    main()

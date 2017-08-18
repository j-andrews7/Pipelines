#!/usr/bin/env python3
"""
Given a bed file of binned genome positions and columns for each sample with a value indicating 
overlap of the SE with that bin, plot a heatmap-like figure for the entire genome. Playing with bin 
size may help the aesthetics of the image. The row labels are left off by default, but they can be
included so the user can know the start/end of each chromosome.

Usage: python3 plot_se_bins.py <input.bed> 

Args:
    input.bed = The input file with the data to be plotted.
"""

import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

sns.set()

####-Functions-####


def MakeLabel(pos_list):
    """
    Given a list of [chrom, start, end], return it as a formatted string (chr:start-end).
    """

    # Ensure proper data.
    if len(pos_list) != 3:
        print("Wrong number of position data being passed. Unable to create row label. Fix it.")
        sys.exit()

    label = pos_list[0] + ":" + pos_list[1] + "-" + pos_list[2]  # Create the label.

    return label


def ParseFile(input_file):
    """
    Parse the input file to get a list of row labels and a list of dicts, 
    each with column header + column data.
    """

    row_labels = []
    data_dict = {}

    # Open input file.
    with open(input_file) as f:
        header = f.readline().strip().split()
        samples = [x.split("_")[0] for x in header[3:]]  # Get the samples from the header.

        # Initialize entry in dict for each column.
        for samp in samples:
            data_dict[samp] = []

        for line in f:
            line = line.strip().split()
            label = MakeLabel(line[0:3])  # Make label from the first three elements of the line.

            i = 0  # Used as a counter to keep index for adding data to appropriate dicts.

            # Add data to appropriate dict.
            for item in line[3:]:  # Ignore lines with no data.
                if item == ".":
                    break
                elif i == 0:
                    row_labels.append(label)  # Append the label to the list holding them.
                samp = samples[i]  # Grab the corresponding sample.
                values = data_dict[samp]  # Get the list of values for the sample.
                values.append(float(item))  # Add the new data point.
                data_dict[samp] = values  # Reassign updated values list to the sample in the dict.
                i += 1  # Increment counter by one.

    return row_labels, data_dict


####-Main-####

def main():
    """
    Plot the data using the row_labels as, you guessed it, row labels and the data_dict as columns.
    The plot will be saved as input_file.png with the previous extension removed.
    """

    input_file = sys.argv[1]

    row_labels, data_dict = ParseFile(input_file)  # Get the data needed.

    df = pd.DataFrame(data_dict, index=row_labels)  # Create the dataframe.
    # EDIT THIS TO CHANGE FIGURE SIZE.
    plt.figure(figsize=(8, 11), dpi=1200)
    # Set colors [-5 to 5]. Can use html hex codes, recognized html colors, or rgb triplets.
    colors = ['#8c510a', "#bf812d", "#f5f5f5",
              "#f5f5f5", "#80cdc1", "#01665e"]
    cmap = ListedColormap(colors, name="cmap", N=6)  # Change N if you have a greater range.

    # Set colors for over/under the bound limits.
    cmap.set_over("#003c30")
    cmap.set_under("#543005")
    bounds = [-20, -10, -3, 0, 3, 10, 20]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    # Create the plot without y axis labels. Change 'extend' to 'both' or 'min' to make cbar extend
    # in opposite direction.
    heatmap = sns.heatmap(df, cbar=True, cbar_kws={'extend':'max'} ,cmap=cmap, norm=norm, yticklabels=False)
    plt.xticks(rotation=90)
    plt.title(input_file.split(".")[0])

    out_name = input_file.split(".")[0] + ".pdf"  # EDIT extension to change output format.
    heatmap.figure.savefig(out_name)


# Actually run.
if __name__ == '__main__':
    main()

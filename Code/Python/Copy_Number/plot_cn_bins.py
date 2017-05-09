#!/usr/bin/env python3
"""
Given a bed file of binned genome positions and columns for each sample with a '1' indicating overlap
of the CN with that bin, plot a heatmap-like figure for the entire genome. Playing with bin size may
help the aesthetics of the image. The row labels will almost certainly have to be redone in post,
they are simply included so the user can know the start/end of each chromosome.

Usage: python3 plot_cn_bins.py <input.bed>

Args:
    input.bed = The input file with the data to be plotted.
"""

import sys
import pandas as pd
import matplotlib as mpl
# matplotlib.use('SVG')  # I like PNG ('Agg') or SVG ('SVG') formats.
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

sns.set()

input_file = sys.argv[1]

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
	Parse the input file to get a list of row labels and a list of dicts, each with column header + column data.
	"""

	row_labels = []
	data_dict = {}

	# Open input file.
	with open(input_file) as f:
		header = f.readline().strip().split()
		samples = header[3:]  # Get the samples from the header.

		# Initialize entry in dict for each column.
		for samp in samples:
			data_dict[samp] = []

		for line in f:
			line = line.strip().split()
			label = MakeLabel(line[0:3])  # Make label from the first three elements of the line.
			row_labels.append(label)  # Append the label to the list holding them.

			i = 0  # Used as a counter to keep index for adding data to appropriate dicts.

			# Add data to appropriate dict.
			for item in line[3:]:
				samp = samples[i]  # Grab the corresponding sample.
				values = data_dict[samp]  # Get the list of values for the sample.
				values.append(int(item))  # Add the new data point.
				data_dict[samp] = values  # Reassign updated values list to the sample in the dict.
				i += 1  # Increment counter by one.

	return row_labels, data_dict


####-Main-####

def main(input_file):
	"""
	Plot the data using the row_labels as, you guessed it, row labels and the data_dict as columns.
	The plot will be saved as input_file.png with the previous extension removed.
	"""

	row_labels, data_dict = ParseFile(input_file)  # Get the data needed.

	df = pd.DataFrame(data_dict, index=row_labels)  # Create the dataframe.
	# EDIT figure size here if wanted.
	plt.figure(figsize=(8, 11), dpi=1200)
	colors = ['black','#8c510a', "#bf812d",
			  "#f5f5f5", "#80cdc1", "#01665e"]  
	# Set colors [-1 to 4]. To get a discrete values for each, use the color for the central value twice.
	# Can use html hex codes, recognized html colors, or rgb triplets.

	cmap = ListedColormap(colors, name="cmap", N=6)  # Change N if you have a greater range.
	bounds = [-1, 0, 1, 2, 3, 4, 5]  # Change this if your number ranges differ.
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	heatmap = sns.heatmap(df, cmap=cmap, cbar=True, yticklabels=False, norm=norm)  # Create the plot without a color bar or y axis labels.
	plt.xticks(rotation=90)
	out_name = input_file.split(".")[0] + ".pdf"  # Can just change the extension here to change output format.
	heatmap.figure.savefig(out_name)


# Actually run.
if __name__ == '__main__':
	main(input_file)

#!/usr/bin/env python3
"""
05/17/2016
jared.andrews07@gmail.com
-------------------------

Given a del or amp matrix, sum the data columns and make a histogram from the counts. Outputs the sum as an extra column in a new file.
Important to remember that the histogram's bins upper bounds are exclusive (i.e., a bin from 2-3 will only include values that are 2-2.99, NOT 3). 

Usage: python3 hist_from_cn_matrix.py <matrix.bed> <output.bed> 

Args:
    matrix.bed = Name of the input file to sum.
    output.bed = Positions with sum.
"""

import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

sns.set()

in_file = sys.argv[1]
out_file = sys.argv[2]

out = open(out_file, "w")

sums = []  # List to hold sum for each line.

with open(in_file) as f:
	# Get header and print to output.
	header = f.readline().strip() + "\t" + "SUM"
	print(header, file=out)

	for line in f:
		line = line.strip()
		line_list = line.split()
		pos = line_list[0:3]
		data = [1 if abs(int(x)) > 1 else int(x) for x in line_list[3:]]  # Converts anything greater than 1 to 1 and grabs the data for the line.
		data_sum = abs(sum(data))
		new_line = line + "\t" + str(data_sum)
		print(new_line, file=out)
		sums.append(data_sum)

print("Max value: " + str(max(sums)))
print("Number of bins: " + str((max(sums) + 1)))
bins = list(range(0, (max(sums) + 1), 1))  # Set number of bins from a range of 0 to the max number + 1, using a step of 1.

out.close()

fig_out = out_file.split(".")[0]

plt.figure(figsize=(6,8), dpi=1200)
x = pd.Series(sums, name="Recurrence of Gained CN Across FL/DL Samples")
binwidth=1
histplot = sns.distplot(x, color="#8c510a", bins=bins, kde=False, norm_hist=True)  # Create the histogram. #8c510a for dels. #01665e for amps.
out_name = out_file.split(".")[0] + ".png"
histplot.figure.savefig(out_name)
#!/usr/bin/env python3
"""
Quantile normalizes data using R. But without the idiotic, wonky output of actually using R. 

Usage: python3 quantile_normalize.py <input file> <data start column (1-based)>
"""


#Necessary packages
import sys
#Needed for quantile normalization
import readline
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import argparse
import numpy as np

#R package needed
preprocessCore = importr('preprocessCore')

#Get input file
input_file = sys.argv[1]
data_start = int(sys.argv[2]) - 1
output_file = "QN_" + input_file

def Quantile_Normalize(input_file, data_start):
	"""
	Take an input file, parse each line up to the data_start column and add those position elements
	to a list as a string. Take the elements of each line from data_start to end and add to an array,
	using headers to keep track of where to add each element to array. Quantile normalizes final array
	and returns both the list of positions and quantile normalized numpy array.

	Args:
		input_file = The input file to quantile normalize.
		data_start = Index of column in which actual data to be normalized starts.

	Returns:
		header = Header of output file.
		pos_list = List of positions for each line.
		norm_matrix = Quantile normalized matrix of data.
	"""

	#Open input file
	with open(input_file) as f:

		print("Creating data matrix, may take a few minutes.")

		#Get header and print to output
		header = f.readline().strip()

		#Determine number of samples in file
		samples = header.strip().split("\t")[data_start:]

		#Initialize list to hold all the other lists
		pos_list=[]
		sample_list = []

		#Add appropriate number of lists to master list
		for item in samples:
			sample_list.append([]) 

		# Debug
		chroms = []


		#Iterate through file and store each column in a list
		for line in f:

			#Used to keep track of data index later
			count=0

			line = line.strip().split("\t")
			if line[0] not in chroms:
				print(line[0])
				chroms.append(line[0])
			position = line[0:data_start]
			pos_list.append("\t".join(position))
			data=line[data_start:]

			#Add data to appropriate list
			for entry in data:
				#Add pseudocount
				sample_list[count].append(float(entry) + 0.1)
				count += 1

		print("Converting to R matrix.")

		#Actually do the QN
		matrix = sample_list
		del sample_list
		v = robjects.FloatVector([ element for col in matrix for element in col ])
		m = robjects.r['matrix'](v, ncol = len(matrix), byrow=False)
		print("Performing quantile normalization.")
		Rnormalized_matrix = preprocessCore.normalize_quantiles(m)
		norm_matrix = np.array(Rnormalized_matrix)

		return header, pos_list, norm_matrix

		# m_count = 1
		# m = robjects.r['matrix'](0.0, ncol=len(matrix), nrow=len(matrix[0]))
		# for samp in matrix:
		# 	i_count = 1
		# 	for entry in samp:
		# 		m.rx[i_count, m_count] = entry
		# 		i_count += 1
		# 	m_count += 1

		# print("Performing quantile normalization.")

		# Rnormalized_matrix = preprocessCore.normalize_quantiles(m)
		# norm_matrix = np.array(Rnormalized_matrix)

		# return header, pos_list, norm_matrix


####-MAIN-####

def main(input_file, output_file, data_start):
	"""
	Quantile normalize input file using data from data_start columns to end of line.
	Print a quantile normalize copy of the input file.

	Args:
		input_file = Path to file to quantile normalize.
		data_start = Column in which data to quantile normalize begins in the input file.
		output_file = Name of file to write to output.
	"""

	# Quantile normalize the data.
	header, pos_list, normalized_matrix = Quantile_Normalize(input_file, data_start)
	
	print("Creating output file...")


	output = open(output_file, "w")

	print(header, file=output)

	# Iterate through each position and print QN'd data
	count = 0
	for thing in pos_list:
		thing_index = count
		count += 1

		# Print actual data
		norm_data = normalized_matrix[thing_index]
		out_data = "\t".join("{0:.2f}".format(piece) for piece in norm_data)
			
		print(thing + "\t" + out_data, file=output)

	output.close()

# Run if ran as script rather than imported.
if __name__ == '__main__':
	main(input_file, output_file, data_start)

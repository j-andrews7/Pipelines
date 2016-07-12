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
		temp_out = The name of the temporary output file containing the quantile normalized matrix.
	"""

	temp_out = input_file + ".tmp"
	# temp_out = '"' + temp_out + '"'

	# Open input file
	with open(input_file) as f:

		# input_file = '"' + input_file + '"' # Formatting for input to R function.
		print("Getting header and creating position list from " + input_file)

		# Get header and print to output
		header = f.readline().strip()
		header_l = header.split()

		# Determine number of samples in file
		samples = header.strip().split("\t")[data_start:]
		ncols = len(header_l)
		print("ncols = " + str(ncols))

		# Initialize list to hold all the other lists
		pos_list=[]

		chroms = []  # For debugging.

		# Iterate through file and store each column in a list
		for line in f:
			line = line.strip().split("\t")
			if line[0] not in chroms:
				print("Pulling " + line[0])
				chroms.append(line[0])
			position = line[0:data_start]
			pos_list.append("\t".join(position))

		data_start = data_start + 1  # R uses 1-based columns.

		qn = robjects.r('''
		function(input_file,temp_out,data_start,ncols) {
			require('preprocessCore')
			print("Reading file in as matrix.")
			all <- data.matrix(read.table(input_file,sep="\t",header=TRUE));
			print("Adding pseudocounts.")
			all <- all[,data_start:ncols] + 0.1
			print("Quantile normalizing.")
			all=normalize.quantiles(all);
			print("Printing to temporary output.")
			write.table(all,temp_out,sep="\t",row.names=FALSE);
		}
		''')

		print("R quantile normalization function defined.")
		print("Reading file in as R matrix and quantile normalizing.")

		qn(input_file,temp_out,data_start,ncols)  # Run the R function.

		return header, pos_list, temp_out


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
	header, pos_list, temp_out = Quantile_Normalize(input_file, data_start)
	
	print("Creating output file...")


	output = open(output_file, "w")

	print(header, file=output)

	# Iterate through each position and print QN'd data
	with open(temp_out.strip('"')) as f:
		f.readline() # Skip screwed up header.

		count = 0

		for line in f:
			line = line.strip().split("\t")
			out_data = "\t".join("{0:.2f}".format(float(piece)) for piece in line)

			
			pos = pos_list[count]  # Use a count to keep track of position list index.
			count += 1

			# Print actual data	
			print(pos + "\t" + out_data, file=output)

	output.close()

# Run if ran as script rather than imported.
if __name__ == '__main__':
	main(input_file, output_file, data_start)

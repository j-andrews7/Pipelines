#!/usr/bin/env python3
"""
Used to filter a file of super enhancers based on numerous criteria. Utilizes log2 microarray expression data for associated genes and ChIP data for each loci.
User can set the -min and -max number of genes they require to be associated with a loci for it to pass. 
Loci that don't have the wanted number of associated genes will not be included in the resulting output file. 

A key file is used to demarcate sets of samples to compare to each other. Each sample set should contain one sample per line of the key file, delimited by a tab.
Samples should be named without the mark (e.g. CC012314, not CC012314_K27AC). It is not necessary for there to be equal numbers of samples in each set.

Columns of input files must be tab-delimited.

Column layout of expression file must be <GENE_SYM> <data columns> 
Note: A header is required to identify sample columns. Must contain the sample name, other info (mark) doesn't matter.

Column layout of key file must be <dataset1 sample name> <dataset2 sample name> 
Note: A header is required to name the two datasets. These will be used to create some column headers in the output files, so be descriptive. 

Usage: python3 compare_expression.py -exp <microarray expression file> -key <key file> -o <output file> [OPTIONS]

Args:
    (required) -exp <expression file> = Name of expression data file to process
    (required) -key <key file> = Name of key file in format as described above
    (required) -o <output file> = Name of output file to be created
    (optional) -fc <Fold change cutoff for expression data, default=2.0>
"""

#Necessary packages
import sys
import argparse
import numpy as np
from scipy import stats
import math
import random


####-Functions-####
def Parse_Key(key_file):
	"""
	Parse the key file to determine the two sample sets to compare to each other.

	Parameters:
		key_file = The key file to be processed.

	Returns:
		set1_name = String to be used as the name of the first dataset.
		set2_name = String to be used as the name of the second dataset.
		samp_list1 = A list containing the names of all samples of the first data set.
		samp_list2 = A list containing the names of all samples of the second data set.
	"""

	#Initialize the two lists
	samp_list1 = []
	samp_list2 = []

	#Open the key file
	with open(key_file) as f:

		#Get header and store dataset names
		header = f.readline()
		header = header.strip().split("\t")

		#Check if the header has the correct number of elements, and if not, yell at user and quit the script. If yes, store the column headers as the data set names.
		if len(header) == 2:
			set1_name = header[0]
			set2_name = header[1]
		else:
			print("Key file header has wrong number of headers. Check key file - must only have two columns.")
			sys.exit()

		#Iterate through file and add samples in each column to the appropriate list
		for line in f:

			#Split line by tab
			line = line.rstrip().split("\t")

			#Check if the line has multiple elements
			if len(line) == 2:

				#Check if any elements are empty
				if line[0] == "":
					samp_list2.append(line[1])
				elif line[1] == "":
					samp_list1.append(line[0])
				else:
					samp_list1.append(line[0])
					samp_list2.append(line[1])

			#If only one element, then it must go to the first sample set
			else:
				samp_list1.append(line[0])

	return (set1_name, set2_name, samp_list1, samp_list2)



def Calc_FC_Genes(dataset1, dataset2):
	"""
	Calculates the log2 fold change between two lists of values. Assumes data is already in log2.

	Parameters:
		(required) dataset1 = List of values for 1st dataset.
		(required) dataset2 = List of values for 2nd dataset.

	Returns:
		FC_val = Log2 fold change value between the two datasets.
	"""

	#Average each dataset
	dataset1_avg = (sum(dataset1)/len(dataset1))
	dataset2_avg = (sum(dataset2)/len(dataset2))

	#Get FC value
	FC_val = dataset1_avg - dataset2_avg

	#Return the foldchange value
	return round(FC_val,4)


# TO DO - Change to take a list of lists, each pertaining to a sample set. Return a list of lists.
def Get_Data_Cols(set1, set2, header):
	"""
	For two datasets, determines which columns of a header correspond to each sample and sticks the index of the columns into lists.

	Parameters:
		(required) set1 = List containing all sample names in first dataset.
		(required) set2 = List containing all sample names in second dataset.
		(required) header = The header of a file, split as appropriate into a list.

	Returns:
		dataset1_cols = List containing the indexes of columns that contain data for samples in dataset1.
		dataset2_cols = List containing the indexes of columns that contain data for samples in dataset2.
	"""

	#Initialize lists to hold the indices of the data columns for each dataset.
	dataset1_cols = []
	dataset2_cols = []

	#Assign data columns to the correct dataset by index.
	for item in set1:

		#If found, add the index of the column to the set1_cols list
		for entry in header:
			if entry.startswith(item):
				col_ind = header.index(entry)
				dataset1_cols.append(col_ind)

		#If no data found, tell user
		if len(dataset1_cols) == 0:
			print("Data not found for " + item + ". If you think data should be found, check the key file to be sure the sample name is correct.")

	if len(dataset1_cols) == 0:
			print("No data found for first set of samples. Check expression and key file.")
			sys.exit()

	for item in set2:

		#If found, add the index of the column to the set1_cols list
		for entry in header:
			if entry.startswith(item):
				col_ind = header.index(entry)
				dataset2_cols.append(col_ind)

		#If data not found, tell user.
		if len(dataset2_cols) == 0:
			print("Data not found for " + item + ". If you think data should be found, check the key file to be sure the sample name is correct.")

	if len(dataset2_cols) == 0:
			print("No data found for second set of samples. Check expression and key file.")
			sys.exit()

	#Return results
	return dataset1_cols, dataset2_cols



def Get_Diff_Genes(exp_file, FC_cutoff, name1, name2, set1, set2):
	"""
	Finds genes from the microarray data that have differential expression between set1 and set2 based on the FC_cutoff.

	Parameters:
		(required) exp_file = Path to the expression data file.
		(required) FC_cutoff = Float corresponding to the fold change cutoff to be used as the definition of differential expression.
		(required) name1 = Name of the first dataset.
		(required) name2 = Name of the second dataset.
		(required) set1 = List of the samples that are part of set1.
		(required) set2 = List of the samples that are part of set2.

	Returns:
		gene_list = A list of tuples, each containing a gene that's differentially expressed and the log2 fold change value for it.
	"""

	print("Identifying differentially expressed genes between " + name1 + " and " + name2 + " datasets using a fold change cutoff of " + str(FC_cutoff))

	#Initialize gene list
	gene_list = []

	#Initialize lists to determine which columns belong to which dataset
	set1_cols = []
	set2_cols = []

	#Convert the FC_cutoff value to log2.
	FC_log = math.log2(FC_cutoff)

	#Open expression file and get header
	with open(exp_file) as f:

		header=f.readline().strip().split("\t")

		#Check if the samples of each set are present in the header or not. If not, tell the user that their data is missing.
		set1_cols, set2_cols = Get_Data_Cols(set1, set2, header)

		#Iterate through file and compare the expression of each gene between the two sample sets.
		for line in f:

			#Initialize lists to hold the data for each set.
			data1 = []
			data2 = []

			#Split the line by tab
			line = line.strip().split("\t")

			#Store gene symbol
			gene = line[0]

			#Add the appropriate data to each data list
			for entry in set1_cols:
				data1.append(float(line[entry]))

			for entry in set2_cols:
				data2.append(float(line[entry]))

			#Get the log2 fold change values between each dataset.
			diff_exp = Calc_FC_Genes(data1, data2)

			#Add to the gene_list if the differential expression value meets the user specified FC cutoff
			if diff_exp >= FC_log or diff_exp <= -FC_log:
				gene_list.append((gene,diff_exp))

	#Feedback is good
	print(str(len(gene_list)) + " genes meeting the fold change cutoff of " + str(FC_cutoff))

	#Return results
	return gene_list


####-Variables-####
#Create argument parser
parser = argparse.ArgumentParser(usage=__doc__)

#Add accepted arguments
parser.add_argument("-exp", "--expression", dest = "exp_file", required=True)
parser.add_argument("-key", "--key", dest = "key_file", required=True)
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-fc", "--foldchange", dest ="fc", type=float, default=2)

#Grab arguments from the parser
args = parser.parse_args()


#Output so user can double check options
print( "Expression file: {}\nKey file: {}\nOutput file: {}\nExpression FC cutoff: {} \n \n".format(
        args.exp_file,
        args.key_file,
        args.output_file,
        args.fc,
        ))

#Assign arguments to more straightforward variables
exp_file = args.exp_file
key_file = args.key_file
output_file = args.output_file
loci_fc = args.fc



#####--MAIN--#####

#Parse key file to get sample lists.
dataset1_name, dataset2_name, sample_list1, sample_list2 = Parse_Key(key_file)

#Grab all genes differentially expressed between the two datasets from the expression file.
diff_genes = Get_Diff_Genes(exp_file, gene_fc, dataset1_name, dataset2_name, sample_list1, sample_list2)

# TO DO - Hard code group combinations to get all pairwise comparisons
#!/usr/bin/env python3
"""
Used to filter a file of super enhancers based on numerous criteria. Utilizes log2 microarray expression data for associated genes and ChIP data for each loci.
User can set the -min and -max number of genes they require to be associated with a loci for it to pass. 
Loci that don't have the wanted number of associated genes will not be included in the resulting output file. 

The user should use alternative methods to get genes associated with each loci and place them in the last column of the loci
file (see get_genes_in_range.py script). 

A key file is used to demarcate two sets of samples to compare to each other. Each sample set should contain one sample per line of the key file, delimited by a tab.
The first column will be used as one sample set, the second column as the other sample set. Samples should be named without the mark (e.g. CC012314, not CC012314_K27AC).
It is not necessary for there to be equal numbers of samples in each set.

Fold change calculations are run on the ChIP data to determine which loci differ between the two sample sets. The user may adjust the FC 
cutoff with the -sfc option (default is magnitude of fold change 2). The ChIP data of the loci for each sample set will be averaged, converted to log2, and 
then compared to each other to calculate fold change. A Welch's t-test will also be performed for each sample set to determine loci that differ significantly.
The p-value for this test can be changed using the -p option (default 0.05). 
Those that do not reach the FC or p-value cutoffs will not be printed to the output file.

Additionally, a tag cutoff can be used to further filter results (-t to use, False by default.) The -tcut option allows user to set the cutoff. 

The -u option (--unique) will only print the loci that reach the tag cutoff in one dataset and not the other. 

A "dense" output file is also created (output_name + "dense.txt") containing only the data for the two datasets.

Columns of input files must be tab-delimited.

Column layout of loci file must be <loci_ID, chr, start, and end in first four columns> <data columns> <associated genes, delimited by ";">.
Note: A header is required to identify sample columns. Must contain the sample name, other info doesn't matter.

Column layout of expression file must be <GENE_SYM> <data columns> 
Note: A header is required to identify sample columns. Must contain the sample name, other info (mark) doesn't matter.

Column layout of key file must be <dataset1 sample name> <dataset2 sample name> 
Note: A header is required to name the two datasets. These will be used to create some column headers in the output files, so be descriptive. 

Usage: python3 Analyze_Loci.py -i <ChIP data file> -exp <microarray expression file> -key <key file> -o <output file> [OPTIONS]

Args:
    (required) -i <input file> = Name of input file (ChIP data) to process 
    (required) -exp <expression file> = Name of expression data file to process
    (required) -key <key file> = Name of key file in format as described above
    (required) -o <output file> = Name of output file to be created
    (optional) -min <minimum number of genes that can be associated with each loci for it to pass, default=2>
    (optional) -max <maximum number of genes that can be associated with each loci for it to pass, default=7>
    (optional) -p <p-value for Welch's t-test to consider significant, default=0.05>
    (optional) -lfc <Fold change cutoff for ChIP data from loci, default=2.0>
    (optional) -gfc <Fold change cutoff for ChIP data from loci, default=2.0>
    (optional) -t <If used, will use tcut value as a tag cutoff for the ChIP data, requiring at least one sample in the key file
    			to have more tags for the loci to be included in the result, default=False>
    (optional) -tcut <Value to define the tag cutoff to be used for -t and -u options, default = 1500>
    (optional) -u <If used, will only output loci that meet the tag cutoff in one dataset or another>
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



def Calc_FC_Loci(dataset1, dataset2):
	"""
	Calculates the log2 fold change between two lists of values. Changes data to log2.

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
	FC_val = math.log2(dataset1_avg) - math.log2(dataset2_avg)

	#Return the foldchange value
	return FC_val



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


def Check_Gene_FC(gene, genelist):
	"""
	Checks if gene is located within the genelist containing all of the differentially expressed genes. If so, returns the value.

	Parameters:
		(required) gene = Gene symbol to check.
		(required) genelist = List of differentially expressed genes between the two samples with the corresponding log2 fold change value.

	Returns:
		gene_FC = Log2 fold change value for the gene.
	"""

	#Iterate through the gene list and check if gene is in the list
	for entry in genelist:
		if gene == entry[0]:

			#Grab FC value
			gene_FC = entry[1]

			#Return the results
			return gene_FC

		#Returns fold change of 1 to account for loci with no associated genes (NA value)
		if gene == "NA":

			gene_FC = 1

			return gene_FC

	return None


def Check_Tags(cutoff, data1_l, data2_l, unique):
	"""
	Checks if one of the samples in the datasets meets the tag cutoff. If the unique option is true, makes sure only 
	one of the datasets meets the cutoff.

	Parameters:
		(required) cutoff = Value to be used as the cutoff.
		(required) data1_l = List of tags for first data set.
		(required) data2_l = List of tags for second data set.
		(required) unique = Boolean to determine whether to check for only one of the datasets meeting the cutoff.

	Returns:
		result = Boolean to determine whether the cutoff was met or not.
	"""

	#Bools to determine if each dataset meets the cutoff
	data1_meetcut = False
	data2_meetcut = False

	#Check if each list has a value that meets the cutoff
	for samp in data1_l:

		if samp >= cutoff:
			data1_meetcut = True
			break

	for samp in data2_l:

		if samp >= cutoff:
			data2_meetcut = True
			break

	#Check if only one dataset meeting cutoff is wanted or if it doesn't matter
	if unique:
		
		if (data1_meetcut == True and data2_meetcut == False) or (data2_meetcut == True and data1_meetcut == False):
			result = True
		else:
			result = False

		return result

	else:
		if data1_meetcut or data2_meetcut:
			result = True
		else:
			result = False

		return result


#This function is way too big, but I can't be assed to go break it up as I should. 
def Loci_Filter(input_file, p_val, loci_fc_cutoff, min_genes, max_genes, diff_geneslist, 
	set1_n, set2_n, samp_l1, samp_l2, output_file, tag=False, tag_cut = 1500, unique=False):
	"""
	Filters through the super loci file, only printing those that meet all criteria (associated genes, significance, magnitude) to the output file(s).

	Parameters:
		(required) input_file = File containing the loci data and associated genes.
		(required) p_val = The p-value to use as a cutoff for significance for Welch's t-test between sample sets. 
		(required) loci_fc_cutoff = Value to be used as the fold change cutoff between the loci data sets. 
		(required) min_genes = Minimum number of genes that must be associated with the loci.
		(required) max_genes = Maximum number of genes that can be associated with the loci.
		(required) diff_geneslist = List containing tuples with the differentially expressed genes and the log2 fold change value.
		(required) set1_n = Name of the first dataset.
		(required) set2_n = Name of the second dataset.
		(required) samp_l1 = List containing all sample names of the first dataset.
		(required) samp_l2 = List containing all sample names of the second dataset.
		(required) output_file = Name of the output file.
		(optional) tag = Boolean determining whether to use a tag cutoff to filter the loci.
		(optional) tag_cut = Tag cutoff to use, default=1500.
		(optional) unique = Boolean to determine whether to only report loci that reach the tag cutoff in one dataset or the other.
	"""

	#Initialize lists to determine which columns are for each dataset
	data1_cols = []
	data2_cols = []

	#Intialize counts to hold the number of loci that pass each filter.
	gene_num_count = 0
	gene_exp_count = 0
	loci_mag_count = 0
	loci_sig_count = 0
	loci_tag_count = 0
	total_loci_count = 0

	#Get log2 se_fc_cutoff
	log_fc_loci = math.log2(loci_fc_cutoff)

	#Let user know what's going on
	print("Filtering loci now. This may take a few minutes. Take a break, champ, you earned it.")

	#Important feedback list
	insults = ["ya dingus.", "come on, man, you had one job. One!", "you cretinous inbred.", "you big, beautiful mess.", "you fantastically incapable data jockey.", "you incredibly dense screwup."]

	#Open the se file and get the header
	with open(input_file) as f:

		#Strip newline character, split by tab
		header = f.readline().strip().split("\t")

		#Assign data columns to the correct dataset by index.
		data1_cols, data2_cols = Get_Data_Cols(samp_l1, samp_l2, header)

		#Open the output file(s) and print the header
		output = open(output_file, "w")

		#Print to normal output file after added wanted columns
		header.extend(["ASSOC_GENES", ("LOG2_FC_SIG_GENES_" + set1_n + "v" + set2_n), ("LOG2_FC_LOCI_"+ set1_n + "v" + set2_n), "LOCI_PVAL"])
		if unique:
			header.extend(["ASSOC_GENES", ("LOG2_FC_SIG_GENES_" + set1_n + "v" + set2_n), ("LOG2_FC_LOCI_"+ set1_n + "v" + set2_n), "LOCI_PVAL", "UNIQUE_TO"])
		print(*header[0:], sep="\t", file=output)

		#If so, create dense file name and open it.
		output_base = output_file.split(".")[0]
		output_dense_file = output_base + "_dense.txt"
		output_dense = open(output_dense_file, "w")

		#Get correct column headers for dense output
		header_dense = [header[0], header[1], header[2], header[3]]
		for item in data1_cols:
			header_dense.append(header[item])
		for item in data2_cols:
			header_dense.append(header[item])

		#Add new column headers
		header_dense.extend(["ASSOC_GENES", ("LOG2_FC_SIG_GENES_" + set1_n + "v" + set2_n), ("LOG2_FC_LOCI_"+ set1_n + "v" + set2_n), "LOCI_PVAL"])
		if unique:
			header_dense.extend(["ASSOC_GENES", ("LOG2_FC_SIG_GENES_" + set1_n + "v" + set2_n), ("LOG2_FC_LOCI_"+ set1_n + "v" + set2_n), "LOCI_PVAL", "UNIQUE_TO"])

		#Print to the dense file
		print(*header_dense[0:], sep="\t", file=output_dense)

		#Iterate through file line by line
		for line in f:

			#Increase total count 
			total_loci_count += 1

			#Split the line by tab to get data list
			data = line.strip().split("\t")

			#Initialize lists to hold data for the two datasets
			datalist1 = []
			datalist2 = []

			#Make a list for the genes and split by ";".
			genes = data[-1].split(";")

			#Take care of not real gene (NA value)
			if "NA" in genes:
				genes.remove("NA")


			#First check if there are the appropriate number of genes associated with the locus, and if not, go to next line.
			if len(genes) >= min_genes and len(genes) <= max_genes:

				#Increase gene_num_count
				gene_num_count += 1

				#Set a boolean to determine whether a differential gene is associated with the locus and a list to hold said associated genes.
				differ_exp = False
				locus_diff_genes = []

				#Check that at least one gene is differentially expressed
				for item in genes:

					gene_foldchange = Check_Gene_FC(item, diff_geneslist)

					#If the gene is found in the differential gene list, then set the differ_exp bool to True and add the gene and FC value to the SE_diff_genes list
					if gene_foldchange is not None:
						differ_exp = True
						locus_diff_genes.append((item, gene_foldchange))

				#If a differentially expressed gene was found, proceed. Otherwise go to next line.
				if differ_exp:

					#Increase gene_exp_count
					gene_exp_count += 1

					#Add the appropriate data to each data list
					for entry in data1_cols:
						datalist1.append(float(data[entry]))

					for entry in data2_cols:
						datalist2.append(float(data[entry]))

					#Check if only loci that reach the tag cutoff are wanted and skip as appropriate
					if tag or unique:
						if Check_Tags(tag_cut, datalist1, datalist2, unique):

							#Increase loci_tag_count
							loci_tag_count += 1

							#Get the log2 fold change values between each dataset.
							diff_levels = Calc_FC_Loci(datalist1, datalist2)

							#Continue if the FC meets the criteria, but otherwise 
							if diff_levels >= log_fc_loci or diff_levels <= -log_fc_loci:
								
								#Increase the loci_mag_count
								loci_mag_count += 1

								#Run Welch's t-test assuming population variances might not be equal
								locus_tval, locus_pval = stats.ttest_ind(datalist1, datalist2, equal_var=False)

								#Check if the p-value meets the cutoff for significance
								if locus_pval <= p_val:

									#Increase loci_sig_count
									loci_sig_count += 1

									#Print results to output file
									print(*data[0:], sep="\t", end="\t", file=output)
									print(*locus_diff_genes[0:], sep=";", end="\t", file=output)
									print("{0:.4f}".format(diff_levels), "{0:.4f}".format(locus_pval), sep="\t", file=output)

									#This is so, so terrible.
									print(*data[0:4], sep="\t", end="\t", file=output_dense)
									print(*datalist1[0:], sep="\t", end="\t", file=output_dense)
									print(*datalist2[0:], sep="\t", end="\t", file=output_dense)
									print(data[-1], end="\t", file=output_dense)
									print(*locus_diff_genes[0:], sep=";", end="\t", file=output_dense)
									print("{0:.4f}".format(diff_levels), "{0:.4f}".format(locus_pval), sep="\t", file=output_dense)

					else:

						#Get the log2 fold change values between each dataset.
						diff_levels = Calc_FC_Loci(datalist1, datalist2)

						#Continue if the FC meets the criteria, but otherwise 
						if diff_levels >= log_fc_loci or diff_levels <= -log_fc_loci:
							
							#Increase the loci_mag_count
							loci_mag_count += 1

							#Run Welch's t-test assuming population variances might not be equal
							locus_tval, locus_pval = stats.ttest_ind(datalist1, datalist2, equal_var=False)

							#Check if the p-value meets the cutoff for significance
							if locus_pval <= p_val:

								#Increase loci_sig_count
								loci_sig_count += 1

								#Print results to output file
								print(*data[0:], sep="\t", end="\t", file=output)
								print(*locus_diff_genes[0:], sep=";", end="\t", file=output)
								print("{0:.4f}".format(diff_levels), "{0:.4f}".format(locus_pval), sep="\t", file=output)

								#This is so, so terrible.
								print(*data[0:4], sep="\t", end="\t", file=output_dense)
								print(*datalist1[0:], sep="\t", end="\t", file=output_dense)
								print(*datalist2[0:], sep="\t", end="\t", file=output_dense)
								print(data[-1], end="\t", file=output_dense)
								print(*locus_diff_genes[0:], sep=";", end="\t", file=output_dense)
								print("{0:.4f}".format(diff_levels), "{0:.4f}".format(locus_pval), sep="\t", file=output_dense)
			
		#Close the output files	
		output_dense.close()
		output.close()

	#Print last bit of feedback for user.
	print("Of total " + str(total_loci_count) + " provided loci:")
	print(str(gene_num_count) + " loci passed the associated gene filters of " + str(min_genes) + " minimum and " + str(max_genes) + " maximum.\n")
	print(str(gene_exp_count) + " of those then passed the differential expression filter of " + str(gene_fc) + " of at least one associated gene.\n")
	if tag or unique:
		print(str(loci_tag_count) + " of those passed the tag cutoff of " + str(tag_cutoff) + ".\n")
	print(str(loci_mag_count) + " of those then passed the magnitude filter of " + str(loci_fc_cutoff) + " fold change between " + set1_n + " and " + set2_n + " sample averages.\n")
	print(str(loci_sig_count) + " of those then passed the significance filter of " + str(p_val) + " between " + set1_n + " and " + set2_n + "\n")


####-Variables-####
#Create argument parser
parser = argparse.ArgumentParser(usage=__doc__)

#Add accepted arguments
parser.add_argument("-i", "--input", dest = "input_file", required=True)
parser.add_argument("-exp", "--expression", dest = "exp_file", required=True)
parser.add_argument("-key", "--key", dest = "key_file", required=True)
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-min", "--min", dest ="gene_min", type=int, default=2)
parser.add_argument("-max", "--max", dest ="gene_max", type=int, default=7)
parser.add_argument("-p", "--pval", dest ="p_val", type=float, default=0.05)
parser.add_argument("-lfc", "--locifoldchange", dest ="loci_fc", type=float, default=2)
parser.add_argument("-gfc", "--genefoldchange", dest ="gene_fc", type=float, default=2)
parser.add_argument("-t", "--tag", action="store_true")
parser.add_argument("-tcut", "--tagcutoff", dest ="tag_cutoff", type=float, default=1500)
parser.add_argument("-u", "--unique", action="store_true")

#Grab arguments from the parser
args = parser.parse_args()


#Output so user can double check options
print( "Input file: {}\nExpression file: {}\nKey file: {}\nOutput file: {}\nGene minimum: {}\nGene maximum: {}\nLoci p-value cutoff: {}\nLoci FC cutoff: {}\nExpression FC cutoff: {} \n \n".format(
        args.input_file,
        args.exp_file,
        args.key_file,
        args.output_file,
        args.gene_min,
        args.gene_max,
        args.p_val,
        args.loci_fc,
        args.gene_fc,
        ))

#Assign arguments to more straightforward variables
input_file = args.input_file
exp_file = args.exp_file
key_file = args.key_file
output_file = args.output_file
gene_min = args.gene_min
gene_max = args.gene_max
p_val = args.p_val
loci_fc = args.loci_fc
gene_fc = args.gene_fc
tag = args.tag
tag_cutoff = args.tag_cutoff
unique = args.unique


#####--MAIN--#####

#Parse key file to get sample lists.
dataset1_name, dataset2_name, sample_list1, sample_list2 = Parse_Key(key_file)

#Grab all genes differentially expressed between the two datasets from the expression file.
diff_genes = Get_Diff_Genes(exp_file, gene_fc, dataset1_name, dataset2_name, sample_list1, sample_list2)

#Filter the SEs and print to output
Loci_Filter(input_file, p_val, loci_fc, gene_min, gene_max, 
	diff_genes, dataset1_name, dataset2_name, sample_list1, sample_list2, output_file, tag, tag_cutoff, unique)


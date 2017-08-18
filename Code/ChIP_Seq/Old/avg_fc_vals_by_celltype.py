#!/usr/bin/env python3
"""
Calculates the average value for samples of identical cell types and then performs all possible fold change comparisons. 
Can create gzipped bedgraph files from the FC data columns.

Usage: python3 avg_fc_vals_by_celltype.py -i <input.txt> -o <output.txt> -ucsc <Optional>

Args:
    (required) -i input.txt = Name of input file to process 
    (required) -o output.txt = Name of output file to be created
    (optional) -ucsc Will create gzipped bedgraph files for each FC column if included.

NOTE:Samples must be included in the script dictionary with the appropriate cell type for the script to work properly. Check code for more info.
"""

import sys
import math
from subprocess import call
import argparse
import gzip

parser = argparse.ArgumentParser(usage=__doc__)

####-Variables-####

#List to hold which columns end up with average data. Used for fold change logic.
filled_cols = []

#Initialize cell type dictionary. Can add additional samples as needed.
sample_data = {
	"CB011514_X": "CB",
	"CB011514_H3AC": "CB",
	"CB011514_K27AC": "CB",
	"CB011514_K4ME3": "CB",
	"CB011514_K4ME1": "CB",
	"CB011514_K27ME3": "CB",
	"CB012214_X": "CB",
	"CB012214_H3AC": "CB",
	"CB012214_K27AC": "CB",
	"CB012214_K4ME3": "CB",
	"CB012214_K4ME1": "CB",
	"CB012214_K27ME3": "CB",
	"CB012314_X": "CB",
	"CB012314_H3AC": "CB",
	"CB012314_K27AC": "CB",
	"CB012314_K4ME3": "CB",
	"CB012314_K4ME1": "CB",
	"CB012314_K27ME3": "CB",
	"CB020514_X": "CB",
	"CB020514_H3AC": "CB",
	"CB020514_K27AC": "CB",
	"CB020514_K4ME3": "CB",
	"CB020514_K4ME1": "CB",
	"CB020514_K27ME3": "CB",
	"CB021314_X": "CB",
	"CB021314_H3AC": "CB",
	"CB021314_K27AC": "CB",
	"CB021314_K4ME3": "CB",
	"CB021314_K4ME1": "CB",
	"CB021314_K27ME3": "CB",
	"CC011514_X": "CC",
	"CC011514_H3AC": "CC",
	"CC011514_K27AC": "CC",
	"CC011514_K4ME3": "CC",
	"CC011514_K4ME1": "CC",
	"CC011514_K27ME3": "CC",
	"CC012214_X": "CC",
	"CC012214_H3AC": "CC",
	"CC012214_K27AC": "CC",
	"CC012214_K4ME3": "CC",
	"CC012214_K4ME1": "CC",
	"CC012214_K27ME3": "CC",
	"CC012314_X": "CC",
	"CC012314_H3AC": "CC",
	"CC012314_K27AC": "CC",
	"CC012314_K4ME3": "CC",
	"CC012314_K4ME1": "CC",
	"CC012314_K27ME3": "CC",
	"CC020514_X": "CC",
	"CC020514_H3AC": "CC",
	"CC020514_K27AC": "CC",
	"CC020514_K4ME3": "CC",
	"CC020514_K4ME1": "CC",
	"CC020514_K27ME3": "CC",
	"CC021314_X": "CC",
	"CC021314_H3AC": "CC",
	"CC021314_K27AC": "CC",
	"CC021314_K4ME3": "CC",
	"CC021314_K4ME1": "CC",
	"CC021314_K27ME3": "CC",
	"TS081414_MEMORY_X" : "MEMORY",
	"TS081414_MEMORY_H3AC" : "MEMORY",
	"TS081414_MEMORY_K27AC" : "MEMORY",
	"TS081414_MEMORY_K4ME3" : "MEMORY",
	"TS081414_MEMORY_K4ME1" : "MEMORY",
	"TS081414_MEMORY_K27ME3" : "MEMORY",
	"TS102214_MEMORY_X" : "MEMORY",
	"TS102214_MEMORY_H3AC" : "MEMORY",
	"TS102214_MEMORY_K27AC" : "MEMORY",
	"TS102214_MEMORY_K4ME3" : "MEMORY",
	"TS102214_MEMORY_K4ME1" : "MEMORY",
	"TS102214_MEMORY_K27ME3" : "MEMORY",
	"TS081414_NAIVE_X" : "NAIVE",
	"TS081414_NAIVE_H3AC" : "NAIVE",
	"TS081414_NAIVE_K27AC" : "NAIVE",
	"TS081414_NAIVE_K4ME3" : "NAIVE",
	"TS081414_NAIVE_K4ME1" : "NAIVE",
	"TS081414_NAIVE_K27ME3" : "NAIVE",
	"TS102214_NAIVE_X" : "NAIVE",
	"TS102214_NAIVE_H3AC" : "NAIVE",
	"TS102214_NAIVE_K27AC" : "NAIVE",
	"TS102214_NAIVE_K4ME3" : "NAIVE",
	"TS102214_NAIVE_K4ME1" : "NAIVE",
	"TS102214_NAIVE_K27ME3" : "NAIVE",
	"TS081414_NAIVECD5P_X" : "NAIVECD5P",
	"TS081414_NAIVECD5P_H3AC" : "NAIVECD5P",
	"TS081414_NAIVECD5P_K27AC" : "NAIVECD5P",
	"TS081414_NAIVECD5P_K4ME3" : "NAIVECD5P",
	"TS081414_NAIVECD5P_K4ME1" : "NAIVECD5P",
	"TS081414_NAIVECD5P_K27ME3" : "NAIVECD5P",
	"VAKM_X" : "ACTB",
	"VAKM_H3AC" : "ACTB",
	"VAKM_K27AC" : "ACTB",
	"VAKM_K4ME3" : "ACTB",
	"VAKM_K4ME1" : "ACTB",
	"VAKM_K27ME3" : "ACTB",
	"VANN_X" : "ACTB",
	"VANN_H3AC" : "ACTB",
	"VANN_K27AC" : "ACTB",
	"VANN_K4ME3" : "ACTB",
	"VANN_K4ME1" : "ACTB",
	"VANN_K27ME3" : "ACTB",
	"VAQQQ_X" : "ACTB",
	"VAQQQ_H3AC" : "ACTB",
	"VAQQQ_K27AC" : "ACTB",
	"VAQQQ_K4ME3" : "ACTB",
	"VAQQQ_K4ME1" : "ACTB",
	"VAQQQ_K27ME3" : "ACTB",
	"VGA3_X" : "ACTB",
	"VGA3_H3AC" : "ACTB",
	"VGA3_K27AC" : "ACTB",
	"VGA3_K4ME3" : "ACTB",
	"VGA3_K4ME1" : "ACTB",
	"VGA3_K27ME3" : "ACTB",
	"VGA5_X" : "ACTB",
	"VGA5_H3AC" : "ACTB",
	"VGA5_K27AC" : "ACTB",
	"VGA5_K4ME3" : "ACTB",
	"VGA5_K4ME1" : "ACTB",
	"VGA5_K27ME3" : "ACTB",
	"VGA8_X" : "ACTB",
	"VGA8_H3AC" : "ACTB",
	"VGA8_K27AC" : "ACTB",
	"VGA8_K4ME3" : "ACTB",
	"VGA8_K4ME1" : "ACTB",
	"VGA8_K27ME3" : "ACTB",
	"VGA9_X" : "ACTB",
	"VGA9_H3AC" : "ACTB",
	"VGA9_K27AC" : "ACTB",
	"VGA9_K4ME3" : "ACTB",
	"VGA9_K4ME1" : "ACTB",
	"VGA9_K27ME3" : "ACTB",
	"VGA12_X" : "ACTB",
	"VGA12_H3AC" : "ACTB",
	"VGA12_K27AC" : "ACTB",
	"VGA12_K4ME3" : "ACTB",
	"VGA12_K4ME1" : "ACTB",
	"VGA12_K27ME3" : "ACTB",
	"VGR1_X" : "RESTB",
	"VGR1_H3AC" : "RESTB",
	"VGR1_K27AC" : "RESTB",
	"VGR1_K4ME3" : "RESTB",
	"VGR1_K4ME1" : "RESTB",
	"VGR1_K27ME3" : "RESTB",
	"VGR2_X" : "RESTB",
	"VGR2_H3AC" : "RESTB",
	"VGR2_K27AC" : "RESTB",
	"VGR2_K4ME3" : "RESTB",
	"VGR2_K4ME1" : "RESTB",
	"VGR2_K27ME3" : "RESTB",
	"VGR4_X" : "RESTB",
	"VGR4_H3AC" : "RESTB",
	"VGR4_K27AC" : "RESTB",
	"VGR4_K4ME3" : "RESTB",
	"VGR4_K4ME1" : "RESTB",
	"VGR4_K27ME3" : "RESTB",
	"VGR7_X" : "RESTB",
	"VGR7_H3AC" : "RESTB",
	"VGR7_K27AC" : "RESTB",
	"VGR7_K4ME3" : "RESTB",
	"VGR7_K4ME1" : "RESTB",
	"VGR7_K27ME3" : "RESTB",
	"VGR10_X" : "RESTB",
	"VGR10_H3AC" : "RESTB",
	"VGR10_K27AC" : "RESTB",
	"VGR10_K4ME3" : "RESTB",
	"VGR10_K4ME1" : "RESTB",
	"VGR10_K27ME3" : "RESTB",
	"VGR11_X" : "RESTB",
	"VGR11_H3AC" : "RESTB",
	"VGR11_K27AC" : "RESTB",
	"VGR11_K4ME3" : "RESTB",
	"VGR11_K4ME1" : "RESTB",
	"VGR11_K27ME3" : "RESTB",
	"VRKM_X" : "RESTB",
	"VRKM_H3AC" : "RESTB",
	"VRKM_K27AC" : "RESTB",
	"VRKM_K4ME3" : "RESTB",
	"VRKM_K4ME1" : "RESTB",
	"VRKM_K27ME3" : "RESTB"
}



#####--FUNCTIONS--#####

def Calc_Avg(inputfile, outputfile):
	"""
	Calculates the average QN'd, RPM'd values for each cell type as defined in the sample_data dictionary and prints the results to an output file.

	Args:
		inputfile: Name of the input file.
		outputfile: Name of the output file.
	"""

	#Lists to hold sample values of each cell type. Can add additional ones and modify code as needed.
	cb_vals = []
	cc_vals = []
	memory_vals = []
	naive_vals = []
	naivecd5p_vals = []
	actb_vals = []
	restb_vals = []

	chrom_correct = {"1":"1", "12":"2", "16":"3", "17":"4", "18":"5","19":"6", "20":"7", "21":"8", "22":"9", "2":"10",
	"3":"11", "4":"12", "5":"13", "6":"14", "7":"15", "8":"16", "9":"17", "10":"18", "11":"19", "13":"20", "14":"21", "15":"22", "23":"X"}

	#List to hold all the above lists to allow for easy interation through newly added cell types, etc.
	master_list = [cb_vals, cc_vals, memory_vals, naive_vals, naivecd5p_vals, actb_vals, restb_vals]

	#Dictionary to hold cell type for each column
	col_types = dict()

	#Open input file
	with open(inputfile) as f:

		#Open output file, "w" to make it writable
		output_file = open(outputfile, "w")

		#Create a dict for each cell type in the header to associate values in each column with a specific cell type
		header = f.readline()
		header = header.strip().split()

		#boolean to determine which columns have data just from the first line
		check_for_data = True

		#Remove the quotation marks around each column header that R adds for whatever reason. Can use quote=FALSE when writing the table out from R to get rid of this actually.
		for item in header:
			item_index = header.index(item)
			item = item.replace('"', "")
			header[item_index] = item
			
			#Check if item is in the cell type dictionary
			if item_index > 3:
				if item in sample_data:
					#If it is, add the index of the item and the sample type to the col_types dictionary
					col_types.update({item_index:sample_data[item]})
				else:
					print("Column header not found in cell type dictionary, check script!")

		#Print header to output file
		print(*header[0:4], sep="\t", end="\t",file=output_file)
		print("CB_avg","CC_avg","MEMORY_avg","NAIVE_avg","NAIVECD5P_avg","ACTB_avg","RESTB_avg", sep="\t",file=output_file)

		#Iterate through each line in input file
		for line in f:

			#The float function apparently can't handle the small e for scientific notation
			line = line.replace("e","E")

			#Strip newline character and split line by white space
			line = line.strip()
			line = line.split()

			#Add "chr" back to the chromosome element and correct for R changing the chromosome to a random number. 
			line[0] = "chr" + chrom_correct[line[0]]

			#Convert start, stop, peak_ID of line to ints because R output is sometimes silly and will write chromosomal positions in scientific notation
			for i in line[1:4]:
				i_index = line.index(i)
				line[i_index] = int(float(i))

			#Print chr, start, stop, & peak_ID values.
			print(*line[0:4], sep="\t", end="",file=output_file)

			#Iterate through line and add the data to the appropriate list based on cell type
			for data in line:

				#Get index of the value
				data_index = line.index(data)

				#For elements that contain data, add to the appropriate list based on the col_types dictionary
				if data_index in col_types:

					if col_types[data_index] == "CB":
						cb_vals.append(float(data))
					elif col_types[data_index] == "CC":
						cc_vals.append(float(data))
					elif col_types[data_index] == "MEMORY":
						memory_vals.append(float(data))
					elif col_types[data_index] == "NAIVE":
						naive_vals.append(float(data))
					elif col_types[data_index] == "NAIVECD5P":
						naivecd5p_vals.append(float(data))
					elif col_types[data_index] == "ACTB":
						actb_vals.append(float(data))
					elif col_types[data_index] == "RESTB":
						restb_vals.append(float(data))

			#Iterate through each list of data values for each cell type and average them
			for val_list in master_list:

				total = 0
				avg = float('nan')

				#Get running total of values for that cell type in the line
				if len(val_list) > 0:
					for val in val_list:
						total += float(val)

					#Calculate average based on total values and # of values in the list
					avg = (total/len(val_list))

					#Clear the list
					del val_list[:]

				#Checks which columns have valid averages and stores them in a list for fold change usage
				if check_for_data:
					if math.isnan(avg):
				 		filled_cols.append(False)
					else:
				 		filled_cols.append(True)
				
				#Print the average value
				print("\t" + "{0:.4f}".format(float(avg)), end="", file=output_file)

			#Change to false so it doesn't bother trying to figure out which columns have data in subsequent iterations
			check_for_data = False

			#New line
			print("", file=output_file)

		#Close output file
		output_file.close()


def Calc_FC(values):
	"""
	Calculates fold change between two values and returns the fold change value.

	Args:
		values: Tuple of the values to get fold change for. (value/denominator)

	Returns:
		FC_val = The fold change value.
	"""

	#Check if the values are both valid. Return nan if not.
	if math.isnan(values[0]) and math.isnan(values[1]):
		FC_val = float('nan')
	else:
		FC_val = ((values[0])/(values[1]))

	#Return result
	return FC_val


def Get_FC(inputfile, outputfile):
	"""
	Gets the fold change between each cell type as defined in the sample_data dictionary for each peak and prints the results to an output file.
	Used to iterate back through the file after all averages have been calculated. Potentially the laziest function ever written. 

	Args:
		inputfile: Name of the input file.
		outputfile: Name of the output file.
	"""

	#Lists to hold avg values of each cell type per line. Can add additional ones and modify code as needed. NaN by default.
	cb_val = float('nan')
	cc_val = float('nan')
	memory_val = float('nan')
	naive_val = float('nan')
	naivecd5p_val = float('nan')
	actb_val = float('nan')
	restb_val = float('nan')



	#List to hold all the fold change data for each line for easy printing. 
	FC_data = []

	#Open input file
	with open(inputfile) as f:

		#Open output file, "w" to make it writable
		output_file = open(outputfile, "w")

		#Create a dict for each cell type in the header to associate values in each column with a specific cell type
		header = f.readline()
		header = header.strip().split()

		#boolean to determine which columns have data just from the first line
		check_for_data = True

		#Remove the quotation marks around each column header that R adds for whatever reason. Can use quote=FALSE when writing the table out from R to get rid of this actually.
		for item in header:
			item_index = header.index(item)
			item = item.replace('"', "")
			header[item_index] = item

		#Print header to output file
		print(*header[0:], sep="\t", end="\t",file=output_file)
		print("CB_CC","MEMORY_CC","NAIVE_CC","NAIVECD5P_CC","ACTB_CC","RESTB_CC","MEMORY_CB","NAIVE_CB","NAIVECD5P_CB","ACTB_CB","RESTB_CB", 
			"NAIVE_MEMORY", "NAIVECD5P_MEMORY", "ACTB_MEMORY", "RESTB_MEMORY", "NAIVECD5P_NAIVE","ACTB_NAIVE","RESTB_NAIVE","ACTB_NAIVECD5P",
			"RESTB_NAIVECD5P", "ACTB_RESTB",sep="\t",file=output_file)

		#Iterate through each line in input file
		for line in f:

			#Strip newline character and split line by white space
			line = line.strip()
			line = line.split()

			#Convert start, stop, peak_ID of line to ints because R output is sometimes silly and will write chromosomal positions in scientific notation
			for i in line[1:4]:
				i_index = line.index(i)
				line[i_index] = int(float(i))

			#Print all values in the line
			print(*line[0:], sep="\t", end="",file=output_file)

			#Iterate through line and add the data to the appropriate variable based on index
			for data in line:

				#Get index of the value
				data_index = line.index(data)

				#Assign each element of the line to its appropriate variable
				if data_index == 4:
					cb_val = float(data)
				elif data_index == 5:
					cc_val = float(data)
				elif data_index == 6:
					memory_val = float(data)
				elif data_index == 7:
					naive_val = float(data)
				elif data_index == 8:
					naivecd5p_val = float(data)
				elif data_index == 9:
					actb_val = float(data)
				elif data_index == 10:
					restb_val = float(data)

			#Lists to hold tuples of values for which fold change will be calculated
			cb_cc = (cb_val, cc_val)
			memory_cc = (memory_val, cc_val)
			naive_cc = (naive_val, cc_val)
			naivecd5p_cc = (naivecd5p_val, cc_val)
			actb_cc = (actb_val, cc_val)
			restb_cc = (restb_val, cc_val)
			memory_cb = (memory_val, cb_val)
			naive_cb = (naive_val, cb_val)
			naivecd5p_cb = (naivecd5p_val, cb_val)
			actb_cb = (actb_val, cb_val)
			restb_cb = (restb_val, cb_val)
			naive_memory = (naive_val, memory_val)
			naivecd5p_memory = (naivecd5p_val, memory_val)
			actb_memory = (actb_val, memory_val)
			restb_memory = (restb_val,memory_val)
			naivecd5p_naive = (naivecd5p_val, naive_val)
			actb_naive = (actb_val, naive_val)
			restb_naive = (restb_val, naive_val)
			actb_naivecd5p = (actb_val, naivecd5p_val)
			restb_naivecd5p = (restb_val, naivecd5p_val)
			actb_restb = (actb_val, restb_val)

			#List to hold all the above tuples for easy iteration
			tuple_list = [cb_cc, memory_cc, naive_cc, naivecd5p_cc, actb_cc, restb_cc, memory_cb, naive_cb, naivecd5p_cb, 
			actb_cb, restb_cb, naive_memory, naivecd5p_memory, actb_memory, restb_memory, naivecd5p_naive, actb_naive, 
			restb_naive,actb_naivecd5p, restb_naivecd5p, actb_restb]

			#Iterate through each tuple of avg values for each cell type and calculate the fold change for them
			for tup in tuple_list:
				FC_data.append(Calc_FC(tup))

			
			#Print all values in the line
			[print("\t" + "{0:.4f}".format(val), end="",file=output_file) for val in FC_data]
				
			#Clear the fold change data list for the next line
			del FC_data[:]

			#New line
			print("", file=output_file)

		#Close output file
		output_file.close()



def Tidy_Up(temp_output1):
	"""Saves the output in the true, user-defined output file. Then removes the temporary file.

	Args: 
		temp_output: Name of the temporary output file.
		final_output: Name of the final output file. 
	"""
	#Go ahead and call a command line sort after closing the output_file. This can be dangerous if used improperly. Read the docs on subprocess.call()
	#If really needed.
	#call("(head -n 1 " + temp_output2 + " && tail -n +2 " + temp_output2 + " | sort -k1.4,1.5 -k2 -V)" + " > " +  final_output, shell=True)
	#The above is no longer needed due to changes to Calc_Avg function and the fact that R just renamed the column values rather than 
	#mess up the order. 
	call("rm " + temp_output1, shell=True)


def Make_UCSC(column_index):
	"""
	Creates a gzipped bedgraph file from the specified data column of the output file. Bases the trackline and file name on the column header.

	Parameters:
		column_index = The index of the column to be used as data for the 4th column of the resulting output file
	"""

	#Open the input (original output) file
	with open(args.output_file) as f:

		#Get the header and split into a list
		header=f.readline().strip().split()

		#Get the appropriate column
		basename=header[column_index]

		#Tell user what's going on
		print("Checking " + basename +" column.")

		#Set the output name
		output_name="FC_" + mark + "_" + basename + ".bedgraph.gz"

		#Check first line to see if the column has data, and if not, don't bother making a file from it
		first_line = f.readline().strip().split()
		check_col = first_line[column_index]
		first_line_data = first_line[0] + "\t" + first_line[1] + "\t" + first_line[2] + "\t" + first_line[column_index] + "\n"

		if check_col == "nan":
			print("Column contains empty values, moving to next column.")
			return

		#Tell user what's going on
		print("Writing to " + output_name)

		#Open the bedgraph output file
		output = gzip.open(output_name, "wt")

		gzip_header = ("track type=bedGraph name=FC_" + mark + "_" + basename + " description=FC_" + mark + "_" + basename +
			" visibility=2 color=200,0,0 altColor=0,100,200 priority=10 autoScale=off alwaysZero=off gridDefault=on graphType=bar viewLimits=0:3\n")

		output.write(gzip_header)
		output.write(first_line_data)

		#Iterate through file and print to the output file
		for line in f:
			line = line.strip().split()

			#Get data
			keep = line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[column_index] + "\n"
			output.write(keep)

		output.close()


####-Variables-####

#-i <input.txt> -o <output.txt> -ucsc 
parser.add_argument("-i", "--input", dest = "input_file", required=True)
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-ucsc", "--ucsc", action="store_true")

args = parser.parse_args()

#Output so user can double check options
print( "Input file: {}\nOutput file: {} \n \n".format(
        args.input_file,
        args.output_file
        ))


#Creates a temp file for calculating average. Passed to fold change function.
output_file_tmp = args.output_file + ".tmp"
mark=args.input_file.split("_")[1]


#####--MAIN--#####

#Call our functions defined above. First calculate the averages for each cell type.
print("Getting average values for each cell type.")
Calc_Avg(args.input_file, output_file_tmp)
print("Averages calculated.")

#Calculate fold change between averages for cell types of interest.
print("Calculate fold change for comparisons of interest.")
Get_FC(output_file_tmp, args.output_file)
print("Fold change calculated.")

#Check if user wants to create UCSC tracks from each data column
if args.ucsc_bool:

	#Iterate through the output file the appropriate number of times
	for i in range(11,32):

		Make_UCSC(i)

#Remove temporary file.
print("Removing temporary files.")
Tidy_Up(output_file_tmp)
print("Complete.")
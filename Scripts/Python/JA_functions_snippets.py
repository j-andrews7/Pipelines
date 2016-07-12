"""
Up to date as of 01/12/2016
Documented code snippets for simple things I forget and useful functions that I'd rather not write again.
"""


#delete all elements of a list
mylist[:] = []

#Add multiple elements to a list at once
mylist.extend(multiple other lists, etc)

#Go ahead and call a command line sort after closing the output_file. This can be dangerous if used improperly. Read the docs on subprocess.call()
#If really needed. Can call other programs and shell commands in this way.
from subprocess import call
call("(head -n 1 " + output_name_tmp + " && tail -n +2 " + output_name_tmp + " | sort -k1.4,1.5 -k2 -V)" + " > " +  output_name, shell=True)

#More secure way
# Sort and remove temp file.
with open(out_file, "w") as fout:
	cmd = ["sort", "-V", "-k1,1", "-k2,2", tmp_out]
	call(cmd, stdout=fout)

cmd = ["rm", tmp_out]
call(cmd)


#Check if variable is NaN
import math
math.isnan(x)

#Reduce output to four decimal places
print("{0:.4f}".format(val))

#Return r length subsequences of elements from the input iterable. Yields combination tuples. Sort iterable before using if a specific order is needed. 
#Combinations are emitted in lexicographic sort order.
import itertools
itertools.combinations(iterable, r)

#Unpack tuple
my_i, my_card = select_choice()

#Exit python script
import sys
sys.exit()

#Randomly select item from list
import random
random.choice(insults)

#####--FUNCTIONS--#####
#These were copied from various scripts verbatim, but can be easily tweaked to meet other needs. Most have the required modules (if any) listed above them.
#Some may suck, others may be quite good. All work.

#Calculates average for data columns of the same cell type. 
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

			#Add "chr" back to the chromosome element
			line[0] = "chr" + line[0]

			#Convert chr23 back to chrX because R is satan incarnate and converts them all to integers
			if "chr23" in line[0]:
				line[0] = "chrX"

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

#Calculates fold change for a tuple of 2 values.
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

#Gets fold change for all data columns in a bed file. Done in a very lazy way, using itertools would be more "pythonic", but getting everything in the 
#proper column would require some thinking. 
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

#Calls shell commands to remove tmp files used by script and sort the final output file. Saves a bit of time and headache.
def Tidy_Up(temp_output1, temp_output2, final_output):
	"""
	Sorts the temporary output file and saves the output in the true, user-defined output file. Then removes the temporary file.

	Args: 
		temp_output: Name of the temporary output file.
		final_output: Name of the final output file. 
	"""
	#Go ahead and call a command line sort after closing the output_file. This can be dangerous if used improperly. Read the docs on subprocess.call()
	#If really needed.
	call("(head -n 1 " + temp_output2 + " && tail -n +2 " + temp_output2 + " | sort -k1.4,1.5 -k2 -V)" + " > " +  final_output, shell=True)
	call("rm " + temp_output1, shell=True)
	call("rm " + temp_output2, shell=True)

#Plots the average methylation level for the input file as a histogram
#Requires:
import matplotlib as mpl #Imported for plotting 
mpl.use('Agg') #Necessary for saving the plot as an image rather than displaying it
import matplotlib.pyplot as plt #Imported for plotting
def plotMethLevel(fileName):

	dataList = [] #Holds the methylation values for plotting
	lineSplit = [] #Initializes empty list to hold values in each line
	methData = 0 #Holds the methylation data for a line

	#Determine filepath to save plot to based on the input fileName
	basename = os.path.splitext(os.path.basename(fileName))[0] 
	plotname = basename + "_distribution.png"

	#Open input file
	with open(fileName) as cFile:

		#Loop through file
		for line in cFile:

			#Remove newline
			line = line.rstrip()

			#Split line and assign contents to a list and assign the appropriate elements of list to the call variables
			lineSplit = line.split()
			methData = float(lineSplit[5])

            #Add the average methylation of the CGI to the dataList for plotting
			dataList.append(methData)

    ########Plotting########

	plt.figure() #Make one figure
	plt.hist(dataList,10) #Plot the calcualted methylation values, grouping them every 0.1. 
	plt.title("Methylation Distribution Across CGIs") #Set title
	plt.xlabel("Average Methylation Scores") #Set x-axis label
	plt.ylabel("Occurrences") #Set y-axis label
	plt.savefig(plotname) #Save the plot to the appropriate file name

#Gets average and std for a list of values
import math
def average_and_stdev(input_list):
       
    my_list=input_list #Declare a variable and assigns the argument list to it
    my_sum=sum(my_list) #Declare a variable and sum each item in my_list and assign the value to it
    mean=my_sum/len(my_list) #Declare a variable and assign the mean for my_sum to it by dividing my_sum by the number of elements in my_list
    deviations=0 #Declare a variable and assign a value of 0 to it

    for j in my_list: #Iterate through each item in my_list
        deviations=deviations+(j-mean)**2 #Calculates the variance for each item in my_list and adds it to the deviations variable
    stdev=math.sqrt(deviations/(len(my_list)-1)) #Calculates the standard deviation for the line by dividing the deviations variable 
    											 #by the number of elements in the line - 1 and taking the squareroot of the result
    return mean,stdev; #Returns the mean and standard deviation for the line (gene)

#Returns reverse complement of a DNA sequence
def reverse_complement(sequence):
    #This function takes the reverse complement of a sequence
    rcSeq = ""

    #define complement dictionary - No need for actual dictionary 
    compDict = {"A": "T", "T":"A", "G":"C", "C":"G", "N":"N"}

    #Declare variable to hold reverse complement and assign to it the reversed sequence
    revSeq = sequence.upper()[::-1]

    #Create the complement of the reverse sequence by comparing each character of the sequence to the compDict
    for char in revSeq:
        rcSeq+=compDict[char]

    #return the reverse complement
    return rcSeq

#Used for command line options and such
import argparse
parser = argparse.ArgumentParser()

####-Variables-####

#-se <superenhancers.bed> -g <genes.gtf> -size <wing size from center of SE to use> -o <output.bed>
parser.add_argument("-se", "--superenhancer", dest = "SE_file", help="Super enhancer file name")
parser.add_argument("-g", "--gene", dest = "gene_file", help="Gene list file name")
parser.add_argument("-size", "--wingsize", dest ="wing_size", help="Size of wing to use", type=int, default=125000)
parser.add_argument("-o", "--output", dest = "output_file", help="Name of output file")

args = parser.parse_args()

print( "Super enhancer file {} Gene file {} Wing size {} Output {} ".format(
        args.SE_file,
        args.gene_file,
        args.wing_size,
        args.output_file,
        ))


#Grab SE IDs and positions from SE coords file
def SE_Grabber(SE_coords):
	"""
	Grabs SE and chromosomal positions from a text file.

		Args: 
			SE_coords = SE_coords file in txt format

		Returns:
			SE_dict = Dictionary in following format {SE_ID:(chrom, (start, end))}
	"""

	#Initialize SE dictionary
	SE_dict = {}

	#Open provided SE_coords file
	with open(SE_coords) as f:

		#Iterate through file line by line
		for line in f:

			#Split line by tab
			line = line.strip().split()

			SE_ID = line[0]
			chrom = line[1]
			start = line[2]
			stop = line[3]


			#Add the SE to the SE dictrionary
			SE_dict[SE_ID] = (chrom, (start, stop))

	#Return the SE dictionary
	return SE_dict

def Add_Positions(SE_ID_Vals_file, SE_Pos_Dict, Output):
	"""
	Adds chromosomal position to the SE_ID_Vals_file for each SE_ID and writes output to the output file.

		Args: 
			SE_ID_Vals_file = SE_ID_vals file in txt format
			SE_Pos_Dict = Dictionary containing the chromosomal positions of each SE in the SE_coords file.
			Output = Name of output file to write to.
	"""

	#Open SE_IDS_Vals_file
	with open(SE_ID_Vals_file) as f:


		#Get header
		header = f.readline().strip().split()

		#Open output file
		output_file = open(Output, "w")

		#Print header
		print("SEID", "CHR", "START", "STOP", *header[1:], sep="\t", file=output_file)

		#Iterate through line by line
		for line in f:

			line = line.strip().split()
			SEID = line[0]
			data = line[1:]

			if SEID in SE_Pos_Dict:
				positions = SE_Pos_Dict[SEID]

				chrom = positions[0]
				start_stop = positions[1]
				start = start_stop[0]
				stop = start_stop[1]

				#Print to output
				print(SEID,chrom,start,stop,*data[0:], sep="\t", file=output_file)

			else:
				print(SEID + " not found in super enhancer dictionary! Check files!")

		#Close output file
		output_file.close()


def Get_Midpoint(start_pos, stop_pos):
	"""
	From a start and stop position, finds the midpoint between the two.

		Args: 
			start_pos = Start position of the loci.
			stop_pos = Stop position of the loci.

		Returns:
			midpoint = Integer midpoint between the start and stop positions
	"""
	#Add two points together
	total = int(start_pos) + int(stop_pos)

	#Find midpoint using integer division to round
	midpoint = total//2

	#Return result
	return midpoint


def Add_Wings_Mid(position, wing_length):
	"""
	From a single position, add wings of either side of a specified length and return the two resulting positions as a tuple.

		Args: 
			position = Position to add the wings to.
			wing_length = The length of the wings to add to each side of the position.

		Returns:
			wing_positions = A tuple containing the last position of each wing.
	"""
	wing_length = int(wing_length)

	#Add and subtract wing length from each position
	wing_start = position - wing_length
	wing_stop = position + wing_length

	#Create tuple
	wing_positions = (wing_start, wing_stop)

	#Return result
	return wing_positions

def Add_Wings_Edge(start, stop, wing_length):
	"""
	From start and stop positions, add wings to both of a specified length and return the two resulting positions as a tuple.

		Parameters: 
			start = Start position to which the wing will be added.
			stop = Stop position to which the wing will be added.
			wing_length = The length of the wings to add to each side of the position.

		Returns:
			wing_positions = A tuple containing the last position of each wing.
	"""
	wing_length = int(wing_length)

	#Add and subtract wing length from each position
	wing_start = start - wing_length
	wing_stop = stop + wing_length

	#Create tuple
	wing_positions = (wing_start, wing_stop)

	#Return result
	return wing_positions


def Overlap(start1, end1, start2, end2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    return end1 >= start2 and end2 >= start1


def Check_Num(s):
	""" Check if a number is in the string, returns True if so, otherwise False."""
    return any(i.isdigit() for i in s)


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

		#Set the output name
		output_name=basename+".bedgraph.gz"

		#Tell user what's going on
		print("Writing to " + output_name)

		#Open the bedgraph output file
		output = gzip.open(output_name, "wt")

		gzip_header = ("track type=bedGraph name=QN_" + basename + " description=QN_" + basename + 
			" visibility=2 color=200,0,0 altColor=0,100,200 priority=10 autoScale=off alwaysZero=off gridDefault=on graphType=bar viewLimits=0:50\n")

		output.write(gzip_header)

		#Iterate through file and print to the output file
		for line in f:
			line = line.strip().split()

			#Get data
			keep = line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[column_index] + "\n"
			output.write(keep)

		output.close()

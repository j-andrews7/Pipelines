#!/usr/bin/env python3
"""
Grab SEIDs from a SE bed positional file and add them to the appropriate column of the tags file for bins that overlap the SE. Ignores lines that don't overlap 
Also condenses the bins and sums the RPM data for each sample column.

Usage: python3 get_SEID_condense.py <SE position input file> <tags file> <output file>
"""

#Necessary packages
import sys

#Get input file
SE_input_file = sys.argv[1]
tags_input_file = sys.argv[2]
output_file = sys.argv[3]
out_condensed = "condensed_" + output_file 


####-Functions-####

def Overlap(start1, end1, start2, end2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    return end1 >= start2 and end2 >= start1

#Grab SE IDs and positions from SE coords file
def SE_Grabber(SE_coords):
	"""
	Grabs SE and chromosomal positions from a text file.

		Parameters: 
			SE_coords = SE_coords file in txt format

		Returns:
			SE_dict = Dictionary in following format {SE_ID:(chrom, (start, end))}
	"""

	print("Creating SE_ID dictionary.")

	#Initialize SE dictionary
	SE_dict = {}

	#Open provided SE_coords file
	with open(SE_coords) as f:

		#Iterate through file line by line
		for line in f:

			#Split line by tab
			line = line.strip().split()

			SE_ID = line[3]
			chrom = line[0]
			start = line[1]
			stop = line[2]


			#Add the SE to the SE dictionary
			SE_dict[SE_ID] = (chrom, (start, stop))

	#Return the SE dictionary
	return SE_dict

def ADD_SEID(SE_dictionary, tags_input, out_file):
	"""
	For each line of the tags input file, adds the SE_ID to the fourth column. 

		Parameters:
			SE_dictionary = Dictionary containing (chrom, (start, stop)) for each SE.
			tags_input = The input file with the tags data.
			out_file = File to write output to.
	"""

	#Open output file
	output = open(out_file, "w")

	print("Finding overlap with SE positions from tag file positions.")

	#Open input file and iterate through line by line
	with open(tags_input) as f:

		#Get header and print to output file
		header = f.readline().strip().split("\t")
		header.insert(3, "SE_ID")
		print(*header[0:], sep="\t", file=output)

		#Iterate through file line by line
		for line in f:

			line = line.strip().split() 

			SE_found = False

			#Assign line elements to variables
			locus_chrom = line[0]
			locus_start = int(line[1])
			locus_stop = int(line[2])

			#Iterate through each gene to see if it overlaps the winged locus positions
			for SE in SE_dictionary:
				SE_position = SE_dictionary[SE]
				SE_chrom = SE_position[0]

				#Check if it's even on the same chromosome as the locus, if not, skip to the next SE
				if SE_chrom == locus_chrom:
					SE_start_stop = SE_position[1]
					SE_start = int(SE_start_stop[0])
					SE_stop = int(SE_start_stop[1])

					#Check for overlap of the two ranges and add the gene to the gene_list if this returns True.
					if Overlap(locus_start, locus_stop, SE_start, SE_stop):
						SE_found = True
						line.insert(3, SE) 
						print(*line[0:], sep="\t", file=output)
						break

	output.close()


def Condense_By_SEID(full_output, cond_output):
	"""
	Condenses all lines by SE_ID, summing the RPM data for each sample in the process. 

		Parameters:
			full_output = Uncondensed input file with SE_IDs in 4th column followed by data.
			cond_output = File to output condensed data to.
	"""

	#Open output file, "w" to make it writable
	output_file = open(cond_output, "w")

	print("Condensing by SE_ID.")

	#Open input file
	with open(full_output) as f:

		#Get header and write it to the output file
		header = f.readline().strip().split()
		del header[4]
		print(*header[0:], sep="\t",file=output_file)

		#Initiate variable to determine if SE is "new" or needs to be appended to
		new_SE = True

		#Initialize a list to hold RPM values for each SE
		RPM_list = []

		#Initialize a list to hold SEs
		SE_list = []

		#Iterate through each line in input file
		for line in f:

			#Strip newline character and split line by white space
			line = line.strip().split()

			#Get SE_ID for the line
			SE_ID = line[3]

			#Check if the SE_ID is already in the SE list and adjust the new SE variable if necessary
			if SE_ID in SE_list:
				new_SE = False
			else:
				new_SE = True

			#If the SE is new, add it to the SE list and get the start position
			if new_SE:

				#Check if this is the first SE
				if len(SE_list) is 0:
					SE_list.append(SE_ID)
					chrom = line[0]
					start = int(line[1])
					curr_SE = SE_ID
					end = int(line[2])

					#Slice line list to get the data elements
					RPM_list = line[4:]

				#If not the first SE, adjust accordingly and print the previous SE before setting up the new one
				else:			
					
					#Print results for the old SE to the output file
					print(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + str(curr_SE), end="", file=output_file)

					#Print each item in the RPM_list
					for item in RPM_list:
						print("\t" + "{0:.4f}".format(float(item)), end="", file=output_file)

					#Print newline
					print("", file=output_file)

					#Set all variables to new SE
					SE_list.append(SE_ID)
					chrom = line[0]
					start = int(line[1])
					curr_SE = SE_ID
					end = int(line[2]) 

					#Slice line list to get the data elements
					RPM_list = line[5:]

			#If the SE is not new, add the values of all RPMs to the RPM_list
			else:

				#Throw RPM values for current line into a list
				curr_RPM_vals = line[5:]

				#Get end position, this won't be printed or used except for the last SE in the file
				end = int(line[2])

				#Iterate through all values in the RPM_list and add the appropriate RPM value from the current line to it
				for i in range(0,len(RPM_list)):
					RPM_list[i] = float(RPM_list[i]) + float(curr_RPM_vals[i])

		#Print results for the old SE to the output file
		print(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + str(curr_SE), end="", file=output_file)

		#Print each item in the RPM_list
		for item in RPM_list:
			print("\t" + "{0:.4f}".format(float(item)), end="", file=output_file)

		#Print newline
		print("", file=output_file)

		#Close output file
		output_file.close()


####-Main-####

#Open SE_input file and make dictionary for each SEID
SE_di = SE_Grabber(SE_input_file)

#Add SE_ID
ADD_SEID(SE_di, tags_input_file, output_file)

#Condense by SE_ID
Condense_By_SEID(output_file, out_condensed)

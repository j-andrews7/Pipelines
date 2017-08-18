#!/usr/bin/env python3
"""For specific super enhancer IDs, get the chromosomal positions from a different file and stick in columns 2-4 of original file.

Usage: python3 get_SE_coords.py <SE_ID_vals.txt> <SE_coords.txt> <output.txt>

Args:
    SE_ID_vals.txt = Name of input file with SE_IDs in first column followed by tags from samples
    SE_coords.txt = Name of input file with SE_IDs in first column followed by the chromosomal positions of the SE
    output.txt = Name of output file to be created
"""

import sys

####-Functions-####

#Grab SE IDs and positions from SE coords file
def SE_Grabber(SE_coords):
	"""Grabs SE and chromosomal positions from a text file.

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
	"""Adds chromosomal position to the SE_ID_Vals_file for each SE_ID and writes output to the output file.

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


####-Variables-####

#Store file names
SE_ID_file = sys.argv[1]
SE_coords_file = sys.argv[2] 
output_name = sys.argv[3]


####-MAIN-####
#Run functions defined above
SE_dictionary = SE_Grabber(SE_coords_file)
Add_Positions(SE_ID_file, SE_dictionary, output_name)

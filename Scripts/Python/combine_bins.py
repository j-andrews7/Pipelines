#!/usr/bin/env python3
"""For directory of .bin files, combines each bin file into a single, sorted text file

Usage: python3 combine_bins.py <directory> <output.txt>

Args:
    directory = Path to directory of bin files to combine 
    output.txt = Name of output file to be created
"""
import sys
import operator
import os, os.path

#Store each filename in a list
directory = sys.argv[1]
file_list = os.listdir(directory)

#Store output file name as variable
outputname = sys.argv[2]

#Initialize dictionary to hold each unique bin and the values for it
bins_dict = {}

#Initialize list for eventual output header
header = ["chr", "bin"]

#If bins need to be gotten from first .bin file
get_bins = True


#Iterate through each file
for item in file_list:


	#Get complete path of file
	fullpath = directory + item

	#Get sample name and store in a list
	sample_name = item.split("_")[0] + "_" + item.split("_")[1]
	header.append(sample_name)

	#If first file has not yet been opened
	if get_bins:

		#Feedback is good
		print("Creating bin IDs from " + sample_name + "...")

		#Open first file
		with open(fullpath) as f:

			#Iterate through each line in input file
			for line in f:

				#Skips lines with non-standard chromosomes
				if "_g" in line or "random" in line or "M" in line:
					next(f)

				#If line is not for an X chromosome, store the chromosome as an integer. Necessary for later sorting.
				else:
					#Remove newline character at end of line and split by whitespace
					line = line.strip().split()

					#Create unique ID for the bin in the line
					bin_ID = (line[0],line[1])

					#Get tags for bin
					value = line[2]

					#Add new key/value pair to the dictionary
					bins_dict[bin_ID] = value

			print("Number of unique bins: " + str(len(bins_dict.keys())))

			#Change variable so we know the temp output file now has each unique bin
			get_bins = False

	#Open rest of files, get last data element in line and append to temp output file
	else:

		#More feedback
		print("Get values for " + sample_name + " for each unique bin")

		#Open file
		with open(fullpath) as fi:

			#Iterate through each line in input file
			for line in fi:

				#Skips lines with non-standard chromosomes and doesn't bother with Y chromosome
				if "_g" in line or "random" in line or "M" in line:
					next(fi)

				#If line is not for an X chromosome, store the chromosome as an integer. Necessary for later sorting.
				else:
					#Remove newline character at end of line and split by whitespace
					line = line.strip().split()

					#Get key and value from the line to be used in comparison to the dictionary
					ID = (line[0],line[1])
					val = line[2]

					#If the bin already exists in the bins_dict as a key (which it should), add the value to the value already associated with the key, separated by tab
					if ID in bins_dict:
						bins_dict[ID] = bins_dict[ID] + "\t" + val

					#Yells at you if input files have different numbers of bins and breaks out of loop
					else:
						print("Bins not uniform across files - check input files!")
						break


###############-OUTPUT-####################

#Additional feedback
print("Printing output to " + outputname)

#Open output file
outfile = open(outputname, "w")

#Print header
header = "\t".join(header)

print(header, file=outfile) 

#Print results
for key in sorted(bins_dict, key=lambda x: (x[0], x[1])):

	#Split key back up into chrom/bin and convert to strings
	chrom = "chr" + str(key[0])
	bin_val = str(key[1])

	#Get list of values for the key
	values = bins_dict[key]

	#Print results
	print(chrom + "\t" + bin_val + "\t" + values, file=outfile)


#Close output file
outfile.close()

print("Done!")


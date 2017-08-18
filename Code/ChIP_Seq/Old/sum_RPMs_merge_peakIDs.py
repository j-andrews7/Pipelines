#!/usr/bin/env python3
"""Merge all lines with the same peak_ID, summing the RPMs for each sample in the peak and getting chromosome position accordingly

Usage: python3 sum_RPMs_merge_peakIDs.py <input.bed> <output.bed>

Args:
    input.bed = Name of input file to process 
    output.bed = Name of output file to be created

NOTE:PEAK_ID must be the 5th column, followed by all the RPM values. Check script for more info.
"""

import sys

#Store file names
input_file = sys.argv[1]
output_name = sys.argv[2]

print("""\

                                       ._ o o
                                       \_`-)|_
                                    ,""       \ 
                                  ,"  ## |   @ @.      ________________
                                ," ##   ,-\__    `.   /                |
                              ,"       /     `--._;) <  WOW, SCIENCE   |  
                            ,"     ## /               \________________|
                          ,"   ##    /


                    """)

print("GIRAFFE PERFORMING INCREDIBLY COMPLEX FUNCTIONS, PLEASE WAIT SIR/MADAM")





#Open output file, "w" to make it writable
output_file = open(output_name, "w")

#Open input file
with open(input_file) as f:

	#Get header and write it to the output file
	header = f.readline().strip().split()

	#Remove the bin column
	del header[3]
	
	print(*header[0:], sep="\t", file=output_file)

	#Initiate variables to hold chr, start/end, and ID of current peak
	start = 0
	end = 0
	chrom = ""
	curr_peak = ""

	#Initiate variable to determine if peak is "new" or needs to be appended to
	new_peak = True

	#Initialize a list to hold RPM values for each peak
	RPM_list = []

	#Initialize a list to hold peaks
	peak_list = []

	#Iterate through each line in input file
	for line in f:

		#Strip newline character and split line by white space
		line = line.strip()
		line = line.split()

		#Get peak_ID for the line
		peak_ID = line[4]

		#Check if the peak_ID is already in the peak list and adjust the new peak variable if necessary
		if peak_ID in peak_list:
			new_peak = False
		else:
			new_peak = True

		#If the peak is new, add it to the peak list and get the start position
		if new_peak:

			#Check if this is the first peak
			if len(peak_list) is 0:
				peak_list.append(peak_ID)
				chrom = line[0]
				start = int(line[1])
				curr_peak = peak_ID
				end = int(line[2])

				#Slice line list to get the data elements
				RPM_list = line[5:]

			#If not the first peak, adjust accordingly and print the previous peak before setting up the new one
			else:			

				#Print results for the old peak to the output file
				print(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + str(curr_peak), end="", file=output_file)

				#Print each item in the RPM_list
				for item in RPM_list:
					print("\t" + "{0:.4f}".format(float(item)), end="", file=output_file)

				#Print newline
				print("", file=output_file)

				#Set all variables to new peak
				peak_list.append(peak_ID)
				chrom = line[0]
				start = int(line[1])
				curr_peak = peak_ID
				end = int(line[2]) 

				#Slice line list to get the data elements
				RPM_list = line[5:]

		#If the peak is not new, add the values of all RPMs to the RPM_list
		else:

			#Throw RPM values for current line into a list
			curr_RPM_vals = line[5:]

			#Get end position, this won't be printed or used except for the last peak in the file
			end = int(line[2])

			#Iterate through all values in the RPM_list and add the appropriate RPM value from the current line to it
			for i in range(0,len(RPM_list)):
				RPM_list[i] = float(RPM_list[i]) + float(curr_RPM_vals[i])

	#Print results for the old peak to the output file
	print(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + str(curr_peak), end="", file=output_file)

	#Print each item in the RPM_list
	for item in RPM_list:
		print("\t" + "{0:.4f}".format(float(item)), end="", file=output_file)

	#Print newline
	print("", file=output_file)

	#Close output file
	output_file.close()

print("COMPLETE")
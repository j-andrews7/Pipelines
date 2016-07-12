#!/usr/bin/env python3
"""Merge all lines with the same transcript ID, smashing exons together.

Usage: python3 merge_exons_gtf.py <input.gtf> <output.gtf>

Args:
    input.gtf = Name of input file to process 
    output.gtf = Name of output file to be created
"""

import sys

#Store file names
input_file = sys.argv[1]
output_name = sys.argv[2]

print("""\

                                       ._ o o
                                       \_`-)|_
                                    ,""       \ 
                                  ,"  ## |   @ @.      _____________________
                                ," ##   ,-\__    `.   /                     |
                              ,"       /     `--._;) <  WOW, EXON MERGING   |  
                            ,"     ## /               \_____________________|
                          ,"   ##    /


                    """)

print("GIRAFFE PERFORMING INCREDIBLY COMPLEX FUNCTIONS, PLEASE WAIT SIR/MADAM")


#Open output file, "w" to make it writable
output_file = open(output_name, "w")

#Open input file
with open(input_file) as f:

	#Initiate variables to hold chr, start/end, and ID of current transcript
	start = 0
	end = 0
	chrom = ""
	curr_transcript = ""
	curr_line_end = ""

	#Initiate variable to determine if transcript is "new" or needs to be appended to
	new_transcript = True

	#Initialize a list to hold transcript
	transcript_list = []

	#Iterate through each line in input file
	for line in f:

		#Strip newline character and split line by white space
		line = line.strip()
		line = line.split("\t")

		#Get transcript ID for the line
		line_end = line[8]
		transcript_ID = line[8].split("; ")[1]

		#Check if the transcript_ID is already in the transcript_list and adjust the new_transcript variable if necessary
		if transcript_ID in transcript_list:
			new_transcript = False
		else:
			new_transcript = True

		#If the transcript is new, add it to the transcript_list and get the start position
		if new_transcript:

			#Check if this is the first transcript
			if len(transcript_list) is 0:
				transcript_list.append(transcript_ID)
				chrom = line[0]
				line_beg=[line[1], line[2], line[5], line[6], line[7]]
				start = int(line[3])
				curr_transcript = transcript_ID
				curr_line_end = line_end
				end = int(line[4])

			#If not the first peak, adjust accordingly and print the previous transcript before setting up the new one
			else:			


				#Print results for the old transcript to the output file
				print(chrom + "\t" + line_beg[0] + "\t" + line_beg[1] + "\t" + str(start) + "\t" + str(end) + "\t" + 
					line_beg[2] + "\t" + line_beg[3] + "\t" + line_beg[4] + "\t" + curr_line_end, file=output_file)

				#Set all variables to new transcript
				transcript_list.append(transcript_ID)
				chrom = line[0]
				line_beg=[line[1], line[2], line[5], line[6], line[7]]
				start = int(line[3])
				curr_transcript = transcript_ID
				curr_line_end = line_end
				end = int(line[4]) 


		#If the transcript is not new, change the end variable as appropriate
		else:

			#Get end position, this won't be printed or used except for the last peak in the file
			end = int(line[4])

	#Print results for the old transcript to the output file
	print(chrom + "\t" + line_beg[0] + "\t" + line_beg[1] + "\t" + str(start) + "\t" + str(end) + "\t" + 
	line_beg[2] + "\t" + line_beg[3] + "\t" + line_beg[4] + "\t" + curr_line_end, file=output_file)

	#Close output file
	output_file.close()

print("COMPLETE")
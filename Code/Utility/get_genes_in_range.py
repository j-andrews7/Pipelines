#!/usr/bin/env python3
"""
For the specified gene and locus lists, determine which genes are located within a specific range from the middle or edges of each locus. 
The gene symbols are then added to the last column of the locus table and printed to the output file.

Usage: python3 get_genes_in_range.py -i <input.bed> -g <genes.gtf> -o <output.bed> [OPTIONS]

Args:
    (required) -i <input.bed> = Name of locus list file to process.
    (required) -g <genes.gtf or genes.bed> = Name of gene list file to intersect with the loci list.
    (required) -o <output.bed> = Name of output file to be created.
    (optional) -C = Column with chromosome - one-based (default=1).
	(optional) -S = Column with start position - one-based (default=2).
	(optional) -E = Column with end position - one-based (default=3).
    (optional) -size <wing size> = An integer to define the distance from the center or edges of the loci to look for gene overlap. 125 kb by default.
    (optional) -H = Use if the input file has a header. If left out, it will be assumed there is no header.
    (optional) -m = Will add wings from the midpoint of the loci rather than the edges. 
    (optional) -p = Will not print loci that overlap the TSS (+/- psize) of any of the transcripts in the gtf file. Good for typical RE searches.
    (optional) -twsize = Size of wing to add to each edge of the TSS of each transcript to set the "promoter" position. Default is 2000 bp.
    (optional) -noovlp = Will not include loci that directly overlap with a gene in the annotation file. Good for SE searches.
"""

import sys
#Used for command line options and such
import argparse
parser = argparse.ArgumentParser(usage=__doc__)


####-Functions-####

def Get_Midpoint(start_pos, stop_pos):
	"""
	From a start and stop position, finds the midpoint between the two.

		Parameters: 
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

		Parameters: 
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
    """
    Does the range (start1, end1) overlap with (start2, end2)?
    """
    return int(end1) >= int(start2) and int(end2) >= int(start1)


#Grab gene names, chromosomal positions, and potentially TSSs from gene list file
def Gene_Grabber(gene_f, get_TSS=False, TSS_wing=2000):
	"""
	Grabs gene names, chromosomal positions, and potentially the TSSs from a gtf or bed file.

		Parameters: 
			gene_f = Gene file in gtf or bed format

		Returns:
			gene_dict = Dictionary in following format {gene_name:(chrom, (start, end))}
			(optional) TSS_dict = Dictionary containing 
	"""

	#Initialize gene dictionary
	gene_dict = {}
	TSS_dict = {}

	# Get file extension.
	ext = gene_f.split(".")[-1]

	if ext.upper() == "BED":

		#Check if the TSS positions are wanted
		if get_TSS:

			# Can't use TSS with bed file. Well, can, but the file I wanted to use didn't have strand info, so I'm not bothering.
			print("Cannot exclude TSS positions with a BED file. Use a GTF if this option is needed.")
			print("Well, you can, but I'm guessing your file probably doesn't have strand info.")
			sys.exit()

		else:

			#Open provided gene file
			with open(gene_f) as f:

				#Iterate through file line by line, ignoring those without "gene" in the 3rd column
				for line in f:

					#Split line by tab
					line = line.strip().split("\t")

					chrom = line[0]
					start = int(line[1])
					stop = int(line[2])
					gene_symb = line[3]
					

					#Add the gene to the gene dictrionary
					gene_dict[gene_symb] = (chrom, (start, stop))

			#Return the gene dictionary
			return gene_dict


	else:

		#Check if the TSS positions are wanted
		if get_TSS:

			#Open provided gene file
			with open(gene_f) as f:

				#Iterate through file line by line, ignoring those without "gene" in the 3rd column
				for line in f:

					#Split line by tab
					line = line.strip().split("\t")

					#Check if "gene" is in the 3rd column
					if line[2].lower() == "gene":

						chrom = line[0]
						start = int(line[3])-1
						stop = int(line[4])

						#Get the gene symbol from the gtf file. 
						line_info = line[8]
						gene_name = line_info.split(";")[4]
						gene_symb_quotes = gene_name.strip().split()[1]
						gene_symb = gene_symb_quotes.strip('"')

						#Add the gene to the gene dictrionary
						gene_dict[gene_symb] = (chrom, (start, stop))


					#Check if "transcript" is in the 3rd column
					elif line[2].lower() == "transcript":

						chrom = line[0]

						#Check strand and use the appropriate start for the TSS
						strand = line[6]

						if strand == "+":
							TSS = int(line[3])
						else:
							TSS = int(line[4])

						#Add wings to each side of TSS
						TSS_up, TSS_down = Add_Wings_Mid(TSS, TSS_wing)
						

						#Get the transcript name from the gtf file. 
						line_info = line[8]
						transcript_name = line_info.split(";")[7]
						trans_name_quotes = transcript_name.strip().split()[1]
						trans_name = trans_name_quotes.strip('"')

						#Add the transcript to the TSS dictionary
						TSS_dict[trans_name] = (chrom, (TSS_up, TSS_down))

			#Return the gene dictionary
			return gene_dict, TSS_dict

		else:

			#Open provided gene file
			with open(gene_f) as f:

				#Iterate through file line by line, ignoring those without "gene" in the 3rd column
				for line in f:

					#Split line by tab
					line = line.strip().split("\t")

					#Check if "gene" is in the 3rd column
					if line[2].lower() == "gene":

						chrom = line[0]
						start = int(line[3]) - 1
						stop = int(line[4])

						#Get the gene symbol from the gtf file. 
						line_info = line[8]
						gene_name = line_info.split(";")[4]
						gene_symb_quotes = gene_name.strip().split()[1]
						gene_symb = gene_symb_quotes.strip('"')

						#Add the gene to the gene dictrionary
						gene_dict[gene_symb] = (chrom, (start, stop))

			#Return the gene dictionary
			return gene_dict


def Find_Genes_In_Wings(inp_file, gene_file, wing_size, output, midpoint=False, excl_proms=False, TSS_wing=2000, no_overlap=False, header=False):
	"""
	Strings together above functions to find genes that overlap the wing sizes from the midpoint or edges of each locus.

		Parameters: 
			(required) inp_file = Input file with locus positions and ID.
			(required) gene_file = Gene annotation file in gtf format.
			(required) wing_size = The length of the wings to add to each side of the midpoint or edges of the locus.
			(required) output = Name of the output file.
			(optional) midpoint = Determines whether the wings are added to the midpoint or the edges of the locus.
			(optional) header = Determines whether the input file has a header or not.
			(optional) excl_proms = Determines whether loci that overlap promoters should be excluded from output or not.
			(optional) TSS_wing = Size of the wings to add to each side of the TSS to define the promoter of each transcript.
			(optional) no_overlap = Determines if loci that overlap the gene are included in output or not. Default=False (included).
	"""

	#Get the positions for each gene in the gtf annotation file
	print("Pulling genes from annotation file. Enjoy doge in the meantime.")
	doge ="""\ 				─────────▄──────────────▄
				────────▌▒█───────────▄▀▒▌
				────────▌▒▒▀▄───────▄▀▒▒▒▐
				───────▐▄▀▒▒▀▀▀▀▄▄▄▀▒▒▒▒▒▐
				─────▄▄▀▒▒▒▒▒▒▒▒▒▒▒█▒▒▄█▒▐
				───▄▀▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▀██▀▒▌
				──▐▒▒▒▄▄▄▒▒▒▒▒▒▒▒▒▒▒▒▒▀▄▒▒▌
				──▌▒▒▐▄█▀▒▒▒▒▄▀█▄▒▒▒▒▒▒▒█▒▐
				─▐▒▒▒▒▒▒▒▒▒▒▒▌██▀▒▒▒▒▒▒▒▒▀▄▌
				─▌▒▀▄██▄▒▒▒▒▒▒▒▒▒▒▒░░░░▒▒▒▒▌
				─▌▀▐▄█▄█▌▄▒▀▒▒▒▒▒▒░░░░░░▒▒▒▐
				▐▒▀▐▀▐▀▒▒▄▄▒▄▒▒▒▒▒░░░░░░▒▒▒▒▌
				▐▒▒▒▀▀▄▄▒▒▒▄▒▒▒▒▒▒░░░░░░▒▒▒▐
				─▌▒▒▒▒▒▒▀▀▀▒▒▒▒▒▒▒▒░░░░▒▒▒▒▌
				─▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▐
				──▀▄▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▄▒▒▒▒▌
				────▀▄▒▒▒▒▒▒▒▒▒▒▄▄▄▀▒▒▒▒▄▀
				───▐▀▒▀▄▄▄▄▄▄▀▀▀▒▒▒▒▒▄▄▀
				──▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▀
				"""

	print(doge)

	#Check if promoter dictionary is needed, and get it if it is. If not, only return gene_dictionary
	if excl_proms:
		gene_dictionary, prom_dictionary = Gene_Grabber(gene_file, excl_proms, TSS_wing)
		
	else:
		gene_dictionary = Gene_Grabber(gene_file)

	#Open the output file
	output_file = open(output, "w")

	#Open the locus file
	with open(inp_file) as f:

		if header:
			header_line = f.readline().strip() + "\t" + "Gene_Symb"
			print(header_line, file=output_file)

		print("Comparing loci positions to gene positions. This will take a while for large datasets.")

		#Used to count the number with no genes that overlap
		no_ovlp_count=0
		ovlp_count=0

		genes_ovlp = 0
		loci_excl_count = 0

		#Iterate through file line by line
		for line in f:

			#Create list to hold genes for that locus ID
			gene_list = []

			line = line.strip().split()

			#Assign line elements to variables
			locus_chrom = line[chrom_col]
			locus_start = int(line[start_col])
			locus_stop = int(line[end_col])

			#Check if the midpoint or edges of the locus are to be used
			if midpoint:
				#Get the midpoint of the locus
				locus_midpoint = Get_Midpoint(locus_start, locus_stop)

				#Add wings to the midpoint
				winged_positions = Add_Wings_Mid(locus_midpoint, wing_size)
			else:
				winged_positions = Add_Wings_Edge(locus_start, locus_stop, wing_size)

			wing_start, wing_stop = winged_positions

			#Iterate through each gene to see if it overlaps the winged locus positions
			for gene in gene_dictionary:
				gene_position = gene_dictionary[gene]
				gene_chrom = gene_position[0]

				#Check if it's even on the same chromosome as the locus, if not, skip to the next gene
				if gene_chrom == locus_chrom:
					gene_start_stop = gene_position[1]
					gene_start = gene_start_stop[0]
					gene_stop = gene_start_stop[1]


					#Check that the original start/stop of the locus isn't wider than with the wings added, and if so, just use the original start/stop
					if int(locus_start) < int(wing_start) and int(locus_stop) > int(wing_stop):
						wing_start = locus_start
						wing_stop = locus_stop

					#Check for overlap of the two ranges and add the gene to the gene_list if this returns True.
					if Overlap(wing_start, wing_stop, gene_start, gene_stop):
						gene_list.append(gene)

			#Make sure all TSSs associated with the gene don't overlap the locus if excl_promoter==True.
			if excl_proms:

				ovlp=False

				#Iterate through the prom_dictionary to check each TSS associated with the gene
				for item in gene_list:
					for prom in prom_dictionary:
						if item in prom:
							prom_position = prom_dictionary[prom]
							prom_start_stop = prom_position[1]
							prom_start = prom_start_stop[0]
							prom_stop = prom_start_stop[1]

							#Check for overlap between the promoter and the locus, and if it occurs, just move to the next line
							ovlp = Overlap(locus_start, locus_stop, prom_start, prom_stop)

							if ovlp:
								
								break

					if ovlp:
						break

				if not ovlp:
					if len(gene_list) > 0:
						ovlp_count+=1
						genes_ovlp+=len(gene_list)
						line.append(";".join(gene_list))
						print(*line[0:], sep="\t", file=output_file)
					else:
						no_ovlp_count+=1
						line.append("NA")
						print(*line[0:], sep="\t", file=output_file)
				else:
					loci_excl_count += 1

			elif no_ovlp:

				overlap=False

				#Iterate through the gene_list to check for overlap
				for item in gene_list:

					gene_position = gene_dictionary[item]
					gene_chrom = gene_position[0]

					#Check if it's even on the same chromosome as the locus, if not, skip to the next gene
					if gene_chrom == locus_chrom:
						gene_start_stop = gene_position[1]
						gene_start = gene_start_stop[0]
						gene_stop = gene_start_stop[1]

						#Check for overlap between the promoter and the locus, and if it occurs, just move to the next line
						overlap = Overlap(locus_start, locus_stop, gene_start, gene_stop)

					if overlap:		
						break

				if not overlap:
					if len(gene_list) > 0:
						ovlp_count+=1
						genes_ovlp+=len(gene_list)
						line.append(";".join(gene_list))
						print(*line[0:], sep="\t", file=output_file)
					else:
						no_ovlp_count+=1
						line.append("NA")
						print(*line[0:], sep="\t", file=output_file)
				else:
					loci_excl_count += 1

			else:

				#Print to the output file
				#Check if there are genes in the gene_list, if not it will print NA
				if len(gene_list) > 0:
					ovlp_count+=1
					genes_ovlp+=len(gene_list)
					line.append(";".join(gene_list))
					print(*line[0:], sep="\t", file=output_file)
				else:
					no_ovlp_count+=1
					line.append("NA")
					print(*line[0:], sep="\t", file=output_file)

	if excl_proms:
		print("Number of loci excluded for overlap with a TSS: " + str(loci_excl_count))
	if no_overlap:
		print("Number of loci excluded for overlap with a gene: " + str(loci_excl_count))

	print("Number of loci with no overlapping genes: " + str(no_ovlp_count))
	print("Number of loci with overlapping genes: " + str(ovlp_count))

	avg_genes = genes_ovlp/ovlp_count
	print("Average overlapping genes per loci with overlap: " + str(avg_genes))

	output_file.close()



####-Variables-####

#Create arguments and options
parser.add_argument("-i", "--input", dest = "input_file", required=True)
parser.add_argument("-g", "--gene", dest = "gene_file", required=True)
parser.add_argument("-size", "--wingsize", dest ="wing_size", type=int, default=125000)
parser.add_argument("-C", "--Chrom", dest ="chrom", type=int, default=1)
parser.add_argument("-S", "--Start", dest ="start", type=int, default=2)
parser.add_argument("-E", "--End", dest ="end", type=int, default=3)
parser.add_argument("-o", "--output", dest = "output_file", required=True)
parser.add_argument("-m", "--midpoint", action="store_true")
parser.add_argument("-H", "--Header", action="store_true")
parser.add_argument("-p", "--promoter", action="store_true")
parser.add_argument("-twsize", "--tsswingsize", dest ="tss_wing", type=int, default=2000)
parser.add_argument("-noovlp", "--nooverlap", action="store_true")

args = parser.parse_args()

#Have user enter column numbers of various data
chrom_col = args.chrom - 1
start_col = args.start - 1
end_col = args.end - 1


#Make sure all the columns are appropriate, if not, yell at user and exit script.
if chrom_col < 0 or start_col < 0 or end_col < 0 or chrom_col == start_col or chrom_col == end_col or start_col == end_col:
	print("Column input should be 1 based and different numbers must be entered for each!")
	sys.exit()

#Output so user can double check options
print( "Input file {}\n Gene file {}\n Wing size {}\n Output {}\n\n".format(
        args.input_file,
        args.gene_file,
        args.wing_size,
        args.output_file
        ))

#Easier to use argument variables
inp_file = args.input_file
g_file = args.gene_file
wing_size = args.wing_size
out_file = args.output_file
header = args.Header
mid = args.midpoint
prom = args.promoter
tss_wing_size = args.tss_wing
no_ovlp = args.nooverlap



####-Main-####
#Run the main function
Find_Genes_In_Wings(inp_file, g_file, wing_size, out_file, mid, prom, tss_wing_size, no_ovlp, header)


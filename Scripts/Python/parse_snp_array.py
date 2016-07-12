#!/usr/bin/env python3

"""
Parse SNP array call table columns into individual VCF files using the SNP array annotation file.

Usage: python3 parse_snp_array.py <snp_annotations.csv> <snp_array_table.txt>
"""

#EDITED 01/07/2016 to only pull het calls (2), as there's no good way to determine whether allele A or allele B is the variant.

import sys

####-Functions-####
def get_ref_var_alleles(snp_annotations):
	""" Return position and ref/var alleles for each variant in the snp_annotations.

	Args:
		snp_annotations = SNP annotation file file.

	Returns:
		rs_pos = Dict containing variant rsid and it's position.
		rs_alleles = Dict containing ref allele and variant allele for a specific 
	"""
	print("Getting variant positions and alleles.")

	# Store chromosome and variant position.
	rs_pos = {}
	rs_alleles = {}

	# Open provided file
	with open(snp_annotations) as f:

		# Iterate through file line by line.
		for line in f:
			if line.startswith("#") or line.startswith("P"):
				continue
			line = line.strip().split(",")
			rs_id = line[1].strip('"')
			chrom = line[2].strip('"')
			pos = line[3].strip('"')
			allele_a = line[8].strip('"')
			allele_b = line[9].strip('"')
			rs_pos[rs_id] = (chrom,pos)
			rs_alleles[rs_id] = (allele_a,allele_b)

	return rs_pos, rs_alleles
	

def check_call(call_number):
	"""Check if line needs to be printed or if variant is not het."""
	if int(call_number) == 2:
		return True
	else:
		return False


def get_calls(snp_call_table):
	"""Get all calls from snp table for each column. Build dictionary for each and return a list of all of them.
	Returns a list of samples based on the header as well.
	"""
	print("Collecting calls for each sample.")
	# Hold all dicts.
	master_list = []

	# Open file
	with open(snp_call_table) as f:
		header = f.readline()
		header = header.strip().split("\t")
		samples = header[1:-3]

		# Initialize dictionaries.
		for entry in samples:
			new_dict = {}
			master_list.append(new_dict)

		# Iterate through lines.
		for line in f:
			line = line.strip().split("\t")
			rs_id = line[-3]
			calls = line[1:-3]

			# Add call to appropriate dict.
			for i in range(0,len(master_list)):
				master_list[i][rs_id] = calls[i]

	return samples, master_list
	


####-Main-####
# Grab input files
snp_anno = sys.argv[1]
snp_table = sys.argv[2]

# Get variant positions and alleles.
rs_positions, rs_alleles = get_ref_var_alleles(snp_anno)
# Get sample names and calls for each.
samples, calls = get_calls(snp_table)

# Print to output
print("Print to output files.")
for entry in samples:
	idx = samples.index(entry)
	sample_calls = calls[idx]
	output = open(entry+"_variants.vcf","w")
	count = 0

	# Print necessary header lines.
	print("##fileformat=VCFv4.2", file=output)
	print("##source=parse_snp_array.py", file=output)
	print("##reference=hg19", file=output)
	print('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">', file=output)
	print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=output)
	print("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT" + "\t" + snp_table + "_" + entry, file=output)

	for item in sample_calls:
		call = sample_calls[item]
		# Check the call to see if variant is there (2 or 3)
		if check_call(call):
			# Get position and other info if it is.
			try:
				var_position = rs_positions[item]
				chrom = "chr"+var_position[0]
				pos = var_position[1]
				var_alleles = rs_alleles[item]
				ref_allele = "."
				alt_allele = "."

				# ref_allele = var_alleles[0]
				# alt_allele = var_alleles[1]
			except KeyError as e:
				print("rs_id ",e," not found in annotation file, skipping.")
				count += 1
				continue

			# if int(call) == 2:
			# 	freq = "AF=0.5"
			# 	genotype = "0/1"
			# else:
			# 	freq = "AF=1.0"
			# 	genotype = "1/1"
			freq = "AF=0.5"
			genotype = "0/1"


			# Actually print info for variant to output.
			print(chrom,pos,item,ref_allele,alt_allele,".",".",freq,"GT",genotype,sep="\t",file=output)

	print("Number skipped: " + str(count))

	# Close output
	output.close()


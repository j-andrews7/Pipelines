#!/usr/bin/env python3
"""
Determines which genomic features are unique to a group of samples. Lists within script specify groups - modify as needed.

Usage: python3 unique_SEs_FLDL_subgroups.py <input.bed> 
"""

import sys

input_file = sys.argv[1]


#Modify these lists to change groupings.
# DL_ABC = ['DL188','DL237','DL252','DL273','DL551']
# DL_GCB = ['DL3A538','DL135']
DL = ['DL188','DL237','DL252','DL273','DL551','DL3A538','DL135']
# FL_1 = ['FL3A145','FL202','FL238','FL301','FL303','FL313']
# FL_2 = ['FL153','FL174','FL255']
FL = ['FL3A145','FL202','FL238','FL301','FL303','FL313','FL153','FL174','FL255']
# CLL = ['CLL140','CLL144','CLL347','CLL378','CLL396','CLL412','CLL580','CLL583','CLL589','CLL597','CLL600','CLL612','CLL618','CLL633','CLL634','CLL668','CLL674','CLL679']
CC_CB = ['CC011514','CC012214','CC012314','CC020514','CC021314','CB011514','CB012214','CB012314','CB20514','CB021314']
TS = ['TS012612B','TS050411','TS072111A','TS072611A']
# VGA = ['VGA3','VGA5','VGA8']
# VGR = ['VGR1','VGR2','VGR4','VGR6','VGR7','VGR10']

#Normal = ['CC011514','CC012214','CC012314','CC020514','CC021314','CB011514','CB012214','CB012314','CB20514','CB021314','VGA3','VGA5','VGA8','VGR1','VG2','VGR4','VGR6','VGR7','VGR10']
#Tumor = ['DL188','DL237','DL252','DL273','DL551','DL3A538','DL135','FL3A145','FL202','FL238','FL301','FL303','FL313','FL153','FL174','FL255','CLL140','CLL144','CLL347','CLL378','CLL396','CLL412','CLL580','CLL583','CLL589','CLL597','CLL600','CLL612','CLL618','CLL633','CLL634','CLL668','CLL674','CLL679']

# DL_ABC_out = []
# DL_GCB_out = []
DL_out = []
# FL_1_out = []
# FL_2_out = []
FL_out = []
#CLL_out = []
CC_CB_out = []
# VGA_out = []
# VGR_out = []
TS_out = []

#Normal_out = []
#Tumor_out = []

with open(input_file) as f:

	for line in f:

		unique = ''

		line = line.strip()
		data = line.split("\t")
		samples = data[3].split(";")

		test = samples[0]

		#Find proper group
		# if test in DL_ABC:
		# 	DL_ABC_out.append(line)
		# elif test in DL_GCB:
		# 	DL_GCB_out.append(line)
		if test in DL:
			DL_out.append(line)
		# elif test in FL_1:
		# 	FL_1_out.append(line)
		# elif test in FL_2:
		# 	FL_2_out.append(line)
		elif test in FL:
			FL_out.append(line)
		# elif test in CLL:
		# 	CLL_out.append(line)
		elif test in CC_CB:
			CC_CB_out.append(line)
		# elif test in CB:
		# 	CB_out.append(line)
		elif test in TS:
			TS_out.append(line)
		# elif test in VGA:
		# 	VGA_out.append(line)
		# elif test in VGR:
		# 	VGR_out.append(line)
		#elif test in Normal:
		#	Normal_out.append(line)
		#elif test in Tumor:
		#	Tumor_out.append(line)
		else:
			print(test + " not found in any group. Check groups in script to be sure it's included!")

# DL_ABC_file = open("Unique_DL_ABC_SEs.bed", "w")
# for entry in DL_ABC_out:
# 	print(entry,file=DL_ABC_file)
# DL_ABC_file.close()

# DL_GCB_file = open("Unique_DL_GCB_SEs.bed", "w")
# for entry in DL_GCB_out:
# 	print(entry,file=DL_GCB_file)
# DL_GCB_file.close()

DL_file = open("Unique_DL_SEs.bed", "w")
for entry in DL_out:
	print(entry,file=DL_file)
DL_file.close()

# FL_1_file = open("Unique_FL_1_SEs.bed", "w")
# for entry in FL_1_out:
# 	print(entry,file=FL_1_file)
# FL_1_file.close()

# FL_2_file = open("Unique_FL_2_SEs.bed", "w")
# for entry in FL_2_out:
# 	print(entry,file=FL_2_file)
# FL_2_file.close()

FL_file = open("Unique_FL_SEs.bed", "w")
for entry in FL_out:
	print(entry,file=FL_file)
FL_file.close()

# CLL_file = open("Unique_CLL_SEs.bed", "w")
# for entry in CLL_out:
# 	print(entry,file=CLL_file)
# CLL_file.close()

CC_CB_file = open("Unique_CC_CB_SEs.bed", "w")
for entry in CC_CB_out:
	print(entry,file=CC_CB_file)
CC_CB_file.close()

TS_file = open("Unique_TS_SEs.bed", "w")
for entry in TS_out:
	print(entry,file=TS_file)
TS_file.close()

# VGA_file = open("Unique_VGA_SEs.bed", "w")
# for entry in VGA_out:
# 	print(entry,file=VGA_file)
# VGA_file.close()

# VGR_file = open("Unique_VGR_SEs.bed", "w")
# for entry in VGR_out:
# 	print(entry,file=VGR_file)
# VGR_file.close()

# Normal_file = open("Unique_Normal_SEs.bed", "w")
# for entry in Normal_out:
# 	print(entry,file=Normal_file)
# Normal_file.close()

# Tumor_file = open("Unique_Tumor_SEs.bed", "w")
# for entry in Tumor_out:
# 	print(entry,file=Tumor_file)
# Tumor_file.close()



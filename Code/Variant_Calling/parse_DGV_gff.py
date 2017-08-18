#!/usr/bin/env python3
"""Parse output from DGV gold standard variants GFF file.

Usage: python3 parse_DGV_gff.py > output 
"""

import sys

for line in sys.stdin:

	line = line.strip().split("\t")

	data = line[8]
	select = data.split(";")[3]
	cnv_type = select.split("=")[1]

	out = "\t".join(line[0:8]) + "\t" + cnv_type + "\n"

	sys.stdout.write(out)
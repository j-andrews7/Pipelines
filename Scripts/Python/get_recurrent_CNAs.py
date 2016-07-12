#!/usr/bin/env python3
"""For a merged, annotated CNV file, pull those that are recurrent. Assumes samples in 14th column and ',' as delimiter.  

Usage: $​ ​python3 get_recurrent_CNAs.py <stdin> > output.bed
"""

import sys

for line in sys.stdin:

	l = line.strip().split("\t")
	data = l[3]
	samples = data.split(",")
	
	if len(samples) >= 2:
		sys.stdout.write(line)

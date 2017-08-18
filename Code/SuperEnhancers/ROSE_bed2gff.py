#!/usr/bin/env python3
"""
For a .bed file after calling peaks, merging, cutting, and sorting, convert
this file to gff format for use with ROSE.

Usage: python3 ROSE_bed2gff.py <input.bed>

Args:
    input.bed = Name of input file to process.
"""

import sys

# Store file names
input_file = sys.argv[1]
base = ".".join(input_file.split(".")[0:-1])
output_file = base + ".gff"

with open(input_file) as f:

    # Open output file, "w" to make it writable
    output_f = open(output_file, "w")

    # Initiate variable to hold what next ID should be
    new_ID = 1

    for line in f:
        line = line.strip().split()

        # Handle if bed file already has an ID that should be retained.
        if len(line) > 3:
            if line[3] is not None and line[3] is not ".":
                ID = line[3]
            else:
                ID = new_ID
        else:
            ID = new_ID

        print(line[0], ID, "RE", str(int(line[1]) + 1), line[2],
              ".", ".", ".", ID, sep="\t", file=output_f)

        # Increase ID by one for next line
        new_ID += 1

    output_f.close()

print("COMPLETE")

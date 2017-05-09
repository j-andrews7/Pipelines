#!/usr/bin/env python
"""
For a CNV bed file like so:

chr1    228694  356530  AOKACTB 9   3

Swap the value in column 6 to column 5, make the column 4 ID more sensical, and add rgb colors to
column 9 based on the new column 5 value (light to darker blue for amps, light to darker red for dels).

Final output will be like:

chr1    12862229    12867954    DL3A193_amp1 3  .   .  .    0,0,255

Usage: cnv_bed_color.py -i <input.bed> -o <output.bed>

Args:
    -i (str): Annotated CNV bed file to process.
    -o (str): Name of output file.
"""
import argparse
import sys


def colorize_line(line, amp_id, del_id):
    """
    Parses a line, removes crap, adds rgb colors to the 9th column based on values in 5th column.

    Returns:
        new_line (str): The new line with colors added and the columns swapped around as necessary.
        amp_id (int): ID to be used for next amplification.
        del_id (int): ID to be used for next deletion.
    """

    line = line.split()
    cn = int(line[5])
    line[4] = str(cn)

    if cn == 3:
        cn_id = line[3] + "_amp" + str(amp_id)
        amp_id += 1
        rgb = "0,0,255"
    elif cn > 3:
        cn_id = line[3] + "_amp" + str(amp_id)
        amp_id += 1
        rgb = "0,0,127"
    elif cn == 1:
        cn_id = line[3] + "_del" + str(del_id)
        del_id += 1
        rgb = "255,0,0"
    elif cn == 0:
        cn_id = line[3] + "_del" + str(del_id)
        del_id += 1
        rgb = "127,0,0"
    else:
        print("Remove normal copy number entries from file. Exiting.")
        sys.exit()

    line[3] = cn_id
    line.append(line[1])
    line.append(line[2])
    line.append(rgb)
    line[5] = "."

    new_line = "\t".join(line)

    return new_line, amp_id, del_id


def main(inp, out):

    # Starting IDs for dels and amps, incremented by line as necessary.
    amp_id = 1
    del_id = 1

    out_file = open(out, "w")
    with open(inp) as f:
        print("Parsing and colorizing.")

        for line in f:
            line = line.strip()
            new_line, amp_id, del_id = colorize_line(line, amp_id, del_id)
            print(new_line, file=out_file)

    out_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument("-i", "--input", dest="input_file", required=True)
    parser.add_argument("-o", "--output", dest="output_file", required=True)

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    main(input_file, output_file)

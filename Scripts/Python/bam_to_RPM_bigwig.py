#!/usr/bin/env python3
"""
Updated 06/28/2016
Author: jared.andrews07@gmail.com

Convert a bam file to a RPM normalized bigwig file. 

Usage: python3 bam_to_RPM_bigwig.py <input.bam> <output.bw> 
"""

import sys
import argparse
from pybedtools.contrib.bigwig import bam_to_bigwig
"""
Module to help create scaled bigWig files from BAM
"""
import pybedtools
import os
import subprocess

inp_fi = sys.argv[1]
out_fi = sys.argv[2]


def mapped_read_count(bam, force=True):
    """
    Scale is cached in a bam.scale file containing the number of mapped reads.
    Use force=True to override caching.
    """
    scale_fn = bam + '.scale'
    if os.path.exists(scale_fn) and not force:
        for line in open(scale_fn):
            if line.startswith('#'):
                continue
            readcount = float(line.strip())
            return readcount

    cmds = ['samtools',
            'view',
            '-c',
            '-F', '0x4',
            bam]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if stderr:
        raise ValueError('samtools says: %s' % stderr)

    readcount = float(stdout)

    # write to file so the next time you need the lib size you can access
    # it quickly
    if not os.path.exists(scale_fn):
        fout = open(scale_fn, 'w')
        fout.write(str(readcount) + '\n')
        fout.close()
    return readcount

def bedgraph_to_bigwig(bedgraph, genome, output):
    genome_file = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    cmds = [
        'bedGraphToBigWig',
        bedgraph.fn,
        genome_file,
        output]
    os.system(' '.join(cmds))
    return output

def wig_to_bigwig(wig, genome, output):
    genome_file = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    cmds = [
        'wigToBigWig',
        wig.fn,
        genome_file,
        output]
    os.system(' '.join(cmds))
    return output

def bam_to_bigwig(bam, genome, output, scale=False):
    """
    Given a BAM file `bam` and assembly `genome`, create a bigWig file scaled
    such that the values represent scaled reads -- that is, reads per million
    mapped reads.

    (Disable this scaling step with scale=False; in this case values will
    indicate number of reads)

    Assumes that `bedGraphToBigWig` from UCSC tools is installed; see
    http://genome.ucsc.edu/goldenPath/help/bigWig.html for more details on the
    format.
    """
    genome_file = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    kwargs = dict(bg=True, split=True, g=genome_file)
    if scale:
        readcount = mapped_read_count(bam)
        _scale = 1 / (readcount / 1e6)
        kwargs['scale'] = _scale
    x = pybedtools.BedTool(bam).genome_coverage(**kwargs)
    y = x.sort()
    cmds = [
        'bedGraphToBigWig',
        y.fn,
        genome_file,
        output]
    os.system(' '.join(cmds))



bam_to_bigwig(bam=inp_fi, genome='hg19', output=out_fi, scale=True)
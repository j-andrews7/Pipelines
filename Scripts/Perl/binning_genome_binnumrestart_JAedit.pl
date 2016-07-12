#! /usr/bin/perl -w

use strict;
use warnings; 

if (@ARGV != 2) {
    die "\nUsage: perl binning_genome_binnumrestart.pl Chrom_lengths.txt <output_filename.bed> \n\n";
}

print "How many nucleotides per bin?\n";
# 
# ##Variables##
my $binlength = <STDIN>;
my $binstop = $binlength;
my @data;
my $chroms = shift(@ARGV);
my $outputfile = shift(@ARGV);
my $chrom;
my $length;
my $start = 0;
my $stop;
my $bin = 0;

#opens the file with chromosome lengths or quits the program if it can't
open (IN, "<$chroms") or die ("Couldn't open file: $chroms\n");

open(OUT, ">$outputfile") or die ("Couldn't open file: $outputfile");

# reads in each line of the chromosome file
while (my $line = <IN>){
	# removes newline characters
	chomp($line);
	
	#puts both chromosome name and length into a small array.
	@data = split('\t', $line);
	
	#designates $chrom as the name of the chromosome
	#and length as the total length of the chrom in bp
	$chrom = $data[0];
	$length = $data[1];
	
	#while loop that executes as long as the starting nucleotide position is 
	#less than the total length of the chromosome
	while ($start < $length){
	
	#designates stop site as 199 nucleotides from start 
	$stop = $start + $binstop;
	
	print OUT "$chrom\t$start\t$stop\t$bin\n";
	$start = $stop;
	$bin++;
	}

$bin = 0;
$start =0;	
}

close IN;
close OUT;
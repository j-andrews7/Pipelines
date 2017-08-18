#! /usr/bin/perl -w

use strict;
use warnings; 

# if (@ARGV = 0) {
#     die "\nUsage: perl Calculating_FDR_Tags.pl <inputfiles> <outputfilename.txt>\n\n";
# }


my $output = pop(@ARGV);
my $lines=0;
my $total;
my $size = 0;
my $current = 0;
my $FDR = 0;
my $tags = 0;

my $tagave = 0;
my $FDRave = 0;
my $FDRmed = 0;
my $tagmed = 0;
my @FDR;
my @tags;
my $intline = 0;

# opens an outfile 
open (OUT, ">$output") or die ("Couldn't open outfile\n");

print OUT "Filename\tGenome Coverage\tNumber of peaks\tAverage Tag\tMedian Tag\tAverage FDR\tMedian FDR\n";

foreach my $file (@ARGV)
{
	# Opens the file or ends the program
	open (IN, "<$file") or die ("Couldn't open file: $file\n");
	
	# Takes all characters before the first "_" and saves it as $name. May only work if there 
	# are 2 underscores in the file name....
	my ($name) = $file =~ /^(\S+_.+)_/;

	
	# Iterates through your input file
	while (my $line = <IN>) {
	
		# removes endline characters
		chomp($line);	
		
		if ($line =~ m/^chr\S+/){
# 			print OUT "$line\n";
			# splits each line at tabs and saves tab-delimited variables in array @values
			my (@values) = split("\t", $line);
		
			###Count chromosome coverage###
			$current = $values[2] - $values[1] +1;
			$size = $size+ $current;
			
			###count lines###
			$lines++;
		
			##Calculate average FDR###
			$FDR = $FDR + $values[8];
			push(@FDR, $values[8]);
		
			###Calculate Average Tags###
			$tags = $tags + $values[5];
			push(@tags, $values[5]);
	
		}
	
	}

$FDRave = $FDR/$lines;
$tagave = $tags/$lines;
@FDR = sort { $a <=> $b }(@FDR);
@tags = sort { $a <=> $b }(@tags);
$intline = int($lines/2);

#prints the new rows into the outfile
print OUT "$file\t$size\t$lines\t$tagave\t$tags[$intline]\t$FDRave\t$FDR[$intline]\n";

$lines=0;
$size = 0;
$current = 0;
$FDR = 0;
$tags = 0;
$tagave = 0;
$FDRave = 0;
$FDRmed = 0;
$tagmed = 0;
$intline = 0;
undef @FDR;
undef @tags;

# closes both the in an out files before the next file is opened. 
close IN;
}


close OUT;
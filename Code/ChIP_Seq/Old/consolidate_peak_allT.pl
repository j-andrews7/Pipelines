#!/usr/bin/perl -w
use strict;

# declare the $usage
my $usage = "consolidate_peak_allT.pl <input> <output>
input file has to be ordered by chr, start, end
columns are peakid id chr start end length summit tags minustenlogp fe fdr delimited by tab
\n";

# terminate program if there are too many or too few command line inputs 
die $usage unless @ARGV == 2;

# store the command line inputs to the @ARGV array
my ($input, $output) = @ARGV;

# open input file for reading and output file for writing
open IN, $input || die "can not open $input\n";

open OUT, ">$output" || die "can not open $output\n";


my %ids;

# store all the samples ids in the %ids hash
while(<IN>)
{
	chomp;

	next unless /^\d/;
	
	#PKM_PEAKID	PKM_ID	PKM_CH	PKM_START	PKM_STOP	PKM_LENGTH	PKM_SUMMIT	PKM_TAGS	PKM_MINUSTENLOGP	PKM_FOLDENRICHMENTNUMBER	PKM_FDR

	my ($peakid, $id, $ch, $st, $ed, $length, $summit, $tags, $p, $fe, $fdr) = split;
	
	$ids{$id}++;
}

close IN;


# sort all the sample ids and store them in the @ids array
my @ids = sort(keys(%ids));

my @zero;

# create an array of all 0's whose length is the same as the number of samples
for(my $i=0; $i<scalar(@ids); $i++)
{
	push @zero, 0;
}

open IN, $input || die "can not open $input\n";

my %ct; # to store the number of peaks from each sample id in a consolidated region

my %fextags; #to store the fextags for each sample id in a consolidated region

my $prev_st = 0;

my $prev_ed = 0;

my $prev_ch = '0';

# print the header line
print OUT join("\t", 'CH', 'ST', 'ED', @ids, @ids, 'TOTAL'),"\n";

# main algorithm
while(<IN>)
{
	next unless /^\d/;	# skip the headline
	
	# remove the new line character at the end of the line
	chomp;
	
	# chop the input line and store each element in an variable
	my ($peakid, $id, $ch, $st, $ed, $length, $summit, $tags, $p, $fe, $fdr) = split;
			
	# if this is a new chromosome
	if($ch ne $prev_ch)
	{
		# if this is not the first line
		if( $prev_ch ne '0' )
		{
			print OUT join("\t", $prev_ch, $prev_st, $prev_ed, @ct{@ids}, @fextags{@ids}, samplenumber(@ct{@ids})), "\n";
		}
		
		# update new consolidated segment info
		$prev_ch = $ch;
		
		$prev_st = $st;
		
		$prev_ed = $ed;
		
		# reset and record counts, fextags, etc
		@ct{@ids} = @zero;
		
		@fextags{@ids} = @zero;
		
		$ct{$id}++;
		
		$fextags{$id}+=$tags; # taking sum of all peaks in consolidated region for a given sample
		
		# $fextags{$id}=$fextags{$id}>$qn?$fextags{$id}:$qn; # taking the max
	}
	# if this is still the same chromosome
	else
	{
		# if the current segment overlaps with the consolidated segment
		if( $st <= $prev_ed and $st >= $prev_st )
		{	
			# update consolidated segment info
			$prev_ed = $ed>$prev_ed?$ed:$prev_ed;
			
			# record counts, fextags, etc
			$ct{$id}++;
			
			$fextags{$id}+=$tags;
		}
		# if the current segment does not overlap with the consolidated segment
		else
		{
			print OUT join("\t", $prev_ch, $prev_st, $prev_ed, @ct{@ids}, @fextags{@ids}, samplenumber(@ct{@ids})), "\n";
			
			# update new consolidated segment info
			$prev_ch = $ch;
		
			$prev_st = $st;
			
			$prev_ed = $ed;
			
			# reset and record counts, fextags, etc
			@ct{@ids} = @zero;
			
			@fextags{@ids} = @zero;
		
			$ct{$id}++;
			
			$fextags{$id}+=$tags;
		}
	}
}

print OUT join("\t", $prev_ch, $prev_st, $prev_ed, @ct{@ids}, @fextags{@ids}, samplenumber(@ct{@ids})), "\n";

# subroutine to calculate the number of distinct samples in a consolidated segment
sub samplenumber
{
	my @array = @_;
	
	my $sum;
	
	foreach my $val (@array)
	{
		$sum ++ if $val > 0;
	}
	
	return $sum;
}
#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $usage = "bin_whole_genome.pl <input> <output>\n";

# die $usage unless @ARGV == 2;


my $input;
my $output;

foreach my $file (@ARGV) {
	$input = $file;
	my ($newname) = $file =~ /^(\S+).bed/;
	my ($fileholder) = $file =~ /^(\S+).wig/;
	$output = "$fileholder".".bin";

my %ch_size;

read_ch_size( 'hg19.chrom_sizes_trimmed', \%ch_size);

my %binned;

#my ($input, $output) = @ARGV;

open OUT, ">$output" || die "can not open output\n";

open IN, $input || die "can not open $input\n";
	
my $ch;
	
while(<IN>)
{
	next if /^track/;
	
	#variableStep chrom=chr1 span=10
		
	#print "$_\n";
		
	if(/^var/)
	{
		($ch) = /=chr(\S+)/;
			
		#last if $ch==10;
	}
	else
	{
		s/\s+$//;
			
		my ($po, $depth) = split;
			
		$binned{$ch}[int($po/200)] += $depth;
	}
}
	
close IN;
	

foreach my $ch (sort(keys %binned))
{
	for(my $i=0; $i<=int($ch_size{$ch}/200); $i++)
	{
		
		$binned{$ch}[$i] = 0 unless defined($binned{$ch}[$i]);
		print OUT join("\t", $ch, $i-1, $binned{$ch}[$i]), "\n";
	}
}


sub read_ch_size
{
	my ( $ch_size_file, $ch_size_ref) = @_;
	
	open CH, $ch_size_file || die "can not open $ch_size_file\n";
	
	while(<CH>)
	{
		s/\s+$//;
		
		my ($ch, $size) = split;
		
		$ch =~ s/^chr//;
		
		$ch_size_ref->{$ch} = $size;
	}
	
	close CH;
}
		
}	
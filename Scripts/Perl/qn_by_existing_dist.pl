#!/usr/bin/perl -w
use strict;

my $usage = "qn_by_existing_dist.pl <bin file> <dist file> <output>\n";
# take a sample bin file and a previously qned file and use is as an existing distribution to qn the bin file

die $usage unless @ARGV == 3;

my ($binfile, $distfile, $output) = @ARGV;

my %dist;

my @key;

open DIST, $distfile || die "can not open $distfile\n";

while(<DIST>)
{
	chomp;
	
	my ($ch, $bin, $tag) = split;
	
	$dist{$ch.'_'.$bin} = $tag;
	
	push @key, $ch.'_'.$bin;
}

my @sort_dist = sort{$a<=>$b}(values(%dist));

my %samplebin;

open BIN, $binfile || die "can not open $binfile\n";

while(<BIN>)
{
	chomp;
	
	my ($ch, $bin, $tag) = split;
	
	$ch = 23 if $ch eq 'X';
	
	$ch = 24 if $ch eq 'Y';
	
	$samplebin{$ch.'_'.$bin} = $tag if exists($dist{$ch.'_'.$bin});
}

foreach my $key(@key)
{
	$samplebin{$key} = 0 unless exists($samplebin{$key});
}

my %qn_samplebin;

foreach my $key (sort{$samplebin{$a}<=>$samplebin{$b}}(keys %samplebin))
{
	$qn_samplebin{$key} = shift @sort_dist;
}

open OUT, ">$output" || die "can not open $output\n";

foreach my $key(@key)
{
	print OUT join("\t", split(/_/,$key), $qn_samplebin{$key}), "\n";
} 




#!/usr/bin/perl -w
use strict;

my $usage = "qn.pl <input> <output>\n";

die $usage unless @ARGV == 2;

my ($input, $output) = @ARGV;

open IN, $input || die "can not open $input\n";

open OUT, ">$output" || die "cannot open $output\n";

my @tags;

my @sort_tags;

my @key;

while(<IN>)
{
	chomp;
	
	my ($ch, $bin, @data) = split;
	
	foreach(my $i=0; $i<@data; $i++)
	{
		$tags[$i]{$ch.'_'.$bin} = $data[$i];
	}
	
	push @key, $ch.'_'.$bin;
}

foreach(my $i=0; $i<@tags; $i++)
{
	@{$sort_tags[$i]} = sort{$tags[$i]{$a}<=>$tags[$i]{$b}}(keys %{$tags[$i]});
}

my @qn_tags;

for(my $i=0; $i<@{$sort_tags[0]}; $i++)
{
	my $qn_value;
	
	foreach(my $j=0; $j<@tags; $j++)
	{
		$qn_value += $tags[$j]{$sort_tags[$j][$i]};
	}
	
	$qn_value = sprintf("%.2f", $qn_value/scalar(@tags));
	
	foreach(my $j=0; $j<@tags; $j++)
	{
		$qn_tags[$j]{$sort_tags[$j][$i]} = $qn_value;
	}
}

foreach my $key (@key)
{
	print OUT join("\t", split(/_/,$key));
	
	foreach(my $j=0; $j<@tags; $j++)
	{
		print OUT "\t$qn_tags[$j]{$key}";
	}
	
	print OUT "\n";
} 

#! /usr/bin/perl -w

use strict;
use warnings; 


if (@ARGV != 2) {
    die "\nUsage: perl Add_XID.pl <mergefile.bed> <outfile.bed>\n\n";
}

#pulls files from command line
my $mergefile = shift(@ARGV);
my $outfile = shift(@ARGV);

print "ID name:\n";
my $ID = <STDIN>;

print "Column #:\n";
my $col = <STDIN>;

print "Do you want a header? (y/n)(0/1)\n";
my $ask= <STDIN>;

chomp($ID, $col, $ask);

my $newcol= $col-1;
my $newline;
my @data;
my $count = $ask;
my $me;

# Opens the file or ends the program
open (IN, "<$mergefile") or die ("Couldn't open file: $mergefile\n");
	
#Opens the file or ends the program
open (OUT, ">$outfile") or die ("Couldn't open file: $outfile\n");

while (my $line = <IN>){
	chomp($line);
	@data = split("\t", $line);
	
	if ($count==0){
		splice(@data, $newcol, 0, $ID);
		$newline = join ("\t", @data);
		print OUT "$newline\n";
	
	}else{
        $me = $count;
		splice(@data, $newcol, 0, $me);
		$newline = join ("\t", @data);
		print OUT "$newline\n";
	}
	$count++;	
}
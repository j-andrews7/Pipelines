#! /usr/bin/perl -w

use strict;
use warnings; 


if (@ARGV != 3) {
    die "\nUsage: perl condense_by_bin.pl <compare_file> <File_to_be_condensed> <output_filename>\n\n";
}

#pulls files from command line
my $XIDfile = shift(@ARGV);
my $tags_file = shift(@ARGV);
my $outfile = shift(@ARGV);
my %binhash;
my %IDhash;
my $count = 0;

###Variables###
my @data;

	# Opens the file or ends the program
	open (IN, "<$XIDfile") or die ("Couldn't open file: $XIDfile\n");

 	
while (my $line = <IN>){
	
	chomp($line);
	@data = split("\t", $line);
	
	#Change the first two values in $data to look at different column values in compare_file
 	$binhash{"$data[0]_$data[3]"} = "$data[1]_$data[2]";

 	#Save the peak ID of each bin in a hash for output
 	$IDhash{"$data[0]_$data[3]"} = "$data[4]";
}

# closes the input XID file
close IN;


#Opens the output file or ends the program
open (IN2, "<$tags_file") or die ("Couldn't open file: $tags_file\n");


#Opens the file or ends the program
open (OUT, ">$outfile") or die ("Couldn't open file: $outfile\n");

# reads the file line by line
while (my $line = <IN2>){
	
    #removes the newline character
	chomp($line);

	# saves the values from all columns in an array
	@data = split("\t", $line);

	#Prints header, inserts PEAK_ID column - Added 05/12/2015
	if ($count < 1){
		splice @data, 4, 0, 'PEAK_ID';
		print OUT join("\t", @data), "\n";
		$count++;
	}
   
    #change values in $data to look at different values in to_be_condensed file.
 	my $currtag = "$data[0]_$data[3]";
 	
	if (exists $binhash{$currtag}){
		splice @data, 4, 0, "$IDhash{$currtag}";
		print OUT join("\t", @data), "\n";
	}else{next;}
}

close IN2;
close OUT;

##Uncomment to print out the hash table...just so you can see what it looks like. 
# foreach my $key (keys(%binhash)){
# print "$key -> $binhash{$key}\n";
# }
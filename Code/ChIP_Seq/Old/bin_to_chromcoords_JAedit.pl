#!/usr/bin/perl -w
use strict;


if (@ARGV <2) {die "Usage: perl $0 <Binned_chromosomes3.txt> <any number of binned files with a header. Use * to run multiple files in a folder. ie \"*RPM.txt\">\n\nBinned_chromosomes3.txt has 200 bp bins starting at 1-200. bins restart at 0 for each chromosome. Can use any input chromsome bin file that has format \"chr# start stop bin#\""};

print "\nIn which column do your samples start?\n";
my $sample = <STDIN>;
chomp $sample;
$sample--;

print "Bin # column:\n";
my $bin = <STDIN>;
chomp $bin;
$bin--;

print "chr # column:\n";
my $chr = <STDIN>;
chomp $chr;
$chr--;

print "\nExcellent. Now leave me be. I have work to do....\n\n";

#open in the binned chromosomes file
my $chroms = shift(@ARGV);
open(IN, "<$chroms") or die ("Could not open file: $chroms\n");

##Variables##
my %beg;
my %end;

#iterate through chromsomes file and store a hash of the beginning and end positions for each chr/bin combinations
while (my $line = <IN>) {
    chomp($line);
    my @data = split("\t", $line);
    my $position = "$data[0]_$data[3]";
    $beg{$position} = $data[1];
    $end{$position} = $data[2];
    #print "$position\n";
}

close IN;

print "I have hashed your bins. Jeez...you mammals have a lot of genetic information....\n";

my $many =1;
foreach my $file (@ARGV){
    open(IN, "<$file") or die ("Could not open file: $file\n");
    open(OUT, ">chrcoords_$file") or die ("Could not open file: coords_$file\n");
    my $count=0;
    while (my $line = <IN>) {
        chomp($line);
        
        if ($count==0) {
            my @data = split("\t", $line);
            my $number = @data;
            $number--;
            my $newline = join("\t", $data[$chr], "START", "STOP", $data[$bin], @data[$sample..$number]);
            print OUT "$newline\n";
            $count=1;
            next;
        }
        my @data = split("\t", $line);
        my $number = @data;
        $number--;
        my $start = "$data[$chr]_$data[$bin]";
        my $end2 = "$data[$chr]_$data[$bin]";
        my $last = join("\t", @data[$sample..$number]);
        print OUT "$data[$chr]\t$beg{$start}\t$end{$end2}\t$data[$bin]\t$last\n";
    }

    print "\nI have completed $many files!\n";
    $many++;
    
close IN;
close OUT;
}

print "\nPlease check the new file header to insure that it still matches your new file output...I'm gonna go take a nap...\n";
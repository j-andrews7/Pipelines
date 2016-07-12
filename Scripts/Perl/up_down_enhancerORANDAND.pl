#!/usr/bin/perl -w
use strict;

if (@ARGV != 2){ die "\nUsage: perl $0 <promoter_file> <enhancerfile>\n\n"}

my $file =shift(@ARGV);
open(IN, "<$file") or die "Could not open file: $file\n";

my %prom;
while (my $line=<IN>) {
    chomp($line);
    my @data = split("\t", $line);
    $prom{$data[14]} = 0;
}

close IN;

#MMPID   CHR     MMPIDST MMPIDED 5KMMPIDST       5KMMPIDED       2KBTSSST        2KBTSSED        TSS     SAMPLE_ID       FC_X    FC_H3   FC_K27  ABS_K4  GENE    EXPRESSION



my $file2 =shift(@ARGV);
open(IN, "<$file2") or die "Could not open file: $file2\n";

open(UP, ">up_$file2") or die "Could not open file: up_$file\n";
open(DOWN, ">down_$file2") or die "Could not open file: down_$file\n";

while (my $line=<IN>){ 
    chomp($line);
    ## current problem: can have the same enhancer considered up and down.
    my @data = split("\t", $line);
    if (($data[10]>=2 && $data[11]>=1.5 && $data[12]>=1.5) || ($data[10]>=2 && $data[11]>=1.5 && $data[12]==0) || ($data[10]>=2 && $data[11]==0 && $data[12]>=1.5) || ($data[10]>=2 && $data[11]==0 && $data[12]==0)){ #  || $data[13]>0
        if (exists $prom{$data[14]}) {
            print UP "$line\n";
        }
        

    }if (($data[10]<= 0.5 && $data[11]<= 0.6 && $data[12]<= 0.6 && $data[10] != 0 && $data[11] != 0 && $data[12] != 0)|| ($data[10]<= 0.5 && $data[11]<= 0.5 && $data[12]==0 && $data[10] != 0 && $data[11] != 0) || ($data[10]<= 0.5 && $data[11]== 0 && $data[12]==0 && $data[10] != 0) || ($data[10]<= 0.5 && $data[11]==0 && $data[12]<= 0.6 && $data[10] != 0 && $data[12] != 0)){

        if (exists $prom{$data[14]}) {
        print DOWN "$line\n";
        }
        
        

    }
    
}# $data[10] != 0 && $data[11] != 0 && $data[12] != 0

close IN;
close UP;
close DOWN;

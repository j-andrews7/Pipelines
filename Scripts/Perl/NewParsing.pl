#!/usr/bin/perl -w
use strict;

#This code takes the following as input:
#1) The massive data table that basically has one line for each possibl MMPID/gene interaction as well all the MMPID histone FCs and gene expression FC
#2) The limitBasic table of genes...basically any genes to which you would like to assign a promoter
#3) an additional list of genes that Jackie sent me at some point. Easy enough to cut down to just one reference gene file


if (@ARGV != 3){ die "\nUsage: perl $0 <DL_circuits_file> <Limitbasic_file> <Jackie's file>\n\n"}

my $file =shift(@ARGV);
open(IN, "<$file") or die "Could not open file: $file\n";

my $limit =shift(@ARGV);
open(LIM, "<$limit") or die "Could not open file: $limit\n";


open(PROM, ">promoter_$file") or die "Could not open file: promoter_$file\n";

open(ENHAN, ">enhancer_$file") or die "Could not open file: enhancer_$file\n";


my %lim;

while (my $line=<LIM>) {
    chomp($line);
    my @data = split("\t", $line);
    $lim{$data[4]} =0;
}

close LIM;

my $Jackie =shift(@ARGV);
open(JACK, "<$Jackie") or die "Could not open file: $Jackie\n";

while (my $line=<JACK>) {
    chomp($line);
    my @data = split("\t", $line);
    $lim{$data[1]} =0;
}

close JACK;

my %prom;
my %promout;
my $count = 0;
my $overlap =0;
while (my $line=<IN>) {
    chomp($line);
    my @data = split("\t", $line);

    

    
    if ($count==0) {
        $count++;
        next;
    }
    #next if ($data[15] == 0);
    next if(!defined $lim{$data[14]});
        
    my $ID = "$data[9]_$data[14]";
    my $mid = ((($data[3]-$data[2])/2) + $data[2]);
    my $start = abs($data[2] - $data[8]);
    my $stop = abs($data[3] - $data[8]);
    my $distance = abs($data[8]-$mid);

    #if ($data[2] <= $data[8] && $data[3] >= $data[8]) {
    #    if ($data[16] eq "+" && $distance > 0) {
    #    print "$line\n";
    #    $overlap++;
    #    }elsif ($data[16] eq "-" && $distance < 0) {
    #    print "$line\n";
    #    $overlap++;
    #    }
    #    
    #
    #}
    #
    ##This currently lets there be only one promoter for each gene per sample and decides by which MMPID is closer to the TSS
    ##Conflicting promoters are redirected to the enhancer file
    
    
    #$data[2] and $data[3] are start and end of MMPID $data[6] and $data[7] are 2KBTSS
    #if($x > $y && $x < $z)
    if (($distance <= 2000) or ($data[2] <= $data[6] and $data[3] >= $data[7]) or ($data[2] >= $data[6] and $data[2] <= $data[7]) or ($data[3] >= $data[6] and $data[3] <= $data[7]) or ($data[2] >= $data[6] and $data[3] <=$data[7])){
        if ((!defined $prom{$ID})){       
            $promout{$ID} = "$line\t$distance";
            #print PROM "$line\n";
            $prom{$ID} = $distance;
            next;#code
        }elsif (defined $prom{$ID}){
            if ($prom{$ID} <= $distance) {
                print ENHAN "$line\n";
                next;
            }else{
                print ENHAN "$promout{$ID}\n";
                $promout{$ID} = "$line\t$distance";
                #print PROM "$line\n";
                $prom{$ID} = $distance;
                next;
            }
        } 
    }else{print ENHAN "$line\n";}
        

    
    
}

foreach my $key (keys %promout){
    print PROM "$promout{$key}\n";
}

close IN;
close PROM;
close ENHAN;

open(PROM, "<promoter_$file") or die "Could not open file: promoter_$file\n";

my %prom2; 
while (my $line=<PROM>) {
    chomp($line);
    my @data = split("\t", $line);
    $prom2{$data[0]}=0;
}

close PROM;

open(ENHAN, "<enhancer_$file") or die "Could not open file: enhancer_$file\n";
open(OUT, ">noprom_enhancer_$file") or die "Could not open file: noprom_enhancer_$file\n";

while (my $line=<ENHAN>) {
    chomp($line);
    my @data = split("\t", $line);
    
    if (exists ($prom2{$data[0]})) {
        next;
    }else{print OUT "$line\n";}
    
}

close ENHAN;
close OUT;








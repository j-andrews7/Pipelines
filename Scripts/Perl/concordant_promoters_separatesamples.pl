#!/usr/bin/perl -w
use strict;

if (@ARGV != 2){ die "\nUsage: perl $0 <promoter_file> <up or down enhancerfile>\n\n"}

my $file =shift(@ARGV);
open(IN, "<$file") or die "Could not open file: $file\n";

##Has sample_Gene as keys and FAIRE\tH3AC as values
my %prom;
#hash with each of the samples
my %samples;
my %consamples;

my @all;
while (my $line=<IN>) {
    chomp($line);
    my @data = split("\t", $line);
    my $ID = "$data[9]_$data[14]";
    my $info = join("\t", @data[10..11]);
    $prom{$ID} = $info;

}

close IN;

#MMPID   CHR     MMPIDST MMPIDED 5KMMPIDST       5KMMPIDED       2KBTSSST        2KBTSSED        TSS     SAMPLE_ID       FC_X    FC_H3   FC_K27  ABS_K4  GENE    EXPRESSION


my $file2 =shift(@ARGV);
open(IN, "<$file2") or die "Could not open file: $file2\n";



my @promdata;
my $lines =0;
my %limit;
my $promoters=0;
my $concordant =0;
my $yessample=0;
my $yespromoter=0;


while (my $line=<IN>) {
    chomp($line);
    my @data = split("\t", $line);
    my $ID = "$data[9]_$data[14]";
    


    #next if defined($limit{$ID});
    $limit{$ID} = 0;

        
    if( exists $samples{$data[9]} ) {
        push @{ $samples{$data[9]} }, $data[15];
    } 
    else {
        $samples{$data[9]} = [$data[15]];
    }

    
    #push(@{$samples{$data[9]}}, $data[15]); 
    
    if (exists $prom{$ID}) {
        @promdata = split("\t", $prom{$ID})
    }else {print "Could not find promoter info for $ID\n"; next;}
    
    ##One should be able to be above 1.25 and the other not...
    # ($promdata[0]>1.25 || $promdata[1]>1.25)
    #($promdata[0]>1.25 && $promdata[1]>1.25) ||($promdata[0]>1.25 && $promdata[1]==0) || ($promdata[0]==0 && $promdata[1]>=1.25)
        
        
        
        if ($promdata[0]>=1.5 || $promdata[1]>=1.25){
        
           if( exists $consamples{$data[9]} ) {
            push @{ $consamples{$data[9]} }, $data[15];
            } 
            else {
                $consamples{$data[9]} = [$data[15]];
            }
            #push(@{$consamples{$data[9]}}, $data[15]);
            #print OUT "$line\t$promdata[0]\t$promdata[1]\n";
        }
        
        
        
        ##($promdata[0]< 0.8 && $promdata[1] < 0.8 && $promdata[0] != 0 && $promdata[1] != 0) || ($promdata[0]< 0.8 && $promdata[1]==0 && $promdata[0] != 0) || ($promdata[0]== 0 && $promdata[1] < 0.8 && $promdata[1] != 0)
        #if (($promdata[0]<= 0.66 && $promdata[0] != 0) || ($promdata[1]<=0.8 && $promdata[1] != 0)){
        #
        #        if( exists $consamples{$data[9]} ) {
        #        push @{ $consamples{$data[9]} }, $data[15];
        #        } else {
        #            $consamples{$data[9]} = [$data[15]];
        #        }
        #   #push(@{$consamples{$data[9]}}, $data[15]);
        #   #print OUT "$line\t$promdata[0]\t$promdata[1]\n";
        #    
        #}
        
        
        
$lines++;
}

close IN;
open(OUT, ">concordant_$file2") or die "Could not open file: concordant_$file\n";

my $longest=0;
my $largest;

foreach my $sample (keys %samples){
    my $now = scalar(@{$samples{$sample}});
    if ($now > $longest) {
        $longest = $now;
        $largest= $sample;
    }else{next;}
    
}

###print out the header for the outputfile. 
foreach my $sample (keys %samples){
    print OUT "$sample\t";
}
foreach my $sample (keys %consamples){
    print OUT "con_$sample\t";
}

print OUT "END\n";

#print scalar(%samples)."\t".scalar(%consamples)."\n";

for (my $i = 0; $i<$longest; $i++){
    
    foreach my $sample (keys %samples){
  
        if (defined(@{$samples{$sample}}[$i])) {
            print OUT "@{$samples{$sample}}[$i]\t";
            $promoters++;
        }else{print OUT "\000\t"};
        
        
    }
    foreach my $sample (keys %consamples){
        if (defined(@{$consamples{$sample}}[$i])) {
            print OUT "@{$consamples{$sample}}[$i]\t";
                        $concordant++;
        }else{print OUT "\000\t"};
    }
    
    
print OUT "END\n";

    
}
close OUT;

#foreach my $key (keys %consamples){
#    print "$key => @{$consamples{$key}}\n";
#}


$promoters = keys %prom;


print "not really Total # promoters: $promoters
Number of concordant promoters: $concordant
Total number of lines: line: $lines\n";














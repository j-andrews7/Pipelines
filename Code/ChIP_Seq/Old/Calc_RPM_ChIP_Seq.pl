#!/usr/bin/perl -w
use strict;

if (@ARGV != 2){
    die ("Usage: perl $0 <input bedfile with samples in 11th column> <output filename>");
}

my $file = shift(@ARGV);
open(IN, "<$file") or die ("Could not open file: $file");


my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die ("Could not open file: $outfile");

my $count = 0;
my @avdata;
my $average;
my $samp;
my %readshash;
my @samplenames;
my $total; #needs to be the total number of reads for each sample...
#### samtools view -c -F 4 file.bam This should find only mapped reads. Used mapped reads #s from each sample summary.
my $curr =0;
##Making hash with # of mapped reads for each file##

%readshash = (
BMT267_K27AC => 15938476,
CLL580_K27AC => 17769286,
CLL583_K27AC => 18610398,
CLL589_K27AC => 11359039,
CLL597_K27AC => 24075469,
CLL600_K27AC => 19504552,
CLL612_K27AC => 23844378,
CLL618_K27AC => 20676525,
CLL633_K27AC => 18141975,
CLL634_K27AC => 19189274,
CLL668_K27AC => 21331994,
CLL674_K27AC => 22124490,
CLL679_K27AC => 17312598,
MEC1_K27AC => 15936861,
TS081414_MEMORY_K27AC => 21049934,
TS081414_NAIVECD5P_K27AC => 17286175,
TS081414_NAIVE_K27AC => 16894383,
TS102214_MEMORY_K27AC => 5730045,
TS102214_NAIVE_K27AC => 20880958,
CB011514_K27AC => 15109336,
CB012214_K27AC => 18832183,
CB012314_K27AC => 15629200,
CB020514_K27AC => 22881954,
FL313_K27AC => 22362508,
HBL1_K27AC => 17942813,
HELA_K27AC => 25468005,
JURKAT_K27AC => 18097464,
K562_K27AC => 21517615,
KARPAS422_K27AC => 24897141,
LINE293T_K27AC => 23198705,
LN182_K27AC => 22651979,
OCILY3_K27AC => 19721935,
OCILY7_K27AC => 19588692,
PBB135_K27AC => 23563446,
PBB145_K27AC => 27039011,
PBB166_K27AC => 22740857,
PBB174_K27AC => 16121045,
PBB252_K27AC => 22465026,
PBB267_K27AC => 23129004,
PBB301_K27AC => 17944348,
PBB303_K27AC => 28255813,
PFEIFFER_K27AC => 19354105,
TS012612A_K27AC => 20957912,
TS012612B_K27AC => 29603736,
TS020212A_K27AC => 15303410,
TS020212B_K27AC => 24519931,
TS040611B_K27AC => 22458983,
TS050411_K27AC => 15105264,
TS060911A_K27AC => 21944247,
TS061611B_K27AC => 20372971,
TS072111A_K27AC => 20863821,
TS072611A_K27AC => 16191030,
TS072611B_K27AC => 22248278,
TS081111A_K27AC => 19087473,
TS102711A_K27AC => 26562874,
TS111711_K27AC => 21414537,
TS120310_K27AC => 20442629,
TS121511A_K27AC => 19406363,
TS121511B_K27AC => 14571779,
U2932_K27AC => 23616383,
VGA3_K27AC => 19917478,
CB021314_K27AC => 18154675,
CC011514_K27AC => 13507870,
CC012214_K27AC => 19190371,
CC012314_K27AC => 22256884,
CC020514_K27AC => 20426933,
VGA5_K27AC => 23132703,
VGA8_K27AC => 22320178,
VGR1_K27AC => 17518402,
VGR2_K27AC => 22671954,
VGR10_K27AC => 22094318,
VGR4_K27AC => 18306542,
VGR6_K27AC => 14736433,
VGR7_K27AC => 21634454,
WSUNHL_K27AC => 22276689,
CC021314_K27AC => 21310686,
CLL140_K27AC => 20816590,
CLL144_K27AC => 24541627,
CLL347_K27AC => 13952523,
CLL378_K27AC => 12617269,
CLL396_K27AC => 14590302,
CLL412_K27AC => 12193907,
DB_K27AC => 23187020,
DL135_K27AC => 24497448,
DL188_K27AC => 15911075,
DL3A538_K27AC => 12281752,
DL237_K27AC => 16060190,
DL252_K27AC => 18780144,
DL273_K27AC => 21131044,
DL551_K27AC => 14623405,
DOHH2_K27AC => 18410519,
FARAGE_K27AC => 22635226,
FL3A145_K27AC => 21175580,
FL153_K27AC => 20277502,
FL174_K27AC => 20577216,
FL202_K27AC => 19830791,
FL238_K27AC => 22447267,
FL255_K27AC => 21984828,
FL301_K27AC => 22982759,
FL303_K27AC => 21206622
);





while (my $line = <IN>) {
    chomp($line);

    if ($count >= 1) {
      
        my @data = split ("\t", $line);
        
        my $hi = @data;
        my @samples = @data[4..($hi-1)]; #JARED: Change from 10 to whatever column the data starts in 
        
        
        undef(@avdata);
        $curr =0;
        foreach $samp (@samples){
            $samp = $samp+1; ####added Jan 23 2015
            if (exists $readshash{$samplenames[$curr]}){
                $total = $readshash{$samplenames[$curr]};
            }else{
                print "Could not find data for sample: $samplenames[$curr]\n";
            }
            
            
            my $RPKM = sprintf("%.3f", ($samp/($total/1000000)));
            push(@avdata, $RPKM);
            $curr++;
        }
    
        unshift(@avdata, @data[0..3]); #JARED: This adds the data from the non/normalized columns back on, so just change it to 0..last-column-non-normalized-data. 
    
        my $newline = join("\t", @avdata);
        print OUT "$newline\n";
    }else{
        print OUT "$line\n";
        my @data = split("\t", $line);
        my $hi = @data;
        @samplenames = @data[4..($hi-1)]; #JARED: Change from 10 to whatever column the data starts in 
    }
$count++;
}
    
close IN;
close OUT;
                                                    
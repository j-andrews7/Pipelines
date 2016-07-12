#!/usr/bin/perl -w
use strict;

if (@ARGV <= 1){ die "\nUsage: perl $0 <output filename> <input files>\n\n*NOTE*: you can input any number of input files using the * wildcard feature.\nPlease check the content of this script to make sure that your filenames match the current regex.\n\n"};

#assigns the first file on the command line as the outputfile
my $outputfile = shift(@ARGV);
#assigns the next file as the first input file. 
my $file =shift(@ARGV);
# opens the second file or prints an error message
open(IN, "<$file") or die "Could not open file: $file\n";

### Quick summary of what this program does. ###
# Takes the first file from the list and creates an index (hash) from it. Hashes have identifiers (keys) with associated data (values)
# in this case, the keys are stored in the variable $ID and are the chromosome and position of each snp. This means that you can search the hash
# by chromosome and position ($ID) and it will output the associated value. For the first file, the values are the tab-separated CHR POSITION and LOG2RATIO for
# each position for each $ID. In each subsequent file, the program will identify each SNP's chr and position, search the hash for the SNP and add the new LOG2ratio
# separated by a tab to the pre-existing value. An array called @header will keep track of the order in which files are added. @header will eventuall.y be printed out as column headers

# initializing a number of variables
my %snphash; # the % denotes that %snphash is a hash
my $filename; # This variable will hold the number/name of the file, ie 20248

#saves the file number as "filename" from the first file. 
if ($file =~ /(\S+)_accepted_hits.sorted_FPKM.txt/) { # regular expression that matches each filename. It's a built-in function of perl that the first string in parentheses in a regex can be referred to as $1 (second thing in parenthesis is $2 etc..)
    $filename = $1; # saves file name/number as $filename
}

## Opens the first file and creates a hash that has the chromosome and chromosome position as keys
## and "CHR POSITION LOG2RATIO" as values separated by tabs. 

##opens the first file
while (my $line = <IN>) { 
    
    #removes newline character from line
    chomp($line);
    # these files have \r\n windows line endings. This second line will take care of those endings
    $line =~ s/\r|\n//g;
    
    # splits the line at tabs and saves each tab-separated value in an array called @data.
    my @data = split("\t", $line);
    
    #ignores any line starting with a "#" or Probeset, which is the line with the current column headers
    next if ($line =~ /^\#/ || $line =~ /^gene_id/);
    

    
    # $ID is the chr and position separated by an underscore and is used to initialize a hash key for each position
    # These keys can later be search for within the hash to add new tag numbers for each new file
    my $ID = "$data[0]"."_"."$data[1]";
    my $value = $data[3];
    
    #saves the key/value pairs in the hash. 
    $snphash{$ID} = $value;
}

# Creates a header for the output file that has the descriptive columns and the data from the first file. 
my @header = ("chr","bin","$filename");

#closes the first input file. 
close IN;


# iterates through each additional file at the command line. 
foreach my $file2 (@ARGV){
    
    # saves the sample name as $filename
    if ($file2 =~ /(\S+)_accepted_hits.sorted_FPKM.txt/) {
        $filename = $1;
    }
    #adds filename to the header array 
    push(@header, $filename);
    
    # opens the next file
    open(IN, "<$file2") or die "Could not open file: $file2\n";
    
    #initializes a count hash 
    my %count;
   
   #reads through the file
    while (my $line = <IN>) {
        
        #ignores any line starting with a "#" or geneid
        next if ($line =~ /^\#/ || $line =~ /^gene_id/);
        
        # removes newline characters and splits the data at tabs into the array @data
        chomp($line);
        $line =~ s/\r|\n//g;
        my @data = split("\t", $line);
        
        
        my $ID = "$data[0]"."_"."$data[1]";
        
        # asks, if the current chr_bin location already exists in the hash made from the first file. If it does
        #the tag # from the current file is appended with a tab to the tag value that is already stored for the chr_bin key
        if (exists $snphash{$ID}) {
            #adding new tags with tab (.= concatenates rather than changing the value)
            $snphash{$ID} .= "\t$data[3]";
            
            #keeping count
            $count{$ID}++;
            
            # good idea to have something like this when using hashes, since a duplicate value will overwrite (or in this case append) extra values. 
            if ($count{$ID} == 2) {
                print "$ID appears twice!\n";
            }
        # prints an error if that chr_bin combination cannot be found in the current file
        }else{print "Could not find position $ID for sample $filename\n";}
        
        
    }
    # closes the current file
    close IN;

}

# opens the output file you want to write to.
open(OUT, ">$outputfile") or die "Could not open file: $outputfile\n";

#prints out the header array joined by tabs. 
my $finalheader = join("\t", @header);
print OUT "$finalheader\n";

# prints out all the chr_bin and tag key/value pairs. 
foreach my $key (keys %snphash){
        print OUT "$snphash{$key}\n";
}

close OUT;
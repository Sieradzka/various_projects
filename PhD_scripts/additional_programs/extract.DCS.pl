#!/usr/bin/env perl
use strict;
use warnings;
#----------------------------------------------------------------------------------------------------------------------
our $title = "E(eV)\n#Ang (deg)	Elastic	Rotational\n"; # a title of each file
my $beginning = '^(\d+\.?\d+)\s*$'; # the first line from which the text will be extracted
my $end1 = ""; # the final line of extracting 
my $end2 = ""; # the alternative final line of extracting
my $end3 = '^\s*$'; # the empty line
my $file = "DCS.dat"; # the file from which the lines will be extracted
my $output_file = "DCS"; # a file where the extracted lines go; the name of file will be changed by adding number at the end
my $directory = "";
#----------------------------------------------------------------------------------------------------------------------
print "beginning=$beginning, end=$end3, file=$file, output_file=$output_file\n";
#	system("mkdir $directory");
	parse($beginning, $end3, $file, $output_file);
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
sub parse { 	# the subroutine extracts lines between two patterns: $start and $stop1 or $stop2
		# $start - the first line from which the text will be extracted
		# $stop1 - the final line of extracting
		# $stop2 - the alternative final line of extracting
		# $filename - the file from which the lines will be extracted
		# $output_file - a file where the extracted lines go
  my ($start, $stop1, $filename, $output_file) = (@_);
  my $output;
  my $output_name = "";
  my $counter=1;
  my $found=0;
  my $numb=0.0;

  open(INPUT, $filename) or die "Can't open input file.\n";
  while (<INPUT>) {
#	chomp($_);
    if($_ =~ /$start/) {$numb = $1;}
#print "$start1\n";
#print "$stop1\n"; 
    # Find block of lines to extract                                                           
    if( /$start/ ... /$stop1/ ) {

        # Start of block                                                                       
        if( /$start/ ) {
    	    $output_name="${output_file}.$numb"; # To extract each line-set into separate files.
            open($output,'>'."$output_name") or die $!;
	    print $output '#'.$numb.$title; 
        }
        # End of block                                                                         
        elsif ( /$stop1/ ) {
            close($output);
            #$counter++;
            #$found = 0;
        }
        # Middle of block                                                                      
        else{
            #if($found == 0) { # This is to comment first line
            #    print $output '#'.$numb."\n".$title; 
            #    $found=1;
            #}
            #else {
		$_ =~ s/^\t//;
                print $output $_;
            #}
        }

    }
    # Find block of lines to extract                                                           

  }
  close(INPUT);
}

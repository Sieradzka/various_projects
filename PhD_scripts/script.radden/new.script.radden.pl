#!/usr/bin/env perl
use strict;
use warnings;
#----------------------------------------------------------------------------------------------------------------------
my $beginning = "Calculating densities for density matrix:"; # the first line from which the text will be extracted
my $end1 = "Integrated charge density for density matrix"; # the final line of extracting 
my $end2 = "Integrated charge density for density matrix"; # the alternative final line of extracting
my $end3 = '^\s*$'; # the empty line
my $file = "radden.thymine.out"; # the file from which the lines will be extracted
my $output_file = "density.matrix."; # a file where the extracted lines go; the name of file will be changed by adding number at the end
my $directory = "plots_thymine";
#----------------------------------------------------------------------------------------------------------------------
	system("mkdir $directory");
	parse($beginning, $end1, $end2, $file, $output_file);
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
sub parse { 	# the subroutine extracts lines between two patterns: $start and $stop1 or $stop2
		# $start - the first line from which the text will be extracted
		# $stop1 - the final line of extracting
		# $stop2 - the alternative final line of extracting
		# $filename - the file from which the lines will be extracted
		# $output_file - a file where the extracted lines go
  my ($start, $stop1, $stop2, $filename, $output_file) = (@_);
  my $output;
  my $output_name = "";
  my $counter=1;
  my $found=0;

  open(INPUT, $filename) or die "Can't open input file.\n";
  while (<INPUT>) {

    # Find block of lines to extract                                                           
    if( /^.*($start)/ ... /.*($stop1|$stop2).*/ ) {

        # Start of block                                                                       
        if( /$start/ ) {
    	    $output_name=sprintf("${output_file}%02d",$counter); # To extract each line-set into separate files.
            open($output,'>'."$directory/$output_name") or die $!;
        }
        # End of block                                                                         
        elsif ( /$stop1/||/$stop2/ ) {
            close($output);
            $counter++;
            $found = 0;
        }
        # Middle of block                                                                      
        else{
            if($found == 0) { # This is to comment first line
                print $output '#'.$_; 
                $found=1;
            }
            else {
                print $output $_;
            }
        }

    }
    # Find block of lines to extract                                                           

  }
  close(INPUT);
}


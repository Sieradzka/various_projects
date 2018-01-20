#!/usr/bin/env perl
use strict;
use warnings;
#----------------------------------------------------------------------------------------------------------------------
my $beginning = "EIGENPHASE SUMS"; # the first line from which the text will be extracted
my $end1 = "Symmetry:"; # the final line of extracting 
my $end2 = "<---------done:free_scattering_mod:free_scattering"; # the alternative final line of extracting
my $end3 = '^\s*$'; # the empty line
my $file = "log_file.0"; # the file from which the lines will be extracted
my $output_file = "BBBBBB"; # a file where the extracted lines go; the name of file will be changed by adding number at the end
my $symmetry  = 0;
#----------------------------------------------------------------------------------------------------------------------
#                         system("awk '/$start/ {flag=1;next} (/$end1/||/$end2/){flag=0} flag {print}' $file");
#			  system("sed -n -e '/$start/,/($end1||$end2)/p' $file | sed '/$start/d' >$output_file");
#----------------------------------------------------------------------------------------------------------------------
=begin comment
  if (open(INPUT, $file)) {
        while (my $line = <INPUT>) {
		if ($line =~ m/\s*Symmetry: (\d+)\s*/) {
			$symmetry = $1;
                }
	}
        close(INPUT);
  }
  else {die "  Error: can't open $file !\n";}
=end comment
=cut
#----------------------------------------------------------------------------------------------------------------------
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
    	    $output_name=sprintf("${output_file}_%02d.out",$counter); # To extract each line-set into separate files.
            open($output,'>'.$output_name) or die $!;
        }
        # End of block                                                                         
        elsif ( /$stop1/||/$stop2/ ) {
            close($output);
            $counter++;
            $found = 0;
        }
        # Middle of block                                                                      
        else{
#            if($found == 0) {				# The first line has a space which is removing
#                print $output (split(/ /))[1]; 	# Comment this line if you want to remove first line from extracting
#                $found=1;
#            }
#            else {
                print $output $_;
#            }
        }

    }
    # Find block of lines to extract                                                           

  }
  close(INPUT);
}
########################################################################################################################
########################################################################################################################
=begin comment
sub parse { # the subroutine extracts lines between two patterns: $start and $stop1 or $stop2
  my ($start, $stop1, $stop2) = (@_);
  my $printable = 0;
  my $on = 0;
  my $counter=1;

  while (<INPUT>) {
    if ($_ =~ /^.*($start).*/ ... /.*($stop1|$stop2).*/) {
	print OUTPUT $_;
	$printable =1;
	$on = 1;
        $counter++;
    } elsif($on){
	if ($_ =~ /($stop1|$stop2)/) {
	last;
#	 print OUTPUT $_;
	}
    } elsif ($_ =~ /^>/) {
	$printable = 0;
    } elsif ($printable) {
	print OUTPUT $_;
    }
  }
}
=end comment
=cut
########################################################################################################################
########################################################################################################################
#                map{print MYFILE2 "@$_\n"}@exponent;
#		 $line =~ s/\|((D|Q))M//;

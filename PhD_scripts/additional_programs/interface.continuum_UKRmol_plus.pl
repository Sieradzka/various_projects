#!/usr/bin/env perl
use strict;
use warnings;
#----------------------------------------------------------------------------------------------------------------------
  my $no = 6; # a number of l
  my @file = glob("fcoefl*");
  my @exponent = @_;
  my $line = "";
  my $output_file = "continuum.out"; # test.out
#-----------------------------
#=begin comment
  my $i = 0;
  foreach (@file) {
  	my $k = 1;
  	$exponent[$i][0] = "exponents(:,$i) =";
  	if (open(MYFILE, $_)) {
        	while ($line = <MYFILE>) {
			if ($line =~ m/\s*(0\.\d+)/) {
                        	$exponent[$i][$k] = sprintf("%02.6f,",$1);
				$k++;
                	}
        	}
        close(MYFILE);
  	} 
	$i++;
 } 
#  map{print  "@$_\n"}@exponent;
  open (MYFILE2, '>'.$output_file) or die $!;
  map{print MYFILE2 "@$_\n"}@exponent;
  close (MYFILE2);
#=end comment
#=cut
########################################################################################################################
########################################################################################################################
=begin comment
#  my @string = @_;
  my $file = "";
for (my $i = 0; $i le $no; $i++){
  my $k=1;
  $file = glob("fcoefl${i}_EM*");
  $exponent[$i][0] = "exponents(:,$i) =";
  if (open(MYFILE, $file)) {
        while ($line = <MYFILE>) {
		if ($line =~ m/\s*(0\.\d+)/) {
                        $exponent[$i][$k] = sprintf("%02.6f,",$1);
			$k++;
                }
        }
        		#@string[$i] = join(",", @exponent);
        		#print "@string[$i] \n";
        close(MYFILE);
  }
#  else {
#        die "  Error: can't open $file !\n";
#  }
}
  map{print  "@$_\n"}@exponent;
=end comment
=cut
########################################################################################################################
########################################################################################################################

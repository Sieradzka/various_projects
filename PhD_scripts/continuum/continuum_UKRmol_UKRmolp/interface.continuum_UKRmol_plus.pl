#script to rewrite continuum exponents from fcoefl* or swmol3.continuum.r* file to integrals input
#!/usr/bin/env perl
use strict;
use warnings;
#----------------------------------------------------------------------------------------------------------------------
  my $no = 6; # a number of l
  my @file = glob("$ARGV[0]"); #file containing the continuum exponents e.g. "swmol3.continuum.r13.l4" or "fcoefl*"
  my @exponent = @_;
  my $line = "";
  my $output_file = "$ARGV[1]"; #"continuum.out"; # test.out
#-----------------------------
#=begin comment
  if ($file[0] =~ /swmol3.continuum.r.*/){
	my @no_exp = @_;
	my $i = 0;
	my $j = 0;
	my $k = 1;
	if (open(MYFILE, $file[0])) {
        	while ($line = <MYFILE>) {
			if ($line =~ m/\s*jco\s*=[\s*d+\s*,]*/) {
				$line =~ s/\s*jco\s*=\s*//;
				@no_exp = split(/,/,$line);
			}
			if ($line =~ m/\s*(0\.\d\d+)/) {
                        	$exponent[$i][$k] = sprintf("%02.6f,",$1);
				if ($k < $no_exp[$i]) {$k++;}
				else {
					$exponent[$i][0] = "exponents(:,$i) =";
					$k = 1;
					$i++;
				}
                	}
        	}
        close(MYFILE);
  	} 
  	open (MYFILE2, '>'.$output_file) or die $!;
  	map{print MYFILE2 "@$_\n"}@exponent;
  	close (MYFILE2);  
  }
  else{
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
#  	map{print  "@$_\n"}@exponent;
  	open (MYFILE2, '>'.$output_file) or die $!;
  	map{print MYFILE2 "@$_\n"}@exponent;
  	close (MYFILE2);
  }
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

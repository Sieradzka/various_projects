#Purpose of this script is to find transmision moments of first symetry and calculate polarizability
#!/usr/bin/env perl
use strict;
use warnings;
#----------------------------------------------------------------------------------------------------------------------
  my $prefix = "/padata/gamma/users/asieradzka/";
  my $suffix = "pyridine_2water/pyridine_2water/cc-pVDZ.SE.14frozen.12active.40virtual.1states.r15.L6.old.double_prec/geom1"; #"$ARGV[0]"; 
  my $path = "$prefix/$suffix/outputs/target.molpro.out";
  my $line = "";
  my @geom = @_; #
  my @geom_matrix = @_; #to save 
  my $i = 0;

 #--------------Searching transmision moments in prop file---------------

  if (open(MYFILE, "$path")) {
        while ($line = <MYFILE>) {
		if ($line =~ m/Variable\s*Last\s*Current\s*Next.*/) {
  			$i = 0; # will save here the number of rows in '@geom_matrix' matrix


		}
		if ($line =~ m/([A-Z]).* \/ BOHR\s*([\+\-]?\d+\.\d+)\s*([\+\-]?\d+\.\d+)\s*([\+\-]?\d+\.\d+)\s*/) {

 #                   		$line =~ s/([\+\-]?\d+\.\d+\s*)//;
				$geom[$i][0] = $1;
				$geom[$i][1] = $3;

				$i++;
			}
		#}

        }
        close(MYFILE);
  }
  else {die "  Error: can't open target.molpro.out !\n";}
#  map{print "@$_\n"}@geom;

  for (my $k = 0; $k < $i; $k++) {
	my $l = $k % 3;
	my $m = $k / 3;
	$geom_matrix[$m][0] = $geom[$k][0].",";
	$geom_matrix[$m][$l+1] = $geom[$k][1].",";
  }

#--------------Saving '@geom_matrix' to geom_opt.out file---------------
        open (MYFILE2, ">$prefix/$suffix/geom_opt.out");
 	map{print MYFILE2 "@$_\n"}@geom_matrix;
        close (MYFILE2);
#----------------------------------------------------------------------------------------------------

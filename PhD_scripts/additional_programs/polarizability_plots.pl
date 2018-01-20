#Purpose of this script is to find transmision moments of first symetry and calculate polarizability
#!/usr/bin/env perl
use strict;
use warnings;
#----------------------------------------------------------------------------------------------------------------------
  my $prefix = "/home/physastro/asieradzka/aga/pyridine_water/polarizability.SEP";
  my $suffix = "$ARGV[0]"; #"pyridine/cc-pVDZ.pol.6frozen.15active.20virtual.4states.r15.L6.old.quadru_prec/geom1/changing.prop.out"; 
  my $path = "$prefix/$suffix/changing.prop.out";
  my $line = "";
  my @energy = @_; #the energy of state
  my @trans_matrix = @_; #to save transmision moments and symmetry of state
  my @sorted_trans = @_; #to sort energy from lowest to highest

 #--------------Searching transmision moments in prop file---------------
  my $i_e = 0; # will save here the number of rows in '@energy' matrix
  my $i = 0; # will save here the number of rows in '@trans_matrix' matrix
  if (open(MYFILE, "$path")) {
        while ($line = <MYFILE>) {
                if ($line =~ m/State No.\s*([\**\d+]).*([\+\-]?\d+\.\d+e[\+\-]\d+)\s*/) {
                $energy[$i_e][0] = $1;
                $energy[$i_e][1] = $2;
		$i_e++;
                }
                if ($line =~ m/^1\s*(\d+)*\s*([\+\-]?\d+\.\d+)*.*/) {
#$trans_matrix[$i][0] - the number of state - first component
#$trans_matrix[$i][1] - the number of symmetry - first component
#$trans_matrix[$i][2] - the number of state - second component
#$trans_matrix[$i][3] - the number of symmetry - second component 
#$trans_matrix[$i][4] - 12
#$trans_matrix[$i][5] - the number of transmision moment: 1 -dipole moment, 2 -quadruple moment
#$trans_matrix[$i][6] - the number of 'm' quantum number
#$trans_matrix[$i][7] - the value of transmision moment
                    $line =~ s/\s*(\d+)//;
                    for (my $k = 0; $k < 7; $k++) { 
                    	$line =~ s/\s*([\+\-]?\d+)//;
                        $trans_matrix[$i][$k] = $1;
		    } 
                    $line =~ s/\s*([\+\-]?\d+\.\d+e[\+\-]\d+)\s*//;
                    $trans_matrix[$i][7] = $1;
                    $i++;
                }
        }
        close(MYFILE);
  }
  else {die "  Error: can't open changing.prop.out !\n";}

  my $k_used = 0;
  for (my $k = 0; $k < $i; $k++) {
	if (($trans_matrix[$k][2] == 1) && ($trans_matrix[$k][3] == 0) && ($trans_matrix[$k][5] == 1) ){	#($trans_matrix[$k][1] == 0) && && ($trans_matrix[$k][6] == 1)
		for (my $l = 0; $l < 8; $l++) {
		$sorted_trans[$k_used][$l] = $trans_matrix[$k][$l];
		}
	my $state = $sorted_trans[$k_used][0];
	if ($state == 1){$sorted_trans[$k_used][8] = 0;}
	else{
		$sorted_trans[$k_used][8] = ($sorted_trans[$k_used][7])**2/$energy[$state-1][1]*27.21*2 + $sorted_trans[$k_used - 1][8];
	}
	$k_used++;
	}
  }
#--------------Saving '@trans_matrix' to polariz.out file---------------
                open (MYFILE2, ">$prefix/$suffix/polariz.out");
# 		map{print MYFILE2 "@$_\n"}@energy;
 		map{print MYFILE2 "@$_\n"}@sorted_trans;
#		$check == 1 ? print MYFILE2 "\n" : print MYFILE2 "" ;
               close (MYFILE2);
#----------------------------------------------------------------------------------------------------

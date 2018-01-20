#Purpose of this script is interfacing with MOLPRO output and produce output compatible with DENPROP output i.e. prop.out
#sub interface_molpro {
#  my ($path) = $_;
  my ($r_par) = @_;
  my $prefix = "pyridine_plus/cc-pVDZ.CAS.16frozen.8active.0virtual.33states.r-1CAS.geom3";
  my $suffix = "geom1/outputs/target.molpro.out"; 
  my $path = "$prefix/$suffix";
  my $multiplet = ""; #the name of multiple e.g. "Singlet", "Triplet"
  my $no_multiplet = $_; #the name of multiple e.g. 1, 3
  my $symmetry = 0; #to point out the symmetry of a state
  my $states = 0; #the number of a state
  my $energy = 0; #the energy of state
  my @energy_matrix = @_; #to save energy and symmetry of state
  my @sorted_energy = @_; #to sort energy from lowest to highest
  my $i = 0; # will save here the number of rows in '@energy_matrix' matrix
  my @dipole_moment = @_; #matrix where all information about dipole moment will be saved
  my @sorted_dipole_moment = @_; # sorted "@dipole_moment" matrix
  my $atoms = 0; #the number of atoms

 #--------------Searching geometry in MOLPRO output---------------
  system("sed -n -e '/geometry={/,/}/p' $path |sed '/ }/d' >$prefix/geom1/molpro.prop.out");
 #--------------Counting lines in molpro.prop.out file---------------
  if (open(MYFILE2, "molpro.prop.out")) {
        while ($line = <MYFILE2>) {
		$atoms++;
        }
        close(MYFILE2);
  }
  else {
        die "  Error: can't open molpro.prop.out !\n";
  }
 #--------------Searching and 'grepping' energies and dipole moments in MOLPRO output---------------
 #There are two places in MOLPRO output where it can be found: the first for optimised orbitals, the second for all.
 #This part of code is to find this second place (where are energies and dipole moments for all orbitals) and to save it in temporary file.

  system("awk '/Warning: orbitals not optimized!/,0' $path >temporary.molpro");
 #--------------Searching energies in temporary file---------------
  if (open(MYFILE, "temporary.molpro")) {
	$no_multiplet = -1; #to distinguish singlet states between triplet ones; in MOLPRO output singlet states are before the triplets, there is no direct way to distinguish the states.
        while ($line = <MYFILE>) {
                if ($line =~ m/\s*!MCSCF STATE\s*1\.1\s*Energy\s*([\+\-]?\d+\.\d+)/) {
                        $no_multiplet = $no_multiplet + 2;
                }
                if ($line =~ m/\s*!MCSCF STATE\s*(\d+)\.(\d+)\s*Energy\s*([\+\-]?\d+\.\d+)/) {
                        $energy_matrix[$i][0] = $1; #the number of state    $states = $1;
                        $energy_matrix[$i][1] = $2; #the number of symmetry $symmetry = $2;
                        $energy_matrix[$i][2] = sprintf("%-002.15e", $3); #the value of energy    $energy = $3/1000;
                        if ($no_multiplet == 1){
                                $multiplet = "Singlet";
                        }
                        else {
                                $multiplet = "Triplet";
                        }
                        $energy_matrix[$i][3] = $no_multiplet; #the name of multiple e.g. 1, 3
                        $energy_matrix[$i][4] = $multiplet; #the name of multiple e.g. "Singlet", "Triplet"
                        $i++;
                }
        }
        close(MYFILE);
  }
  else {
        die "  Error: can't open temporary.molpro !\n";
  }
#--------------Sorting and ordering matrix to match to prop.out file---------------
  @sorted_energy = sort {$a->[2] <=> $b->[2]}@energy_matrix; #sorting matrix according to energy (from lowest to highest)
  @energy_matrix = 0;
  for (my $k = 0; $k < $i; $k++) {
        $energy_matrix[$k][0] = sprintf("%2d",5);
        $energy_matrix[$k][1] = sprintf("%2d",$k+1);
        $energy_matrix[$k][2] = sprintf("%2d",0);
        $energy_matrix[$k][3] = sprintf("%2d",0);
        $energy_matrix[$k][4] = sprintf("%2d",$sorted_energy[$k][1] - 1); #the number of symmetry
        $energy_matrix[$k][5] = sprintf("%2d",$sorted_energy[$k][3]); #the name of multiple e.g. 1, 3
        $energy_matrix[$k][6] = sprintf("%2d",0);
        $energy_matrix[$k][7] = sprintf("%2d",0);
        $energy_matrix[$k][8] = $sorted_energy[$k][2]; #the value of energy
        $energy_matrix[$k][9] = "State No.";
        $energy_matrix[$k][10] = sprintf("%2d",$k+1);
        $energy_matrix[$k][11] = $sorted_energy[$k][4]; #the name of multiple e.g. "Singlet", "Triplet"
  }

#--------------Searching dipole and quadrupole moments in temporary file---------------
  my $j = 0; # number of all dipole moments
  my %DQ = ('D', 1,
	    'Q', 2,
	     );
  my %XYZ = ('X',   1, #ok
	     'Y',  -1,
	     'Z',   0,
	     'XX',  2,
	     'YY', -2,
	     'ZZ',  0,
	     'XY',  0,
	     'XZ',  1,
	     'YZ', -1,
	     );
  if (open(MYFILE, "temporary.molpro")) {
	$no_multiplet = -1; #to distinguish singlet states between triplet ones; in MOLPRO output singlet states are before the triplets, there is no direct way to distinguish the states.
        while ($line = <MYFILE>) {
               if ($line =~ m/\s*!MCSCF (trans|expec)\s*<(\d+)\.(\d+)\|(D|Q)M((X.*|Y.*|Z.*))\|(\d+)\.(\d+)>/) {
                     if ($line =~ m/\s*!MCSCF (trans|expec)\s*<1\..\|(D|Q)M.*\|1\.1>/) {
                            $no_multiplet = ($no_multiplet + 2)%4; #in MOLPRO output Singlet and Triplet states appearing periodically.
                    }
                    $dipole_moment[$j][0] = 1; #in prop.out file first value is always equal 1
                    for (my $i = 11; $i < 15; $i++) { # writing numbers of states and symmetry to columns which will not appear in output file
                    $line =~ s/\s*(\d+)//; 
                    $dipole_moment[$j][$i] = sprintf("%2d",$1); # $i = 11 - no. of first state; $i = 12 - no. of first state symmetry; $i = 13 - no. of second state; $i = 14 - no. of second state symmetry;
                    }
                    $dipole_moment[$j][5] = $atoms; # the number of atoms in molecule + 1
		    $line =~ s/\|((D|Q))M//;
                    $dipole_moment[$j][6] = sprintf("%2d",$DQ{$1}); # gives 1 if dipole, 2 if quadrupole
		    $line =~ s/((X.*|Y.*|Z.*))\|//;
                    $dipole_moment[$j][7] = sprintf("%2d",$XYZ{$1}); 
                    $line =~ s/\s*([\+\-]?\d+\.\d+)\s*//;
                    $dipole_moment[$j][8] = $1; # gives energy
                    $dipole_moment[$j][10] = $no_multiplet ; # gives number of spin
		    for (my $k = 0; $k < $i; $k++) { # changes number of state in MOLPRO notation to DENPROP notation
		    	if ($dipole_moment[$j][11] == $sorted_energy[$k][0] && $dipole_moment[$j][12] == $sorted_energy[$k][1] && $dipole_moment[$j][10] == $sorted_energy[$k][3]) {
		    	$dipole_moment[$j][1] = sprintf("%2d",$k+1);
			}
		    }
                    $dipole_moment[$j][2] = sprintf("%2d",$dipole_moment[$j][12]-1); # number of symmetry
		    for (my $k = 0; $k < $i; $k++) { # changes number of state in MOLPRO notation to DENPROP notation
		    	if ($dipole_moment[$j][13] == $sorted_energy[$k][0] && $dipole_moment[$j][14] == $sorted_energy[$k][1] && $dipole_moment[$j][10] == $sorted_energy[$k][3]) {
		    	$dipole_moment[$j][3] = sprintf("%2d",$k+1);
			}
		    }
                    $dipole_moment[$j][4] = sprintf("%2d",$dipole_moment[$j][14]-1); # number of symmetry
		    $j++;
        	}
	}
     	@sorted_dipole_moment = sort {$a->[12] <=> $b->[12]|| $a->[10] <=> $b->[10] || $a->[14] <=> $b->[14]|| $a->[11] <=> $b->[11] || $a->[13] <=> $b->[13]|| $a->[6] <=> $b->[6] || $a->[7] <=> $b->[7]  }@dipole_moment;
#     	@sorted_dipole_moment_triplet = sort {$a->[3] <=> $b->[3]|| $a->[1] <=> $b->[1] || $a->[0] <=> $b->[0] || $a->[2] <=> $b->[2] }@dipole_moment_triplet ;
  	close(MYFILE);
  }
  else {
	die "  Error: can't open temporary.molpro !\n";
  }

#--------------Saving '@energy_matrix' to molpro.prop.out file---------------
                open (MYFILE2, ">>$prefix/geom1/molpro.prop.out");
                    map{print MYFILE2 "@$_\n"}@energy_matrix;
		    for (my $k = 0; $k < $j; $k++) {
		    	for (my $l = 0; $l < 9; $l++) {
				if ($sorted_dipole_moment[$k][7]!= -2){
				#if (!($sorted_dipole_moment[$k][6] == 2 && $sorted_dipole_moment[$k][7] == 0)){
				print MYFILE2 "$sorted_dipole_moment[$k][$l] ";
				}
			}
				if ($sorted_dipole_moment[$k][7]!= -2){
				#if (!($sorted_dipole_moment[$k][6] == 2 && $sorted_dipole_moment[$k][7] == 0)){
				print MYFILE2 "\n";
				}
#			print MYFILE2 "\n";
		    }	
		
#                    map{print MYFILE2 "@$_\n\n\n"}@sorted_dipole_moment;
               close (MYFILE2);
#----------------------------------------------------------------------------------------------------
system ("rm temporary.molpro");
#return 1;
#}

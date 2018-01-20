# To do:
#  * idtarg in outerres input 
#  * to use run for previous geometry in initial_orbitals_occupation_guess 

use dirfile;

# Settings of swmol3 input for different symmetries
our %symmetries = (
  'D2h',   [ 3, "'X', 'Y', 'Z'" ],
  'C2v',   [ 2, "'X', 'Y'" ],
  'C2h',   [ 2, "'YZ', 'X'" ], # [ 2, "'XY', 'Z'" ], 
  'D2',    [ 2, "'XZ', 'YZ'" ],
  'C2',    [ 1, "'XY'" ],
  'Cs',    [ 1, "'Z'" ], # ZM change to X from original 'Z'
  'Ci',    [ 1, "'XYZ'" ],
  'C1',    [ 0, "''" ],
);

# Order of irreducible representations in Sweden codes
our %irred_repr = (
  'D2h',   [ "Ag", "B3u", "B2u", "B1g", "B1u", "B2g", "B3g", "Au" ], # incorrect in documentation of swmol3 (in output 2px orb. are in 2 = B3u and 3d1+ are in 6 = B2g)
  'C2v',   [ "A1", "B1", "B2", "A2" ],
  'C2h',   [ "Ag", "Bu", "Au", "Bg" ],
  'D2',    [ "A", "B3", "B2", "B1" ],
  'C2',    [ "A", "B" ],
  'Cs',    [ "Ap", "App" ], # p stands for prime ', there is a problem on linux machines with '
  'Ci',    [ "Ag", "Au" ],
  'C1',    [ "A" ],
);

# Group multiplication table as in swmol3 output
# but minus 1 to have the first IR 0 as used in congen and scatci
# For smaller groups than D2h only part of this table (e.g. 4x4 for C2v or 2x2 for C2) is used
our @group_table = (           # our @group_table = (
  [ 0, 1, 2, 3, 4, 5, 6, 7 ],  #   [ 1, 2, 3, 4, 5, 6, 7, 8 ],
  [ 1, 0, 3, 2, 5, 4, 7, 6 ],  #   [ 2, 1, 4, 3, 6, 5, 8, 7 ],
  [ 2, 3, 0, 1, 6, 7, 4, 5 ],  #   [ 3, 4, 1, 2, 7, 8, 5, 6 ],
  [ 3, 2, 1, 0, 7, 6, 5, 4 ],  #   [ 4, 3, 2, 1, 8, 7, 6, 5 ],
  [ 4, 5, 6, 7, 0, 1, 2, 3 ],  #   [ 5, 6, 7, 8, 1, 2, 3, 4 ],
  [ 5, 4, 7, 6, 1, 0, 3, 2 ],  #   [ 6, 5, 8, 7, 2, 1, 4, 3 ],
  [ 6, 7, 4, 5, 2, 3, 0, 1 ],  #   [ 7, 8, 5, 6, 3, 4, 1, 2 ],
  [ 7, 6, 5, 4, 3, 2, 1, 0 ]   #   [ 8, 7, 6, 5, 4, 3, 2, 1 ]
);                             # );

our %spin_multiplicity = (
  'singlet', 1,
  'doublet', 2,
  'triplet', 3,
  'quartet', 4,
  'quintet', 5,
  'sextet', 6,
  'septet', 7,
  'octet', 8,
  'nonet', 9,
);

# Relative atomic masses - used to determine the center of mass
# from http://www.nist.gov/pml/data/comp.cfm
# only the most abundant isotopes
our %mass = (
  'H' ,   1.007825,
  'D' ,   2.014102,
  'T' ,   3.016049,
  'He',   4.002603, # He4 100.0% (almost)
  'Li',   7.016005, # Li7  92.4%
  'Be',   9.012182,
  'B' ,  11.009305, # B11  80.1%
#  'C' ,  12.000000, # C12  98.9%
  'C' ,  12.0110,   # Molpro
  'N' ,  14.003074, # N14  99.6%
  'O' ,  15.9994,   # Molpro
#  'O' ,  15.994915, # O16  99.8%
  'F' ,  18.998403,
  'Ne',  19.992440, # Ne20 90.5%
  'Na',  22.989769,
  'Mg',  23.985042, # Mg24 79.0%
  'Al',  26.981539,
  'Si',  27.976927, # Si28 92.2%
  'P' ,  30.973762,
  'S' ,  31.972071, # S32  95.0%
  'Cl',  34.968853  # Cl35 75.8%
);

# --------------- subroutines for geometries ----------------

sub generate_geometries {
  my ($r_par, $r_geometry) = @_;

  my @atoms = @{$r_par->{'model'}->{'atoms'}};
  my $na = scalar @atoms;                # number of atoms
  my $dir = $r_par->{'dirs'}->{'model'}; # model directory

  my $r_geom;           # auxiliary reference to one geometry
  my $dir_geom;         # geometry directory
  my @geometries = ();  # array of all geometries - returned at the end
  my $ng = 0;           # number of geometries

  open(GEOM, ">$dir${bs}geometries"); # File with geometries

  $r_par->{'data'}->{'geom_labels'} = "  $r_geometry->{'geometry_labels'}";
  print GEOM '#geom', $r_par->{'data'}->{'geom_labels'}, "\n";
  foreach my $r_given_geom ( @{$r_geometry->{'geometries'}} ) {

    # directory
    $ng++;
    $dir_geom = "$dir${bs}geom$ng";
    &make_dir($dir_geom);

    # general settings for geometry
    $r_geom = {
      'geometry', "  $r_given_geom->{'description'}", # string to use in output files
      'dir', $dir_geom,                               # directory                              \
      'symmetry', $r_par->{'model'}->{'symmetry'}, # symmetry (this is not used right now)  -- of a given geometry
      'natype', $r_given_geom->{'natype'},            # number of atoms needed for QChem codes /
      'atoms', $r_given_geom->{'atoms'}
    };
    if ($r_geometry->{'correct_cm'} == 1) { &make_cm_correction($r_geom); }
    &coord_to_strings($r_geom);
    push @geometries, $r_geom; 

    # Write information about geometry into the string (used also later in files for potential curves etc.)
    # and in the file 'geometries'
    print GEOM sprintf("%5d", $ng), $r_geom->{'geometry'}, "\n";
  }
  &print_info("Number of geometries to run: $ng\n", $r_par);

  close(GEOM);
  return @geometries;
}

sub deg2rad {
  my ($angle) = @_;
  return $angle * 3.14159265358979323846 / 180;
}

sub make_cm_correction {
  my ($r_geom) = @_;
  my $n = scalar @{$r_geom->{'atoms'}};

  # determine cm
  my $mT = 0.0;             # total mass
  my @xT = (0.0, 0.0, 0.0); # coordinates of CM
  for (my $i = 0; $i < $n; $i++) { # read all atoms
     my $ma = $mass{$r_geom->{'atoms'}->[$i][0]};
     $mT += $ma;
     for (my $ir = 1; $ir <= 3; $ir++) { $xT[$ir - 1] += $ma * $r_geom->{'atoms'}->[$i][$ir]; }
  }
  foreach (@xT) { $_ = $_/$mT; }

  # shift all atoms
  for (my $i = 0; $i < $n; $i++) {
     for (my $ir = 1; $ir <= 3; $ir++) { $r_geom->{'atoms'}->[$i][$ir] -= $xT[$ir - 1]; }
  }
  return 1;
}

sub coord_to_strings {
  my ($r_geom) = @_;
  my $n = scalar @{$r_geom->{'atoms'}};
  for (my $i = 0; $i < $n; $i++) {
     foreach (@{$r_geom->{'atoms'}->[$i]}[1..3]) { $_ = sprintf("%13.9f", $_); }
  }
}

# ============ end of subroutines for GEOMETRIES ============

# ================= subroutines for INPUTS ==================

# ------------------------- swmol3 --------------------------

sub make_swmol3_input {
  my ($r_par, $r_str) = @_;
  my $r_geom = $r_par->{'data'}->{'geom'};
  my $natype = $r_geom->{'natype'};
  my $natoms = $natype;
  if ($r_par->{'data'}->{'task'} eq "scattering") { $natype++; }

  # Substitutions to input
  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "BASIS",    $r_par->{'model'}->{'basis'});
  &replace_in_template($r_str, "R_UNIT",   $r_par->{'model'}->{'r_unit'});
  &replace_in_template($r_str, "NATYPE",   $natype);
  &replace_in_template($r_str, "NSYMOP",   $symmetries{$r_geom->{'symmetry'}}->[0]);
  &replace_in_template($r_str, "SYMOP",    $symmetries{$r_geom->{'symmetry'}}->[1]);

  # Adding basis sets for each atom
  my $str_atoms = "";
  for (my $i = 0; $i < $natoms; $i++) {
    $str_atoms .= &get_basis_for_swmol3_input($r_par->{'dirs'}->{'basis'}, $r_par->{'model'}->{'basis'}, @{$r_geom->{'atoms'}->[$i]});
  }
  if ($r_par->{'data'}->{'task'} eq "scattering") {
    $str_atoms .= &get_basis_for_swmol3_input($r_par->{'dirs'}->{'basis'}, "r$r_par->{'model'}->{'radius'}.L$r_par->{'model'}->{'L_basis'}_$r_par->{'model'}->{'cont_suffix'}", "continuum",  0.0, 0.0, 0.0);
  }
  &replace_in_template($r_str, "ATOMS", $str_atoms);
  return 1;
}

sub get_basis_for_swmol3_input {
  my ($dir, $basis, $atom, $xxx, $yyy, $zzz) = @_;
  foreach ($xxx, $yyy, $zzz) { $_ = sprintf("%13.9f", $_); }
  my $str = "";
  &read_file("$dir${bs}swmol3.$atom.$basis", \$str);
  &replace_in_template(\$str, "XXX", $xxx);
  &replace_in_template(\$str, "YYY", $yyy);
  &replace_in_template(\$str, "ZZZ", $zzz);
  return $str;
}

# ------------------------ gaustail -------------------------

sub make_gaustail_input {
  my ($r_par, $r_str) = @_;

  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "RMAT", "$r_par->{'model'}->{'radius'}.0");
  return 1;
}

# -------------------------- sword --------------------------

sub make_sword_input {
  my ($r_par, $r_str) = @_;

  if ($r_par->{'data'}->{'task'} eq "target" || $r_par->{'run'}->{'bound'} == 1) {
    &replace_in_template($r_str, "ZTAIL", ".false.");
  }
  else {
    &replace_in_template($r_str, "ZTAIL", ".true.");
  }
  return 1;
}

# -------------------------- swscf --------------------------

sub make_swscf_input {
  my ($r_par, $r_str) = @_;

  # Substitutions to input
  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "IOCC",  join(",", @{$r_par->{'data'}->{'orbitals'}->{'occupied'}}));
  &replace_in_template($r_str, "IFOCK", join(",", @{$r_par->{'data'}->{'orbitals'}->{'ifock'}}));

  return 1;
}

# ------------------------- molpro --------------------------

sub make_molpro_input {
  my ($r_par, $r_str) = @_;
  my $r_geom = $r_par->{'data'}->{'geom'};
  my $natoms = scalar @{$r_geom->{'atoms'}};
  my $sym_op = "SYMMETRY,".$symmetries{$r_geom->{'symmetry'}}->[1]; # Symmetry operations
  if ($r_geom->{'symmetry'} eq "C1") { $sym_op = "NOSYM"; }
  $sym_op =~ s/\'//g; # In MOLPRO input there are no ''

  # Substitutions to input
  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "BASISNAME",$r_par->{'model'}->{'basis'});
  if ($r_par->{'model'}->{'r_unit'} == 1) { &replace_in_template($r_str, "ANGSTROM", "ANGSTROM"); }
  else                                       { &replace_in_template($r_str, "ANGSTROM", ""); }
  &replace_in_template($r_str, "SYMOP", $sym_op);

  # Adding geometries and basis sets for each atom
  my $str_basis = "";
  my $str_geometry = "";
  my %used_basis = ();
  for (my $i = 0; $i < $natoms; $i++) {
    my $str = "";
    my ($atom, $xxx, $yyy, $zzz) = @{$r_geom->{'atoms'}->[$i]};
    if (!$used_basis{$atom}) {
      &read_file("$r_par->{'dirs'}->{'basis'}${bs}molpro.$atom.$r_par->{'model'}->{'basis'}", \$str);
      $str_basis .= $str;
      $used_basis{$atom} = 1;
    }
    foreach ($xxx, $yyy, $zzz) { $_ = sprintf("%13.9f", $_); }
    $str_geometry .= "$atom,, $xxx, $yyy, $zzz\n";
  }
  $str_basis = "{\n$str_basis\n}";

  if ($r_par->{'run'}->{'molpro_basis'} == 1) { $str_basis = $r_par->{'model'}->{'basis'}; }

  &replace_in_template($r_str, "BASIS", $str_basis);
  &replace_in_template($r_str, "GEOMETRY", $str_geometry);
######################## NEW PART OF MOLPRO #################################
#=begin comment
  &replace_in_template($r_str, "ADDITIONAL", $r_par->{'model'}->{'molpro_setting'});
  if ($r_par->{'model'}->{'model'} =~ /(SE|SEP|pol)/) { #|| $r_par->{'data'}->{'scf_ok'} == 0
    &replace_in_template($r_str, "METHOD", "HF");
    &replace_in_template($r_str, "HFORBPRINT", "200");
    &replace_in_template($r_str, "CASSCF", "");
    &replace_in_template($r_str, "CASWF", "");
    &replace_in_template($r_str, "ALLWF", "");
    &replace_in_template($r_str, "PRINT", "");
  }
  else {
    # Hartree-Fock wave function (must be set if target has charge)
    my $nelec = $r_par->{'model'}->{'nelectrons'};
#    &replace_in_template($r_str, "HFWF", "wf,nelec=$nelec;");

    # Adding closed and occupied orbitals and wf
#    &replace_in_template($r_str, "CLOSED", join(",", @{$r_par->{'data'}->{'orbitals'}->{'frozen'}}));
#    &replace_in_template($r_str, "OCC", join(",", @{$r_par->{'data'}->{'orbitals'}->{'target'}}));
    # Adding all 'wf' to MOLPRO input for which we want to optimise orbitals
    # and also for all target states needed
    my $casscf_states = "";
    my $all_states = "";
    my $occupied = "";
    my $m_closed = "";
    my $m_occ = "";
    for (my $i = 1; $i <= $r_par->{'data'}->{'nir'}; $i++) {
      $m_closed .= "m_closed($i) ";
      $m_occ  .= "m_occ($i) ";
      $occupied[$i - 1] = $r_par->{'model'}->{'frozen_orbs'}->[$i - 1] + $r_par->{'model'}->{'active_orbs'}->[$i - 1];
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
        my $spin = $spin_multiplicity{$statespin} - 1;
        my $nstates = $r_par->{'model'}->{'ncasscf_states'}->{$statespin}->[$i - 1];
        if ($nstates > 0) {
          $casscf_states .= "wf,nelec=$nelec,sym=$i,spin=$spin; state,$nstates;\n";
        }
        my $nstates = $r_par->{'model'}->{'ntarget_states'}->{$statespin}->[$i - 1];
        if ($nstates > 0) {
          $all_states .= "wf,nelec=$nelec,sym=$i,spin=$spin; state,$nstates;\n";
        }
      }
    }
    my $casscf_begin = "";
    my $casscf_end = "";
    my $nactive_electrons = $nelec - 2 * $r_par->{'model'}->{'nfrozen'};
    $casscf_begin = "\{multi; closed, $m_closed; occ, $m_occ;\n";
    $casscf_end = "natorb,print=200; dm; tran,qm;\}";
    $casscf_states = $casscf_begin.$casscf_states."start,2100.2; ".$casscf_end;
    $all_states = $casscf_begin.$all_states."dont,orbital; ".$casscf_end;
    my %input_CASSCF = (
      'M_NAME',   "m_name = CASSCF($nactive_electrons,$r_par->{'model'}->{'nactive'})",
      'M_CLOSED', "m_closed = [".join(",", @{$r_par->{'model'}->{'frozen_orbs'}})."]",
      'M_OCC',    "m_occ = [".join(",", @occupied)."]",
      'PRINT',    "{matrop; load,dm,den,2140.2; write,dm;}",
    );
      &replace_in_template($r_str, "METHOD", "CASSCF");
      &replace_in_template($r_str, "HFORBPRINT", "-1");
      &replace_in_template($r_str, "CASSCF", join("\n", $input_CASSCF{'M_NAME'}, $input_CASSCF{'M_CLOSED'}, $input_CASSCF{'M_OCC'}));
      &replace_in_template($r_str, "CASWF", $casscf_states);
      &replace_in_template($r_str, "ALLWF", $all_states);
      &replace_in_template($r_str, "PRINT", $input_CASSCF{'PRINT'});
  }
#    print "Key: $_ and Value: $input_CASSCF{$_}\n" foreach (keys%input_CASSCF);
#=end comment
#=cut
  return 1;
}

# ------------------------- mpoutrd -------------------------

sub make_mpoutrd_input {
  my ($r_par, $r_str) = @_;

  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "MOLPROOUTPUT", "$r_par->{'dirs'}->{'outputs'}${bs}target.molpro.out");
  &replace_in_template($r_str, "SWEDMOSINPUT", "$r_par->{'dirs'}->{'inputs'}${bs}target.swedmos.inp");
  &replace_in_template($r_str, "NOBT", join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_all'}}));
  &replace_in_template($r_str, "ZTAIL", ".false.");
  return 1;
}

# ------------------------- integrals -------------------------

sub make_integrals_input {
  my ($r_par, $r_str) = @_;

  &replace_in_template($r_str, "RMATR",    $r_par->{'model'}->{'radius'}.".0");
  &replace_in_template($r_str, "NSYMOP",   $symmetries{$r_geom->{'symmetry'}}->[0]);
  &replace_in_template($r_str, "SYMOP",    $symmetries{$r_geom->{'symmetry'}}->[1]);
  &replace_in_template($r_str, "MOLECULE", lc($r_par->{'model'}->{'molecule'}));
#  if ($r_par->{'data'}->{'task'} eq "target") {
#      &replace_in_template($r_str, "NOB", join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_all'}}));
#  }
#  else {
    &replace_in_template($r_str, "NOB", join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}}));
    &replace_in_template($r_str, "DELTHRES", $r_par->{'model'}->{'delthres'});
    #my $nir = $r_par->{'data'}->{'nir'};
    #&replace_in_template($r_str, "DELTHRES", join(",", @{$r_par->{'model'}->{'delthres'}}[0..$nir-1]));
    &replace_in_template($r_str, "MAXL", $r_par->{'model'}->{'L_basis'});
#  }
	if ($r_par->{'run'}->{'scattering'}) {
		&replace_in_template($r_str, "!", "");
    		$str_continuum .= &get_basis_for_integrals_input($r_par->{'dirs'}->{'basis'}, "r$r_par->{'model'}->{'radius'}_$r_par->{'model'}->{'cont_suffix'}", "continuum");
    		&replace_in_template($r_str, "CONTINUUM", $str_continuum);
   		#&replace_in_template($r_str, "EXPONENTS", &get_exponents_from_swmol3_continuum_basis($r_par));
	}
	else {
		&replace_in_template($r_str, "!", "!");
    		&replace_in_template($r_str, "CONTINUUM", "");
   		#&replace_in_template($r_str, "EXPONENTS", "");
	}
  return 1;
}

sub get_basis_for_integrals_input {
  my ($dir, $basis, $continuum) = @_;
  my $str = "";
  &read_file("$dir${bs}integrals.$continuum.$basis", \$str);
  return $str;
}

sub get_exponents_from_swmol3_continuum_basis {
  my ($r_par) = @_;
  my $str_exp = "";
  my $str_basis = &get_basis_for_swmol3_input($r_par->{'dirs'}->{'basis'}, "r$r_par->{'model'}->{'radius'}.L$r_par->{'model'}->{'L_basis'}", "continuum",  0.0, 0.0, 0.0);
  my @nexp = ();
  if ($str_basis =~ m/jco = ((\d+,)+)\s/s) {
    @nexp = split(',', $1);
  }
  else { die "Error: failed to read 'jco' from continuum basis file !\n"; }
  my @all_exp = split(/\s+1\s+1\s+/s, $str_basis);  # the first element is &atom namelist
  foreach (@all_exp) { $_ =~ s/\s+1\.\s*$// ; }     # but other elements are only exponents
  my $ne = 0;
  for (my $n = 0; $n < scalar @nexp; $n++) {
    $str_exp .= "exponents(:,$n) = ";
    for (my $i = 1; $i <= $nexp[$n]; $i++) {
      $ne++;
      $str_exp .= "$all_exp[$ne], ";
    }
    if ($n < (scalar @nexp) - 1) { $str_exp .= "\n  "; }
  }
  return $str_exp;
}

# ------------------------- swedmos -------------------------

sub make_swedmos_input {
  my ($r_par, $r_str) = @_;

  my $which_orbitals = 'target_used';
  if ($r_par->{'data'}->{'task'} eq "target") { $which_orbitals = 'target_all'; }

  my $ind = "";
  for (my $i = 1; $i <= $r_par->{'data'}->{'nir'}; $i++) {
    $ind .= "1,$i,1,$r_par->{'data'}->{'orbitals'}->{$which_orbitals}->[$i - 1], ";
  }

  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "IND", $ind);
  &replace_in_template($r_str, "DELTHRES", $r_par->{'model'}->{'delthres'});
#  my $nir = $r_par->{'data'}->{'nir'};
#  &replace_in_template($r_str, "DELTHRES", join(",", @{$r_par->{'model'}->{'delthres'}}[0..$nir]));

  return 1;
}

# ------------------------ swtrmo ---------------------------

sub make_swtrmo_input {
  my ($r_par, $r_str) = @_;

  # For scattering calculation, we specify all used orbitals, target + continuum - deleted
  # For target, we specify only target molecular orbitals, nothing was deleted
  my $which_orbitals = 'all_used';
  if ($r_par->{'data'}->{'task'} eq "target") { $which_orbitals = 'target_all'; }

  my $iget = "";
  for (my $i = 1; $i <= $r_par->{'data'}->{'nir'}; $i++) {
    $iget .= "1,$i,1,$r_par->{'data'}->{'orbitals'}->{$which_orbitals}->[$i - 1], ";
  }

  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "IGET", $iget);
  return 1;
}

# ------------------------ congen ---------------------------

sub make_congen_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym    = $r_par->{'data'}->{$task}->{'symmetry'};                      # number 0,1,2,...
  my $stateir     = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym]; # string A1,B1,...
  my $statespin   = $r_par->{'data'}->{$task}->{'spin'};                          # string singlet,doublet,...
  my $statespinno = $spin_multiplicity{lc($statespin)};                           # number 1,2,3,...
  my $nelectrons  = $r_par->{'model'}->{'nelectrons'};                         # number of electrons, first target
  my $nactive_electrons = $nelectrons - 2 * $r_par->{'model'}->{'nfrozen'};    # number of target active electrons

  # &state namelist
  #================

  my %state_namelist = (
    'MOLECULE', $r_par->{'model'}->{'molecule'},           # For scattering: "e + " will be added
    'SPIN',     $statespin,
    'SYMMETRY', $stateir,
    'MEGUL',    "70",                                         # For target: e.g. 713, 1 for singlet, 3 for symmetry
    'ISCAT',    1,                                            # For scattering: changed to 2
    'SYMTYP',   2,
    'QNTOT',    "$statespinno,$statesym,0",
    'NELECT',   $nelectrons,                                  # For scattering: +1 will be added
    'NOB',      "",
    'NOB0',     "",
    'NREFO',    0,
    'REFORB',   "",
    'LNDO',     $r_par->{'model'}->{'lndo'}
  );

  &print_info("\nSearching for reference orbitals for $statespin $stateir state.\n", $r_par);
  my $nrefo  = 0;
  my @reforb = ();

  # Special settings for target
  if ($task eq "target") {

    $state_namelist{'MEGUL'} = "7".$statespinno.$statesym;
    $state_namelist{'NOB'} = ($r_par->{'run'}->{'ukrmolplus'} == 0 ? join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_all'}}) : join(",", @{$r_par->{'data'}->{'orbitals'}->{'all_used'}})),
    $state_namelist{'NOB0'} = ($r_par->{'run'}->{'ukrmolplus'} == 0 ? $state_namelist{'NOB'} : join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}})),

    # reference orbitals depend on target state spin and symmetry
    # only frozen + active orbitals are used to find reference orbital
    &print_info("Putting $nelectrons electrons into (".join(',', @{$r_par->{'data'}->{'orbitals'}->{'target'}}).") orbitals\n", $r_par);
    ($nrefo, @reforb) = &get_reference_orbitals($r_par->{'data'}->{'orbitals'}->{'target'}, $nelectrons, $statesym, $statespinno, $r_par);
    if ($r_par->{'model'}->{'model'} =~ /(pol)/ && $statesym == 0){ 
	$r_par->{'data'}->{'target'}->{'pol_reforb'} = join(",", @reforb);
    }
    if ($nrefo == 0) {
      die   "Increase your active space or change the target settings:\n
             put 0 for target state: $statespin $stateir($statesym) !\n";
    }
  }

  # Special settings for scattering calculation
  else {

    $nelectrons++;
    $state_namelist{'MOLECULE'} = "e + $state_namelist{'MOLECULE'}";
    $state_namelist{'ISCAT'} = 2;
    $state_namelist{'NOB'}  = join(",", @{$r_par->{'data'}->{'orbitals'}->{'all_used'}});
    $state_namelist{'NOB0'} = join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}});

    # reference orbitals depend on target state spin and symmetry
    # now all available orbitals are used
    &print_info("Putting $nelectrons electrons into (".join(',', @{$r_par->{'data'}->{'orbitals'}->{'target_used'}}).") orbitals\n", $r_par);
    my @reference = @{$r_par->{'data'}->{'orbitals'}->{'target'}};
    $reference[$statesym] = $r_par->{'data'}->{'orbitals'}->{'target'}->[$statesym] + 1;
    ($nrefo, @reforb) = &get_reference_orbitals(\@reference, $nelectrons, $statesym, $statespinno, $r_par); #$r_par->{'data'}->{'orbitals'}->{'all_used'}
    if ($nrefo == 0) {
      die   "This should never happen for N+1 problem. Let's go doing something else :->\n";
    }
  }

  $state_namelist{'NELECT'} = $nelectrons;  
  $state_namelist{'NREFO'} = $nrefo;
  $state_namelist{'REFORB'} = join(",", @reforb);
  $state_namelist{'REFORB'} =~ s/(\d+,\d+,\d+,\d+,\d+,)/\1 /g; # Add spaces for better reading
  &print_info("Reference orbitals found:\n$state_namelist{'REFORB'}\n\n", $r_par);

  &replace_all_in_template($r_str, \%state_namelist);

  # &wfngrp namelist(s)
  #====================

  my %wfngrp_namelist = (
    'GNAME',    "",
    'QNTAR',    "-1,0,0", # Default for target, no constraints
    'NELECG',   $r_par->{'model'}->{'nelectrons'},         # For scattering: +1
    'NREFOG',   $state_namelist{'NREFO'},
    'REFORG',   $state_namelist{'REFORB'},
    'NDPROD',   0,
    'NELECP',   "",
    'NSHLP',    "",
    'PQN',      "",
    'MSHL',     "",
  );

  # Read &wfngrp namelist template  
  my $wfngrp_template = "";
  if(!&read_file("$r_par->{'dirs'}->{'templates'}${bs}congen.wfngrp.inp", \$wfngrp_template)) {
    &print_info("Warning: no template file for congen.wfngrp.inp !\n", $r_par);
  }
  my $wfngrp = $wfngrp_template;

  # For target we add one &wfngrp namelist
  # it is assumed that information about frozen orbitals is 
  #   in the array $r_par->{'data'}->{'orbitals'}->{'frozen'}
  # and about active (= valence for SE and SEP) orbitals are
  #   in the array $r_par->{'data'}->{'orbitals'}->{'active'}
  if ($task eq "target") {
    if ($r_par->{'model'}->{'model'} =~ /(pol)/) {
	    if ($statesym ==0){
		    $wfngrp_namelist{'GNAME'} = "$state_namelist{'MOLECULE'} - $state_namelist{'SPIN'} $state_namelist{'SYMMETRY'}";
		    &add_frozen_shells($r_par, \%wfngrp_namelist);
		    &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
		    &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
		    $$r_str .= $wfngrp;
 	    }
	    $wfngrp_namelist{'GNAME'}  = "Ground state^-1 x virtual^1";
	    &add_frozen_shells($r_par, \%wfngrp_namelist);
	    &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons - 1);
	    &add_virtual_shells($r_par, \%wfngrp_namelist, 1);
	    $wfngrp = $wfngrp_template;
	    &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
	    $$r_str .= $wfngrp;
    }
    else {
	    $wfngrp_namelist{'GNAME'} = "$state_namelist{'MOLECULE'} - $state_namelist{'SPIN'} $state_namelist{'SYMMETRY'}";
	    &add_frozen_shells($r_par, \%wfngrp_namelist);
	    &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
	    &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
	    $$r_str .= $wfngrp;
    }	  
  }
  # For scattering calculation the number of &wfngrp namelists
  # depends on model
  else {
    $wfngrp_namelist{'NELECG'}++;

    # SE and SEP models
    if ($r_par->{'model'}->{'model'} =~ /(SE|SEP)/) {

      # Check whether settings are consistent
      if ($r_par->{'model'}->{'ntarget_states_used'} > 1) {
        &print_info("Warning: ntarget_states_used is $r_par->{'model'}->{'ntarget_states_used'}\n", $r_par);
        &print_info("         which is inconsistent with $r_par->{'model'}->{'model'} model !\n", $r_par);
        &print_info("         Only one (the lowest) state will be used !\n", $r_par);
      }

      # Determin spin and symmetry of the ground state
      my ($ground_state_spin, $ground_state_sym) = split(/\./, $r_par->{'data'}->{'target'}->{'ground_state'});

      # In both SE and SEP we use 
      # ground state x continuum^1
      $wfngrp_namelist{'GNAME'} = "Ground state x continuum^1";
      $wfngrp_namelist{'QNTAR'} = "$ground_state_spin,$ground_state_sym,0";
      &add_frozen_shells($r_par, \%wfngrp_namelist);
      &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
      &add_continuum_shells($r_par, \%wfngrp_namelist, $ground_state_sym);
      &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
      $$r_str .= $wfngrp;
      
      # ground state x virtual^1 (if any)
      if ($r_par->{'data'}->{'orbitals'}->{'virtual'}->[$statesym] > 0) {
        $wfngrp_namelist{'GNAME'}  = "Ground state x virtual^1";
        $wfngrp_namelist{'QNTAR'} = "-1,0,0";
        &add_frozen_shells($r_par, \%wfngrp_namelist);
        &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
        &add_virtual_shells($r_par, \%wfngrp_namelist, 1);
        $wfngrp = $wfngrp_template;
        &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
        $$r_str .= $wfngrp;
      }

      # and for SEP only
      if ($r_par->{'model'}->{'model'} eq "SEP") {
        # core valence^-1 x virtual^2
        if ($r_par->{'data'}->{'orbitals'}->{'virtual'}->[$statesym] > 0) {
          $wfngrp_namelist{'GNAME'}  = "Core valence^-1 x virtual^2";
          &add_frozen_shells($r_par, \%wfngrp_namelist);
          &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons - 1);
          &add_virtual_shells($r_par, \%wfngrp_namelist, 2);
          $wfngrp = $wfngrp_template;
          &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
          $$r_str .= $wfngrp;
        }
      }

    }
    # CAS models
    else {
      # First loop over all target states spin symmetries 
      # and get Target state x continuum^1 wfngrp
      for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
        foreach my $targetstatespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
          if ($r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$spin_multiplicity{lc($targetstatespin)}.".$i"}) { # skip if target states of a given spin-symmetry are not used
            my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i];
            $wfngrp_namelist{'GNAME'} = "$targetstatespin $str_ir x continuum^1";
            $wfngrp_namelist{'QNTAR'} = "$spin_multiplicity{lc($targetstatespin)},$i,0";
            &add_frozen_shells($r_par, \%wfngrp_namelist);
            &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
            &add_continuum_shells($r_par, \%wfngrp_namelist, $i);
            $wfngrp = $wfngrp_template;
            &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
          }
        } # foreach spin state (multiplicity)
      } # for each IRs
      
      # then add (core+cas)^N+1
      $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc+1";
      $wfngrp_namelist{'QNTAR'} = "-1,0,0";
      &add_frozen_shells($r_par, \%wfngrp_namelist);
      &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons + 1);
      $wfngrp = $wfngrp_template;
      &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
      $$r_str .= $wfngrp;

      # and finally (core+cas)^N-h x virtual^1+h if there are any virtual orbitals
      if ($r_par->{'model'}->{'nvirtual'} > 0) {
        if ($r_par->{'model'}->{'model'} =~ /CAS-A/) {
          for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
            foreach my $targetstatespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
              if ($r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$spin_multiplicity{lc($targetstatespin)}.".$i"} &&                          # skip if target states of a given spin-symmetry are not used
                  $r_par->{'data'}->{'orbitals'}->{'virtual'}->[get_orbitals_symmetry($i, $r_par->{'data'}->{'scattering'}->{'symmetry'})] > 0) { #   or if there are no virtual orbitals of required symmetry
                my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i];
                $wfngrp_namelist{'GNAME'} = "$targetstatespin $str_ir x virtual^1";
                $wfngrp_namelist{'QNTAR'} = "$spin_multiplicity{lc($targetstatespin)},$i,0";
                &add_frozen_shells($r_par, \%wfngrp_namelist);
                &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
                &add_virtual_shell_for_given_target_symmetry($r_par, \%wfngrp_namelist, $i);
                $wfngrp = $wfngrp_template;
                &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                $$r_str .= $wfngrp;
              }
            } # foreach spin state (multiplicity)
          } # for each IRs
        }
        else { # CAS or CAS-B or CAS-C
          $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc virtual^1";
          &add_frozen_shells($r_par, \%wfngrp_namelist);
          &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
          &add_virtual_shells($r_par, \%wfngrp_namelist, 1);
          $wfngrp = $wfngrp_template;
          &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
          $$r_str .= $wfngrp;
        }

        # for CAS-C add (core+cas)^N-1 x virtual^2 yet
        if ($r_par->{'model'}->{'model'} =~ /CAS-C/) {
          $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc-1 virtual^2";
          &add_frozen_shells($r_par, \%wfngrp_namelist);
          &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons - 1);
          &add_virtual_shells($r_par, \%wfngrp_namelist, 2);
          $wfngrp = $wfngrp_template;
          &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
          $$r_str .= $wfngrp;
        }
      }

    }
    
  } # end of scattering calculation setting
  
  return 1;
}

sub add_frozen_shells {
  my ($r_par, $r_wfngrp_namelist) = @_;
  my ($nelecp, $nshlp) = (0, 0);

  # clean from previous &wfngrp
  $r_wfngrp_namelist->{'NDPROD'} = 0;
  $r_wfngrp_namelist->{'NELECP'} = "";
  $r_wfngrp_namelist->{'NSHLP'}  = "";
  $r_wfngrp_namelist->{'PQN'}    = "";
  $r_wfngrp_namelist->{'MSHL'}   = "";
  
  if ($r_par->{'model'}->{'nfrozen'} > 0) {
    $r_wfngrp_namelist->{'NDPROD'}++;
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      if ($r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i] > 0) {
        $nelecp += 2 * $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i];
        $nshlp++;
        $r_wfngrp_namelist->{'PQN'}  .= "0,1,$r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i], ";
        $r_wfngrp_namelist->{'MSHL'} .= "  $i,   ";
      }
    }
    $r_wfngrp_namelist->{'NELECP'} .= "$nelecp,";
    $r_wfngrp_namelist->{'NSHLP'}  .= "$nshlp,";
  }
  
  return 1;
}

sub add_active_shells {
  my ($r_par, $r_wfngrp_namelist, $nel) = @_;

  if ($r_par->{'model'}->{'nactive'} > 0) {
    $r_wfngrp_namelist->{'NDPROD'}++;
    my $nshlp = 0;
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      if ($r_par->{'data'}->{'orbitals'}->{'active'}->[$i] > 0) {
        $nshlp++;
        my $nfrozen = $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i];
        $r_wfngrp_namelist->{'PQN'}  .= "0,".($nfrozen + 1).",".($nfrozen + $r_par->{'data'}->{'orbitals'}->{'active'}->[$i]).", ";
        $r_wfngrp_namelist->{'MSHL'} .= "  $i,   ";
      }
    }
    $r_wfngrp_namelist->{'NELECP'} .= "$nel,";
    $r_wfngrp_namelist->{'NSHLP'}  .= "$nshlp,";
  }
    
  return 1;
}

sub add_continuum_shells {
  my ($r_par, $r_wfngrp_namelist, $target_state_sym) = @_;

  $r_wfngrp_namelist->{'NDPROD'}++;
  my $cont_sym = get_orbitals_symmetry($target_state_sym, $r_par->{'data'}->{'scattering'}->{'symmetry'});
  my $ntarget = $r_par->{'data'}->{'orbitals'}->{'target_used'}->[$cont_sym];
  $r_wfngrp_namelist->{'PQN'}  .= "0,".($ntarget + 1).",".($ntarget + 2).", ";
  $r_wfngrp_namelist->{'MSHL'} .= "  $cont_sym,   ";
  $r_wfngrp_namelist->{'NELECP'} .= "1,";
  $r_wfngrp_namelist->{'NSHLP'}  .= "1,";
  
  return 1;
}

sub add_virtual_shells {
  my ($r_par, $r_wfngrp_namelist, $nel) = @_;

  $r_wfngrp_namelist->{'NDPROD'}++;
  my $nshlp = 0;
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    if ($r_par->{'data'}->{'orbitals'}->{'virtual'}->[$i] > 0) {
      $nshlp++;
      my $ntarget = $r_par->{'data'}->{'orbitals'}->{'target'}->[$i];
      $r_wfngrp_namelist->{'PQN'}  .= "0,".($ntarget + 1).",".($ntarget + $r_par->{'data'}->{'orbitals'}->{'virtual'}->[$i]).", ";
      $r_wfngrp_namelist->{'MSHL'} .= "  $i,   ";
    }
  }
  $r_wfngrp_namelist->{'NELECP'} .= "$nel,";
  $r_wfngrp_namelist->{'NSHLP'}  .= "$nshlp,";
  
  return 1;
}

sub add_virtual_shell_for_given_target_symmetry {
  my ($r_par, $r_wfngrp_namelist, $target_state_sym) = @_;

  $r_wfngrp_namelist->{'NDPROD'}++;
  my $virt_sym = get_orbitals_symmetry($target_state_sym, $r_par->{'data'}->{'scattering'}->{'symmetry'});
  my $ntarget = $r_par->{'data'}->{'orbitals'}->{'target'}->[$virt_sym];
  if ($r_par->{'data'}->{'orbitals'}->{'virtual'}->[$virt_sym] > 1) {
    $r_wfngrp_namelist->{'PQN'}  .= "0,".($ntarget + 1).",".($ntarget + 2).", ";
  }
  else {
    $r_wfngrp_namelist->{'PQN'}  .= "0,".($ntarget + 1).",".($ntarget + 1).", ";
  }
  $r_wfngrp_namelist->{'MSHL'} .= "  $virt_sym,   ";
  $r_wfngrp_namelist->{'NELECP'} .= "1,";
  $r_wfngrp_namelist->{'NSHLP'}  .= "1,";
  
  return 1;
}

# Auxiliary subroutine to determine symmetry of orbitals
# from group multiplication table in such a way that
# target state symmetry  x  orbitals symmetry = scattering state symmetry
sub get_orbitals_symmetry {
  my ($target_sym, $scat_sym) = @_;
  my $orb_sym = 0;
  IR: for (my $ir = 0; $ir < scalar @group_table; $ir++) {
    if ($group_table[$target_sym]->[$ir] == $scat_sym) { $orb_sym = $ir; last IR; }
  }
  return $orb_sym;
}

# subroutine which returns reference orbitals of a CSF with a given total spin and symmetry
sub get_reference_orbitals {
  my ($r_norbitals, $nelectrons, $state_sym, $state_spin, $r_par) = @_;
  my @el_pos  = (); # current positions of all electrons
  my @orb_sym = (); # array of symmetries of orbitals
  
  # first fill up all electrons starting from the first orbital of the symmetry 0
  my $nel = $nelectrons;  # number of remaining electrons to put in orbitals
  my $norb = 0;           # number of all orbitals
  for (my $i = 0; $i < scalar @$r_norbitals; $i++) { # loop over all IRs
    for (my $j = 0; $j < $r_norbitals->[$i]; $j++) { # loop over all orbitals in IR $i
      if    ($nel >  1) { $el_pos[$nelectrons - $nel] = $norb + 1; $el_pos[$nelectrons - $nel + 1] = -$norb - 1; }
      elsif ($nel == 1) { $el_pos[$nelectrons - $nel] = $norb + 1; }
      $orb_sym[$norb] = $i;
      $nel -= 2;
      $norb++;
    }
  }

  my $found = 0;
  while ($found == 0) {
    # first check the configuration, if symmmetry and spin are correct then stop
#    &print_info("@el_pos\n", $r_par);
    my ($sym, $spin) = &get_csf_symmetry_and_spin(\@orb_sym, \@el_pos); 
#    &print_info("Total sym and spin is $sym and $spin\n", $r_par);
    if ($sym == $state_sym && $spin == $state_spin) { 
      $found = 1;
    }
    else {
      # find electron to move
      my $move_el = $nelectrons;
      my $n = scalar @orb_sym;
      if ($el_pos[$move_el - 1] == -$n) {                              # if the last electron is in the last orbital with spin down
        while (abs($el_pos[$move_el - 1]) == $n && $move_el > 0) {     # then look for an electron (in reverse order) untill there is one which can be moved 
          if    ($nelectrons % 2 == 0 && $move_el % 2 == 1) { $n--; }  # along the way we have to change the highest possible occupied orbital
          elsif ($nelectrons % 2 == 1 && $move_el % 2 == 0) { $n--; }  # bacause there can be only two electrons in each orbital
          $move_el--;
        }
      }
      if ($move_el == 0) { # if no electron can be moved there is no configuration for a given state spin and symmetry
        $found = -1;
      }
      else {               # otherwise move the electron 
        $move_el--;        # shift because the first element has an index 0
        if ($el_pos[$move_el] > 0) {                            # if the electron has spin up then
          $el_pos[$move_el] = -$el_pos[$move_el];               # change spin
        }
        else {                                                  # if the electron has spin down then
          $el_pos[$move_el] = -$el_pos[$move_el] + 1;           # move it
        }
        for (my $i = $move_el + 1; $i < $nelectrons; $i++ ) {   # the rest of electrons place to the next available orbitals
          if ($el_pos[$i - 1] < 0) { $el_pos[$i] = -$el_pos[$i - 1] + 1; }
          else                     { $el_pos[$i] = -$el_pos[$i - 1]; }
        }
      } # if ($move_el == 0)
    } # if ($sym == $state_sym && $spin == $state_spin)
  } # while ($found == 0)
  if ($found == -1) {
 #this condition is to get reference orbitals for irreducible representations other that first one. It should be done better in the future!
    if ($r_par->{'model'}->{'model'} =~ /pol/ && $state_sym != 0) { 
	    my @pol_ref = @_;
	    @pol_ref = split(',',$r_par->{'data'}->{'target'}->{'pol_reforb'});
	    @pol_ref[2] += 1;
	    @pol_ref[$state_sym*5 + 2] -= 1;
	    @pol_ref[$state_sym*5 + 3] += 1;
	    my $nrefo = scalar @pol_ref/5;
	    return ($nrefo,@pol_ref);
    }
    else {
	    &print_info("\n  Warning: no configuration found for $nelectrons electrons\n",$r_par);
	    &print_info("             for a state with (sym, spin) = ($state_sym, $state_spin) !!!\n",$r_par);
	    &print_info("             Orbitals used: @{$r_norbitals}\n\n",$r_par);
	    return (0, "");
    }
  }
  else {
    return &get_reforb_quintets_for_csf(\@orb_sym, \@el_pos);
  }
}

# subroutine which returns total spin and symmetry from @orb_sym and @el_pos
sub get_csf_symmetry_and_spin {
  my ($r_orb_sym, $r_el_pos) = @_;
  my ($sym, $spin) = (0, 0);
  
  for (my $i = 0; $i < scalar @$r_el_pos; $i++) { # loop over all electrons
    if ($r_el_pos->[$i] > 0) { $spin++; }         # add spin up
    else                     { $spin--; }         # add spin down
    # multiply the current total state symmetry by symmetry of the orbital in which the electron is
    $sym = $group_table[$sym]->[$r_orb_sym->[abs($r_el_pos->[$i]) - 1]];
  }
  if ($spin < 0) { $spin = -$spin; }
  $spin++; # to have 1 for singlet etc. 

  return ($sym, $spin);
}

# subroutine which returns nrefo and quintets in reforb 
# from arrays of orbitals symmetries @orb_sym and positions of electrons @el_pos
sub get_reforb_quintets_for_csf {
  my ($r_orb_sym, $r_el_pos) = @_;
  my $nrefo = 0;
  my @reforb = ();

  my $nelectrons = scalar @$r_el_pos;
  my $nel = 1;                                         # number of electrons for a quintet
  my $start_orb = 1;                                   # starting orbital to be used in a quintet
  my ($prev_orb, $prev_orb_sym);                       # number and symmetry of the previous orbital
  my $curr_orb = $r_el_pos->[0];                       # number of the current orbital
  my $curr_orb_sym = $r_orb_sym->[abs($curr_orb) - 1]; # symmetry of the current orbital
  for (my $i = 1; $i < $nelectrons; $i++) {            # loop over all electrons but the first

    $prev_orb = $curr_orb;
    $prev_orb_sym = $curr_orb_sym;
    $curr_orb = $r_el_pos->[$i];
    $curr_orb_sym = $r_orb_sym->[abs($curr_orb) - 1];
    
    if ($curr_orb_sym == $prev_orb_sym                 # if the same symmetry of the current orbital and previous one
       && ($curr_orb == -$prev_orb                     # and the electron is at the same orbital as the previous one
           || (abs($curr_orb) == abs($prev_orb) + 1    # or it is in the next orbital
               && $nel % 2 == 0) )) {                  #    and at the same time there is even number of electrons in the quintet
      $nel++;                                          # then just add one electron to the quintet
    }
    else {                                             # otherwise we have to save the quintet and start a new one
      $nrefo++;
      push @reforb, $prev_orb_sym, $start_orb, $nel, ($nel % 2 == 1 && $prev_orb < 0) ? 1 : 0, 0;
      if ($prev_orb_sym == $curr_orb_sym) { $start_orb += ($nel + 1) / 2; } # if the same symmetry then get the next available orbital
      else                                { $start_orb = 1; }               # if not then start from the first
      $nel = 1;
    }
  }
  $nrefo++;
  push @reforb, $curr_orb_sym, $start_orb, $nel, ($nel % 2 == 1 && $curr_orb < 0) ? 1 : 0, 0;

  return ($nrefo, @reforb);
}


# ------------------------ scatci ---------------------------

sub make_scatci_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};

  my %input_namelist = (
    'MOLECULE', $r_par->{'model'}->{'molecule'},
    'UKRMOLP',  ($r_par->{'run'}->{'ukrmolplus'} == 1 ? ".true." : ".false."),
    'SPIN',     $statespin,
    'SYMMETRY', $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym],
  );

  # For each target csf we calculate given number of eigenvalues
  # from which later in make_denprop_input select states with lowest energies
  if ($task eq "target") {

    # First update auxiliary variables used later for denprop and scattering.scatci inputs
    push @{$r_par->{'data'}->{'target'}->{'spinsym_order'}}, "$spin_multiplicity{lc($statespin)}.$statesym";
    my $nciset = scalar @{$r_par->{'data'}->{'target'}->{'spinsym_order'}};

    $input_namelist{'MEGUL'} = "7".$spin_multiplicity{lc($statespin)}.$statesym;

    if ($r_par->{'model'}->{'model'} =~ /pol/){
      $input_namelist{'NSTAT'} = $r_par->{'data'}->{'target'}->{'pol_NOCSF'}->[$statesym];
      $input_namelist{'NCISET'} = ($statesym == 0 ? 1 : 0),	
    }
    else { 
      $input_namelist{'NSTAT'} = $r_par->{'model'}->{'ntarget_states'}->{$statespin}->[$statesym];
      $input_namelist{'NCISET'} = $nciset;
    }	   
  }
  # For scattering calculation we need all eigenvalues
  else {

    my @mcont = map($group_table[$statesym]->[$_], @{$r_par->{'data'}->{'target'}->{'mcont'}});
    my $notgt = "";
    for (my $i = 0; $i < scalar @mcont; $i++) {
      $notgt .= "$r_par->{'data'}->{'orbitals'}->{'cont_used'}->[$mcont[$i]],";
    }
    $input_namelist{'MEGUL'}  = "70";
    $input_namelist{'NTGSYM'} = $r_par->{'data'}->{'target'}->{'ntgt'};
    $input_namelist{'MCONT'}  = join(",", @mcont);
    $input_namelist{'NOTGT'}  = $notgt;
    $input_namelist{'NUMTGT'} = join(",", @{$r_par->{'data'}->{'target'}->{'ntgtl'}});
    $input_namelist{'ICIDG'} = ($r_par->{'run'}->{'parallel_diag'} == 1 ? 0 : 1),
  }
  &replace_all_in_template($r_str, \%input_namelist);
  return 1;
}
# ------------------------ hamdiag ---------------------------

# ----------------------- gausprop --------------------------

sub make_gausprop_input {
  my ($r_par, $r_str) = @_;

  # The same as in swtrmo for target
  my $iget = "";
  for (my $i = 1; $i <= $r_par->{'data'}->{'nir'}; $i++) {
    $iget .= "1,$i,1,$r_par->{'data'}->{'orbitals'}->{'target_all'}->[$i - 1], ";
  }
  
  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "IGET", $iget);
  return 1;
}

# ----------------------- denprop ---------------------------

sub make_denprop_input {
  my ($r_par, $r_str) = @_;
  my $r_states = $r_par->{'data'}->{'target'}->{'states'};

  # According to number of target calculations (different spin-symmetry states) set NTGTF as "1,2,3,..."
  my $ntgtf = "";
  for (my $i = 1; $i <= $r_par->{'data'}->{'target'}->{'ntgt'}; $i++) {
    $ntgtf .= "$i,";
  }

  # Get number of target states used in each spin-symmetry from the number of states 
  # as required in model -> ntarget_states_used (only states with lowest energies are used)
  
  # First order states by energy and check whether there is enough states as required if not use only calculated ones
  my @ordered_states = sort { $r_states->{$a} <=> $r_states->{$b} } keys %{$r_states};
  $r_par->{'data'}->{'target'}->{'ordered_states'} = \@ordered_states;
  $r_par->{'data'}->{'target'}->{'ground_state'} = $ordered_states[0];
  my $r_ordered = $r_par->{'data'}->{'target'}->{'ordered_states'};
  my $nstates = scalar @$r_ordered;
  if ($nstates < $r_par->{'model'}->{'ntarget_states_used'}) {
    &print_info("  Warning: There is less target states ($nstates) in scatci output(s) than required ($r_par->{'model'}->{'ntarget_states_used'}) !!!\n", $r_par);
    &print_info("           Increase numbers of states in model->ntarget_states.\n", $r_par);
    &print_info("           Only those states which were found will be used !\n", $r_par);
    $r_par->{'model'}->{'ntarget_states_used'} = $nstates;
  }

  # Now get number of states of given spin and symmetry within the lowest-lying states
  my $tgt_state_used = $r_par->{'model'}->{'ntarget_states_used'};
  if ($r_par->{'model'}->{'model'} =~ /pol/){$tgt_state_used += 1;}
  $r_par->{'data'}->{'target'}->{'used_tgt_states'} = {};
  for (my $i = 0; $i < $tgt_state_used; $i++) {
    my ($spin, $sym, $number) = split(/\./, $r_ordered->[$i]);
    $r_par->{'data'}->{'target'}->{'used_tgt_states'}->{"$spin.$sym"} = $number;
  }
  $r_par->{'data'}->{'target'}->{'ntgt'} = scalar (keys %{$r_par->{'data'}->{'target'}->{'used_tgt_states'}});
  
  # Finally set variables for scatci input 
  # and also auxiliary array 'mcont' used later in scattering.scatci inputs
  my $nftsor = "";
  my $ntgtf = "";
  $r_par->{'data'}->{'target'}->{'ntgtl'} = [];
  my $r_spinsym_order = $r_par->{'data'}->{'target'}->{'spinsym_order'};
  for (my $i = 0; $i < scalar @{$r_spinsym_order}; $i++) {
    if ($r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$r_spinsym_order->[$i]}) {
      my ($spin, $sym) = split(/\./, $r_spinsym_order->[$i]);
      $nftsor .= "7$spin$sym,";
      $ntgtf .= ($i+1).",";
      push @{$r_par->{'data'}->{'target'}->{'ntgtl'}}, $r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$r_spinsym_order->[$i]};
      push @{$r_par->{'data'}->{'target'}->{'mcont'}}, $sym;
    }
  }
  if ($r_par->{'model'}->{'model'} =~ /pol/){
      @{$r_par->{'data'}->{'target'}->{'ntgtl'}} = @{$r_par->{'data'}->{'target'}->{'pol_NOCSF'}};
  }
  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "IPOL",     $r_par->{'model'}->{'max_multipole'});
  &replace_in_template($r_str, "NTGT",     $r_par->{'data'}->{'target'}->{'ntgt'});
  &replace_in_template($r_str, "NFTSOR",   $nftsor);
  &replace_in_template($r_str, "NTGTF",    $ntgtf);
  &replace_in_template($r_str, "NTGTL",    join(",", @{$r_par->{'data'}->{'target'}->{'ntgtl'}}));
  &replace_in_template($r_str, "UKRMOLP",  ($r_par->{'run'}->{'ukrmolplus'} == 1 ? ".true." : ".false."));
  return 1;
}

# ---------------------- outerres ---------------------------

sub make_outerres_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};

  &replace_in_template($r_str, "MOLECULE", "e + ".$r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "SPIN",     $statespin);
  &replace_in_template($r_str, "SYMMETRY", $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym]);
  &replace_in_template($r_str, "MGVN",     $statesym);
  &replace_in_template($r_str, "STOT",     $spin_multiplicity{lc($statespin)});
  &replace_in_template($r_str, "NTARG",    $r_par->{'model'}->{'ntarget_states_used'});
  &replace_in_template($r_str, "IDTARG",   join(",", @{$r_par->{'data'}->{'target'}->{'idtarg'}}));
  &replace_in_template($r_str, "RMATR",    $r_par->{'model'}->{'radius'}.".0");
  &replace_in_template($r_str, "RAF",      $r_par->{'model'}->{'raf'});
  &replace_in_template($r_str, "ISMAX",    $r_par->{'model'}->{'max_multipole'});
  &replace_in_template($r_str, "NERANG",   scalar split(/\s*,\s*/, $r_par->{'model'}->{'nescat'}));
  &replace_in_template($r_str, "NESCAT",   $r_par->{'model'}->{'nescat'});
  &replace_in_template($r_str, "EINC",     $r_par->{'model'}->{'einc'});
  &replace_in_template($r_str, "IEUNIT",   $r_par->{'model'}->{'e_unit'});
  &replace_in_template($r_str, "IXSN",     $r_par->{'model'}->{'x_unit'});
  &replace_in_template($r_str, "MAXI",     $r_par->{'model'}->{'maxi'});
  &replace_in_template($r_str, "MAXF",     $r_par->{'model'}->{'maxf'});
  &replace_in_template($r_str, "UKRMOLP",  ($r_par->{'run'}->{'ukrmolplus'} == 1 ? ".true." : ".false."));
  &replace_in_template($r_str, "LCONT",    $r_par->{'data'}->{'scattering'}->{'cont_csf'});
  return 1;
}

# ---------------------- timedel ---------------------------
sub make_timedel_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};
  my @einc = ();
  my $einit = $r_par->{'model'}->{'einit'};
  my $efinal = $r_par->{'model'}->{'efinal'};

  @einc = split(',', $r_par->{'model'}->{'einc'});
  if ($einit eq "") {$einit = $einc[0];}
  if ($efinal eq "") {$efinal = $einc[1] * $r_par->{'model'}->{'nescat'};}

  &replace_in_template($r_str, "EINIT",    $einit);
  &replace_in_template($r_str, "EFINAL",   $efinal);
  &replace_in_template($r_str, "MOLECULE", "e + ".$r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "SPIN",     $statespin);
  &replace_in_template($r_str, "SYMMETRY", $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym]);
  &replace_in_template($r_str, "MGVN",     $statesym);
  &replace_in_template($r_str, "STOT",     $spin_multiplicity{lc($statespin)});
  &replace_in_template($r_str, "NERANG",   scalar split(/\s*,\s*/, $r_par->{'model'}->{'nescat'}));
  &replace_in_template($r_str, "NESCAT",   $r_par->{'model'}->{'nescat'});
  &replace_in_template($r_str, "EINC",     $r_par->{'model'}->{'einc'});
  &replace_in_template($r_str, "IEUNIT",   $r_par->{'model'}->{'e_unit'});
  &replace_in_template($r_str, "ISMAX",    $r_par->{'model'}->{'max_multipole'});
  &replace_in_template($r_str, "RAF",      $r_par->{'model'}->{'raf'});
  return 1;
}

# ---------------------- time-delay ------------------------

#sub make_time-delay_input {
#  my ($r_par, $r_str) = @_;
#  &replace_in_template($r_str, "SYMMETRY", $irred_repr{$r_par->{'settings'}->{'symmetry'}}->[$statesym]);
#  return 1;
#}

# ----------------------- radden ---------------------------

sub make_radden_input {
  my ($r_par, $r_str) = @_;
  my @orbitals = ($r_par->{'model'}->{'dm_list'}..($r_par->{'model'}->{'dm_list'} + $r_par->{'model'}->{'req_dm'}));

  &replace_in_template($r_str, "R_START",  	$r_par->{'model'}->{'radius_start'});
  &replace_in_template($r_str, "R_FINISH", 	$r_par->{'model'}->{'radius_finish'});
  &replace_in_template($r_str, "NO_ORBIT",	$r_par->{'model'}->{'req_dm'});
  &replace_in_template($r_str, "ILIST",		($r_par->{'model'}->{'req_dm'} < 0 ? ".false." : ".true."));
  &replace_in_template($r_str, "!",  		($r_par->{'model'}->{'req_dm'} < 0 ? "!" : ""));
  &replace_in_template($r_str, "WHICH_ORBIT",	join(",", @orbitals));
  &replace_in_template($r_str, "WHICHDM",  	($r_par->{'run'}->{'ukrmolplus'} == 1 ? 1 : 2));
  &replace_in_template($r_str, "UKRMOLP",  	($r_par->{'run'}->{'ukrmolplus'} == 1 ? ".false." : ".true."));
  &replace_in_template($r_str, "MOLECULE", 	lc($r_par->{'model'}->{'molecule'}));
  return 1;
}

# ============= end of subroutines for INPUTS ===============

# ================= subroutines for OUTPUTS =================

# ------------------------- sword --------------------------

# From sword output we can determine number of basis functions in each IR
# For target it is the number of all target orbitals
#   saved in $r_par->{'data'}->{'orbitals'}->{'target_all'}
# For scattering calculation it is the number of all target and continuum basis functions
#   saved in $r_par->{'data'}->{'orbitals'}->{'all'}
sub read_sword_output {
  my ($r_par) = @_;
  my $orbitals = 'target_all';
  if ($r_par->{'data'}->{'task'} eq "scattering") { $orbitals = 'all'; }

  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
    while (my $line = <OUTPUT>) {
      chomp($line);
      if ($line =~ m/^\s*BASIS FUNCTIONS PER SYMMETRY\s*$/i) {
        for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
          my $line2 = <OUTPUT>;
          if ($line2 =~ m/^\s*(\d+)\s+(\d+)\s*$/i) {
            $r_par->{'data'}->{'orbitals'}->{$orbitals}->[$i] = $2; 
          }
          else { die "  Something wrong in sword output!!!\n"}
        }
        &print_info("  Number of basis functions in each IR read successfully:\n", $r_par);
        &print_info("      All orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{$orbitals}})."\n", $r_par);
      }
    }
    close(OUTPUT);
  }
  else {
    die "Error: can't open output of sword !\n";
  }
  
  # For scattering calculation we also determine number of continuum orbitals
  if ($r_par->{'data'}->{'task'} eq "scattering") {
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) { 
      $r_par->{'data'}->{'orbitals'}->{'continuum'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'all'}->[$i] 
                                                          - $r_par->{'data'}->{'orbitals'}->{'target_all'}->[$i];
    }
    &print_info("  Continuum orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'continuum'}})."\n", $r_par);
  }
  
  return 1;
}

# ------------------------- molpro --------------------------

# From molpro the energies of the target molecular orbitals are read in
#   and saved in $r_par->{'data'}->{'orbitals'}->{'energies'}
# the number of occupied orbitals is determined from molpro output
# then the number of frozen, active and virtual orbitals 
#   in each IR is determined and the script continues with the next program
# WARNING: at the moment this routine does not handle open shell cases and does not check
#          whether molpro experienced problems with convergence!!!
sub read_molpro_output {
  my ($r_par) = @_;
  my $sym = 0;
  my $nir = $r_par->{'data'}->{'nir'};                    # number of IR
  my $r_e = $r_par->{'data'}->{'orbitals'}->{'energies'};
  my $model = $r_par->{'model'}->{'model'};
  my $header = "NATURAL ORBITALS"; # we look for CASSCF orbitals
  if ($r_par->{'data'}->{'scf_ok'} == 0 || $model =~ m/^SE|pol/) { $header = "ELECTRON ORBITALS"; } # we look for HF orbitals #$model =~ m/SE|pol/
  &print_info("Looking for $header in MOLPRO output...\n", $r_par);

  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
    my $look_for_mos = 0;   # auxiliary variable to decide whether to look for orbital energies
    my $in_mos_section = 0;

    $r_par->{'data'}->{'target'}->{'noccupied'} = 0;

    while ($line = <OUTPUT>) {
      chomp($line);
      if ($line =~ m/\s*Final occupancy:((\s*\d+){$nir})/) { 
         for (my $i = 0; $i < $nir; $i++) { 
            $line =~ s/\s*(\d+)//;
            $r_par->{'data'}->{'orbitals'}->{'occupied'}->[$i] = $1;
         }
         &print_info("Occupancy from the Molpro output: ".join(", ", @{$r_par->{'data'}->{'orbitals'}->{'occupied'}})."\n");
      }
      if ($line =~ m/\s*.RHF STATE\s*(\d+)\.(\d+)\s*Energy\s*([\+\-]?\d+\.\d+)/) {
        $r_par->{'data'}->{'target'}->{'hf_energy'} = $3;
        $r_par->{'data'}->{'target'}->{'hf_symmetry'} = $2 - 1;
        &print_info("  Hartree-Fock energy = $3\n");
      }
      if ($r_par->{'run'}->{'ukrmolplus'} == 1){
      	if ($line =~ m/^.*NUMBER OF CONTRACTIONS:\s*(\d+)\s*\(\s*((\s*\d+[A-Z].*\s*.\s*){$nir})/) { #look for NOB in MOLPRO for INTEGRALS input
          for (my $i = 0; $i < $nir; $i++) {
            $line =~ s/\s*(\d+)[A-Z]//;
            $r_par->{'data'}->{'orbitals'}->{'target_all'}->[$i] = $1;
          }
          &print_info("Total number of orbitals for each symmetry: ".join(", ", @{$r_par->{'data'}->{'orbitals'}->{'target_all'}})."\n");
        }
      }
      if ($line =~ /$header/) { #the data that follow contain the orbital data
        $in_mos_section = 1;
      }
      if ($line =~ m/DATASETS/) { $in_mos_section = 0; } #string 'DATASETS' mark the end of the orbital data

      if ($in_mos_section == 1) { #look for the orbital energies if we are in the section containg the orbital coefficients.
        if ($line =~ /^\s*$/) { #we only look for the info on the orbital energies on the line following a blank line.
          $look_for_mos = 1;
        }
        if ($look_for_mos == 1 && $line =~ m/\s*(\d+)\.(\d)\s*(\d)\s*([\+\-]?\d+.\d+) .*/) {
          $iorb = $1;
          $sym = $2;
          $occ = $3;
          $r_e->{"$sym.$iorb"} = $4;
          $val = $r_e->{"$sym.$iorb"};
          if ($occ != 0) { #occupancy of a particular orbital is only YES/NO
             $r_par->{'data'}->{'target'}->{'noccupied'}++;
          }
          $look_for_mos = 0;
        }
      }
    }
    close(OUTPUT);
  }
  else {
    die "  Error: can't open output of molpro !\n";
  }

  # Checking occupied orbitals
  my @sorted_orb = sort { $r_e->{$a} <=> $r_e->{$b} } keys %$r_e;
  my $r_occ = [];
  &get_number_of_mos_automatically($r_occ, $nir, \@sorted_orb, 
                                   1, $r_par->{'data'}->{'target'}->{'noccupied'});

  my $scf_occ_ok = 1;

  if ($scf_occ_ok == 0) {
    &print_info("  Orbital occupation from get_number_of_mos_automatically:\n", $r_par);
    &print_info("     iocc  = ".join(",", @{$r_occ})."\n", $r_par);
    &print_info("Number of occupied orbitals is inconsistent with the HF configuration: something is wrong with the Molpro output or with parsing of the Molpro output.\n", $r_par);
    die "Stop while reading orbitals...";
  }
  else {
    $r_par->{'data'}->{'scf_ok'} = 1;
    &print_info("  which is correct (hopefully).\n", $r_par);

    $r_par->{'data'}->{'target'}->{'nfull_occ'} = int($r_par->{'model'}->{'nelectrons'} / 2); # number of fully occupied orbitals
    $r_par->{'data'}->{'target'}->{'nhalf_occ'} = $r_par->{'model'}->{'nelectrons'} % 2;      # half occupied orbital

    # Number of various orbitals in each IR
    # Frozen orbitals
    my $starting_orb = 1;
    &get_number_of_mos('frozen', $r_par, \@sorted_orb, $starting_orb);
    &print_info("  Frozen orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'frozen'}})."\n", $r_par);

    # Active orbitals
    $starting_orb += $r_par->{'model'}->{'nfrozen'};
    &get_number_of_mos('active', $r_par, \@sorted_orb, $starting_orb);
    &print_info("  Active orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'active'}})."\n", $r_par);

    # Target orbitals = frozen + active
    for (my $i = 0; $i < $nir; $i++) { 
      $r_par->{'data'}->{'orbitals'}->{'target'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i] 
                                                       + $r_par->{'data'}->{'orbitals'}->{'active'}->[$i];
    }
    &print_info(" Target  orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'target'}})."\n", $r_par);


    # Virtual orbitals
    $starting_orb += $r_par->{'model'}->{'nactive'};
    &get_number_of_mos('virtual', $r_par, \@sorted_orb, $starting_orb);
    &print_info("  Virtual orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'virtual'}})."\n", $r_par);
    
    # Target used orbitals = target + virtual
    for (my $i = 0; $i < $nir; $i++) { 
      $r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'target'}->[$i] 
                                                            + $r_par->{'data'}->{'orbitals'}->{'virtual'}->[$i];
    }
    &print_info("  Used    orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}})."\n", $r_par);
  }
  return 1;
}

# ------------------------- integrals --------------------------

sub read_integrals_output {
  my ($r_par) = @_;
  my $nir = $r_par->{'data'}->{'nir'};                    # number of IR

#  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
#    while ($line = <OUTPUT>) {
#      chomp($line);
#      if ($line =~ m/\s*Number of orbitals per symmetry to keep:\s*((\s*\d+))/) { 
#         &print_info("Number of orbitals per symmetry to keep: ", $r_par);
#         for (my $i = 0; $i < $nir; $i++) { 
#            $line =~ s/\s*(\d+)//;
#            $r_par->{'data'}->{'orbitals'}->{'all'}->[$i] = $1;
#            &print_info("$1, ", $r_par);
#         }
#         &print_info("\n", $r_par);
#      }
#    }
#    close(OUTPUT);
#  }
#  else {
#    die "  Error: can't open output of integrals !\n";
#  }

  my $str = "";
  if (read_file($r_par->{'data'}->{'outputfile'}, \$str)) {
    if ($str =~ m/^.*Number of molecular orbitals in each irreducible representation:\s*((\s*\d+)+)/is) {
      &print_info("Number of molecular orbitals in each irreducible representation: ", $r_par);
      $str =$1;
      &print_info("\n $str ", $r_par);
      my @MO = split(/\s+/, $str);
      for (my $i = 0; $i < $nir; $i++) { 
        $r_par->{'data'}->{'orbitals'}->{'all'}->[$i] = $MO[$i];
      }
      &print_info("\n", $r_par);
    }      
  }
  else {
    die "  Error: can't open output of integrals !\n";
  }

  # For scattering calculation, some orbitals could be deleted.
#  if ($r_par->{'data'}->{'task'} eq "scattering") {

      # We also determine number of used continuum orbitals = all - target_all (= all_used - target)
      for (my $i = 0; $i < $nir; $i++) { 
        $r_par->{'data'}->{'orbitals'}->{'cont_used'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'all'}->[$i] 
                                                            - $r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i];
      }
      &print_info("  Cont used orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'cont_used'}})."\n", $r_par);

      # All used orbitals
      for (my $i = 0; $i < $nir; $i++) { 
        $r_par->{'data'}->{'orbitals'}->{'all_used'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'all'}->[$i];  #$r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i] + $r_par->{'data'}->{'orbitals'}->{'cont_used'}->[$i]; 
      }
      &print_info("  All  used orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'all_used'}})."\n", $r_par);

#  }
  
  return 1;
}

# ------------------------- swscf --------------------------

# From swscf output the energies of target molecular orbitals are read
#   and saved in $r_par->{'data'}->{'orbitals'}->{'energies'}
# according to these energies the number of occupied orbitals
#   in each IR is determined
# if it is not the same as in swscf input file, swscf is called again,
# if it is the same, then the number of frozen, active and virtual orbitals 
#   in each IR is determined and the script continues with the next program
sub read_swscf_output {
  my ($r_par) = @_;
  my $sym = 0;
  my $nir = $r_par->{'data'}->{'nir'};                    # number of IR
  my $r_e = $r_par->{'data'}->{'orbitals'}->{'energies'};

  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
    my $next_vectors_for_symmetry = 0;  # auxiliary variable to handle reading orbital energies
    while ($line = <OUTPUT>) {
      chomp($line);
      if ($line =~ m/^\s*HARTREE FOCK ENERGY\s*([\+\-]?\d+\.\d+)/i) {
        $r_par->{'data'}->{'target'}->{'hf_energy'} = $1;
        &print_info("  Hartree-Fock energy = $1\n", $r_par);
      }
      if ($next_vectors_for_symmetry == 1 || $line =~ m/^\s*VECTORS FOR SYMMETRY\s+(\d+)\s*$/i) {
        if ($next_vectors_for_symmetry == 0) { $sym = $1; }
        else                                 { $next_vectors_for_symmetry = 0; }
        $iorb = 0;
        ENERGY: while ($line = <OUTPUT>) {
          if ($line =~ m/^\s*VECTORS FOR SYMMETRY\s+(\d+)\s*$/i) {
            $sym = $1;
            $next_vectors_for_symmetry = 1;
            last ENERGY;
          }
          chomp($line);
          if ($line =~ m/^\s*[\+\-]?\d+\.\d+/) {       # line with real numbers only -> energies
            while ($line !~ /^\s*$/) {
              $line =~ s/^\s*([\+\-]?\d+\.\d+)//;
              $iorb++;
              $r_e->{"$sym.$iorb"} = $1;
            }
          }
          if ($line =~ m/^\s*\d+\s*[\+\-]?\d+\.\d+/) { # line beginning with integer followed by real number -> energies of deleted orbitals
            while ($line !~ /^\s*$/) {
              $line =~ s/^\s*\d+\s*([\+\-]?\d+\.\d+)//;
              $iorb++;
              $r_e->{"$sym.$iorb"} = $1;
            }
          }
        } # ENERGY loop
      } # if ($line =~ m/^\s*SUB SYM\s* ...)
    }
    close(OUTPUT);
  }
  else {
    die "  Error: can't open output of swscf !\n";
  }

  # Checking occupied orbitals
  my @sorted_orb = sort { $r_e->{$a} <=> $r_e->{$b} } keys %$r_e;
  my $r_occ = [];
  &get_number_of_mos_automatically($r_occ, $nir, \@sorted_orb, 
                                   1, $r_par->{'data'}->{'target'}->{'noccupied'});
  my $scf_occ_ok = 1;
  for (my $i = 0; $i < $nir; $i++) { 
    if ($r_occ->[$i] != $r_par->{'data'}->{'orbitals'}->{'occupied'}->[$i]) {
      $scf_occ_ok = 0;
    }
  }
  &print_info("  Orbital occupation from swscf output:\n", $r_par);
  &print_info("     iocc  = ".join(",", @{$r_occ})."\n", $r_par);
  if ($scf_occ_ok == 0) {
    $r_par->{'data'}->{'orbitals'}->{'occupied'} = $r_occ;

    # For open shells we have to check also ifock
    my $r_ifock = $r_par->{'data'}->{'orbitals'}->{'ifock'};
    if ($r_par->{'data'}->{'target'}->{'nhalf_occ'} == 1) {

      # Symmetry of the last occupied orbital
      $sorted_orb[$r_par->{'data'}->{'target'}->{'noccupied'} - 1] =~ m/(\d+)\.(\d+)/;
      my $sym_half_occ = $1 - 1;

      # Set symmetry of HF target state
      $r_par->{'data'}->{'target'}->{'hf_symmetry'} = $sym_half_occ;

      # Setting new ifock with 2 in the right position
      my $p = 0;
      for (my $i = 0; $i < $nir; $i++) { 
        for (my $j = 0; $j < $r_occ->[$i]; $j++) { 
          if ($i == $sym_half_occ && $j == $r_occ->[$i] - 1) { 
            $r_ifock->[$p] = 2;
          }
          else { 
            $r_ifock->[$p] = 1;
          }
          $p++;
        }
      }
    }
    &print_info("     ifock = ".join(",", @{$r_ifock})."\n", $r_par);
    &print_info("  It is necessary to run swscf again !\n", $r_par);
  }
  else {
    $r_par->{'data'}->{'scf_ok'} = 1;
    &print_info("  which is correct (hopefully).\n", $r_par);

    # Number of various orbitals in each IR
    # Frozen orbitals
    my $starting_orb = 1;
    &get_number_of_mos('frozen', $r_par, \@sorted_orb, $starting_orb);
    &print_info("  Frozen orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'frozen'}})."\n", $r_par);

    # Active orbitals
    $starting_orb += $r_par->{'model'}->{'nfrozen'};
    &get_number_of_mos('active', $r_par, \@sorted_orb, $starting_orb);
    &print_info("  Active orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'active'}})."\n", $r_par);

    # Target orbitals = frozen + active
    for (my $i = 0; $i < $nir; $i++) { 
      $r_par->{'data'}->{'orbitals'}->{'target'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i] 
                                                       + $r_par->{'data'}->{'orbitals'}->{'active'}->[$i];
    }
    &print_info("  Target  orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'target'}})."\n", $r_par);

    # Virtual orbitals
    $starting_orb += $r_par->{'model'}->{'nactive'};
    &get_number_of_mos('virtual', $r_par, \@sorted_orb, $starting_orb);
    &print_info("  Virtual orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'virtual'}})."\n", $r_par);
    
    # Target used orbitals = target + virtual
    for (my $i = 0; $i < $nir; $i++) { 
      $r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'target'}->[$i] 
                                                            + $r_par->{'data'}->{'orbitals'}->{'virtual'}->[$i];
    }
    &print_info("  Used    orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}})."\n", $r_par);
  }

#  &print_info("  Energies of molecular orbitals: in X.Y X stands for IR and Y for number of MO\n", $r_par);
#  my $p = 0;
#  foreach my $key (@sorted_orb) {
#    &print_info("    $key    $r_e->{$key}\n", $r_par);
#    if (10 == $p++) { last; }
#  }
  return 1;
}

# This auxiliary subroutine determines number of orbitals
# in each irreducible representation simply by checking
# whether there are specified orbitals in the settings
# or whether to determine them automatically
sub get_number_of_mos {
  my ($which_orbitals, $r_par, $r_sorted_orb, $starting_orb) = @_;
  my $nir = $r_par->{'data'}->{'nir'};
  
  my @given_orbs = @{$r_par->{'model'}->{$which_orbitals.'_orbs'}}[0..$nir-1];
  my $ngiven_orbs = 0;
  $ngiven_orbs += $_ for @given_orbs;

  my $get_orbitals_automatically = 1;
  if ($ngiven_orbs > 0) {
    if ($ngiven_orbs != $r_par->{'model'}->{'n'.$which_orbitals}) {
      &print_info("\nWarning: Number of $which_orbitals orbitals ".$r_par->{'model'}->{'n'.$which_orbitals}." specified in 'nvirtual'\n", $r_par);
      &print_info("         inconsistent with '${which_orbitals}_orbs' ".join(",", @given_orbs)." !\n", $r_par);
      &print_info("         Number of orbitals will be determined automatically !\n\n", $r_par);
    }
    else {
      $get_orbitals_automatically = 0;
      @{$r_par->{'data'}->{'orbitals'}->{$which_orbitals}} = @given_orbs;
    }
  }
  if ($get_orbitals_automatically == 1) {
    &get_number_of_mos_automatically($r_par->{'data'}->{'orbitals'}->{$which_orbitals}, $nir, $r_sorted_orb, 
                                     $starting_orb, $starting_orb - 1 + $r_par->{'model'}->{'n'.$which_orbitals});
  }
  return 1;
}

# This auxiliary subroutine determines automatically number of orbitals 
# in each irreducible representation to an array specified using reference $r_orb
# starting from the $first_orb to the $last_orb orbital ordered by energy
sub get_number_of_mos_automatically {
  my ($r_orb, $nir, $r_sorted_orb, $first_orb, $last_orb) = @_;
  for (my $i = 0; $i < $nir; $i++) { $r_orb->[$i] = 0; }
  for (my $i = $first_orb - 1; $i < $last_orb; $i++) { 
    $r_sorted_orb->[$i] =~ m/(\d+)\.(\d+)/;
    $r_orb->[$1 - 1]++;
  }
  return 1;
}

# ------------------------ swedmos -------------------------

# From swedmos output we can determine 
# number of deleted orbitals after othogonalization
#   saved in $r_par->{'data'}->{'orbitals'}->{'deleted'}
sub read_swedmos_output {
  my ($r_par) = @_;
  my $nir = $r_par->{'data'}->{'nir'};                    # number of IR
  my $ir = 0;
  my $r_del_orb = [];

  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
    while (my $line = <OUTPUT>) {
      chomp($line);
      if ($line =~ m/^\s*No. of deleted orbitals =\s+(\d+)\s*$/i) {
        $r_del_orb->[$ir] = $1;
        $ir++;
      }
    }
    close(OUTPUT);
  }
  else {
    die "Error: can't open output of swedmos !\n";
  }
  
  # For target calculation, no orbitals should be deleted
  if ($r_par->{'data'}->{'task'} eq "target") {
    if ($ir > 0) {
      &print_info("\nWarning: some target orbitals deleted by swedmos !!!\n", $r_par);
    }
    else {
      &print_info("  No orbitals deleted for target.\n", $r_par);
    }
  }
  # For scattering calculation, some orbitals could be deleted.
  else {
    if ($ir == $nir) {
      $r_par->{'data'}->{'orbitals'}->{'deleted'} = $r_del_orb;
      &print_info("  Deleted orbitals: ".join(",", @{$r_del_orb})."\n", $r_par);

      # We also determine number of used continuum orbitals = continuum - deleted
      for (my $i = 0; $i < $nir; $i++) { 
        $r_par->{'data'}->{'orbitals'}->{'cont_used'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'continuum'}->[$i] 
                                                            - $r_par->{'data'}->{'orbitals'}->{'deleted'}->[$i];
      }
      &print_info("  Cont used orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'cont_used'}})."\n", $r_par);

      # And finally number of all used orbitals = target_used + cont_used
      for (my $i = 0; $i < $nir; $i++) { 
        $r_par->{'data'}->{'orbitals'}->{'all_used'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i] 
                                                           + $r_par->{'data'}->{'orbitals'}->{'cont_used'}->[$i];
      }
      &print_info("  All  used orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'all_used'}})."\n", $r_par);
    }
    else {
      &print_info("\nWarning: in swedmos some target orbitals deleted by swedmos !!!\n", $r_par);
    }
  }
  
  return 1;
}

# ------------------------ congen --------------------------

# In congen output we only check that some configuartions were generated
# if there is no configuration
#   $r_par->{'data'}->{'no_scatci'} is set to 1
# to skip scatci calculation
sub read_congen_output {
  my ($r_par) = @_;

  my $str = "";
  my $statesym    = $r_par->{'data'}->{'target'}->{'symmetry'};                      # number 0,1,2,...
  if (read_file($r_par->{'data'}->{'outputfile'}, \$str)) {
    if ($str =~ m/TOTAL NUMBER OF CSF'S GENERATED IS\s+(\d+)/is) {
      &print_info("  Total number of CSF's generated is $1", $r_par);
      $r_par->{'data'}->{'target'}->{'pol_NOCSF'}->[$statesym]= $1;
      if ($1 == 0) {
        $r_par->{'data'}->{'no_scatci'} = 1;
        &print_info("!!!\n  Other runs for this spin-symmetry will be skipped!\n", $r_par);
      }
      else {
        &print_info(".\n", $r_par);
      }
    }
    elsif ($str =~ m/ERROR IN REFORB DATA/is) {
      $r_par->{'data'}->{'no_scatci'} = 1;
      &print_info("  Error in REFORB in congen input !!!\n  Probably there is no available orbital of required symmetry !\n  Other runs for this spin-symmetry will be skipped!\n", $r_par);
    }
    else {
      die "  Error: incorrect output of congen: file $r_par->{'data'}->{'outputfile'} !\n";
    }
  }
  else {
    die "  Error: can't open output of congen: file $r_par->{'data'}->{'outputfile'} !\n";
  }

  return 1;
}

# ------------------------ scatci --------------------------

# From scatci output we can determine energies of states
#   saved in $r_par->{'data'}->{'target'}->{'states'}
#   and   in $r_par->{'data'}->{$task}->{'states_all_geom'}
sub read_scatci_output {
  my ($r_par) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statespin = $spin_multiplicity{lc($r_par->{'data'}->{$task}->{'spin'})};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $igeom = $r_par->{'data'}->{'igeom'} - 1; # index of geometry for storage

  my $str = "";
  if (read_file($r_par->{'data'}->{'outputfile'}, \$str)) {
    if ($str =~ m/^.*EIGEN-ENERGIES\s+(([\+\-]?\d+\.\d+\s*)+)$/is) {
      $energies_list = $1;
      $energies_list =~ s/\s*$//;
      &print_info("  Found energies:\n", $r_par);
      my @energies = split(/\s+/, $energies_list);
      for (my $i = 1; $i <= scalar @energies; $i++) {
        if ($task eq "target") { $r_par->{'data'}->{'target'}->{'states'}->{"$statespin.$statesym.$i"} = $energies[$i-1]; }
#        &print_info("  $statespin.$statesym.$i -> ".$energies[$i-1]."\n", $r_par);
        if ($igeom > 0) { $r_par->{'data'}->{$task}->{'states_all_geom'}->{"$statespin.$statesym.$i"}->[$igeom] = $energies[$i-1]; }
        else {            $r_par->{'data'}->{$task}->{'states_all_geom'}->{"$statespin.$statesym.$i"} = [];
                          $r_par->{'data'}->{$task}->{'states_all_geom'}->{"$statespin.$statesym.$i"}->[$igeom] = $energies[$i-1];
        }
      }
    }
    else {
      die "  Error: incorrect output of scatci: file $r_par->{'data'}->{'outputfile'} !\n";
    }
    if ($r_par->{'run'}->{'parallel_diag'} == 0 && $task eq 'scattering') { #$r_par->{'run'}->{'ukrmolplus'} == 1 &&
      if ($str =~ m/^.*Scat.*NOCSF=\s*(\d+)\s*/is) { 
#      if ($str =~ m/.*MOCSF =\s*(\d+)\s*dimension final Hamiltonian.*/is) {
        $r_par->{'data'}->{'scattering'}->{'cont_csf'} = $1;
        &print_info("LCONT = $1\n", $r_par);
      }
      else { 
        die "  Error: incorrect output of scatci: file $r_par->{'data'}->{'outputfile'} !\n";
      }
    }
  }
  else {
    die "  Error: can't open output of scatci: file $r_par->{'data'}->{'outputfile'} !\n";
  }
  
  return 1;
}

# ------------------------ hamdiag --------------------------

sub read_hamdiag_output {
  my ($r_par) = @_;
  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
     while (my $line = <OUTPUT>) {
        if ($line =~ m/Only the first\s*(\d+)\s*CI coefficients from each eigenvector will be saved./) {
           $r_par->{'data'}->{'scattering'}->{'cont_csf'} = $1;
           &print_info("LCONT = $1\n", $r_par);
        }
     }
  }
  else {
    die "  Error: can't open output of hamdiag: file $r_par->{'data'}->{'outputfile'} !\n";
  }
 
  close OUTPUT; 

  return 1;
}

# ------------------------ denprop -------------------------

# From denprop output we have to determine order of target states again
# because it can be different from ordering obtained from scatci outputs if there are degenerate states
#   $r_par->{'data'}->{'target'}->{'states'}
sub read_denprop_output {
  my ($r_par) = @_;
  my $istates = 0;
  my %nstates = ();
  my %str_idtarg = (); # used below for setting idtarg input for outerres

  # Read denprop output and print information about target states
  &print_info("  Target states:\n", $r_par);
  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
    while (my $line = <OUTPUT>) {
      chomp($line);
      if ($line =~ m/^\s*\d+\s+(\d+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+([^\s]+)\s+state no.\s+\d+\s+([^\s]+)\s+([^\s]+)/i) {
        my $spin_sym = "$3.$2";
        &print_info("    No. $1 - $5($3)  $6($2)  $4\n", $r_par);
        if ($nstates{$spin_sym}) { $nstates{$spin_sym}++; }
        else                     { $nstates{$spin_sym} = 1; }
        $r_par->{'data'}->{'target'}->{'ordered_states'}->[$istates] = "$spin_sym.$nstates{$spin_sym}";
        $istates++;
        $str_idtarg{"$spin_sym.$nstates{$spin_sym}"} = $istates;
      }
    }
    close(OUTPUT);
  }
  else {
    die "Error: can't open output of denprop !\n";
  }
  
  # Finally set auxiliary array 'id_targ' used later in scattering.outerres inputs
  my $r_spinsym_order = $r_par->{'data'}->{'target'}->{'spinsym_order'};
  my $it = 0;
  for (my $i = 0; $i < scalar @{$r_spinsym_order}; $i++) {
    if ($r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$r_spinsym_order->[$i]}) {
      my ($spin, $sym) = split(/\./, $r_spinsym_order->[$i]);
      for ($is = 1; $is <= $r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$r_spinsym_order->[$i]}; $is++) {
        $r_par->{'data'}->{'target'}->{'idtarg'}->[$it] = $str_idtarg{"$spin.$sym.$is"};
        $it++;
      }
    }
  }
 
  return 1;
}

# ----------------------- outerres -------------------------

# From outerres output we can determine position and width of the resonance
#   saved in $r_par->{'data'}->{'scattering'}->{'resonance_position'}
#   and   in $r_par->{'data'}->{'scattering'}->{'resonance_width'}
sub read_outerres_output {
  my ($r_par) = @_;
  my $state = $spin_multiplicity{lc($r_par->{'data'}->{'scattering'}->{'spin'})}.'.'.$r_par->{'data'}->{'scattering'}->{'symmetry'};
  my $igeom = $r_par->{'data'}->{'igeom'}; # index of geometry
  my ($position, $width) = (1.0e10, 0.0);

  my $str = "";
  my @positions = ();
  my @widths    = ();
  if (read_file($r_par->{'data'}->{'outputfile'}, \$str)) {
    while ($str =~ s/Fitted resonance parameters\s+Positions\s*\/\s*Ryd\s*(.*?)\s*Widths\s*\/\s*Ryd\s*(.*?)\s*Background//is) {
      $str_positions = $1;          $str_widths = $2;
      $str_positions =~ s/D/E/g;    $str_widths =~ s/D/E/g;
      push(@positions, split(/\s+/, $str_positions));
      push(@widths, split(/\s+/, $str_widths));
    }
    for (my $i = 0; $i < scalar @positions; $i++) {
      if ($positions[$i] < $position && $positions[$i] > 0.0) {
        $position = $positions[$i]; $width = $widths[$i];
      }
    }
    if ($position < 1.0e10) {
      &print_info("  Resonance position: $position Ryd\n", $r_par);
      &print_info("  Resonance width:    $width Ryd\n", $r_par);
    }
    else {
      &print_info("  Warning: there is no resonance in output of outerres: file $r_par->{'data'}->{'outputfile'} !\n", $r_par);
      $position = 0.0;
    }
    if ($igeom == 1) { 
      $r_par->{'data'}->{'scattering'}->{'resonance_position'}->{"$state"} = []; 
      $r_par->{'data'}->{'scattering'}->{'resonance_width'}->{"$state"} = []; 
    }
    $r_par->{'data'}->{'scattering'}->{'resonance_position'}->{"$state"}->[$igeom - 1] = 0.5 * $position; # 0.5 is for conversion to a.u.
    $r_par->{'data'}->{'scattering'}->{'resonance_width'}->{"$state"}->[$igeom - 1] = 0.5 * $width;       # 0.5 is for conversion to a.u.
  }
  else {
    die "  Error: can't open output of scatci: file $r_par->{'data'}->{'outputfile'} !\n";
  }
  
  return 1;
}

# ============= end of subroutines for OUTPUTS ==============

# ----------------- auxiliary subroutines -------------------

# Run a specific program including making an input before
#   and reading an output afterwards
sub run_code {
  my ($program, $r_par) = @_;
  my $r_dirs = $r_par->{'dirs'};
  my $task = $r_par->{'data'}->{'task'}; 
  
  # Setting auxiliary variables for input and output files
  $r_par->{'data'}->{'program'} = $program; 
  $r_par->{'data'}->{'inputfile'} = "$r_dirs->{'inputs'}$bs$task.$program.inp";
  $r_par->{'data'}->{'outputfile'} = "$r_dirs->{'outputs'}$bs$task.$program.out";
  
  # If CONGEN, SCATCI, OUTERRES or TIMEDEL/TIME-DELAY are to be run then input and output files contain information about spin and symmetry
  # e.g. target.scatci.triplet.A2, or scattering.congen.doublet.B1.inp
  if ($program =~ /(congen|scatci|hamdiag|outerres|timedel|time-delay)/) {
    my $statespin = $r_par->{'data'}->{$task}->{'spin'};
    my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$r_par->{'data'}->{$task}->{'symmetry'}];
    $r_par->{'data'}->{'inputfile'} = "$r_dirs->{'geom'}$bs$r_dirs->{'inputs'}$bs$task.$program.$statespin.$str_ir.inp";
    $r_par->{'data'}->{'outputfile'} = "$r_dirs->{'geom'}$bs$r_dirs->{'outputs'}$bs$task.$program.$statespin.$str_ir.out";
  }
 
  # If MPOUTRD prepared SWEDMOS input (only if MOLPRO was called), skip creating input
  # also if 'use_templates' option is set to 0, skip creating inputs
  # otherwise use template to create input file
  if ($r_par->{'run'}->{'use_templates'} == 1 && 
      !($task eq "target" && $program eq "swedmos" && $r_par->{'run'}->{'molpro'} == 1)) {
    # Read template for input
    my $str = "";
    if (!&read_file("$r_dirs->{'templates'}$bs$task.$program.inp", \$str)) {
      if(!&read_file("$r_dirs->{'templates'}$bs$program.inp", \$str)) {
        &print_info("Warning: no template file for $program.inp !\n", $r_par);
      }
    }
    # Call make_$program_input to modify input if exists and then save it 
    my $subroutine = "make_${program}_input";
    if (exists &$subroutine) {
      &$subroutine($r_par, \$str);
    }
    &save_file($r_par->{'data'}->{'inputfile'}, \$str);
  }

  # Run the program if it is to be run and print info ...
  if (&run_program_this_time($r_par) == 1) {
    if ($program =~ /(congen|scatci|hamdiag|outerres|timedel|time-delay)/) {
      &print_info("Running $program for $r_par->{'data'}->{$task}->{'spin'} $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$r_par->{'data'}->{$task}->{'symmetry'}] ...", $r_par);
    }
    else {
      &print_info("Running $program ...", $r_par);
    }
    my $dir = "";
    if    ($program eq "molpro")                { $dir = $r_dirs->{'molpro'}; }
    elsif ($program eq "integrals")             { $dir = $r_dirs->{'integrals'}; }
    elsif ($program =~ /(outerres|timedel|time-delay)/) { $dir = $r_dirs->{'bin_out'}; }
    else                                        { $dir = $r_dirs->{'bin_in'}; }
    my $command = "";
#--------------------------- command to run MOLPRO -------------------------------#
    if ($program eq "molpro") {
      my $input  = $r_par->{'data'}->{'inputfile'};
      my $output = $r_par->{'data'}->{'outputfile'};
      if ($sys eq win) { 
        $input  =~ s/$rebs/\//g; # On Windows, MOLPRO is running under cygwin and there should be used '/' in a file path
        $output =~ s/$rebs/\//g;
      };
      $command = "$dir$bs$program$ext_exe --no-xml-output $input --output $output";
    #   $command = "echo molpro";
    }
#--------------------------- command to run INTEGRALS -------------------------------#
    if ($program eq "integrals") {
      system("$cp_cmd $r_par->{'data'}->{'inputfile'} inp");
	if ($r_par->{'run'}->{'computer'} == 0) { # if uqbar	
          $command = "$dir$bs$program$ext_exe <inp >$r_par->{'data'}->{'outputfile'}";
	}
	elsif ($r_par->{'run'}->{'computer'} == 1) { # if impact cluster
          $command = "/padata/beta/users/rmatrix/ompi-1.8.5-intel/bin/mpiexec -n $r_par->{'run'}->{'no_processors'}  $dir$bs$program$ext_exe <inp >$r_par->{'data'}->{'outputfile'}"; #/csoft/intel-12/impi/4.0.3.008/intel64/bin
	}
	elsif ($r_par->{'run'}->{'computer'} == 2) { # if archer
          $command = "echo integrals";
	}
    }
#--------------------------- command to run HAMDIAG -------------------------------#
    elsif ($program eq "hamdiag") {
      system("$cp_cmd $r_dirs->{'templates'}${bs}matrixdat .");
	if ($r_par->{'run'}->{'computer'} == 0) { # if uqbar	
          $command = "time mpirun -np $r_par->{'run'}->{'no_processors'} $dir$bs$program$ext_exe >$r_par->{'data'}->{'outputfile'} ";
	}
	elsif ($r_par->{'run'}->{'computer'} == 1) { # if impact cluster
          $command = "mpiexec -np $r_par->{'run'}->{'no_processors'} $dir$bs$program$ext_exe >$r_par->{'data'}->{'outputfile'} ";
	}
    }
#--------------------------- command to run other programs -------------------------------#    
    else {
      $command = "$dir$bs$program$ext_exe <$r_par->{'data'}->{'inputfile'} >$r_par->{'data'}->{'outputfile'}";
    }
#    &print_info("$command\n", $r_par);
    system("$command");
    &print_info(" done.\n", $r_par);
  }
  else {
    &print_info("Skipping running $program.\n", $r_par);
  }

#--------------------------- command to copy INTEGRALS output to outputs directory -------------------------------#
    if ($program eq "integrals") {
      if ($r_par->{'run'}->{'computer'} == 2) { return 1;}
      system("$cp_cmd log_file.0 $r_par->{'data'}->{'outputfile'}");
    }


  # Call read_$program_output to gather information necessary for later inputs
  my $subroutine = "read_${program}_output";
  if (exists &$subroutine) {
    &print_info("Reading $program output ...\n", $r_par);
    &$subroutine($r_par);
    &print_info(" ... done.\n", $r_par);
  }
  
  return 1;
}

sub run_program_this_time {
  my ($r_par) = @_;
  my $task = $r_par->{'data'}->{'task'}; 
  my $program = $r_par->{'data'}->{'program'}; 

  # If 'only' in the hash array 'run' is not empty, then only those programs will run which are specified
  my $run_this_time = 0;
  if ($r_par->{'run'}->{'only'} eq "") {
    $run_this_time = 1;
  }
  else {
    if ("$task-$program" =~ /$r_par->{'run'}->{'only'}/) {
      $run_this_time = 1;
    }
  }
  return $run_this_time;
}

sub run_system {
  my ($command, $r_par) = @_;
  if (&run_program_this_time($r_par) == 1) {
    system($command);
  }
  return 1;
}

sub run_sub {
  my ($sub_name, $r_par, @args) = @_;
  if (&run_program_this_time($r_par) == 1) {
    &$sub_name(@args);;
  }
  return 1;
}

sub print_info {
  my ($message, $r_par) = @_;
  if ($r_par->{'run'}->{'print_info'} eq "file") {
    open(LOG, ">>$r_par->{'data'}->{'logfile'}");
    print LOG $message;
    close(LOG);
  }
  elsif ($r_par->{'run'}->{'print_info'} eq "screen") {
    print $message;
  }
  return 1;
}

# Auxiliary subroutine which adds the current working directory to all user specified directories
sub add_cwd_to_dirs {
  my ($r_dirs) = @_;
  foreach my $dir (keys %$r_dirs) {
    if ($dir =~ /(bin_in|bin_out|integrals|basis|templates|output|molpro)/) {
      if ($r_dirs->{$dir} eq ".") { 
        $r_dirs->{$dir} = $r_dirs->{'cwd'};
      }
      else { 
        if ($r_dirs->{$dir} !~ /^(\/|[A-Z]:|[a-z]:)/) { # if it is not a full path ( /... on linux, X:... on win)
          $r_dirs->{$dir} = "$r_dirs->{'cwd'}$bs$r_dirs->{$dir}";
        }
      }
    }
  }
  return 1;
}

sub replace_in_template {
  my ($r_str, $what, $with) = @_;
  $with =~ s/,\s*$//; # Get rid of ',' at the end
  $$r_str =~ s/>>>$what<<</$with/gi;
  return 1;
}

sub replace_all_in_template {
  my ($r_str, $r_namelist) = @_;
  foreach my $what (keys %{$r_namelist}) {
    &replace_in_template($r_str, $what, $r_namelist->{$what});
  }
  return 1;
}

sub add_range {
  my ($r_array, $from, $to, $step) = @_;
  my $epsilon = 1e-10;       # used to deal with rounding errors
  my $n = scalar @$r_array;  # starting number of array elements
  my $value = $from;         # initial value
  my %present = ();          # auxiliary
  while ($value <= $to + $epsilon) {
    $r_array->[$n] = $value; 
    $value += $step;
    $n++;
  }
  @$r_array = &delete_duplicates(@$r_array);
  @$r_array = sort (@$r_array);
  return 1;
}

sub delete_duplicates { my %seen; grep !$seen{$_}++, @_ }

# This subroutine makes a simple initial guess for HF calculations using swscf
#  - all electrons will be in totally symmetric molecular orbitals
sub consistent_orbital_setting {
  my ($r_par) = @_;
  my $nir = $r_par->{'data'}->{'nir'};                    # number of IR
  my $consistent_setting = 1;
  for my $which_orbitals ('frozen', 'active', 'virtual') {
    my $norbitals = 0;
    $norbitals += $_ for @{$r_par->{'model'}->{$which_orbitals.'_orbs'}}[0..$nir-1];
    if ($norbitals != $r_par->{'model'}->{'n'.$which_orbitals}) {
      $consistent_setting = 0;
      if ($norbitals > 0) { # something wrong
        &print_info("Warning: Number of $which_orbitals orbitals ".$r_par->{'model'}->{'n'.$which_orbitals}."\n", $r_par);
        &print_info("         inconsistent with settings in '${which_orbitals}_orbs' !\n", $r_par);
        &print_info("         Number of orbitals in each IRs will be determined automatically !\n\n", $r_par);
      }
    }
  }
  return $consistent_setting;
}

# This subroutine makes a simple initial guess for HF calculations using swscf
#  - all electrons will be in totally symmetric molecular orbitals
sub initial_orbitals_occupation_guess {
  my ($r_par) = @_;
  my $nir = $r_par->{'data'}->{'nir'};                    # number of IR

  $r_par->{'data'}->{'target'}->{'nfull_occ'} = int($r_par->{'model'}->{'nelectrons'} / 2); # number of fully occupied orbitals
  $r_par->{'data'}->{'target'}->{'nhalf_occ'} = $r_par->{'model'}->{'nelectrons'} % 2;      # half occupied orbital
  $r_par->{'data'}->{'target'}->{'noccupied'} = $r_par->{'data'}->{'target'}->{'nfull_occ'} + $r_par->{'data'}->{'target'}->{'nhalf_occ'};

  $r_par->{'data'}->{'orbitals'}->{'occupied'}->[0] = $r_par->{'data'}->{'target'}->{'noccupied'};
  $r_par->{'data'}->{'orbitals'}->{'ifock'}->[0] = 1;
  for (my $i = 1; $i < $nir; $i++) { 
    $r_par->{'data'}->{'orbitals'}->{'occupied'}->[$i] = 0;
  }
  for (my $i = 1; $i < $r_par->{'data'}->{'target'}->{'nfull_occ'}; $i++) {
    $r_par->{'data'}->{'orbitals'}->{'ifock'}->[$i] = 1;
  }
  if ($r_par->{'data'}->{'target'}->{'nhalf_occ'} == 1) {
    $r_par->{'data'}->{'orbitals'}->{'ifock'}->[$r_par->{'data'}->{'target'}->{'nfull_occ'}] = 2;
  }
  &print_info("SCF: using initial orbital occupation:\n", $r_par);
  &print_info("     iocc  = ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'occupied'}})."\n", $r_par);
  &print_info("     ifock = ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'ifock'}})."\n", $r_par);

  # Set symmetry of HF target state
  $r_par->{'data'}->{'target'}->{'hf_symmetry'} = 0;
  return 1;
}

# Gather information about eigenphase sums from files eigenph.doublet.A1, ...
# Only those spin symmetries are gathered, which were required in settings->nscat_states
sub gather_eigenphases {
  my ($prefix, $r_par) = @_;
  my @fhs = (); 
  my @lines = ();
  # Open all files to filehandles fh0, ...
  my $fhcounter = 0;
  my $data_description = ();
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] == 1) {
        $fhs[$fhcounter] = "fh$fhcounter";
        my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i];
        open($fhs[$fhcounter], "<", "$prefix.$statespin.$str_ir");
        push @data_description, "$statespin $str_ir";
        $fhcounter++;
      }
    }
  }
  open(OUT, ">", "$prefix.total");
  # Print the first line
  print OUT "#Energy    ";
  for (my $i = 0; $i < $fhcounter; $i++) {
    print OUT "  $data_description[$i]";
  }
  print OUT "\n";
  # Read lines in all files at once and save
  my $in = $fhs[0];
  while ($lines[0] = <$in>) {
    for (my $i = 1; $i < $fhcounter; $i++) {
      my $otherin = $fhs[$i];
      $lines[$i] = <$otherin>;
    }
    if ($lines[0] =~ /^\s*\d\.\d+[ED][-\+]\d+/) {
      for (my $i = 0; $i < $fhcounter; $i++) {
        my ($nothing, $energy, $eigenph) = split(/\s+/, $lines[$i]);
        $energy =~ s/D/E/ig;
        $eigenph =~ s/D/E/ig;
        if ($i == 0) { print OUT "$energy  $eigenph"; }
        else         { print OUT "  $eigenph"; }
      }
      print OUT "\n";
    }
  }
  for (my $i = 0; $i < $fhcounter; $i++) {
    close($fhs[$i]);
  }
}

# Convert files with the cross sections in such a way that gnuplot can plot them 
# Also gather information about total cross sections from these files
# Only those spin symmetries are gathered, which were required in settings->scattering_states
sub convert_and_gather_cross_sections {
  my ($prefix, $r_par) = @_;
  my @lines = ();
  my @data = ();
  my $state = 0;
  my $title = "";
  
  # first convert all files to be readable by gnuplot
  for (my $ir = 0; $ir < $r_par->{'data'}->{'nir'}; $ir++) {
    foreach my $statespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$ir] == 1) {
        my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$ir];
        open(IN, "$prefix.$statespin.$str_ir");
        $state = 0;
        while (my $line = <IN>) {
          chomp($line);
          if ($line =~ /(CROSS SECTIONS.*\s)(\d+)\s*$/) {
            $title = $1;
            my $new_state = $2;
            if ($state < $new_state) {
              if ($state > 0) {
                &print_and_add_cross_sections($state, "$prefix.$statespin.$str_ir.$state", "$title$state", \@lines, \@data);
              }
              for (my $i = 0; $i < scalar @lines; $i++) { $lines[$i] = ""; }
            }
            $state = $new_state;
          }
          elsif ($line =~ /^\s*I\s+E/) {
            if ($lines[0] && $lines[0] ne "") { 
              $line =~ s/^\s*I\s+E\([^\)]*\)\s\s\s//;
              $lines[0] .= $line;
            }
            else {
              $line =~ s/^\s*I\s+/#   /;
              $lines[0] = $line;
            }
          }
          elsif ($line =~ /^\s*(\d+)\s+\d+\.\d+/) {
            my $i = $1; 
            $line =~ s/D/E/gi;
            if ($lines[$i] && $lines[$i] ne "") { 
              $line =~ s/^\s*\d+\s+[^\s]+//;
              $lines[$i] .= $line;
            }
            else {
              $line =~ s/^\s*\d+\s+/  /;
              $lines[$i] = $line;
            }
          }
        } # while (my $line = <IN>)
        if ($state > 0) {
          &print_and_add_cross_sections($state, "$prefix.$statespin.$str_ir.$state", "$title$state", \@lines, \@data);
        } 
        close(IN);
      } # if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$ir] == 1) {
    } # foreach my $statespin
  } # for (my $ir = 0; $ir < $r_par->{'data'}->{'nir'}; $ir++)
  for (my $is = 1; $is <= $state; $is++) {
    open(OUT, ">$prefix.total.$is");
    print OUT "# TOTAL $title$is\n";
    print OUT "$lines[0]";
    for (my $i = 0; $i < scalar @{$data[$is]}; $i++ ) {
      print OUT "  ".join("  ", map(sprintf("%14.6e", $_), @{$data[$is]->[$i]}))."\n";
    }
    close(OUT);
  }
  return 1;
}

sub print_and_add_cross_sections {
  my ($state, $filename, $title, $r_lines, $r_data) = @_;
  open(OUT, ">$filename");
  print OUT "# $title\n";
  for (my $i = 0; $i < scalar @{$r_lines}; $i++ ) {
    print OUT "$r_lines->[$i]\n";
  }
  close(OUT); 
  if (!($r_data->[$state])) { $r_data->[$state] = []; }
  for (my $i = 1; $i < scalar @{$r_lines}; $i++) {
    $r_lines->[$i] =~ s/^\s*//;
    $r_lines->[$i] =~ s/\s*$//;
    my @values = split(/\s+/, $r_lines->[$i]);
    if ($r_data->[$state]->[$i]) {
      for (my $j = 1; $j < scalar @values; $j++) {
        $r_data->[$state]->[$i]->[$j] += $values[$j];
      }
    }
    else {
      $r_data->[$state]->[$i] = [];
      for (my $j = 0; $j < scalar @values; $j++) {
        $r_data->[$state]->[$i]->[$j] = $values[$j];
      }
    }
  }
  return 1;
} 

# Save target energies from $r_par->{'data'}->{'target'}->{'states_all_geom'} into a specified file
# and create gnuplot file for all states
sub save_target_energies {
  my ($filename, $r_par) = @_;
  my $igeom = 0;

  open(GNUPLOT, ">$r_par->{'dirs'}->{'model'}${bs}$filename.gp");
  print GNUPLOT "set title '".$r_par->{'model'}->{'molecule'}.", ".$r_par->{'model'}->{'model'}."'\n\n";
  print GNUPLOT "plot [:] \\\n";

  open(OUTPUT, ">$r_par->{'dirs'}->{'model'}${bs}$filename");
  print OUTPUT "# Target energies / a.u.\n";
  # Labels
  print OUTPUT '#'.$r_par->{'data'}->{'geom_labels'};
  my $n = 1;
  my $gnuplot_str = "";
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
      for (my $j = 0; $j < $r_par->{'model'}->{'ntarget_states'}->{$statespin}->[$i]; $j++) {
        print OUTPUT "  ".sprintf("%16s", "$statespin.$irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i].".($j+1));
        $n++;
        $gnuplot_str .= "  '$filename' u 1:$n t '$statespin.$irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i]".($j+1)."' w l, \\\n";
      }
    }
  }
  print OUTPUT "\n";
  $gnuplot_str =~ s/, \\\n$/\n/g;
  print GNUPLOT $gnuplot_str;
  close(GNUPLOT);

  # Data
  foreach $r_geom (@{$r_par->{'data'}->{'geometries'}}) {
    print OUTPUT ' '.$r_geom->{'geometry'};
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
        for (my $j = 0; $j < $r_par->{'model'}->{'ntarget_states'}->{$statespin}->[$i]; $j++) {
          my $state = "$spin_multiplicity{$statespin}.$i.".($j+1);
          print OUTPUT "  ".sprintf("%16.9f", $r_par->{'data'}->{'target'}->{'states_all_geom'}->{$state}->[$igeom]);
        }
      }
    }
    print OUTPUT "\n";
    $igeom++;
  }
  close(OUTPUT);
  return 1;
}

# Save R-matrix energies from $r_par->{'data'}->{'scattering'}->{'states_all_geom'} into a specified file
sub save_rmatrix_energies {
  my ($filename, $r_par, $number) = @_;
  my $igeom = 0;

  open(OUTPUT, ">$r_par->{'dirs'}->{'model'}${bs}$filename");
  print OUTPUT "# R-matrix energies / a.u.\n";
  # Labels
  print OUTPUT '#'.$r_par->{'data'}->{'geom_labels'};
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] > 0) {
        for (my $j = 0; $j < $number; $j++) {
          print OUTPUT "  ".sprintf("%16s", "$statespin.$irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i].".($j+1));
        }
      }
    }
  }
  print OUTPUT "\n";
  # Data
  foreach $r_geom (@{$r_par->{'data'}->{'geometries'}}) {
    print OUTPUT ' '.$r_geom->{'geometry'};
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
        if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] > 0) {
          for (my $j = 0; $j < $number; $j++) {
            my $state = "$spin_multiplicity{$statespin}.$i.".($j+1);
            print OUTPUT "  ".sprintf("%16.9f", $r_par->{'data'}->{'scattering'}->{'states_all_geom'}->{$state}->[$igeom]);
          }
        }
      }
    }
    print OUTPUT "\n";
    $igeom++;
  }
  close(OUTPUT);
  return 1;
}

# Save positions and widths of the resonance from $r_par->{'data'}->{'scattering'}->{'resonance_position'} 
# and $r_par->{'data'}->{'scattering'}->{'resonance_width'}
sub save_resonance_positions_and_widths {
  my ($filename, $r_par) = @_;
  my $igeom = 0;

  open(OUTPUT, ">$r_par->{'dirs'}->{'model'}${bs}$filename");
  print OUTPUT "# A = anion energy / a.u., E = position / a.u., W = width / a.u.\n";
  # Labels
  print OUTPUT '#'.$r_par->{'data'}->{'geom_labels'};
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] > 0) {
        my $state = "$statespin.$irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i]";
        print OUTPUT "  ".sprintf("%16s", "$state.A")."  ".sprintf("%16s", "$state.E")."  ".sprintf("%16s", "$state.W");
      }
    }
  }
  print OUTPUT "\n";
  # Data
  foreach $r_geom (@{$r_par->{'data'}->{'geometries'}}) {
    print OUTPUT ' '.$r_geom->{'geometry'};
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
        if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] > 0) {
          my $state = "$spin_multiplicity{$statespin}.$i";
          print OUTPUT "  ".sprintf("%16.9f", $r_par->{'data'}->{'target'}->{'states_all_geom'}->{$r_par->{'data'}->{'target'}->{'ground_state'}}->[$igeom] + 
                                              $r_par->{'data'}->{'scattering'}->{'resonance_position'}->{$state}->[$igeom]);
          print OUTPUT "  ".sprintf("%16.9f", $r_par->{'data'}->{'scattering'}->{'resonance_position'}->{$state}->[$igeom]);
          print OUTPUT "  ".sprintf("%16.9f", $r_par->{'data'}->{'scattering'}->{'resonance_width'}->{$state}->[$igeom]);
        }
      }
    }
    print OUTPUT "\n";
    $igeom++;
  }
  close(OUTPUT);
  return 1;
}

# Save density matrices to one file which can be plotted by gnuplot
sub save_density_matrices {
  my ($filename, $r_par) = @_; # $filename - a file where the extracted lines go; the name of file will be changed by adding number at the end
  my $title = "Calculating densities for density matrix:"; # the first line from which the text will be extracted
  my $end = "Integrated charge density for density matrix"; # the final line of extracting 
  my $file_input = "$r_par->{'dirs'}->{'outputs'}${bs}target.radden.out"; # the file from which the lines will be extracted
  my $directory = $filename;
  my $output_name = "";
  my $counter=1;
  my $found=0;

  system("mkdir $directory");
  open(INPUT, $file_input) or die "Can't open input file.\n";
  while (<INPUT>) {

    # Find block of lines to extract                                                           
    if( /^.*($title)/ ... /.*($end).*/ ) {

        # Start of block                                                                       
        if( /$title/ ) {
    	    $output_name=sprintf("${filename}.%d",$counter); # To extract each line-set into separate files.
            open(OUTPUT,'>'."$directory${bs}$output_name") or die $!;
        }
        # End of block                                                                         
        elsif ( /$end/ ) {
            close(OUTPUT);
            $counter++;
            $found = 0;
        }
        # Middle of block                                                                      
        else{
            if($found == 0) { # This is to comment first line
                print OUTPUT '#'.$_; 
                $found=1;
            }
            else {
                print OUTPUT $_;
            }
        }
    }
    # Find block of lines to extract
                                                           
  }
  close(INPUT);

  open(GNUPLOT, ">$directory${bs}$filename.gp");
  print GNUPLOT "#set terminal png enhanced  size  1600, 1200\n";
  print GNUPLOT "#set output '".$r_par->{'model'}->{'molecule'}."_radden.png'\n";
  print GNUPLOT "#set title '".$r_par->{'model'}->{'molecule'}.", ".$r_par->{'model'}->{'model'}."'\n";
  print GNUPLOT "#set xlabel 'Radius [a.u]'\n";
  print GNUPLOT "#set ylabel 'Density [a.u]'\n";
  print GNUPLOT "#set xtics 10,0.5,20\n";
  print GNUPLOT "plot for [i=1:".$r_par->{'model'}->{'req_dm'}."] '".$filename.".'.i using 1:3 smooth unique with lp title 'density matrix '.i\n";
  print GNUPLOT "#set output\n";
  print GNUPLOT "#set terminal x11\n";

  close(GNUPLOT);
  return 1;
}

sub gather_timedel {
  my ($prefix, $r_par) = @_;
  my @data = ();
  my @lines = ();
  my $row = 0;

  for (my $ir = 0; $ir < $r_par->{'data'}->{'nir'}; $ir++) {
    foreach my $statespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$ir] == 1) {
        my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$ir];
	open (IN, "< $prefix.$statespin.$str_ir") or die "Can't open $prefix.$statespin.$str_ir!";
        while (my $line = <IN>) {
	  foreach my $i (0, 1) {
	    $line =~ s/\s*(\d+\.\d+E?[\+\-]?\d+)//; 
	    $lines[$row][$i] = $1 ;
	  }
	$row++;
	}
	close IN or die "Cannot close $file: $!"; 
     	@lines = sort {$a->[0] <=> $b->[0] || $b->[1] <=> $a->[1] }@lines; 
	foreach my $i (0..$row-1) { # take only the higest value of timedel
	  if ( $i%3 == 0 ) {
	    $data[$i/3][0] = $lines[$i][0] * 13.605698066 ; # convert form Rydberg enery to eV
	    $data[$i/3][1] = $lines[$i][1];
	  }
	}
	open (OUT, '>', "$prefix.$statespin.$str_ir") or die "Can't write new file: $!";
	  print OUT "# ENERGY\tTIMEDEL\n ";
          map{print OUT "@$_\n"}@data;
	close OUT;
      } # if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$ir] == 1) {
    } # foreach my $statespin
  } # for (my $ir = 0; $ir < $r_par->{'data'}->{'nir'}; $ir++)
  return 1;
}

#=begin comment
sub change_prop {
  my ($prefix, $r_par) = @_;
  my $find_gs = 0;
  my $gs = 0;
  my $i = 1;

  open (IN, "< $prefix") or die "Can't open $prefix!";
  open (OUT, '>', "$prefix.1") or die "Can't overwrite $prefix!";
  while (my $line = <IN>) {
      chomp($line);
      if ($line =~ m/^5\s+.*(\-\d+\.\d+)D([\+\-]?\d+).*/) {
	if ($find_gs == 0) {
	  $gs = $1."E".$2;
	  print OUT $line."\n";
	  $find_gs = 1;
	}
	else {
	  my $ex = ($1."E".$2-$gs)*27.211;
	  
	  print OUT $line.sprintf("%.2d",$i++)." excited state: $ex eV\n";
	}
      }
      elsif ($line =~ m/^1\s+.*([\+\-]?\d+\.\d+)D([\+\-]?\d+).*/) {
	my $dm = sprintf("%.7f",($1."E".$2)/0.393456);
	print OUT $line. "\t$dm Debye\n";
      }
      else {
	print OUT $line."\n";
      }
  }
	close IN or die "Cannot close $prefix!"; 
	close OUT or die "Cannot close  overwriten $prefix!";

  return 1;
}

#=end comment
#=cut
1;


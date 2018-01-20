# Script to run UKRmol-in and -out calculations for polyatomic molecules
# by Karel Houfek (karel.houfek@gmail.com) from Charles University in Prague
use lib qw ( /padata/gamma/users/asieradzka/SCRIPT/ );
use Cwd;         # for changing the working directory
use ForkManager; # for parallelization
use ukrmolplib;  # my package of subroutines to handle geometries, inputs, outputs etc.
use Mail::Sendmail;

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $start = time;

my %personal = ();
my %geometry = ();
my %model = ();
my %run = ();
my %dirs = ();

foreach my $arg (@ARGV) { 
  eval &read_file($arg);
}

my %data = ( # this hash array will also be passed to subroutines calling programs
             # various subroutines which read outputs update data in this array

  # Information about molecular orbitals (read from scf output if not specfied otherwise):
  'orbitals',   {'energies',    {}, # energies of molecular orbitals as hash array, e.g. "2.1" -> -1.23456 is energy of the 1st state with symmetry 2
                 'target_all',  [], # from target sword:       number of all target MOs in each irreducible representation (IR)
                 'occupied',    [], # from target swscf:       number of occupied target MOs in each IR
                 'frozen',      [], # user specified:          number of frozen target MOs in each IR
                 'active',      [], # user specified:          number of active target MOs in each IR
                 'target',      [], # 'frozen' + 'active':     number of target MOs in each IR used to describe target states
                 'virtual',     [], # user specified:          number of virtual MOs in each IR, determined from energies
                 'target_used', [], # 'target' + 'virtual':    number of target MOs in each IR, which will be used in scattering calculation (IND in swedmos)
                 'all',         [], # from inner sword:        number of all basis functions in each IR
                 'continuum',   [], # 'all' - 'target_all'     number of continuum orbitals
                 'deleted',     [], # from inner swedmos:      number of deleted orbitals by swedmos in each IR
                 'cont_used',   [], # 'continuum' - 'deleted': number of continuum orbitals used for scattering after orthogonalization
                 'all_used',    [], # 'target_used' + 'cont_used'
                 'ifock',       []},# occupancy of orbitals in occ, 1 = full, 2 = half full (see swscf documentation, iocc and ifock)

  # Information about target states
  'target',     {'hf_energy',      0,  # Hartree-Fock energy of target
                 'hf_symmetry',    0,  # Hartree-Fock target state symmetry (needed when it is doublet)
                 'nfull_occ',      0,  # number of fully occupied orbitals in SCF calculation
                 'nhalf_occ',      0,  # number of half  occupied orbitals in SCF calculation
                 'noccupied',      0,  # total number of occupied orbitals in SCF calculation
                 'spin',           "", # current spin multiplicity of the target state
                 'symmetry',       0,  # current symmetry of the target state
                 'states',         {}, # energies of all target states as read from target.scatci output with keys of the form "spin.sym.n" such as "1.0.n" or "3.2.n" where n numbers states of the same spin-symmetry
                 'ordered_states', [], # array contains strings "spin.sym.n" (see above) describing target states ordered by energy
                 'ground_state',   "", # contains string "spin.sym.1" as above of the target ground state
                 'spinsym_order',  [], # order in which different target state spin-symmetries were calculated (each element contains string "spin.symmetry")
                 'used_tgt_states',{}, # number of target state for each spin-symmetry which are used in scattering calculations (e.g. "1.0" -> 2, "1.1" -> 1, ...)
                 # auxiliary variables to use for scatci and denprop inputs (generated automatically according to target states requirements)
                 'ntgt',           0,  # number of target spin-symmetries used in scattering calculations (ntgt -> denprop, ntgsym -> scattering.scatci)
                 'ntgtl',          [], # number of states used in scattering calculations in each symmetry (ntgtl -> denprop, numtgt -> scattering.scatci)
                 'mcont',          [], # symmetry of the continuum orbitals associated with each target spin-symmetry (see MCONT in scattering.scatci), but changed in make_scatci_input according to symmetry of the final state if necessary)
                 # auxiliary variables to use for outerres input
                 'idtarg',         [], # mapping of target states to energy order (idtarg for swinterf in outerres input)
                 # auxiliary arrays for energies etc for all geometries
                 'states_all_geom',{}, # arrays of energies of all target states for all geometries as read from target.scatci output with keys of the form "spin.sym.n" such as "1.0.n" or "3.2.n" where n numbers states of the same spin-symmetry
                                       # e.g. $data{'target'}->{'states_all_geom'}->{'1.0.1'} = [E_geom1, E_geom2, ...]
                 'pol_NOCSF',      [], # number of NOCSF in target congen output, necessary only in case when 'pol' model is chosen, i.e. to calculate polarizability.
                 'pol_reforb',     [], # that array keeps reference orbitals 'reforb' of ground state (for first irreducible representations), it is needed to create &wfngrp namelist  with Ground state^-1 x virtual^1 configurations for other that first irreducible representations (the way how the script does this should be improved in the future). It is necessary only in case when 'pol' model is chosen, i.e. to calculate polarizability. 
                },

  # Information about inner states
  'scattering', {'spin',     "", # current spin multiplicity of the scattering state
                 'symmetry', 0,  # current symmetry of the scattering state
                 'cont_csf', [0,0,0,0,0,0,0,0],  # number of N+1 configurations containing the continuum orbitals (LCONT)
                 'states_all_geom',{}, # arrays of energies of all R-matrix poles for all geometries as read from scattering.scatci output with keys of the form "spin.sym.n" such as "1.0.n" or "3.2.n" where n numbers states of the same spin-symmetry
                                       # e.g. $data{'scattering'}->{'states_all_geom'}->{'1.0.1'} = [E_geom1, E_geom2, ...]
                 'resonance_position',   {}, # arrays of positions of resonance for all geometries as read from scattering.outerres output with keys of the form "spin.sym"
                 'resonance_width',      {}, # arrays of widths    of resonance for all geometries as read from scattering.outerres output with keys of the form "spin.sym"
                },

  # determined automatically (whatever is set will be overwritten)
  'geometries',    [], # array of all geometries
  'geom_labels',   "", # will be labels for geometries used at the first line in files for potential curves etc.
  'task',          "", # will be the current regime, "target" or "scattering"
  'program',       "", # will be the current program to be executed
  'inputfile',     "", # will be the current input file
  'outputfile',    "", # will be the current output file
  'logfile',	   "", # will be the current standard output file
  'igeom',         0,  # index of current geometry to work with
  'geom',          0,  # will be set in the geometry loop as a reference to a given geometry
  'nir',           0,  # number of irreducible representation
  'scf_ok',        0,  # set to 1, if iocc is consistent with result of SCF calculations
  'no_scatci',     0,  # set to 1, if there is no configuration generated during congen run -> skipping scatci

);

# ================== beginning of processing =================

# ----------------------- preliminaries ----------------------

# collection of hash arrays for simple passing of all variables to subroutines
my %parameters = ('personal',   \%personal,
		  'model', 	\%model, 
                  'run',   	\%run, 
                  'dirs',  	\%dirs, 
                  'data',  	\%data);

# Set number of irreducible representations - used as abbreviation
$data{'nir'} = scalar @{$irred_repr{$model{'symmetry'}}};

# Add current working directory to directories specified by user
# (necessary because we change the working directory later)
&add_cwd_to_dirs(\%dirs);

# Set and make the model directory
my $directory = "$model{'molecule'}${bs}$model{'basis'}.$model{'model'}.$model{'nfrozen'}frozen.$model{'nactive'}active.$model{'nvirtual'}virtual.$model{'ntarget_states_used'}states.r$model{'radius'}.L$model{'L_basis'}.$model{'cont_suffix'}.$run{'precision'}$run{'suffix'}";
if ($dirs{'model'} eq "") {
  # example: output/CO/6-311G_dp.SE.15virtual.2frozen.r10
  $dirs{'model'} = "$dirs{'output'}${bs}$directory";
}
else {
  $dirs{'model'} = "$dirs{'output'}${bs}$dirs{'model'}";
}
if ($run{'bound'} == 1) {
  $dirs{'model'} .= '.bound';
}
&make_dir($dirs{'model'}); 

if ($run{'backup'} == 1){
  $dirs{'backup'} = "$personal{'backup_dir'}${bs}$directory";
  &make_dir($dirs{'backup'}); 
}

if ($dirs{'std_out'} eq "") {
  $dirs{'std_out'} = "$dirs{'output'}${bs}$directory${bs}std_out";
}
else {
  $dirs{'std_out'} = "$dirs{'output'}${bs}$dirs{'std_out'}";
}
$data{'logfile'} = $dirs{'std_out'};

&print_info("Current working directory: $dirs{'cwd'}\n\n", \%parameters);
&print_info("Running UK R-matrix codes for e + $model{'molecule'}\n\n", \%parameters);
&print_info("Output in $dirs{'model'}\n\n", \%parameters);
&print_info("Symmetry $model{'symmetry'} with $data{'nir'} IRs for all geometries.\n\n", \%parameters);

# Next we get all geometries using the subroutine generate_geomatries
# which also creates the file 'geometries' and directories for each geometry
# For each geometry it returns hash array { 'geometry', "string of distances and angles as in the file geometries",
#                                           'dir', directory,
#                                           'symmetry', sym,
#                                           'natype', value,
#                                           'atoms', [ ['A',x,y,z for atom 1], ..., ['B',x,y,z for last atom] ] }
# where 'sym' is one of (D2h, C2v, C2, Cs, C2h, D2, Ci) and 'A' or 'B' are symbols of atoms

push( @{$data{'geometries'}}, &generate_geometries(\%parameters, \%geometry));

# ================ main loop over geometries =================

foreach $r_geom (@{$data{'geometries'}}) {
  $data{'igeom'} += 1;
  if ($data{'igeom'} >= $geometry{'start_at_geometry'}) {
    $data{'geom'} = $r_geom;

    # Directories for inputs and outputs
    $dirs{'geom'} = $r_geom->{'dir'};
    chdir($dirs{'geom'});
    $dirs{'inputs'}  = "inputs";      # Directory where all inputs will be stored  - e.g. as target.swmol3.inp
    $dirs{'outputs'} = "outputs";     # Directory where all outputs will be stored - e.g. as target.swmol3.out
    &make_dir($dirs{'inputs'});
    &make_dir($dirs{'outputs'});

    &print_info("\nGeometry \#$data{'igeom'}:\n===========\n", \%parameters);

    # =================== start of target ======================
      $data{'task'} = "target";
      &print_info("\nTarget calculations:\n\n", \%parameters);
      if ($run{'ukrmolplus'} == 0) {
        &run_code("swmol3", \%parameters);
        &run_code("sword", \%parameters);
        # if we have molpro then we use it for everything - including HF 
        if ($run{'molpro'} == 1) {
          &print_info("Using target orbitals generated by molpro.\n", \%parameters);
          &run_code("molpro", \%parameters);
          &run_code("mpoutrd", \%parameters);         
          if ($run{'skip_radden'} != 1) {
            &run_code("radden", \%parameters);
            &save_density_matrices("density.matrices", \%parameters);
          }
        } else {
           &print_info("Using HF target orbitals generated by swscf program.\n", \%parameters);
           &run_code("swfjk", \%parameters);
           # swscf is running several times if necessary to get correct orbital occupation
           &initial_orbitals_occupation_guess(\%parameters);
           $data{"scf_ok"} = 0;
           while ($data{"scf_ok"} == 0) {
             &run_code("swscf", \%parameters);
           }
        }
        &run_code("swedmos", \%parameters);
        &run_system("$cp_cmd fort.21 tgt.orbitals", \%parameters);
        &run_system("$rm_cmd fort.23", \%parameters);
        &run_system("$mv_cmd fort.21 fort.23", \%parameters);
        &run_code("swtrmo", \%parameters);
     } else {
        &run_code("molpro", \%parameters);
        if ($run{'skip_radden'} != 1) {
          &run_code("radden", \%parameters);
          &save_density_matrices("density.matrices", \%parameters);
        }
        if ($run{'computer'} == 2) {
       		&run_code("integrals", \%parameters);
       		&run_system("ssh $personal{'archer_login'} \'mkdir -p $personal{'archer_work'}/$directory/geom1/outputs\'", \%parameters);
       		&run_system("scp $dirs{'cwd'}${bs}$personal{'archer_script'} $dirs{'geom'}${bs}*molden $dirs{'geom'}${bs}inp $personal{'archer_login'}:$personal{'archer_work'}/$directory/geom1/.", \%parameters);
       		&run_system("scp $dirs{'geom'}${bs}outputs/target.molpro.out  $personal{'archer_login'}:$personal{'archer_work'}/$directory/geom1/outputs/.", \%parameters);
       		&run_system("ssh $personal{'archer_login'} \'$cd_cmd $personal{'archer_work'}/$directory/geom1; qsub $personal{'archer_script'}\'", \%parameters);
		&run_system(die "If you can see job number above then Integrals program is submitted on Archer! If not, check if submission script on Archer exist.\n");
        } else {
       		&run_code("integrals", \%parameters);
        }
	  &run_system("$mv_cmd molecular_integrals molecular_integrals_target", \%parameters);
	  foreach my $u (16, 17, 22) { &make_symlink("molecular_integrals_target", "fort.$u"); } 
        }

    # from here we run congen and scatci for each spin state (multiplicity) and symmetry (IR) separately
      for (my $i = 0; $i < $data{'nir'}; $i++) {
        foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'ntarget_states'}}) {
          $data{'target'}->{'spin'} = $statespin;
          if ($model{'ntarget_states'}->{$statespin}->[$i] > 0) { # skip if target states of a given spinsymmetry are not required
            $data{'target'}->{'symmetry'} = $i;
            &run_code("congen", \%parameters);
            if ($data{'no_scatci'} == 1) {
              $data{'no_scatci'} = 0;
            }
            else {
              &run_code("scatci", \%parameters);
            }
          }
        } # foreach spin state (multiplicity)
      } # for each IRs
      if ($run{'ukrmolplus'} == 0) {&run_code("gausprop", \%parameters);}
      &run_code("denprop", \%parameters);
      &run_system("$cp_cmd fort.24 prop.out", \%parameters);

    # ===================== end of target =======================

    # ================== start of scattering ====================
    if ($run{'scattering'} == 1) {
      $data{'task'} = "scattering";
      &print_info("\nScattering calculations with $model{'model'} model:\n\n", \%parameters);
      if ($run{'ukrmolplus'} == 0) {
      &run_code("swmol3", \%parameters);
      if ($run{'bound'} == 0) {
        &run_code("gaustail", \%parameters);
      }
      &run_code("sword", \%parameters);
      &run_system("$cp_cmd tgt.orbitals fort.23", \%parameters);
      &run_code("swedmos", \%parameters);
      &run_system("$rm_cmd fort.23", \%parameters);
      &run_system("$mv_cmd fort.21 fort.23", \%parameters);
      &run_code("swtrmo", \%parameters);
      } 
#---------------- parallel invoking all symmetries -------------------#
    # from here we run congen and scatci for each spin state (multiplicity) and symmetry (IR) separately
      if ($run{'parallel_symm'} > 1) { $psymm = Parallel::ForkManager->new($run{'parallel_symm'}); }
      for (my $i = 0; $i < $data{'nir'}; $i++) {
        foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'scattering_states'}}) {
          if ($run{'parallel_symm'} > 1) { $psymm->start and next; }
          $data{'scattering'}->{'spin'} = $statespin;
        # skip spin-symmetries which are not required
          if ($model{'scattering_states'}->{$statespin}->[$i] == 1) {
            $data{'scattering'}->{'symmetry'} = $i;
            my $spin_sym = "$statespin.$irred_repr{$model{'symmetry'}}->[$i]";
          # each spin-symmetry calculation is done in a local directory
            &make_dir($spin_sym);
            chdir($spin_sym);
            foreach my $u (16, 22, 23, 24) { &make_symlink("..${bs}fort.$u", "fort.$u"); }
            &run_code("congen", \%parameters);
            if ($data{'no_scatci'} == 1) { # skip scatci if there is no configuration generated by congen
              $data{'no_scatci'} = 0;      # reset the flag
            }
            else {
              &run_code("scatci", \%parameters);
              if ($run{'parallel_diag'} == 1) { &run_code("hamdiag", \%parameters); }
              if ($run{'bound'} == 0) { # if bound states are required we do not run outer region codes
                &run_code("outerres", \%parameters);
                if ($run{'save_eigenph'}) { &run_system("$cp_cmd fort.110 ..${bs}eigenph.$spin_sym", \%parameters); }
                if ($run{'save_xsec'})    { &run_system("$cp_cmd fort.100 ..${bs}xsec.$spin_sym", \%parameters); }
                if ($run{'save_Kmatrix'}) { &run_system("$cp_cmd fort.19  ..${bs}K-matrix.$spin_sym", \%parameters); }
                if ($run{'save_Tmatrix'}) { &run_system("$cp_cmd fort.12  ..${bs}T-matrix.$spin_sym", \%parameters); }
                if ($run{'save_outer'})   { &run_system("$cp_cmd fort.10  unit.10", \%parameters);
					    &run_system("$cp_cmd fort.21  unit.21", \%parameters);}
                if ($run{'time-delay'}) { &run_code("time-delay", \%parameters); }
		if ($run{'timedel'}) {
		  &run_code("timedel", \%parameters);
                  &run_system("$cp_cmd fort.9999  ..${bs}timedel.$spin_sym", \%parameters);
		  foreach my $u ('branching_ratios', 'ground_state_widths', 'resonances', 'resonance_X_ratios', 'reson_message') {
		    &run_system("$mv_cmd $u ..${bs}$u.$spin_sym", \%parameters);
		  }
		}
              } # bound
            } # no_scatci
            if ($run{'clean'}) { system("$rm_cmd fort.*"); }
            chdir($dirs{'geom'});
            if ($run{'clean'}) { rmdir($spin_sym); }
          } # run given spin-symmetry
          if ($run{'parallel_symm'} > 1) { $psymm->finish; }
        } # foreach spin state (multiplicity)
      } # for each IRs
      if ($run{'parallel_symm'} > 1) { $psymm->wait_all_children; }
    }
    # =================== end of scattering =====================

    # -------------------- gathering data -----------------------
    if ($run{'gather_data'} == 1) {
      if ($run{'bound'} == 0) {
        &gather_eigenphases("eigenph", \%parameters);
        &convert_and_gather_cross_sections("xsec", \%parameters);
      }
    }
    if ($run{'timedel'}) { &gather_timedel("timedel", \%parameters); }
    if ($model{'model'} eq "CAS") { &change_prop("prop.out", \%parameters); }

    # ------------- cleaning before next geometry ---------------
    foreach my $key (%{$data{'orbitals'}}) {
      $data{'orbitals'}->{$key} = [];
    }
    $data{'orbitals'}->{'energies'} = {};
    $data{'target'}->{'states'}          = {};
    $data{'target'}->{'ordered_states'}  = [];
    $data{'target'}->{'spinsym_order'}   = [];
    $data{'target'}->{'used_tgt_states'} = {};
    $data{'target'}->{'ntgt'}            = 0;
    $data{'target'}->{'ntgtl'}           = [];
    $data{'target'}->{'mcont'}           = [];
    $data{'target'}->{'idtarg'}          = [];

    if ($run{'clean'} == 1) {
        system("$rm_cmd fort.*");
    }
    chdir($dirs{'cwd'});
  } # if ($data{'igeom'} >= $geometry{'start_at_geometry'})  
} # foreach $r_geom

if ($run{'gather_data'} == 1) {
  &save_target_energies("target.energies", \%parameters);
  &save_rmatrix_energies("Rmatrix.energies", \%parameters, 5);
  if ($run{'bound'} == 0) {
    &save_resonance_positions_and_widths("resonance.positions.and.widths", \%parameters);
  }
}

  # ------------- saving in backup directory ---------------
if ($run{'backup'} == 1){
  foreach my $u ('inputs', 'outputs', 'prop.out') { 
    &run_system("$cpdir_cmd $dirs{'geom'}${bs}$u $dirs{'backup'}${bs}.", \%parameters);
  }
  if ($run{'scattering'} == 1) { 
    foreach my $u ('*matrix*', 'xsec*', 'eigenph*') {
      &run_system("$cp_cmd $dirs{'geom'}${bs}$u $dirs{'backup'}${bs}.", \%parameters);
    }
    foreach my $u ('*energies', 'resonance*') {
      &run_system("$cp_cmd $dirs{'model'}${bs}$u $dirs{'backup'}${bs}.", \%parameters);
    }
  }
}

  # ------------- time --------------
my $start_time = sprintf("%02d.%02d.%02d   %02d:%02d:%02d", $mday, $mon + 1, $year + 1900, $hour, $min, $sec);
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $end_time = sprintf("%02d.%02d.%02d   %02d:%02d:%02d", $mday, $mon + 1, $year + 1900, $hour, $min, $sec);
my $duration = time - $start;
&print_info("Job started run at:\t$start_time\n...and stopped at:\t$end_time\n", \%parameters);
&print_info("Execution time: $duration s\n", \%parameters);

  # ------------- sending e-mail --------------
sendmail(
    From    => $personal{'ad_email_from'},
    To      => $personal{'ad_email_to'},
    Subject => 'Your job is done',
    Message => "The calculations are in a directory: $dirs{'model'}.\nThe job started at:\t$start_time.\nFinished at:\t$end_time.",
);

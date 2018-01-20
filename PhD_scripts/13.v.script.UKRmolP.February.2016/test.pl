# ====================== user settings ======================

# --- personal settings ---
%personal = (

  'ad_email_from', 'agnieszka.sieradzka@open.ac.uk', 		# an e-mail address can be given to send information that calculations ended; 
								# from this e-mail address the message will be send;
  'ad_email_to',   'agnieszka.sieradzka@open.ac.uk', 		# the message will be send to this e-mail address.
  'backup_dir',    "/padata/gamma/users/asieradzka/test",  	# path to back up directory where copies of data will be stored
  'archer_login',  "aga\@login.archer.ac.uk",  			# the login to archer machine to run there everything except molpro
  'archer_work',   "~/work",  					# the work directory
  'archer_script', "archer.run.pbs",  			# name of a script to submit on archer; to submit into short queue: "-q short name_of_script"
  
);
# --- geometry ---

# Geometries can be specified either manually in the array 'geometries'
# by copying what is between start copy and end copy for each geometry (first uncomment this block)
# or automatically as shown right below the hash array %geometry

%geometry = (

  'geometry_labels',  "     H-H",         # labels used on the first line of output files
                                          # it should correspond to numbers given in 'geometries'->'description'
  'correct_cm',   0,                      # correct the center of mass to be at the origin

  'geometries', [
#       # start copy
#       { 'description', "    2.40",      # string to use in output files to describe this particular geometry, can be anything
#         'natype', 1,                    # number of atoms needed for QChem codes
#                                         # in 'atoms' below redundant atoms should be specified at the end
#         # !!! specify ALL atoms (even redundant with respect to symmetry elements) 
#         # if they should be moved to have the center of mass of the molecule at the origin
#         'atoms', [ [ "H", 0.0,  0.0, 0.0 ],
#                    [ "H", 0.0,  0.0, 2.4 ] ]
#       },
#       # end copy
  ],

  'start_at_geometry',  0,                # this can be used to continue an interupted run, all geometries are generated but R-matrix codes run only for geometries with index >= than this number
);

# here is an example of automatic generation of geometries
for (my $rrr = 2.4; $rrr <= 2.4; $rrr += 0.1) {
  push(@{$geometry{'geometries'}},
       { 'description', sprintf("%8.5f", $rrr),
         'natype', 2,
         'atoms', [ [ "O", 0.0,			 1.3158162058,                       -0.1043128746 ],
                    [ "H", 1.4771988053,	 1.6990526363,          	      0.8345029966 ],
                    [ "O", 0.0,			-1.3158162058,                       -0.1043128746 ],
                    [ "H", -1.4771988053,	-1.6990526363, 		    	      0.8345029966 ] ]
       }
  );
}

# --- other settings for inputs ---

# this whole hash array will be passed to subroutines calling programs
# so all subroutines can use any information they need
%model = ( 

  # Molecule
  'molecule',   "test",		# Used only for directory and descriptions
		#"hydrogen.peroxide", 
       
  'atoms',        ["O", "H", "O", "H"],	# All atoms must be specified
  'nelectrons',    18,                	# number of target electrons
  'symmetry',      "C2",             	# Options are "D2h", "C2v", "C2", "Cs", "C2h", "D2", "Ci", "C1"
                                      	# Don't forget the restriction that if there is just one C2 axis, then it is Z-axis
                                      	# and for Cs symmetry the reflection plane is XY plane (for C2v planes are XZ and YZ !)

  # Units
  'r_unit',        0,                 	# distance unit: 0 - atomic units, 1 - Angstroms
  'e_unit',        2,                 	# energy unit:   0 - atomic units, 1 - Rydbergs, 2 - eV
  'x_unit',        1,                 	# cross section: 1 - atomic units, 2 - Angstroms^2

  # Model - each model will have its own directory from these settings in the directory of the molecule
  # Important note: Put correct numbers of orbitals even for SE and SEP

  'basis',       #"6-311ppGxx",      	# make sure there are corresponding files in $dir_basis (see below)
		 "cc-pVDZ",   
   
  'model',       "CAS",             	# options are "SE"    = static-exchange, 
                 #"SE",             	#             "SEP"   = SE + polarization, 
                 #"pol",            	#             "CAS"   = complete active space (as a default it is used the model B below), 
                                    	#             "CAS-A" = contracted version of CAS-B, 
                                    	#             "CAS-B" = standard close coupling with (core+cas)^N+1 and (core+cas)^N x (virtual)^1
                                    	#             "CAS-C" = adds (core+cas)^N-1 x (virtual)^2 to CAS-B (to add more polarization)
                                    	#             "pol"   = polarizability; It calculates polarizability of target. Choose number of states in 'ntarget_states'
					#		 	and change 'ntarget_states_used', switch off 'run_scattering', i.e. choose '0'.
  'nfrozen',       "4",            	# number of frozen target orbitals (used for SEP and CAS)      
  'nactive',       "8",             	# number of active target orbitals (used for SEP and CAS)
  'nvirtual',      "0",             	# number of virtual states

  # if the following arrays are empty (zeroes) then orbitals are chosen automatically according to their energies
  # otherwise the script uses chosen orbitals but ONLY IF the total number of chosen orbitals in arrays is consistent
  # with numbers of orbitals specified above in 'nfrozen', 'nactive', 'nvirtual'
  # be careful with these settings, if you choose your active space differently than orbitals are ordered
  # you should specify also which virtual orbitals you want to use
  'frozen_orbs',   [2,2,0,0,0,0,0,0], 	# which orbitals for each symmetry to use as frozen, (also used in MOLPRO model CAS)
  'active_orbs',   [4,4,0,0,0,0,0,0], 	# which orbitals for each symmetry to use as active, (also used in MOLPRO model CAS /MOLPRO input: occ = frozen_orbs + active_orbs/)
#  'virtual_orbs',  [0,0,0,0,0,0,0,0], 	# which orbitals for each symmetry to use as virtual,
#  'frozen_orbs',   [0,0,0,0,0,0,0,0], 	# which orbitals for each symmetry to use as frozen, (also used in MOLPRO model CAS)
#  'active_orbs',   [0,0,0,0,0,0,0,0], 	# which orbitals for each symmetry to use as active, (also used in MOLPRO model CAS /MOLPRO input: occ = frozen_orbs + active_orbs/) 
  'virtual_orbs',  [0,0,0,0,0,0,0,0], 	# which orbitals for each symmetry to use as virtual,

  # Input for MOLPRO
  'molpro_setting',    "", 
  'ncasscf_states',  {'singlet', [2,1,0,0,0,0,0,0],  	# number of states used to optimalised orbitals in MOLPRO
		      'triplet', [1,1,0,0,0,0,0,0]}, 	# Warning! It works only for 'singlet', 'triplet'. 
							# If you want to change it then go to ukrmollib - input molpro and change '0' '1'
						     	# In fact, 'singlet', 'triplet' does not mean here multiplicity but total spin '0' for 'singlet' and '1' for 'triplet'

  'delthres',      "1.0D-07, 1.0D-07, 1.0D-07, 1.0D-07",# deletion threshold in swedmos                                    

  # Target states: for closed-shell target, specify the number of singlets, triplets, ...
  #                for   open-shell target,                       doublets, quartets, ...   
  'ntarget_states', {'singlet', [5,4,0,0,0,0,0,0],  	# number of target states to calculate in each irreducible representation (IR) for CAS
                     'triplet', [4,4,0,0,0,0,0,0]},
#  'ntarget_states', {'singlet', [1,0,0,0,0,0,0,0],  	# number of target states to calculate in each irreducible representation (IR) for SEP
#                     'triplet', [0,0,0,0,0,0,0,0]},
  'ntarget_states_used', "17",         			# number of target states which will be actually used in scattering calculations 
                                      			# (chosen according to their energy from states above)

  # Scattering settings

  # Scattering spin-symmetry: for closed-shell target, specify which of doublets, quartets, ... to use
  #                           for   open-shell target,                  singlets, triplets, ...   
  'scattering_states', {'doublet', [1,1,0,0,0,0,0,0],  	# 1/0 to run/not to run scattering calculation for a given spin-symmetry
			'quartet', [0,0,0,0,0,0,0,0]},
  'L_basis',	   "4",		      	# The maximum of partial waves included into continuum
  'cont_suffix',   "old",	      	# continuum suffix 'old' or 'new, where 'old' indicates standard continuum and 'new' indicated the expanded
  'radius',        "15",              	# R-matrix radius
  'max_multipole', "2",               	# maximum multipole to be retained in expansion of long range potentials
  'raf',           "50.1",            	# radius at which continued fraction method can be used for R-matrix propagation

  'nescat',        "600", 		# number of input scattering energies in each subrange (input for R_SOLVE via namelist &rslvin)
  'einc',          "0.01, 0.025", 	# scattering energies - initial energy, energy increment, there can be more subranges

  'maxi',          "1",               	# the highest initial state for which cross sections are required
  'maxf',          "0",               	# the highest final state for which cross sections are required (zero means all)

  # Other settings specific for R-matrix codes
  'lndo',          150000000,          	# memory control in CONGEN 

  # Input for TIMEDEL
  'einit',	   "0.1",          	# initial energy for timedel, if no value (empty string) is given the initial energy is the same as for R_SOLVE
  'efinal',	   "",          	# final energy for timedel, if no value (empty string) is given the initial energy is the same as for R_SOLVE

  # Input for RADDEN
  'radius_start',  "5",              	# radial interval in a.u. for which the integration will be carried out
  'radius_finish', "20",              	# the last value of radial for which the integration will be calculated
  'req_dm',	   "30",          	# if < 0 (e.g. -1) to all available orbitals/density matrices will be read in and integrated, otherwise the number of orbitals/density matrices to be read in from the density matrix file
  'dm_list',	   "5",          	# first density matrices/orbitals to be read in from the input file. Not used if ilist .eq. 'false' (for non-sequential orbtials, modify input manually)

);


# Settings related to a specific installation of UK R-matrix codes
#   and running and saving options

%run = (

  'suffix',          "_prec.test.old",  		# string added to model directory to distinguish different runs which have the same model and number of orbitals and states
                               		# it can be useful e.g. for different range of geometries or different symmetry runs etc.
  'print_info',      "screen", 		# "screen" or "file" - option whether info messages will be printed on the screen 
					# or into the file of the same name as model directory with the suffix .txt
  'precision',       "double", 		# "double" or "quadru" - option whether run executables with double or quadruple precision. Leave empty if you do not need

  # Running options
  'ukrmolplus',      1,  		# 1 - if you want run UKRmol+, 0 - if UKRmol
  'scattering',      1,  		# run all programs, if you want to run target only, set to 0
  'molpro',          1,  		# calculate natural orbitals using CASSCF with MOLPRO
  'molpro_basis',    1,  		# set 0 if you want MOLPRO to take basis set from files, otherwise it takes basis set according to given basis set name 
  'time-delay',	     1,			# set 1 to run time-delay
  'timedel',	     0,			# set 1 to run timedel 
  'skip_radden',     1,  		# skip calculating radial densities
  'gather_data',     1,  		# gather eigenphase sums, cross sections, energies ...
  'clean',           0,  		# removing fort.* etc
  'use_templates',   1,  		# set this to 0, if you need to modify generated inputs manually and than rerun everything
                        		# but do not change the filenames !

  # Saving options
  'save_eigenph',    1,  		# set 1 to copy fort.110 into a file like eigenph.singlet.Ag (eigenphase sums)
  'save_xsec',       1,  		# set 1 to copy fort.100 into a file like xsec.singlet.Ag (cross sections)  
  'save_Kmatrix',    1,  		# set 1 to copy fort.19  into a file like K-matrix.singlet.Ag
  'save_Tmatrix',    1,  		# set 1 to copy fort.12  into a file like T-matrix.singlet.Ag
  'save_outer',      1,  		# set 1 to copy fort.10 and fort.21  into a file unit.10 and unit.21 respectively
  'backup',	     0,  		# save some important files in backup directory given in %dirs -> backup
                         		# probably usefull only at some computers

  # Parallelization - you can specify a number of processes to run in parallel
  #                   separately for geometries and for symmetries (ireducible repr.)
  # 'parallel_geom' x 'parallel_sym' should be <= number of available CPUs
  # ForkManager package is used for parallelization within perl scripts
#  'parallel_geom',   1,  		# number of geometries to be invoked in the same time
  'parallel_symm',   2,  		# number of symmetries to be invoked in the same time
  'parallel_diag',   0,  		# 1 for run hamdiag - parallel diagonalization of Hamiltionian matrix

  # Using of parallelized codes
  # Probably not reasonable to use with the parallelization above
#  'parallel',	     0,  		# 1 for run hamdiag - parallel diagonalization of Hamiltionian matrix

  # For very special run to get bound states of e + target without running scattering
  #   gaustail is skipped, thus integrals are not limited to the R-matrix sphere and scatci calculates regular eigenstates, not R-matrix poles
  # set this option to 1 if you want to calculate bound states
  'bound',           0, 

  # For 'debugging' only (or maybe it can be useful to rerun outer region only, but don't forget to turn off cleaning)
  # if run_only is empty, then all programs will be executed
  # otherwise it is assumed that outputs from the previous run exist
  #           and only specified programs will be executed
  'only',            "scattering-outerres|scattering-time-delay", #"scattering-congen|scattering-scatci|scattering-outerres|scattering-time-delay", # target-congen|target-scatci|target-denprop| "scattering-outerres", #"target-denprop|scattering-congen|scattering-scatci|scattering-outerres|scattering-time-delay", # scattering-hamdiag| "target-mpoutrd|target-swedmos|target-swtrmo|target-congen|target-scatci|target-gausprop|target-denprop|scattering-swmol3|scattering-gaustail|scattering-sword|scattering-swedmos|scattering-swtrmo|scattering-congen|scattering-scatci|scattering-hamdiag|scattering-outerres"      # e.g. "scattering-congen|scattering-scatci|scattering-outerres" to run only last part of the scastering calculation

  'no_processors',   $ARGV[1],  	# number of processors; it only applies for integrals or hamdiag programs
  'computer', 	     1,  	 	# 0 - uqbar, 1 - impact cluster, 2 - archer

);

# Path to executables can be spicified directly in %dirs or later if you are using several computers
# (see switch($run{'computer'}) below
# Use ${bs} (see dirfile.pm) in relative paths instead of \ or / for portability
# Some directories (indicated by "") are determined later automatically
# If the full path must be used then it is not necessary to use ${bs}

%dirs = (

  'bin_in',    "/padata/gamma/users/asieradzka/bin_UKRmol_$run{'precision'}",
  'bin_out',   "/padata/gamma/users/asieradzka/bin_UKRmol_$run{'precision'}",
  'integrals', "/padata/gamma/users/asieradzka/bin_UKRmol_$run{'precision'}",
  'molpro',    "/padata/gamma/groups/rmatrix/Molpro.2010.1.26/bin/",
  'basis',     "/padata/gamma/users/asieradzka/SCRIPT/basis.sets",         # Directory where basis sets are - templates named 'swmol3.A.$basis' 
									   # or 'molpro.A.$basis' where 'A' stands for an atom
  'templates', "/padata/gamma/users/asieradzka/SCRIPT/input.templates",    # Directory where input templates are - read if 'use_templates' = 1
  'output',    ".",              	# Main directory for output
  'std_out',   "",              	# standard output, if it is empty, it will be set later according to the model used and it will have name: std_out
  'model',     "", 			# if it is empty, it will be set later according to the model used
  'backup',    "", 			# /specify it in personal settings/ additional parameter to save files on backup directory. If it is empty, it will be set later according to the model used
  'geom',      "", 			# will be set in the geometry loop
  'inputs',    "", 			# will be set in the geometry loop
  'outputs',   "", 			# will be set in the geometry loop
  'cwd',       getcwd() 		# current working directory
);

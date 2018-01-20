my $No = 0; # number of orbitals
my $i = 0;
#my $file = "../radden.6-311pgxx.pyridine_water.out"; 
my $file = "radden.cc-pVDZ.pyridine_water.out";
#my $directory = "diffuse";
my $directory = "compact2";
system("mkdir $directory");
  if (open(MYFILE, $file)) {
        while ($line = <MYFILE>) {
		if ($line =~ m/\s*Required number of orbitals:\s*(\d+)/) {
			$No = $1; 
			print "$No \n"; 
		}		
  	}
  	close(MYFILE);
  }
  else {
  die "  Error: can't open RADDEN output !\n";
  }
  $com = "\#";
  for($i = 1; $i <= 9; $i ++){
	system("sed -n -s \"/Calculating densities for density matrix:     $i/,/Integrated charge density for density matrix     $i/p\" $file >$directory/density.matrix.$i");
	system("sed -i \"s/Calculating/$com Calculating/\" \"$directory/density.matrix.$i\"");
	system("sed -i \"s/Radius/$com Radius/\" \"$directory/density.matrix.$i\"");
	system("sed -i \"s/Integrated/$com Integrated/\" \"$directory/density.matrix.$i\"");
  }
  for($i = 10; $i <= $No; $i ++){
	system("sed -n -s \"/Calculating densities for density matrix:    $i/,/Integrated charge density for density matrix    $i/p\" $file >$directory/density.matrix.$i");
	system("sed -i \"s/Calculating/$com Calculating/\" \"$directory/density.matrix.$i\"");
	system("sed -i \"s/Radius/$com Radius/\" \"$directory/density.matrix.$i\"");
	system("sed -i \"s/Integrated/$com Integrated/\" \"$directory/density.matrix.$i\"");
  }


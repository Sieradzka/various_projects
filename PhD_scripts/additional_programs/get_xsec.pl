#!/usr/bin/perl

 my $input = $ARGV[0]; #file containing paths to the cross section files in the R-matrix output format
 my $emin = $ARGV[1]; #energy range start [ev]
 my $emax = $ARGV[2]; #energy range end [ev]

 open(inp,  "$input")   or die "Can’t open $input: $!";

 @files = <inp>; #read in the list of files
 close inp;

 $i = 0;
 foreach(@files) {
	$i = $i + 1;
	$file = $_;
 	open(ff, "$file") or die "Can't open $file: $!"; #open the first cross section file

        #position the file on the line continaing the first cross section value
	$line = <ff>;
        until ($line =~ /^\s*1\s*[0-9]+\.[0-9]+[eED]?[+-]?[0-9]+.*/) { #the first energy
		$line = <ff>; #skip the header stuff
	}

	$e = 0;
        while ($line =~ /^\s*[0-9]+\s*([0-9]+\.[0-9]+[eED]?[+-]?[0-9]+)\s*([0-9]+\.[0-9]+[eED]?[+-]?[0-9]+).*/) {
		$n = $1;
		$val = $2;
                $n =~ s/D/E/g;
		if ($i == 1 && $n >= $emin && $n <= $emax) {
			push @en, "$n";
		}
		if ($n >= $emin && $n <= $emax) {
			$e = $e + 1;
			$cs[$e] = $cs[$e] + $val;
			#print "$n $cs[$e]\n";
		}
		$line = <ff>;
	}

	close f;
    
 }

 $i = 0;
 foreach(@en) {
	$i = $i + 1;
	print "$en[$i-1] $cs[$i]\n";
 }

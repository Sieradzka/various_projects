#!/usr/bin/env perl
use strict;
use warnings;
#----------------------------------------------------------------------------------------------------------------------
my $file = "$ARGV[0]";
my $i = 0; 
#my $line = "";

open my $in,  '<',  $file      or die "Can't read old file: $!";
open my $out, '>', "$file.new" or die "Can't write new file: $!";

while( <$in> )
    {
      if ($_ =~ m/\s*Ene=.*/) {
       $i += 1;
#       print "*******************************\n $_  $i\n";
      }
    s/ Ene=/ Sym=     $i.1 \n Ene=/g;
    print $out $_;
    }
close $out;
close $in;

=begin comment
if (open($in, $file)) {
    while ($line = <$in>) {
      if ($line =~ m/\s*Ene=.*/) {
       $i += 1;
#       print "*******************************\n $line  $i\n";
      }
    $line =~ s/ Ene=/ Sym=     $i.1 \n Ene=/g;
    print $out $line;
    }
    close($in);
  }
  else {
      die "  Error: can't open $file file !\n";
  }
=end comment
=cut


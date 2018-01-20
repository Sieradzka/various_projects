### System specific variables:

$sys = "linux";
if ($sys eq "win") {  # Windows 
  $bs = "\\";
  $rebs = "\\\\";
  $rm_cmd = "del";
  $cp_cmd = "copy";
  $mv_cmd = "ren";
  $cd_cmd = "cd";
  $ext_exe = ".exe";   
  $cpdir_cmd = "xcopy /s /e"; 
}
else {                # Linux
  $bs = "/";
  $rebs = "\/";
  $rm_cmd = "rm";
  $cp_cmd = "cp";
  $mv_cmd = "mv";
  $cd_cmd = "cd";
  $ext_exe = "";   
  $cpdir_cmd = "cp -R";
}

### Commands for symbolic links are platform dependent
sub make_symlink {
    my ($file, $symlink) = @_;
    if ($sys eq "win") { system("mklink $symlink $file"); }
    else               { system("ln -sf $file $symlink"); }
}

### Read a file into a string
sub read_file {
    my ($filename, $r_string) = @_;
    if (open(INPUT, "$filename")) {
      if (scalar @_ > 1) {
        $$r_string = "";
        while (<INPUT>) { $$r_string .= $_; }
        close(INPUT);
        return 1;
      }
      else {
        my $str = "";
        while (<INPUT>) { $str .= $_; }
        close(INPUT);
        return $str;
      }
    }
    else { return 0; }
}

### Save a string into a file
sub save_file {
    my ($filename, $r_string) = @_;
    open(OUTPUT, ">$filename");
    print OUTPUT $$r_string;
    close(OUTPUT);
}

### Get all filenames of a given directory
sub dirfilelist {
  my ($dir) = @_;

  if (!opendir(OP, "$dir")) { print "dirfilelist: Can't open $dir\n"; return (); }
  my @names = grep { -f "$dir/$_" } readdir(OP);   # pouze soubory, ne adresare
  closedir(OP);

  return @names;
}

### Get all subdirectories of a given directory
sub subdirlist {
  my ($dir) = @_;

  if (!opendir(OP, "$dir")) { &ERROR("subdirlist: Can't open $dir\n");return (); }
  my @names = grep { (!/^\.+$/) && -d "$dir$bs$_" } readdir(OP);   # pouze adresare (az na "." a ".."), ne soubory
  closedir(OP);

  return @names;
}

### Make all directories according the given path
sub make_dir {
  my ($dir) = @_;
  if ($dir !~ /$rebs$/) { $dir = $dir.$bs; }
  my $up_dir = "";
  if (!opendir(DIR, $dir)) {
    while ($dir =~ s/(\/?.*?)$rebs//) {
      if (!opendir(SUBDIR, "$up_dir$1")) {
        if (!mkdir("$up_dir$1")) {
          print "Failed to make the directory $up_dir$1 !\n";
          print "Make sure you have permissions to do it.\n";
          die "Ending program ...\n";  
        }
      }
      else {
        closedir(SUBDIR);
      }
      $up_dir .= "$1$bs"; 
    }
  } 
  else {
    closedir(DIR);
  }
}

### Replace a string with a new string in some file (similar to "sed" utility on unix)
sub replace_in_file {
    my ($file, $replace_what, $replace_with, $new_file) = @_;
    if (!(defined $new_file) || $new_file eq "") { $new_file = $file; }
    my $str = "";
    &read_file($file, \$str);
    $str =~ s/$replace_what/$replace_with/sg;
    &save_file($new_file, \$str);
}

### Compare two files in a binary mode
sub compare_files {
  my ($file1, $file2) = @_;
  my $buffer1, $buffer2;
  my $same = 1;

  open INFA, $file1
    or die "\nCan't open $file1 for reading: $!\n";
  open INFB, $file2
    or die "\nCan't open $file2 for reading: $!\n";

  binmode INFA;
  binmode INFB;

  while (
    read (INFA, $buffer1, 65536)	# read in (up to) 64k chunks
    and read (INFB, $buffer2, 65536)
  ) {
    if ($buffer1 ne $buffer2) {
      $same = 0;
    }
  };
  die "Problem comparing: $!\n" if $!;

  close INFA
    or die "Can't close $file1: $!\n";
  close INFB
    or die "Can't close $file2: $!\n";

  return $same;
}

### Secure copy of a file in a binary mode
sub copy_file {
  my ($srcfile, $destfile) = @_;
  my $buffer;

  open INF, $srcfile
    or die "\nCan't open $srcfile for reading: $!\n";
  open OUTF, ">$destfile"
    or die "\nCan't open $destfile for writing: $!\n";

  binmode INF;
  binmode OUTF;

  while (
    read (INF, $buffer, 65536)	# read in (up to) 64k chunks, write
    and print OUTF $buffer	# exit if read or write fails
  ) {};
  die "Problem copying: $!\n" if $!;

  close OUTF
    or die "Can't close $destfile: $!\n";
  close INF
    or die "Can't close $srcfile: $!\n";
}

1;

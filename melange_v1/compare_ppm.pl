#!/usr/bin/perl -w


# Parse command-line arguments
if ($#ARGV != 1) {
  print STDERR "Usage: compare_ppm.pl <ppm file #1> <ppm file #2>\n";
  exit(1);
}

$file1 = $ARGV[0];   $file2 = $ARGV[1];

# Open files
open(FILE1, $file1) or die "Cannot open file $file1\n";
open(FILE2, $file2) or die "Cannot open file $file2\n";

$count=0;
@distance=();
while ($line1 = <FILE1>) {
  chop($line1);  

  $line2 = <FILE2>; 
  if (!defined($line2)) {
    print STDERR "Failed - mismatch found! (1)\n"; exit(1);
  }
  chop($line2);
  if ($line1 eq "P3") {
    if ($line2 eq "P3") { 
      next; 
    } else {
      print STDERR "Failed - mismatch found! (2)\n"; exit(1);
    }
  }

  @tokens_1 = split(/ +/,$line1);
  @tokens_2 = split(/ +/,$line2);
  if ($#tokens_1 != $#tokens_2) {
    print STDERR "$#tokens_1   $#tokens_2\n";
    print STDERR "Failed - mismatch found! (3)\n"; exit(1);
  }

  for ($i=0;$i<=$#tokens_1;$i++) {
#    print STDOUT (abs($tokens_1[$i] - $tokens_2[$i]))." ";
    $distance[$count] = abs($tokens_1[$i] - $tokens_2[$i]);
    $count++;
  }
}

if (<FILE2>) {
    print STDERR "Failed - mismatch found! (5)\n"; exit(1);
}

close(FILE1); close(FILE2);

$error_count=0;
$max_error=0;
$average_error=0;
$average_distance=0;
for ($i=0; $i<=$#distance; $i++) {
#  print STDOUT "$distance[$i] ";
  if ($distance[$i] > 0) {
    $error_count++;
    if ($distance[$i] > $max_error) {
      $max_error = $distance[$i];
    }
    $average_error += $distance[$i];
  }
  $average_distance += $distance[$i];
}

if ($error_count > 0) {
  $average_error /= $error_count;
} else {
  $average_error = 0;
}
$average_distance /= ($#distance+1);

print STDERR "# of pixels that differ   = ".$error_count."\n";
print STDERR "  - max difference        = ".$max_error."\n";
print STDERR "  - average difference    = ".$average_error."\n";

if ($error_count == 0) {
  print STDOUT "EQUAL\n";
} else {
  print STDOUT "DIFFERENT\n";
}

#print STDERR "Success - files match (6)\n"; exit(0);



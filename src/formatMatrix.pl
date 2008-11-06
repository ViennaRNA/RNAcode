#!/usr/bin/perl -w

my $name='BLOSUM62';

print "int BLOSUM62[24][24]={";

while (<>){

  next if /^#/;
  next if /^\s+A\s+R/;

  my @entries=split;

  shift @entries;

  foreach my $i (0..$#entries){
    $entries[$i]=" ".$entries[$i] if not $entries[$i]=~/-/;
  }


  print "{",join(',',@entries),"},\n";


}


print "};\n";

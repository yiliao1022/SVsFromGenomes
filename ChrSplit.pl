#!/usr/bin/perl
use strict;
use warnings;

my $target = $ARGV[0];
my $query = $ARGV[1];

iTools Fatools split -InFa $target -OutDir ./
iTools Fatools split -InFa $query -OutDir ./
system "gunzip ./$target\_cut/*.gz";
system "gunzip ./$query\_cut/*.gz";

system "cd ./$target\_cut; ls *.fa > target.lst";
system "cd ./$query\_cut; ls *.fa > query.lst";

open In1, "./$target\_cut/target.lst" or die "$!";
open In2, "./$query\_cut/query.lst" or die "$!";
my @target;
my @query;

while (<In1>) {
chomp;
push (@target, $_);
}
while (<In2>){
chomp;
push (@query, $_);
}


open OUT1, ">target.sh" or die "$!";
open OUT2, ">query.sh" or die "$!";
foreach my $t (@target) {
   foreach my $q (@query) {
      print OUT1 "~/Softwares/minimap2-2.17_x64-linux/minimap2 -c --cs=long ./$target\_cut/$t ./$query\_cut/$q > ./$target\_cut/$t.$q.long.paf\n";
      print OUT2 "~/Softwares/minimap2-2.17_x64-linux/minimap2 -c --cs=long ./$query\_cut/$q ./$target\_cut/$t > ./$query\_cut/$q.$t.long.paf\n";
  }
}



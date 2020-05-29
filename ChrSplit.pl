#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

#split genomes into chrs and make all-by-all minimap2 alignment commands
#Authors: Yi Liao, Shujun Ou (05/23/2020)
#usage: perl ChrSplit.pl $target.fasta $query.fasta

my $target = $ARGV[0];
my $query = $ARGV[1];
my $target_file = basename($target);
my $query_file = basename($query);

#software paths
my $minimap2 = "/home/oushujun/las/bin/miniconda2/envs/EDTA/bin/minimap2";

# split query and subject into single chr files
split_fa($target) unless -s "$target_file\_cut/$target_file.lst";
split_fa($query) unless -s "$query_file\_cut/$query_file.lst";

# read seq names
open In1, "<./$target_file\_cut/$target_file.lst" or die "$!";
open In2, "<./$query_file\_cut/$query_file.lst" or die "$!";
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

# write out sh lines
open OUT1, ">target_$target_file-$query_file.sh" or die "$!";
open OUT2, ">query_$query_file-$target_file.sh" or die "$!";
foreach my $t (@target) {
   foreach my $q (@query) {
      print OUT1 "/usr/bin/time -v $minimap2 -c --no-kalloc --print-qname --cs=long -t 18 ./$target_file\_cut/$t ./$query_file\_cut/$q -o ./$target_file\_cut/$t.$q.long.paf\n";
      print OUT2 "/usr/bin/time -v $minimap2 -c --no-kalloc --print-qname --cs=long -t 18 ./$query_file\_cut/$q ./$target_file\_cut/$t -o ./$query_file\_cut/$q.$t.long.paf\n";
  }
}


## subroutines
# split genome into one sequence per file, keep seq names with "chr" only
sub split_fa {
	my $genome = $_[0];
	my $genome_file = basename($genome);
	system "mkdir ${genome_file}_cut" unless -e "${genome_file}_cut" && -d "${genome_file}_cut";
	chdir "${genome_file}_cut";
	open FA, "<../$genome" or die $!;
	$/ = "\n>";
	while (<FA>){
		s/>//g;
		my ($id, $seq) = (split /\n/, $_, 2);
		$seq =~ s/\s+//g;
		next unless $id =~ /chr/i;
		open OUT, ">$id.fa" or die $!;
		print OUT ">$id\n$seq\n";
		close OUT;
		}
	close FA;
	system "ls *.fa > $genome_file.lst";
	chdir "..";
	$/ = "\n";
	}

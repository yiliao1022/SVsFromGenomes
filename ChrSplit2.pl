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

my $array = "1..10";
$array = $ARGV[2] if defined $ARGV[2];

#software paths
my $minimap2 = "/home/oushujun/las/bin/miniconda2/envs/EDTA/bin/minimap2";
my $paftools = "/home/oushujun/las/bin/miniconda2/envs/EDTA/bin/paftools.js";

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
my $target_id = (split /\./, $target_file)[0];
my $query_id = (split /\./, $query_file)[0];

my $job = "#!/bin/bash
#SBATCH -N 1
#SBATCH -n 11
#SBATCH --mem=180GB
##SBATCH --mem=220GB
#SBATCH -t 21-00:00:00
#SBATCH -J align_$target_id-$query_id

conda activate EDTA
cd /home/oushujun/jfw/TE/MaizeNAM/SV/LY_pipeline/${query_file}_cut
for i in {$array}; do
# Chr map to Chr
#	/usr/bin/time -v $minimap2 -c -x asm5 --no-kalloc --print-qname --cs=long -t 10 ../${target_file}_cut/${target_id}_chr\$i.fa ./${query_id}_chr\$i.fa -o ${query_id}_chr\$i.mapped-to.${target_id}_chr\$i.paf
#	sort --parallel=36 -k6,6 -k8,8n ${query_id}_chr\$i.mapped-to.${target_id}_chr\$i.paf > ${query_id}_chr\$i.mapped-to.${target_id}_chr\$i.sorted.paf
#	k8 $paftools call ${query_id}_chr\$i.mapped-to.${target_id}_chr\$i.sorted.paf > ${query_id}_chr\$i.mapped-to.${target_id}_chr\$i.sorted.var.txt
#	k8 $paftools view -f maf ${query_id}_chr\$i.mapped-to.${target_id}_chr\$i.sorted.paf > ${query_id}_chr\$i.mapped-to.${target_id}_chr\$i.sorted.maf

# Chr map to genome
	/usr/bin/time -v $minimap2 -c -x asm5 --no-kalloc --print-qname --cs=long -t 10 ../${target_file}_cut/${target_id}_chr\$i.fa ../genomes/${query_file} -o ${query_id}_chr\$i.mapped-to.${target_id}.paf
	sort --parallel=36 -k6,6 -k8,8n ${query_id}_chr\$i.mapped-to.${target_id}.paf > ${query_id}_chr\$i.mapped-to.${target_id}.sorted.paf
	k8 $paftools call ${query_id}_chr\$i.mapped-to.${target_id}.sorted.paf > ${query_id}_chr\$i.mapped-to.${target_id}.sorted.var.txt
	k8 $paftools view -f maf ${query_id}_chr\$i.mapped-to.${target_id}.sorted.paf > ${query_id}_chr\$i.mapped-to.${target_id}\.sorted.maf
done

# check if completed successfully
if [ \$? -eq 0 ]
then
  echo Succeed!
  exit 0
else
  echo Failed!
exit 1
fi

";

# print out job info
open OUT1, ">align_$target_id-$query_id.$array.qsub" or die $!;
print OUT1 $job; 
close OUT1;


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

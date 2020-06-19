#!/usr/bin/env perl 
use strict;
use warnings;
#Shujun Ou (shujun.ou.1@gmail.com)
#06/19/2020

my $usage = "perl combine_SV.pl *sort.out | sort -sV > combined.SV";
my $min_len = 1;
my $SV_type = "DEL"; #INS or DEL

my %SV;
my %ref; #store all variant locations on ref of the input data
while (<>){
	my ($ref_chr, $ref_s, $ref_e, $alt_chr, $alt_s, $alt_e, $len) = (split)[0,1,2,4,5,6,7];
	next if $len < $min_len;
	my ($ref_info, $alt_info) = ("${ref_chr}_${ref_s}_$ref_e", "${alt_chr}_${alt_s}_$alt_e");
	($ref_info, $alt_info) = ($alt_info, $ref_info) if $SV_type eq "DEL";
	my $id = $1 if $ref_info =~ /^(.*)\./; #genome id is separate by .
	$SV{$id}{$alt_info} = $ref_info;
	$ref{$alt_info} = $alt_info;
	}

my @ids = sort{$a cmp $b} (keys %SV);
foreach my $ref_loc (sort{$a cmp $b} keys %ref){
	my $len = $2 - $1 + 1 if $ref_loc =~ /.*_([0-9]+)_([0-9]+)/;
	my ($geno_info, $alt_count, $af, $maf) = ('', 0, 0, 0);
	my $total = @ids;
	my $ref = $ref_loc;
	$ref =~ s/_/\t/g;
	print "$ref\t$SV_type\t";
	foreach my $id (@ids){
		my $alt_loc = $ref_loc;
		($alt_loc = $SV{$id}{$ref_loc} and $alt_count++) if exists $SV{$id}{$ref_loc};
		my $alt_len = $2 - $1 + 1 if $alt_loc =~ /.*_([0-9]+)_([0-9]+)/;
		my $geno = "$alt_loc:$alt_len";
		$len = $alt_len if $alt_len > $len;
#		$geno = ("1/1", "0/0")[$alt_loc eq $ref_loc];
		$geno_info .= "$geno\t";
		}
	$geno_info =~ s/\s+$//;
	$af = sprintf("%.3f", $alt_count/$total);
	$maf = ($af, 1 - $af)[$af > 0.5];
	print "$len\t$alt_count\t$af\t$maf\t$geno_info\n";
	}
		
	
#input format
#Ab10.chr1       1594452 1594453 DEL     B73.chr1        1582907 1582907 1       0.95    114088  52668   
#Ab10.chr1       1647121 1647125 DEL     B73.chr1        1635575 1635575 4       0.95    52668   247401  
#Ab10.chr1       1894526 1894528 DEL     B73.chr1        1882990 1882990 2       0.95    247401  19822   
#Ab10.chr1       1914350 1914372 DEL     B73.chr1        1902817 1902817 22      0.95    19822   6281    
#Ab10.chr1       1920653 1920655 DEL     B73.chr1        1909101 1909101 2       0.95    6281    105018  
#Ab10.chr1       2025673 2025679 DEL     B73.chr1        2014122 2014122 6       0.95    105018  159197  
#Ab10.chr1       2184876 2184877 DEL     B73.chr1        2173319 2173319 1       0.95    159197  175079

#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(&basename &dirname);

####################################################################################
# ChainNetSynnet.pl
#
# Authors: Yi Liao (05/28/2020)
# 
# Copyright (C) Not for commercial use
#
# Chaining Netting and adding synteny infortmation for the output (.axt format) from  
# pairwise whole genome alignment(LAST, LASTZ, MUMmer,minimap2 et al...)
# 
# Prerequisite : Kent's utilities;
#
# Usage:  perl ChainNetSynnet.pl -tseq target.fasta -qseq query.fasta -axt . -linearGap medium -tname A -qname B
#
# Options: -axt  [FILE] Folder where the axt files stored                    [REQUIRED]
#          -tseq [FILE] target genome sequence in fasta format               [REQUIRED]
#          -qseq [FILE] query genome sequence in fasta format                [REQUIRED]
#          -linearGap [FILE] linearGap file for ChainNet                     [REQUIRED]
#          -tname the defined name of the target genome assembly
#          -qname the defined name of the query genome assembly
#          -help  print this information
#
#
####################################################################################
my ($target,$query,$axt,$linearGap,$tname,$qname,$Help);
GetOptions( 'tname=s' => \$tname,
            'qname=s' => \$qname,
            'tseq=s' => \$target,
            'qseq=s' => \$query,
            'axt=s' => \$axt,
            'linearGap=s' => \$linearGap,
            'help' => \$Help
            );



if ($Help){
print <<"END.";
  Usage:  perl ChainNetSynnet.pl -tseq target.fasta -qseq query.fasta -axt /path/to/axt -linearGap medium -tname B73 -qname NC350
  
  Options:  -axt  [FILE] Folder where the axt files stored                    [REQUIRED]
            -tseq [FILE] target genome sequence in fasta format               [REQUIRED]
            -qseq [FILE] query genome sequence in fasta format                [REQUIRED]
            -linearGap [FILE] linearGap file for ChainNet                     [REQUIRED]
            -tname the defined name of the target genome assembly
            -qname the defined name of the query genome assembly
            -help  print this information
END.
exit;
}




my $target_file = basename($target);
my $query_file = basename($query);
my $target_path = dirname($target);
my $query_path = dirname($query);
my $axt_path = dirname($axt);

unless (-e "$target_path/$target_file.sizes") {
        `faSize -detailed $target > $target_path/$target_file.sizes`;
} 
unless (-e "$target_path/$target_file.2bit") {
        `faToTwoBit $target $target_path/$target_file.2bit`;
}
unless (-e "$query_path/$query_file.sizes") {
        `faSize -detailed $query > $query_path/$query_file.sizes`;
}
unless (-e "$query_path/$query_file.2bit") {
        `faToTwoBit $query $query_path/$query_file.2bit`;
}

`for i in $axt/*.axt; do axtChain -linearGap=$linearGap \$i $target_path/$target_file.2bit $query_path/$query_file.2bit \$i.chain;done`;
`chainMergeSort $axt/*.chain > $axt/$tname.$qname.all.chain`;
`chainPreNet $axt/$tname.$qname.all.chain $target_path/$target_file.sizes $query_path/$query_file.sizes $axt/$tname.$qname.all.chain.filter`;
`chainNet -minSpace=1 $axt/$tname.$qname.all.chain.filter $target_path/$target_file.sizes $query_path/$query_file.sizes $axt/$tname.$qname.all.chain.filter.tnet $axt/$tname.$qname.all.chain.filter.qnet`;
`netSyntenic $axt/$tname.$qname.all.chain.filter.tnet $axt/$tname.$qname.all.chain.filter.tnet.synnet`;








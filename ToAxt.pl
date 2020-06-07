##!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

####################################################################################
# ToAxt.pl
#
# Authors: Yi Liao (05/28/2020) Shujun Ou (05/31/2020)
# 
# Copyright (C) Not for commercial use
#
# Prerequisite : Kent's utilities; minimap2; paftools.js, please make sure all these utilities are on your $PATH before run it.
#
# Convert genome alignment file2 to axt format (LAST, LASTZ, MUMmer,minimap2 et al...) and call SVs using paftools.js
#
# Usage:  perl ToAxt.pl -aligner minimap2 -wkfolder . -tname A -qname B
#
# Options: -aligner  [last|lastz|MUMmer|minimap2] aligner for whole genome pairwise alignment      [REQUIRED]
#          -wkfolder [FILE] where the alignment files put                                          [REQUIRED]
#          -tname    [FILE] the name of target genome assembly                                     [REQUIRED]
#          -qname    [FILE] the name of query genome assembly                                      [REQUIRED]
#          -help     print this information
####################################################################################

my ($aligner,$wk,$tname,$qname,$Help,$kenUtil,$TBA,$minimap2);

# program path. Will use the ones in $PATH if unspecified.
$minimap2 = '';
$TBA = '~/las/git_bin/multiz-tba.012109/';
$kenUtil = '~/las/git_bin/kentUtils/linux.x86_64.v385/';

# get input parameters
GetOptions( 'tname=s' => \$tname,
            'qname=s' => \$qname,
            'aligner=s' => \$aligner,
            'wkfolder=s' => \$wk,
            'kenUtil=s' => \$kenUtil,
            'help' => \$Help
           );


if ($Help){
print <<"END.";
  Usage:  perl ToAxt.pl -aligner minimap2 -wkfolder /path/to/align -tname A -qname B 

  Options:  -aligner  [last|lastz|MUMmer|minimap2] aligner for whole genome pairwise alignment      [REQUIRED]
            -wkfolder [path] where the alignment files put                                          [Optional]
	    -kenUtil  [path] where the ken Utility programs put					    [Optional]
            -tname    [FILE] the name of target genome assembly                                     [REQUIRED]
            -qname    [FILE] the name of query genome assembly                                      [REQUIRED]
            -help     print this information                                                        [Optional]
END.
exit;
}


# make output directories
`mkdir $wk/Target_$tname` unless -e "$wk/Target_$tname" && -d "$wk/Target_$tname";
`mkdir $wk/Target_$qname` unless -e "$wk/Target_$qname" && -d "$wk/Target_$qname";


if ($aligner eq "minimap2") {
`for i in \`ls $wk/*.paf|grep -v sorted\`; do
 i=\$(echo \$i|perl -nle 's/.paf//; print \$_');
 sort -k6,6 -k8,8n \$i.paf > \$i.sorted.paf;
 ${minimap2}paftools.js call \$i.sorted.paf > \$i.sorted.var.txt;
 ${minimap2}paftools.js view -f maf \$i.sorted.paf | sed 's/^a /a score=/g; s/_/./g' > \$i.sorted.maf;
# perl -i -nle 's/^a ([0-9]+)/a score=\$1/; s/_/./g; print \$_' \$i.sorted.maf;
 ${TBA}maf_order \$i.sorted.maf $qname $tname > \$i.$qname.$tname.maf;
 ${kenUtil}mafToAxt \$i.sorted.maf $tname $qname  \$i.$tname.$qname.axt;
 ${kenUtil}mafToAxt \$i.$qname.$tname.maf $qname $tname \$i.$qname.$tname.axt;
 mv \$i.$tname.$qname.axt $wk/Target_$tname;
 mv \$i.$qname.$tname.axt $wk/Target_$qname;
 rm \$i.$qname.$tname.maf \$i.sorted.maf;
 rm \$i.sorted.paf;
done`;
} 

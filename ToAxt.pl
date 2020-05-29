##!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

####################################################################################
# ToAxt.pl
#
# Authors: Yi Liao (05/28/2020)
# 
# Copyright (C) Not for commercial use
#
# Prerequisite : Kent's utilities; minimap2; paftools.js, please make sure all these utilities are on you $PATH before run it.
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
my ($aligner,$wk,$tname,$qname,$Help);
GetOptions( 'tname=s' => \$tname,
            'qname=s' => \$qname,
            'aligner=s' => \$aligner,
            'wkfolder=s' => \$wk,
            'help' => \$Help
           );



if ($Help){
print <<"END.";
  Usage:  perl ToAxt.pl -aligner minimap2 -wkfolder /path/to/align -tname A -qname B 

  Options:  -aligner  [last|lastz|MUMmer|minimap2] aligner for whole genome pairwise alignment      [REQUIRED]
            -wkfolder [FILE] where the alignment files put                                          [REQUIRED]
            -tname    [FILE] the name of target genome assembly                                     [REQUIRED]
            -qname    [FILE] the name of query genome assembly                                      [REQUIRED]
            -help     print this information                                                        [Optional]
END.
exit;
}



`mkdir $wk/Target_$tname`;
`mkdir $wk/Target_$qname`;

if ($aligner eq "minimap2") {
`for i in $wk/*.paf; 
 do sort -k6,6 -k8,8n \$i > \$i.sort.paf;
 paftools.js call \$i.sort.paf > \$i.sort.var.txt;
 paftools.js view -f maf \$i.sort.paf  > \$i.sort.maf;
 sed -i 's/^a /a score=/g' \$i.sort.maf;
 sed -i 's/_/./g' \$i.sort.maf;
 maf_order \$i.sort.maf $qname $tname > \$i.$qname.$tname.maf;
 mafToAxt \$i.sort.maf $tname $qname  \$i.$tname.$qname.axt;
 mv \$i.$tname.$qname.axt $wk/Target_$tname;
 mafToAxt \$i.$qname.$tname.maf $qname $tname \$i.$qname.$tname.axt;
 mv \$i.$qname.$tname.axt $wk/Target_$qname;
 rm $wk/*.maf;
 rm $wk/*.sort.paf; done`;
} 

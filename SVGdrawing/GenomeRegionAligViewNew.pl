#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use SVG;

####################################################################################
# SVlocalView.pl
#
# Authors: Yi Liao (06/01/2020) <yiliao1022@gmail.com>
#
# Copyright (C) Not for commercial use
#
# Drawing syntenic map for local regions that exhibit a potential SV (strutral vaiation) based on SVcaller.pl.
#
# Prerequisite : Perl module SVG, svg2xxxx, BEDOPS v2.4.39
#
# perl GenomeRegionAligView.pl -sv 1.txt -target /home/yliao/yliao/Ref/Japonica/O_sativa_japonica_sm21.fasta -query /home/yliao/yliao/Ref/kasalath/kasalath.fasta -tbed /home/yliao/yliao/Ref/Japonica/all.bed -qbed /home/yliao/yliao/2018_tRNA/ref/annotation/Kasalath.IGDBv2.bed -len 2000
#      -sv      1.txt
#      -target  target genome sequence
#      -query   query genome sequence
#      -len     sequence length for extending the flanking region
#      -tbed    target bed file
#      -qbed    query bed file
#      -verbose  verbose
#      -help     help
#
####################################################################################
my ($sv,$wk,$target,$query,$Tbed,$Qbed,$length,$lastz,$tba,$svg2xxx,$bedops,$Help);

GetOptions(
  "sv:s"=>\$sv,
  "wk:s"=>\$wk,
	"target:s"=>\$target,
	"query:s"=>\$query,
	"tbed:s"=>\$Tbed,
	"qbed:s"=>\$Qbed,
	"len:i"=>\$length,
  "lastz:s"=>\$lastz,
  "tba:s"=>\$tba,
  "bedops:s"=>\$bedops,
  "svg2xxx:s"=>\$svg2xxx,
	"help"     =>\$Help
);
###################################################################################################
$wk ||='.';
$length ||=10000;
$lastz ||='';
$tba ||='';
$bedops ||='';
$Tbed ||='NA';
$Qbed ||='NA';

##################################################################################################
if ($Help){
print <<"END.";
  Usage: perl $0 -sv SV.txt -target -wk . /home/yliao/yliao/Ref/Japonica/O_sativa_japonica_sm21.fasta -query /home/yliao/yliao/Ref/kasalath/kasalath.fasta -tbed /home/yliao/yliao/Ref/Japonica/all.bed -qbed /home/yliao/yliao/2018_tRNA/ref/annotation/Kasalath.IGDBv2.bed -len 2000 -tba /path/to/tba -lastz /path/to/lastz -bedops /path/to/bedops -svg2xxx /path/to/svg2xxx
  -sv      structural variations output                                             [REQUIRED]
  -wk      working folder where the SV outout file put [default]                    [OPTIONAL]
  -target  target genome sequence                                                   [REQUIRED]
  -query   query genome sequence                                                    [REQUIRED]
  -len     sequence length for extending the flanking region [default: 10000]       [OPTIONAL]
  -tbed    target bed file for gene annotation                                      [OPTIONAL]
  -qbed    query bed file for gene annotation                                       [OPTIONAL]
  -lastz   path to lastz program                                                    [OPTIONAL]
  -bedops  path to bedops tools                                                     [OPTIONAL]
  -svg2xxx Path to svg2xxx tools                                                    [REQUIRED]
  -help     help
END.
exit;
}
###################################################################################################
open IN, "$sv" or die "$!";
$/="\n";
while (<IN>) {
chomp;
SVGlocal($_,$length,$target,$query,$Tbed,$Qbed);
}
close IN;
##################################################################################################
sub SVGlocal {
my ($line,$extend,$target,$query,$TBED,$QBED) = @_;
my @temp = split (/\t/,$line);
my $T_beg=$temp[1]-$extend;
my $Q_beg=$temp[5]-$extend;
my $T_end=$temp[2]+$extend;
my $Q_end=$temp[6]+$extend;
my $T_len = $T_end-$T_beg;
my $Q_len = $Q_end-$Q_beg;
my @unit1 = split (/\./,$temp[0]);
my @unit2 = split (/\./,$temp[4]);
$temp[0]=$unit1[1];
$temp[4]=$unit2[1];
&SubSeq($target,"$temp[0]",$T_beg,$T_len);
&SubSeq($query,"$temp[4]",$Q_beg,$Q_len);

`${lastz}lastz $temp[0]\_$T_beg\_$T_len $temp[4]\_$Q_beg\_$Q_len --format=maf --nogapped --out=$temp[0].$T_beg.$temp[4].$Q_beg.maf`;
`${tba}single_cov2 $temp[0].$T_beg.$temp[4].$Q_beg.maf R=$temp[0]\_$T_beg\_$T_end > $temp[0].$T_beg.$temp[4].$Q_beg.sing.maf`;
`rm $temp[0]\_$T_beg\_$T_len $temp[4]\_$Q_beg\_$Q_len $temp[0].$T_beg.$temp[4].$Q_beg.maf`;

#####################################################
my $svg= SVG -> new (width=>600, height =>300);
my $rate;
my $rectange_height=4;
my $chr_color="black";
my $beg=100;
my $end=500;

if ($T_len > $Q_len) {
$rate = $T_len/400;
} else {
$rate = $Q_len/400;
}

######Draw skeleton

$svg->line (x1=> $beg, y1 => 80, x2=> $beg + $T_len/$rate, y2=>80, style => { 'fill'=> $chr_color, 'stroke'=> $chr_color,'stroke-width'=>'1.5',});
my $inta = $T_len/$rate/10;
my $i;
for ($i=0;$i<11;$i++) {
	 $svg->line ( x1=> $beg+$inta*$i, y1 => 80, x2=> $beg+$inta*$i, y2 => 74, style => { 'fill'=> $chr_color, 'stroke'=> $chr_color,'stroke-width'=>'1.5',} );
  }

$svg->line (x1=> $beg, y1 => 224, x2=> $beg + $Q_len/$rate, y2=>224, style => { 'fill'=> $chr_color, 'stroke'=> $chr_color,'stroke-width'=>'1.5',});
  my $intb = $Q_len/$rate/10;
  my $j;
  for ($j=0;$j<11;$j++) {
  	 $svg->line ( x1=> $beg+$intb*$j, y1 => 224, x2=> $beg+$intb*$j, y2 => 230, style => { 'fill'=> $chr_color, 'stroke'=> $chr_color,'stroke-width'=>'1.5',} );
    }

$svg->text("x", 78,"y", 60,"-cdata","$unit1[0].$unit1[1]","font-family","arial","font-size",10,"fill","black" );
$svg->text("x", 98,"y", 70,"-cdata","$T_beg","font-family","arial","font-size",10,"fill","black" );
$svg->text("x", 492,"y", 70,"-cdata","$T_end","font-family","arial","font-size",10,"fill","black");

$svg->text("x", 78,"y", 250,"-cdata","$unit2[0].$unit2[1]","font-family","arial","font-size",10,"fill","black" );
$svg->text("x", 98,"y", 240,"-cdata","$Q_beg","font-family","arial","font-size",10,"fill","black" );
$svg->text("x", 492,"y", 240,"-cdata","$Q_end","font-family","arial","font-size",10,"fill","black");

$svg->rectangle ( x=>100, y=>100,width=>$T_len/$rate,height=>$rectange_height,style => { 'fill'=> $chr_color, 'stroke'=> $chr_color,});
$svg->rectangle ( x=>100, y=>200,width=>$Q_len/$rate,height=>$rectange_height,style => { 'fill'=> $chr_color, 'stroke'=> $chr_color,});

######Draw synteny relationship from alignment Maf files
my $aligns = &ReadMaf("$wk/$temp[0].$T_beg.$temp[4].$Q_beg.sing.maf");

foreach my $block (@$aligns) {
  my @temp=split("_",$block);
  my $color="lightsalmon";
  my $tleft=$temp[0]/$rate+100;
  my $tright=$temp[1]/$rate+100;
  my $qright=$temp[2]/$rate+100;
  my $qleft=$temp[3]/$rate+100;
  my $theight=114;
  my $qheight=190;
  my $xv=[$tleft,$tright,$qright,$qleft];
  my $yv=[$theight,$theight,$qheight,$qheight];
  my $points=$svg->get_path( x=>$xv,y=>$yv,-type=>'polyline',-closed=>'true');
  my $tag=$svg->polyline( %$points, style=>{fill=>$color});
}

######Draw Gene features from gff files

if (!$TBED eq "NA")  {
`echo -e "$temp[0]\t$T_beg\t$T_end" | ${bedops}bedextract $TBED - > T$temp[0]\_$T_beg\_$T_end.bed`;
open TBED, "T$temp[0]\_$T_beg\_$T_end.bed" or die "$!";
$/="\n";
while (<TBED>) {
chomp;
my @temp = split (/\t/,$_);
if ($temp[7]=~/CDS/ and $temp[3]=~/\.1/ and $temp[5]=~/\-/) {
$svg->rectangle (x=>($temp[1]-$T_beg)/$rate+100, y=>106,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif ($temp[7]=~/CDS/ and $temp[3]=~/\.1/ and $temp[5]=~/\+/) {
$svg->rectangle (x=>($temp[1]-$T_beg)/$rate+100, y=>92,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif( ($temp[7]=~/three/|$temp[7]=~/five/ )and $temp[3]=~/\.1/ and $temp[5]=~/\+/){
$svg->rectangle (x=>($temp[1]-$T_beg)/$rate+100, y=>92,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif(($temp[7]=~/three/| $temp[7]=~/five/)and $temp[3]=~/\.1/ and $temp[5]=~/\-/){
$svg->rectangle (x=>($temp[1]-$T_beg)/$rate+100, y=>106,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
}elsif ($temp[7]=~/gene/ and $temp[5]=~/\-/) {
$svg->line (x1=>($temp[1]-$T_beg)/$rate+100,y1=>109,x2=>($temp[2]-$T_beg)/$rate+100,y2=>109, style=>{'fill'=>"green",'stroke'=>"green",'stroke-width'=>'1',});
} elsif ($temp[7]=~/gene/ and $temp[5]=~/\+/) {
$svg->line (x1=>($temp[1]-$T_beg)/$rate+100,y1=>95,x2=>($temp[2]-$T_beg)/$rate+100,y2=>95, style=>{'fill'=>"green",'stroke'=>"green",'stroke-width'=>'1',});
}
}
close TBED;
`rm T$temp[0]\_$T_beg\_$T_end.bed`; 
}


if (!$QBED eq "NA") {
`echo -e "$temp[4]\t$Q_beg\t$Q_end" | ${bedops}bedextract $QBED - > Q$temp[4]\_$Q_beg\_$Q_end.bed`;
open QBED, "Q$temp[4]\_$Q_beg\_$Q_end.bed" or die "$!";
while (<QBED>) {
chomp;
my @temp = split (/\t/,$_);
if ($temp[7]=~/CDS/ and $temp[9]=~/(\.|Bra)/ and $temp[5]=~/\-/) {
$svg->rectangle (x=>($temp[1]-$Q_beg)/$rate+100, y=>206,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif ($temp[7]=~/CDS/ and $temp[9]=~/(\.|Bra)/ and $temp[5]=~/\+/) {
$svg->rectangle (x=>($temp[1]-$Q_beg)/$rate+100, y=>192,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif( ($temp[7]=~/three/|$temp[7]=~/five/ )and $temp[9]=~/(\.T01|t1)/ and $temp[5]=~/\+/){
$svg->rectangle (x=>($temp[1]-$Q_beg)/$rate+100, y=>192,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif(($temp[7]=~/three/| $temp[7]=~/five/)and $temp[9]=~/(\.T01|t1)/ and $temp[5]=~/\-/){
$svg->rectangle (x=>($temp[1]-$Q_beg)/$rate+100, y=>206,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
}elsif ($temp[7]=~/gene/ and $temp[5]=~/\-/) {
$svg->line (x1=>($temp[1]-$Q_beg)/$rate+100,y1=>209,x2=>($temp[2]-$Q_beg)/$rate+100,y2=>209, style=>{'fill'=>"green",'stroke'=>"green",'stroke-width'=>'1',});
} elsif ($temp[7]=~/gene/ and $temp[5]=~/\+/) {
$svg->line (x1=>($temp[1]-$Q_beg)/$rate+100,y1=>195,x2=>($temp[2]-$Q_beg)/$rate+100,y2=>195, style=>{'fill'=>"green",'stroke'=>"green",'stroke-width'=>'1',});
}
}
`rm Q$temp[4]\_$Q_beg\_$Q_end.bed`;
close QBED;
}

######################################################
open OUT, ">$temp[0]\_$temp[1].$temp[4]\_$temp[5].svg" or die "Can not open my file";
print OUT $svg->xmlify;
close OUT;
`${svg2xxx}svg2xxx $temp[0]\_$temp[1].$temp[4]\_$temp[5].svg -t pdf`;
}

####################################################
################### Sub Routines ###################
####################################################

sub SubSeq {
my ($seq,$chr,$beg,$end) = @_;
local $/="\n>";
open FA, "$seq" or die "$!";
open Out, ">$chr\_$beg\_$end";
   while (my $sequence = <FA>) {
         chomp $sequence;
         my ($id) = $sequence =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
         print "$id\n";
	 if ($id eq $chr) {
         $sequence =~ s/^>*.+\n//;  # remove FASTA header
         $sequence =~ s/\n//g;  # remove endlines
	 my $seqrange = substr $sequence, $beg -1, $end;
         print Out ">$chr\_$beg\_$end\n$seqrange\n";
         last;
        }
    }
}

############
sub ReadMaf {
local $/="a score";
my $maf_file = shift;
my @align;
open MAF, "$maf_file" || die "Fail open $maf_file";
<MAF>;
while (<MAF>){
if ($_=~/\#/ and !$_=~/eof/) {
        next;
} else {
          my @temp = split (/\n/,$_);
          my @unit1 = split (/\s+/,$temp[1]);
          my @unit2 = split (/\s+/,$temp[2]);
          my $tleft=$unit1[2];
          my $tright=$unit1[2]+$unit1[3];
		  my ($qleft,$qright,$orientation);
      if ($unit2[4]=~/\+/) {
           $qleft=$unit2[2];
           $qright=$unit2[2]+$unit2[3];
           $orientation="plus";
         }elsif($unit2[4]=~/\-/){
           $qleft=$unit2[5]-$unit2[2];
           $qright=$unit2[5]-$unit2[2]-$unit2[3];
	       $orientation="minus";
		 }
  my $block = join ("_",($tleft,$tright,$qright,$qleft,$orientation));
  push (@align,$block);
  }
}
close MAF;
my $ref = \@align;
return $ref;
}


sub ReadGff {
local $/="\n";
my $gff_file = shift;
open GFF, "$gff_file" || die "Fail open $gff_file";
my %gff;
my $key;
while (<GFF>) {
chomp;
if ($_=~/\#/) {
  next;
} else {
my @temp = split (/\t/,$_);
if ($temp[2] eq "gene") {
$temp[8]=~/ID=(.*);Name=/;
$key = $1;
my @features=();
my $ref = \@features;
$gff{$1}=$ref;
}
if ($temp[8]=~/$key/) {
   if ($temp[2]=~/five/) {
   my $five = join ("_",($temp[2],$temp[3],$temp[4]));
   push (@{$gff{$key}},$five);
   }
   if ($temp[2]=~/three/) {
   my $three = join ("_",($temp[2],$temp[3],$temp[4]));
   push (@{$gff{$key}},$three);
   }
   if ($temp[2]=~/CDS/) {
   my $cds = join ("_",($temp[2],$temp[3],$temp[4]));
   push (@{$gff{$key}},$cds);
   }
  }
 }
}
close GFF;
my $reference=\%gff;
return $reference;
}

#!/usr/bin/perl -w
#call: analyseBlocks.pl bedfile genesfile outfile
#check the distance for each element to the adjacent blocks
#output: table with dist_to_left dist_to_right dist_left_to_right length_element length_leftblock length_rightblock

#as there are very little none blocks, they are ignored
    
use strict;
use warnings;
#use Pod::Usage;
use Data::Dumper;


my $bedfile = shift;
my $genesfile = shift;
my $outfile = shift;

open(my $outf, ">>",$outfile);

my $header = "Dist2Left\tDist2Right\tDistL2R\tLenElement\tLenLeftBlock\tLenRightBlock\n";
print $outf $header;

open FA,"<$genesfile" or die "can't open $genesfile\n";
while(<FA>){
    chomp;
    my $line=$_;
    if($line =~ /^#/){next;}
    if($line eq ""){next;}
    my @F = split '\t', $line;
    my @G = split '_', $F[1];
    my $spec = $G[0];
    my $elstart = $F[2];
    my $elend = $F[3];
    my $lnum = $F[5];
    my $rnum = $F[6];
    my $dist2left;
    my $dist2right;
    my $distL2R;
    my $ellen;
    my $llen;
    my $rlen;
    if($lnum eq "None"){
	next;
    }
    if($rnum eq "None"){
	next;
    }
   
    my $grpr = "zcat $bedfile \| grep $spec\_$rnum";
    my $grpl = "zcat $bedfile \| grep $spec\_$lnum";
    my @outl = readpipe("$grpl");
    if(scalar @outl == 0){next;}
    my @outr = readpipe("$grpr");
    if(scalar @outr == 0){next;}
    my @L = split '\t', $outl[0];
    my @R = split '\t', $outr[0];

    $dist2left = abs($L[3]-$elstart);
    $dist2right = abs($elend-$R[2]);
    $distL2R = abs($L[3]-$R[2]);
    $ellen = abs($elstart - $elend);
    $llen = abs($L[3]-$L[2]);
    $rlen = abs($R[3]-$R[2]);
    my $outstr = "$dist2left\t$dist2right\t$distL2R\t$ellen\t$llen\t$rlen\n";
    print $outf $outstr;

}
    

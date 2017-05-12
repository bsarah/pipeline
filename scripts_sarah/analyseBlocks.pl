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
    my $elstart;
    my $elend;
    if($F[2] < $F[3]){
	my $elstart = $F[2];
	my $elend = $F[3];
    }
    else{
	my $elstart = $F[3];
	my $elend = $F[2];
    }
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

    my $lstart;
    my $lend;
    my $rstart;
    my $rend;
    if($L[2] < $L[3]){
	$lstart = $L[2];
	$lend = $L[3];
    }
    else{
	$lstart = $L[3];
	$lend = $L[2];
    }
    if($R[2] < $R[3]){
	$rstart = $R[2];
	$rend = $R[3];
    }
    else{
	$rstart = $R[3];
	$rend = $R[2];
    }
    $dist2left = abs($lend-$lstart);
    $dist2right = abs($lend-$rstart);
    $distL2R = abs($lend-$rstart);
    $ellen = abs($lstart - $lend);
    $llen = abs($lend-$lstart);
    $rlen = abs($rend-$rstart);
    my $outstr = "$dist2left\t$dist2right\t$distL2R\t$ellen\t$llen\t$rlen\n";
    print $outf $outstr;

}
    

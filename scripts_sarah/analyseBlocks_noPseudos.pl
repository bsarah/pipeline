#!/usr/bin/perl -w
#call: analyseBlocks.pl TEMPfile genesfile outfile
#check the distance for each element to the adjacent blocks
#output: table with dist_to_left dist_to_right dist_left_to_right length_element length_leftblock length_rightblock

#as there are very little none blocks, they are ignored
    
use strict;
use warnings;
#use Pod::Usage;
use Data::Dumper;


my $bedfile = shift; #file in TEMP folder!
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
    #skip pseudogenes
    if($F[-2] eq "TRUE" || $F[-2] eq "true" || $F[-2] eq "T" || $F[-2] eq "t" || $F[-2] eq "True" || $F[-2] eq "1"){next;}
    my @G = split '_', $F[1];
    my $spec = $G[0];
    my $elstart;
    my $elend;
    if($F[2] < $F[3]){
	$elstart = $F[2];
	$elend = $F[3];
    }
    else{
	$elstart = $F[3];
	$elend = $F[2];
    }
    my $lnum = $F[5];
    my $rnum = $F[6];
    if($lnum eq "None"){
	next;
    }
    if($rnum eq "None"){
	next;
    }

    my $grpr = "zcat $bedfile \| grep -w $spec\_$rnum";
    my $grpl = "zcat $bedfile \| grep -w $spec\_$lnum";
    my @outl = readpipe("$grpl");
    if(scalar @outl == 0){next;}
    my @outr = readpipe("$grpr");
    if(scalar @outr == 0){next;}
    my @L = split '\t', $outl[0];
    my @R = split '\t', $outr[0];

    my $lstart = $L[2];
    my $llen = $L[3];
    my $rstart = $R[2];
    my $rlen = $R[3];

    #Always take distance from element to block start!
    #If not, there will be a problem with the distances in the reference species!
    #in the reference species, blocks are always directly next to each other
    
    my $dist2left = abs($lstart-$elstart);
    my $dist2right = abs($elend-$rstart);
    my $distL2R = abs($lstart-$rstart);
    my $ellen = abs($elstart - $elend);
#    my $outinfo = ">$lnum\t$rnum\t$lstart\t$elstart\t$elend\t$rstart\n";
    my $outstr = "$dist2left\t$dist2right\t$distL2R\t$ellen\t$llen\t$rlen\n";
#    print $outf $outinfo;
    print $outf $outstr;

}
    

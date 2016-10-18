#!/usr/bin/perl -w
## call writeBED.pl filelist outpath

## input file format:
## chromosome \t species_geneID \t genestart \t geneend \t strand
## \t blockleft \t blockright \t sequences and stuff not needed here

## output file format:
## filename
## list of elements in clus
## chromosome \t species_geneID \t genestart \t geneend \t strand
## \t blockleft \t blockright

use Data::Dumper;
use strict;
use warnings;

my $filename = shift;
my $outpath = shift;

my $outfile="$outpath/allClusters.bed";
open(my $outf,'>>',$outfile);


open FA,"$filename" or die "can't open $filename\n";

while(<FA>){
    my $curfile=$_;
    chomp($curfile);
    my $count=0;
    my @lines = ();
    open CF,"$curfile" or die "can't open $curfile\n";
    while(<CF>){
        my $line=$_;
        chomp($line);
	my @F = split '\t', $line;
	my $chr = $F[0];
	my $spec = $F[1];
	my $start = $F[2];
	my $end = $F[3];
	my $strand = $F[4];
	my $left = $F[5];
	my $right = $F[6];
	my $newline = "$chr\t$spec\t$start\t$end\t$strand\t$left\t$right\n";
	push @lines, $newline;
	$count++;
    }
    my @tmp = split '\/', $curfile;
    my $clusname = $tmp[(scalar @tmp)-1];
    my $newnewline = ">$clusname\t$count\n";
    print $outf $newnewline;
    for(my $i=0;$i<scalar @lines;$i++){
	print $outf $lines[$i];
    }
    $count=0;
    @lines=();
}

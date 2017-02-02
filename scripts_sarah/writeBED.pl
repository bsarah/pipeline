#!/usr/bin/perl -w
## call writeBED.pl filelist outpath outname summary

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
my $outname = shift;
my $summary = shift;

my $outfile="$outpath\/$outname";
open(my $outf,">>$outfile");
open(my $outs,">>$summary");


open FA,"$filename" or die "can't open $filename\n";

my $numclus = 0;
my $sumelems = 0;
my $numsingles = 0;
my $nummultcluselems = 0;

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
    my $clusname = $tmp[-1];
    my $newnewline = ">$clusname\t$count\n";
    print $outf $newnewline;
    for(my $i=0;$i<scalar @lines;$i++){
	print $outf $lines[$i];
    }
    if($count == 1){$numsingles++;}
    else{$nummultcluselems = $nummultcluselems + $count;}
    $sumelems = $sumelems + $count;
    $count=0;
    @lines=();
    $numclus++;
}


my $avelemnum = $sumelems/$numclus;
my $nummultclus = $numclus - $numsingles;
my $avmultcluselemnum = $nummultcluselems/$nummultclus;

print $outs "===============$outname\===============\n";
print $outs "The bed file $outname containing the genes sorted into 
the clusters is located here: $outpath \n";
print $outs "Format bed file (tab separated): 
  >cluster_name(including species and coordinates) number of elements
  for each element one line consisting of: chromosome, species, start, end, 
  strand, 5' blocknumber, 3'blocknumber \n";
print $outs "\n";
print $outs "Number of clusters (excluding None clusters): $numclus \n";
print $outs "Number of singleton clusters thereof: $numsingles \n";
print $outs "Number of multi-element clusters thereof: $nummultclus \n";
print $outs "Average number of elements per cluster (including singletons): $avelemnum \n";
print $outs "Average number of elements per cluster (excluding singletons): $avmultcluselemnum \n";
print $outs "\n";

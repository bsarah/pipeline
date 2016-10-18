#!/usr/bin/perl -w
## call: sortEdli edliList inpath outpath
## edliList contains the list of .edli files that need to be checked
## outfolder specifies the folder where noEdge graphs have to be moved to

use Data::Dumper;
use strict;
use warnings;


my $filename = shift;
my $inpath=shift;
my $outpath= shift;

open FA,"<$filename" or die "can't open $filename\n";
my $count = 0;

while(<FA>){
    chomp;
    my $curfile = $_;
    my $cmd = "wc -l $curfile";
    my @out = readpipe("$cmd");
    chomp(@out);
    my @F = split " ", $out[0];
    if($F[0] < 1){
	$count++;
	my $copycmd = "mv $curfile $outpath";
	readpipe("$copycmd");
    }
}

print "$count \n";

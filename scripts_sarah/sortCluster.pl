#!/usr/bin/perl -w
## call: sortCluster clusList outfolder
## clusList contains the list of cluster files that need to be checked
## outfolder specifies the folder where singleton elements have to be moved to

use Data::Dumper;
use strict;
use warnings;


my $filename = shift;
my $outfolder= shift;

open FA,"<$filename" or die "can't open $filename\n";
my $count = 0;

while(<FA>){
    chomp;
    my $curfile = $_;
    my $cmd = "wc -l $curfile";
    my @out = readpipe("$cmd");
    chomp(@out);
    my @F = split " ", $out[0];
    if($F[0] < 2){
	$count++;
	my $copycmd = "mv $curfile $outfolder";
	readpipe("$copycmd");
    }
}

print "$count \n";

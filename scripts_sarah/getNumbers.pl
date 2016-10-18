#!/usr/bin/perl -w                                                                                                                                                                                                                          
## call: joinClusters cluslist outpath

use Data::Dumper;
use strict;
use warnings;
use File::Basename;


my $filename = shift;
my $outpath = shift;

open FA,"<$filename" or die "can't open $filename\n";

my $count = 0;

while(<FA>){
    chomp;
    my $curfile = $_;

    my $name;
    my $path;
    my $suffix;
    ($prename,$path,$suffix) = fileparse($curfile);
    my @F = split '\.', $prename;
    my $name = $F[0];
    my $newname = "$name-$count\.clus";
    system("mv $curfile $path\/$newname");
    $count++;
}
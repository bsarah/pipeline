#!/usr/bin/perl -w                                                                                                                                                                                                                          
## call: getNumbers cluslist outpath summary

use Data::Dumper;
use strict;
use warnings;
use File::Basename;


my $filename = shift;
my $outpath = shift;
my $summary = shift;

open FA,"<$filename" or die "can't open $filename\n";

open(my $outs,">>$summary");

my $count = 0;
my $nonecount = 0;

while(<FA>){
    chomp;
    my $curfile = $_;
    my $none = "None";
    if($curfile =~ /$none/){$nonecount++;}
    my $name;
    my $path;
    my $suffix;
    my $prename;
    ($prename,$path,$suffix) = fileparse($curfile);
    my @F = split '\.', $prename;
    if(scalar @F >= 2){$name = $F[-2];}
    else{$name = $F[0];}
    my $newname = "$name-$count\.clus";
    system("mv $curfile $path\/$newname");
    $count++;
}


print $outs "===============Cluster information\=============== \n";
print $outs "Number of original clusters: $count \n";
print $outs "Number of 'None'-clusters: $nonecount \n"; 
print $outs "(clusters where no start or end coordinates could be detected)\n";
print $outs "\n";

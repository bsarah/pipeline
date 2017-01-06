#!/usr/bin/perl -w

# countNones.pl nonelcuslist
#count how many elements per species in none cluster and return string, separated with - and =

use Data::Dumper;
use strict;
use warnings;

my $file = shift;


my %spec2num = ();

open FA,"<$file" or die "can't open $file\n";
while(<FA>){
    chomp;
    my $curfile = $_;
    open CF,"<$curfile" or die "can't open $curfile\n";
    while(<CF>){
	chomp;
	my $line = $_;
	my @F = split '\t', $line;
	my @G = split '_', $F[1];
	my $spec = $G[0];
	if(exists $spec2num{$spec}){$spec2num{$spec}++;}
	else{$spec2num{$spec}=1;}
    }
}

my $outstr = "";
foreach my $du (sort keys %spec2num){
    $outstr = "$outstr$du\-$spec2num{$du}\=";
}

print $outstr;

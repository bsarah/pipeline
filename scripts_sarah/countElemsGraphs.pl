#!/usr/bin/perl -w

# countElemsGraphs.pl elemlist summary
#count how many elements per species in the graphs and output a summary file with it



use Data::Dumper;
use strict;
use warnings;

my $file = shift;
my $summary = shift;

open(my $outn,">>",$summary);

my %spec2num = ();

open FA,"<$file" or die "can't open $file\n";
while(<FA>){
    chomp;
    my $curfile = $_;
    my %curspec2num = (); ##report all elemenst and later divide by two as nodes are contained twice in egdelist
    my %curids = ();
    open CF,"<$curfile" or die "can't open $curfile\n";
    while(<CF>){
	chomp;
	my $line = $_;
	my @F = split " ", $line;
	my @H = split '_', $F[0];
	my @G = split '_', $F[1];
	my $spec = $G[-7];
	my $id1 = $G[-6];
	my $spec2 = $H[-7];
	my $id2 = $H[-6];
	if(exists $curids{$id1}){}
	else{
	    if(exists $curspec2num{$spec}){$curspec2num{$spec}++;}
	    else{$curspec2num{$spec}=1;}
	    $curids{$id1} = 1;
	}
	if(exists $curids{$id2}){}
	else{
	    if(exists $curspec2num{$spec2}){$curspec2num{$spec2}++;}
	    else{$curspec2num{$spec2}=1;}
	    $curids{$id2} = 1;
	}
    }
    foreach my $sp (keys %curspec2num){
	my $tmp2 = $curspec2num{$sp};	
	if(exists $spec2num{$sp}){
	    my $tmpnum = $spec2num{$sp};
	    $spec2num{$sp} = $tmpnum + $tmp2;
	}
	else{$spec2num{$sp}=$tmp2;}
    }
}


my @M = split '\/', $file;
my $realfile = $M[(scalar @M) -1];

print $outn "===============Details in $realfile\===============\n";
foreach my $du (sort keys %spec2num){
    print $outn "$du $spec2num{$du}\n";
}
print $outn "\n";

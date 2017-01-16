#!/usr/bin/perl -w

# countElems.pl elemlist pseudoscore summary
#count how many elements per species in the clusters listed in the input file cluster and output a summary file with it

use Data::Dumper;
use strict;
use warnings;

my $file = shift;
my $pseudoscore = shift;
my $summary = shift;

open(my $outn,">>",$summary);

my %spec2num = ();
my %spec2pseu = ();

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
	if($pseudoscore >= 0){
	    my $pseu = $F[scalar @F -1];
	    if($pseu < $pseudoscore){
		if(exists $spec2pseu{$spec}){$spec2pseu{$spec}++;}
		else{$spec2pseu{$spec}=1;}
	    }
	}
	if(exists $spec2num{$spec}){$spec2num{$spec}++;}
	else{$spec2num{$spec}=1;}
    }
}

my $psstr = "";
foreach my $k (keys %spec2pseu){
    $psstr = "$psstr$k\-$spec2pseu{$k}\=";
}
if($psstr eq ""){$psstr = "=";}

print $psstr;

my @H = split '\/', $file;
my $realfile = $H[(scalar @H) -1];

print $outn "===============Details in $realfile\===============\n";
foreach my $du (sort keys %spec2num){
    print $outn "$du $spec2num{$du}\n";
}
print $outn "\n";

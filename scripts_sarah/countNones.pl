#!/usr/bin/perl -w

# countNones.pl nonelcuslist mode pseudoscore
#count how many elements per species in none cluster and return string, separated with - and =

use Data::Dumper;
use strict;
use warnings;

my $file = shift;
my $mode = shift;
my $pseudoscore = shift;


my %species=();
my %pseudos = ();


open FA,"<$file" or die "can't open $file\n";
while(<FA>){
    chomp;
    my $curfile = $_;
    open CF,"<$curfile" or die "can't open $curfile\n";
    while(<CF>){
	chomp;
	my $line = $_;
	my @G = split "\t", $line;
	my $prespec = $G[1];
	my @H = split '_', $prespec;
	my $spec = $H[0];
	if(exists $species{$spec}){$species{$spec} += 1;}
	else{$species{$spec} = 1;}
	if($mode == 0){
	    my $pseudo = $G[-1];
	    if($pseudo >= 0 && $pseudoscore < $pseudo){
		if(exists($pseudos{$spec})){$pseudos{$spec}++;}
		else{$pseudos{$spec}=1;}
	    }
	}
	if($mode == 1){
	    my $state = $G[-2];
	    if($state eq "T" || $state eq "t" || $state eq "True" || $state eq "true" || $state eq "1" || $state eq "TRUE"){
		if(exists($pseudos{$spec})){$pseudos{$spec}++;}
		else{$pseudos{$spec}=1;}
	    }
	}
	
    }
}


my $outstr =  "";

##substract number of pseudogenes from number of singletons as we want to omit
##pseudogenes in the analysis

foreach my $k (sort keys %species) {
    my $newval;
    if(exists($pseudos{$k})){$newval = $species{$k}-$pseudos{$k};}
    else{$newval = $species{$k};}
    $outstr = "$outstr$k\-$newval\=";
}
$outstr = "$outstr\!";
if(scalar (keys %pseudos) == 0){
    $outstr = "$outstr\=";
}
else{
    foreach my $p (sort keys %pseudos) {
	$outstr = "$outstr$p\-$pseudos{$p}\=";
    }
}

print $outstr;

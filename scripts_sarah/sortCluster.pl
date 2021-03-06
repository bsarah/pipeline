#!/usr/bin/perl -w
## call: sortCluster clusList mode pseudoscore outfolder tmpfile
## clusList contains the list of cluster files that need to be checked
## outfolder specifies the folder where singleton elements have to be moved to

use Data::Dumper;
use strict;
use warnings;


my $filename = shift;
my $mode = shift;
my $pseudoscore = shift;
my $outfolder= shift;
my $tmpfile = shift; #write in tmpfile the outstr as it doesn't seem to work properly

open(my $outf, ">>", $tmpfile);

open FA,"<$filename" or die "can't open $filename\n";
my $count = 0;

my %species=();
my %pseudos = ();

while(<FA>){
    chomp;
    my $curfile = $_;
    my $cmd = "wc -l $curfile";
    my @out = readpipe("$cmd");
    chomp(@out);
    my @F = split " ", $out[0];
    if($F[0] < 2){
	##if there is only one cluster in the file, count how many elements 
	##for which species are contained in the singleton clusters
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
	
	$count++;
	my $copycmd = "mv $curfile $outfolder";
	readpipe("$copycmd");
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
foreach my $p (sort keys %pseudos) {
    $outstr = "$outstr$p\-$pseudos{$p}\=";
}

print $outf $outstr;

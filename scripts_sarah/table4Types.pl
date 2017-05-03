#!/usr/bin/perl -w
#call: table4Types.pl typelist table_out
#inputlist has the format species_type number
#output should be a table with types as columns and species as rows


use strict;
use warnings;
#use Pod::Usage;
use Data::Dumper;

my $infile = shift;
my $outfile = shift;

open(my $outf, ">>",$outfile);

my %elements = (); #elements get ids
my %species = (); #species get ids
my %lines = ();
my $specid = 0;
my $elid = 0;

open FA,"<$infile" or die "can't open $infile\n";
while(<FA>){
    chomp;
    my $line=$_;
    my @F = split ' ', $line;
    $lines{$F[0]} = $F[1];
    my @G = split '_', $F[0];
    my $el = $G[1];
    my $spec = $G[0];
    if(exists($species{$spec})){}
    else{$species{$spec} = $specid; $specid++;}
    if(exists($elements{$el})){}
    else{$elements{$el}=$elid;$elid++;}
}

my @types = sort (keys %elements);
my $header = "x";
for(my $i=0;$i< scalar @types; $i++){
    $header = "$header\t$types[$i]";
}
print $outf "$header\n";
foreach my $sp (keys %species){
    my $str = "$sp";
    for(my $t=0; $t< scalar @types;$t++){
	my $key = "$sp\_$types[$t]";
	my $num = 0;
	if(exists($lines{$key})){$num = $lines{$key};}
	$str = "$str\t$num";	
    }
    print $outf "$str\n";
}


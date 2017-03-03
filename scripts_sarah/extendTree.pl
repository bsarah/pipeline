#!/usr/bin/perl -w
#call: extendTree.pl treefile outtree
#return treefile with extended tree

#this program extends the tree file in case that inner node identifier are
#missing. We need dummy identifier for the inner nodes in order to be
#able to draw the tree with iTOL.




use Data::Dumper;
use strict;
use warnings;

my $treefile = shift;
my $outfile = shift;

open(my $outf, ">>", $outfile);

my $tree="";

open TF,"<$treefile" or die "can't open $treefile\n";
while(<TF>){
    chomp;
    $tree = $_;
    last;
}

if($tree eq ""){print "tree format doesn't fit!\n"; exit 1;}

my @T = split "", $tree;
my $cbrac = ")";
my $obrac = "(";
my $com = ",";
my $semco = ";";
my $space = " ";
my $count = 0;
my $newtree = "";
for(my $i  = 1; $i < scalar @T; $i++){
    if($T[$i-1] eq $cbrac){
	#check if current position is != ( ) , ;
	if($T[$i] eq $space){$i++;}
	if($T[$i] eq $obrac || $T[$i] eq $com || $T[$i] eq $semco || $T[$i] eq $cbrac){ #add dummy node
	    my $dumstr = "inode$count";
	    $newtree = "$newtree$dumstr";
	}
    }
    else{
	$newtree = "$newtree$T[$i]";
    }
}

print $outf $newtree;

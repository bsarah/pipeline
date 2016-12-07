#!/usr/bin/perl -w                                                                                                                                                                                                                          
## call: joinClusters cluslist outpath joindList

use Data::Dumper;
use strict;
use warnings;
use File::Basename;


my $filename = shift;
my $outpath = shift;
my $outfile = shift;

open(my $outf, ">>$outfile");

open FA,"<$filename" or die "can't open $filename\n";

my $curfile="";
my $curstart = 0;
my $curend = 0;
my $curnum = 0;

my @counts = ("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z");
my $ct=0;
my $num=0;
#print "$curfile \n";

while(<FA>){
    chomp;
    my $prenewfile = $_;
    my $name;
    my $path;
    my $suffix;
    ($name,$path,$suffix) = fileparse($prenewfile);
    my @F=split '-', $name;
    my $newstart = $F[1];
    my $preend = $F[2];
    my @G = split '\.', $preend; #works also fine, if there is no . in the word
    my $newend = $G[0];
    my $newnum;
    if(scalar @F <= 3){ ##if getNumbers didn't work, new numbers are assigned to the clusters
	$newnum = $num;
	$num++;
    }
    else{
	my @tmpsplit = split '\.', $F[3];
	my @numsplit = split ':', $tmpsplit[0];
	$newnum = $numsplit[(scalar @numsplit)-1]; ##take the last entry which is the current higher number
    }

    if($newstart <= $curend)#join the clusters
    {
	my $joinedfile = "$outpath\/cluster-$curstart-$newend-$curnum\:$newnum\.clus";
	##in case a filename occurs more than once (does this happen?)
	if($curfile eq $joinedfile || $prenewfile eq $joinedfile)
	{
	    $joinedfile = "$outpath\/cluster-$curstart-$newend-$curnum\:$newnum$counts[$ct]\.clus";
	    if($ct >= scalar @counts){$ct = 0;}
	    else{$ct++;}
	}
	system("cat $curfile $prenewfile > $joinedfile");
	system("rm $curfile");
	system("rm $prenewfile");
	$curfile = $joinedfile;
	$curstart = $curstart;
	$curend = $newend;
	my $outstr="join: $curnum-$newnum\n";
	print $outf $outstr;
    }
    else#go on
    {
#	print ("go on\n");
	$curfile = $prenewfile;
	$curstart = $newstart;
	$curend = $newend;
	$curnum = $newnum;
    }


}

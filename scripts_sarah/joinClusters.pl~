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

my $count=0;
my $num=0;
#print "$curfile \n";

while(<FA>){
    chomp;
    my $prenewfile = $_;
    my $name;
    my $path;
    my $suffix;
    ($name,$path,$suffix) = fileparse($prenewfile);

 #   print "$prenewfile\n";
#    my $newfile = substr($prenewfile,7,length($prenewfile)-7);
  #  print "newfile: $newfile\n";
    my @F=split '-', $name;
    my $newstart = $F[1];
    my $preend = $F[2];
    my @G = split '\.', $preend;
    my $newend = $G[0];
    my $newnum;
    if(scalar @F <= 3){
	$newnum = $num;
	$num++;
    }
    else{
	$newnum = $F[3];
    }
#    my @G = split '\.',$preend;
#    my $newnum = $G[0];
 #   print "curfiles: $curfile, $curstart, $curend\n";
 #   print "newfiles: $newfile, $newstart, $newend\n";

    if($newstart <= $curend)#join the clusters
    {
#	print("joining\n");
	my @H= split '-', $curfile;
	my @I=split '\.', $H[2];
	my $curnum = $I[0];
	my $joinedfile = "$outpath\/cluster-$curstart-$newend-joined\.clus";
	if($curfile eq $joinedfile || $prenewfile eq $joinedfile)
	{
	    $joinedfile = "$outpath\/cluster-$curstart-$newend-joined$count\.clus";
	    $count++;
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
	print $outf "nojoin: $newnum\n";
    }


}

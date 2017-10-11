#!/usr/bin/perl -w

#this program simulates a set of species and blocks in order to test the pipeline



use Data::Dumper;
use Getopt::Long qw(GetOptions);
use strict;
use warnings;



#options for simulation
my $species=6; #number of species
my $types = 5; #number of different genes
my $len = 20;  #length of genes
my $elems = 100; #total number of genes per species
my $blocks = 100000; #number of genomic anchors per species
my $breakpoints = 0; #number of translocations within blocks
my $noise = 0; #percentage of blocks to delete per species

my $outpath;


GetOptions(
    'out|o=s' => \$outpath,
    'species=f' => \$species,
    'types=f' => \$types,
    'len=f' => \$len,
    'elems=f' => \$elems,
    'blocks=f' => \$blocks,
    'breakpoints=f' => \$breakpoints,
    'noise=f' => \$noise
    ) or die "error in simulation parameters";



#outfolder
my $genespath = "$outpath\/genes";
my $tmppath = "$outpath\/temp";
my $touchcmd1 = "mkdir $genespath";
my $touchcmd2 = "mkdir $tmppath";
readpipe("$touchcmd1");
readpipe("$touchcmd2");

#store the different genetic elements
my @genes = ();
my @strucs = ();
my @kinds = ();

my @chars = ("A","C","G","T");
my @db = (")",".","(");

for(my $i = 0;$i < $types;$i++){
    my $seq;
    $seq .= $chars[rand @chars] for 1..$len;
    my $struc;
    $struc .= $db[rand @db] for 1..$len;
    my $typename = "type$i";
    push @genes, $seq;
    push @strucs, $struc;
    push @kinds, $typename;
}


#create species and blocks and elements
my @speclist = ();
my $noisenum = $blocks/100 * $noise;
for(my $j=0;$j<$species;$j++){
    my $specname = "spec$j";
    push @speclist, $specname;
    #anchors
    my $outfile = "$tmppath\/$specname\_temp.bed";
    open(my $outf,">>",$outfile);
    my @anchors = ();
    my $c=0;
    while($c < $blocks){
	my $randnum = int(rand($blocks));
	if($randnum < $noisenum){$c++;next;}
	my $line = "CHR\t$specname\_$c\tNA\tNA\tNA\tNA\tNA\n";
	print $outf $line;
	push @anchors, $c;
	$c++;
    }
    #elements
    my $genesfile = "$genespath\/$specname\.bed";
    open(my $outg,">>",$genesfile);
    my $anlen = scalar @anchors;
    my @nowblocks = ();
    my @endblocks = ();
    for(my $a=0;$a<$elems;$a++){
	my $randan = int(rand($anlen));
	my $randstart = $anchors[$randan];
	my $randend = $anchors[$randan+1];
	push @nowblocks, $randstart;
	push @endblocks, $randend;
    }
    ##sort nowblocks and include breakpoints
    my @startsort = sort { $a <=> $b } @nowblocks;
    my @endsort = sort { $a <=> $b } @endblocks;
    ##breakpoints, exchange a certain number of blocks
    for(my $b=0;$b < $breakpoints; $b++){
	my $changea = int(rand($elems));
	my $changeb = int(rand($elems));
	my $tmpst = $startsort[$changea];
	my $tmpen = $endsort[$changea];
	$startsort[$changea] = $startsort[$changeb];
	$endsort[$changea] = $endsort[$changeb];
	$startsort[$changeb] = $tmpst;
	$endsort[$changeb] = $tmpen;
    }
    ##create gene lists
    for(my $k=0;$k<$elems;$k++){
	my $start = $k*$len + 1;
	my $end = ($k+1) * $len;
	my $randseq = int(rand($types));
	my $str = "CHR\t$specname\_$k\t$start\t$end\t+\t$startsort[$k]\t$endsort[$k]\t$strucs[$randseq]\t$genes[$randseq]\t$kinds[$randseq]\tFALSE\tNA\n";
	print $outg $str;
    }
    close $outg;
    close $outf;
    my $gzcmd = "gzip $outfile";
    readpipe("$gzcmd");
}

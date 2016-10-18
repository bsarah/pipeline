#!/usr/bin/perl -w

## perl doSummary.pl summary_outfile cographs_file noncographs_file specieslist_file optstr numclus numjoinclus numsingle numnoedge numgraphs numnoneclus numrealgraphs outpath

use Data::Dumper;
use strict;
use warnings;
use List::Util qw(any none first);
#use Time::localtime;
use POSIX qw(strftime);

my $outfile = shift;
my $cographs = shift;
my $noncographs=shift;
my $specieslist = shift;
#my $optstr = shift;
my $numclus = shift;
my $numjoin = shift;
my $numsingle = shift;
my $numnoedge = shift;
my $numgraphs = shift;
my $numnone = shift;
my $numrealgraphs = shift;
my $outpath = shift;

open(my $outs,">>$outfile");

open CG,"<$cographs" or die "can't open $cographs\n";
open NG,"<$noncographs" or die "can't open $noncographs\n";
open SL,"<$specieslist" or die "can't open $specieslist\n";


my $numncg = 0;
my $numecg = 0;
my $numcg=0;
my $numcli = 0;
my $denscg = 0;
my $maxnumncg = 0;
my $minnumncg = 555;
my $maxnumecg = 0;
my $minnumecg = 555;

while(<CG>){
    chomp;
    my $line=$_;
    my @C = split '\t', $line;
    my $numn = $C[1];
    my $nume = $C[2];
    my $dens = $C[3];
    my $cliquenum=(($numn-1)*$numn)/2;
    if($nume == $cliquenum){
	$numcli++;
    }
    if($numn > $maxnumncg){$maxnumncg = $numn;}
    if($nume > $maxnumecg){$maxnumecg = $nume;}
    if($numn < $minnumncg){$minnumncg = $numn;}
    if($nume < $minnumecg){$minnumecg = $nume;}
    
    $numncg = $numncg + $numn;
    $numecg = $numecg + $nume;
    $denscg = $denscg + $dens;
    $numcg++;
}

if($minnumncg == 555){$minnumncg = 0}
if($minnumecg == 555){$minnumecg = 0}


my $avnumncg = 0;
my $avnumecg = 0;
my $avdenscg = 0;
if($numcg > 0){
    $avnumncg = $numncg/$numcg;
    $avnumecg = $numecg/$numcg;
    $avdenscg = $denscg/$numcg;
}

my $numnng = 0;
my $numeng = 0;
my $numng=0;
my $densng = 0;
my $maxnumnng = 0;
my $minnumnng = 555;
my $maxnumeng = 0;
my $minnumeng = 555;


while(<NG>){
    chomp;
    my $line=$_;
    my @C = split '\t', $line;
    my $numn = $C[1];
    my $nume = $C[2];
    my $dens = $C[3];

    if($numn > $maxnumnng){$maxnumnng = $numn;}
    if($nume > $maxnumeng){$maxnumeng = $nume;}
    if($numn < $minnumnng){$minnumnng = $numn;}
    if($nume < $minnumeng){$minnumeng = $nume;}

    $numnng = $numnng + $numn;
    $numeng = $numeng + $nume;
    $densng = $densng + $dens;
    $numng++;
}

if($minnumnng == 555){$minnumnng = 0}
if($minnumeng == 555){$minnumeng = 0}


my $avnumnng = 0;
my $avnumeng = 0;
my $avdensng = 0;

if($numng > 0){
    $avnumnng = $numnng/$numng;
    $avnumeng = $numeng/$numng;
    $avdensng = $densng/$numng;
}

my @species=();
my @genenum = ();
my $spcount = 0;

while(<SL>){
    chomp;
    my $line = $_;
    my @F = split '\/', $line;
    my $spec = $F[(scalar @F)-1];
    push @species, $spec;
    $spcount++;
    my $cmd = "wc -l $line";
    my @out = readpipe("$cmd");
    my @tmp = split " ", $out[0];
    my $numGenes = $tmp[0];
    push @genenum, $numGenes;
}

my @perm=();
for(my $i=0;$i<scalar @species;$i++){
    my $tmpstr = "$species[$i] ($genenum[$i])\n";
    push @perm, $tmpstr;
}

my $spstr = join(" ",@perm);




my $citation = "Beerinformatics Leipzig";
my $contact = "bsarah at bioinf dot uni-leipzig dot de";

##TODO write readme for each output folder


#my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
#print $outs "File created on $now_string\n";
#print $outs "\n\n";
#print $outs "======Output statistics======\n";
#print $outs "Program call: $optstr \n";
print $outs "Species information:\n";
print $outs "Number of Species: $spcount 
Species (number in brackets specifies the number of genetic elements used for 
this species):
 $spstr \n";
print $outs "Files including the genetic elements are located in 
$outpath\/genes \n";
print $outs "\n";
print $outs "The bed file allClusters.bed containing the genes sorted into 
the clusters is located here: $outpath \n";
print $outs "Format bed file (tab separated): 
  >cluster_name(including species and coordinates) number of elements
  for each element one line consisting of: chromosome, species, start, end, 
  strand, 5' blocknumber, 3'blocknumber \n";
print $outs "\n";
print $outs "Cluster information: \n";
print $outs "Number of original clusters: $numclus \n";
print $outs "Number of clusters where no start or end coordinates could be 
detected: $numnone \n";
my $totclus = $numjoin + $numsingle;
print $outs "Number of clusters after joining nested ones: $totclus \n";
if($totclus < 10){print $outs "Number of clusters is very small. In order to 
include more genes, you can lower the infernal threshold by changing the input
 parameter and restart the prgram. See help page for more information. \n";}
print $outs "Out of this multigene and singleton clusters: $numjoin and $numsingle \n";
print $outs "Cluster files (.clus) are located in $outpath\/cluster with 
subfolders NoneCluster and singletons for clusters without coordinate 
information or only one single element. \n";
print $outs "Cluster file format (tab separated): Chromosome, Species_geneID, 
start, end, strand, 5'blocknum, 3'blocknum, secondary structure, sequence, 
(intron1 start,intron1 end);(intron2 start, intron2 end) or (0,0);\n";
print $outs "\n";
print $outs "Graph information: \n";
print $outs "Total number of graphs: $numgraphs \n";
#my $totnoedge = $numsingle;
print $outs "Number of graphs with edges: $numrealgraphs \n";
print $outs "Number of graphs without edges: $numnoedge \n";
print $outs "In the graphs without edges, no edges were added because of 
either not enough nodes (thus, additionally $numsingle singleton clusters) or all nodes
 belonging to the same species or the similarity thresholds were not fulfiled.
 Use input parameters -s and/or -t to change the similarity thresholds. \n";
print $outs "Number of cographs: $numcg \n";
print $outs "Number of cliques (included in cographs): $numcli \n";
print $outs "Number of noncographs: $numng \n";
print $outs "Cographs: 
  average node number: $avnumncg;
  average edge number: $avnumecg; 
  Max\/min number of nodes: $maxnumncg\/$minnumncg; 
  Max\/min number of edges: $maxnumecg\/$minnumecg; 
  average density:  $avdenscg \n";
print $outs "Noncographs: 
  average node number: $avnumnng; 
  average edge number: $avnumeng; 
  Max\/min number of nodes: $maxnumnng\/$minnumnng; 
  Max\/min number of edges: $maxnumeng\/$minnumeng; 
  average density:  $avdensng \n";
print $outs "Graph files (.edli) are located in $outpath\/graphs 
which contains a file listing graphs without any edges. \n";
print $outs "Graph files are edgelists whereas an edge (a,b) appears as (a,b) 
and (b,a). As graphs are weighted (similarity threshold) the graph files are 
tab separated files containing: node1, node2, sequences similarity, structure 
similarity.\n";
print $outs "If the number of graphs with edges is 0, no sequences fulfilled 
the similarity thresholds for sequence and/or structure. In this way, no graphs
 were drawn or alignments were built. If you want to get graphs with edges, 
similarity thresholds can be lowered using parameters -s and/or -t. The default
 values are 0.9, meaning a similarity of 90%.\n";
if($numrealgraphs > 0){
    print $outs "Graph files used to visualize graphs (.gr) and their 
visualization (.pdf) are located in 
$outpath\/graphs/showGraphs \n";
    print $outs "\n";
    print $outs "Duplication alignments and genetic events information: \n";
    print $outs "The results for the analysis of genetic events are written to
 $outpath\/geneticEvents.txt\
The summary file contains information about how many genetic events were 
counted in a certain combination of species. This can be used to draw a 
phylogenetic tree with genetic events at its nodes. \n";
    print $outs "The files containing the duplications alignments (.aln) for 
each cluster and pairs of species can be found here: 
$outpath\/graphs/alignments \n";
    print $outs "Format of alignment files (.aln): At first the cluster and 
its species are defined. Then, connected nodes get mapped to the same 
one-letter-code (needed for the alignment). Afterwards for each species, its 
genetic elements are sorted by coordinate and depicted by the letter code. The
 letter code is aligned to the ones of each other species in the cluster. '~' 
stand for duplications, '-' for insertions or deletions in the alignment. \n";
}
print $outs "\n";
print $outs "===In case of any errors..===\n";
print $outs "..please first check the output of the program and the file 
errors.txt in $outpath. 
Errors can be caused by file formats that are not 
accepted by the program, names that are not unique or just missing data.\n";
print $outs "\n\n";
print $outs "======Further information======\n";
print $outs "Further information can be found in the README files contained in
 the provided folders or in the corresponding publication.\n";
print $outs "Further analysis can be done ...TODO \n";
print $outs "\n";
print $outs "If you use this data, please cite:\n  $citation \n";
print $outs "If you have any further questions please write to:\n  $contact \n";


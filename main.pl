#!/usr/bin/perl -w


##to write the perl main program use pero getopts package!

#requirements: species.bed files, scripts that will be used later, list of species with exactly the same names as the .bed files
#get path to Output folder (there has to be Cameron's genes folder in there)

#path to altnw: /homes/biertank/bsarah/Documents/drosophilas/tRNAevoGit/tRNAevo/alternativeNW/dist/build/altNW/

#inp: path to Output folder
#inp: cm or gene list
#inp: path to infernal
#inp: path to Genomes
#inp: path to Maf files
#inp: name of ref genome
#inp (optional): similarity threshold for edges in orthology graph for sequence and secondary structure!!!
##note: if one of the similarity values is not important, user should set it to 0
# python3 must be installed, inp: path to python3
# >= perl5.22 needed
# R >= 3.2 needed, libraries qgraph, igraph
# Infernal >= 1.1.1

#options -o -g -m -f -s -t -c -l -q -w -b -y -e -r -i -h -v -a -x -z -n -d
#-output folder
#-genomes folder
#-maf folder
#-refgenome name
#-seq sim number (optional, default 0.9)
#-struc sim number (optional, default 0.9)
#-cm path to file (either -cm or -genelist)
#-genelist path to file (either -cm or -genelist)
#-percentage of lowest scoring maf blocks number (cams program) (optional)
#-w sets the evalue threshold for infernal (default 0.01)
#-b sets the bitscore threshold for infernal (default: unset)
#-python path
#-perl path
#-r path
#-infernal path
#-help
#-version & citation
#-again (run the program again but without creating the blocks (without camerons part)
#-x do NOT check graphs for cographs
#-z do NOT create alignments
#-u -v are sequence and structure thresholds for pseudogenes (thus their range is between -u/-v and -s/-t); default???? 0.7?
#-n is the newick tree file
#-d is the id translation between species names in newick tree and in the maf alignments (not needed if the species names in the tree already got the maf names)

#k61 perl: /opt/localperl/bin/
#

use Data::Dumper;
use strict;
use warnings;
#use List::Util qw(any none first);
use Getopt::Std;
use POSIX qw(strftime);
use File::Basename;
use File::Find;

# declare the perl command line flags/options we want to allow
my %options=();
getopts("hVKo:g:m:f:s:t:c:l:q:w:b:a:y:e:r:i:xzu:v:n:d:", \%options);


##define variables
my $outpath;
my $genomes="";
my $mafs="";
my $refspecies="";
my $seqsim;
my $strucsim;
my $pseqsim;
my $pstrucsim;
my $cmfile="";
my $genefile="";
my $cmoption;
my $perc="";
my $pathtocam; #path to output of cams program (genes folder)
my $pythonpath;
my $evalin=0.01;
my $perlpath;
my $rpath;
my $infernalpath;
my $extravalue;
my $newicktree;
my $ids;


my $dirname = dirname(__FILE__);

#test print
#print "current directory of main.pl: $dirname \n";


#define other information
my $toolname = "FindAFancyAbbrevation";
my $scripts_cam = "$dirname/scripts_cam";
my $scripts_sarah = "$dirname/scripts_sarah";
my $version = "0.1";
my $citation = "Beerinformatics Leipzig";
my $contact = "bsarah at bioinf dot uni-leipzig dot de or anneh at bioinf.uni-leipzig dot de";
my $altnwpath = "$dirname/scripts_sarah";
##test if it work when altNW executable is in scripts folder. (altNW is only for 64 bit!)
#"/scr/gin/bsarah/software/alternativeNW/dist/build/altNW/";


my $checkgraphs =1;
my $createalns = 1;


##Help page
if ($options{h})
{
    my $helpstr = createHelp($toolname,$version,$citation,$contact);
    print $helpstr;
    exit 0;
}
##Contact
if ($options{K})
{
    print "If you use this program, please cite: $citation \n";
    print "If you have further questions, please contact $contact \n";
    exit 0;
}
##Version
if ($options{V})
{
    print "Program version $version \n";
    exit 0;
}


##outpath
my $outpathstr = "";
if ($options{o}){$outpath = $options{o}; $outpathstr = "-o $outpath";}
else{print "No output path given (option -o)!\n"; exit 1;}
##if output folder doesnt exist, create it!
if(-e $outpath){}
else{
    my $cmd42= "mkdir $outpath";
    my @out42 = readpipe("$cmd42");
}

##Genomes
my $genomesstr="";
if ($options{g}){$genomes = $options{g}; $genomesstr = "-e genomes";
		 #check if folder is not empty
		 if(-e $genomes){} else{print "Genomes folder is empty! (option -g)\n"; exit 1;}
}
elsif($options{a}){}
else{print "No path to genomes given! (option -g)\n"; exit 1;}


##Maf
my $mafstr = "";
if ($options{m}){$mafs = $options{m};$mafstr = "-m $mafs";}
elsif($options{a}){}
else{print "No path to maf files given! (option -m)\n"; exit 1;}


##reference species
my $refspeciesstr = "";
if ($options{f}){$refspecies = $options{f};$refspeciesstr="-f $refspecies";
		 #check if file exists
		 my $n1 = "$refspecies\.fa";
		 my $n2 = "$refspecies\.fasta";
		 my $n3 = "$refspecies\.fa\.gz";
		 my $n4 = "$refspecies\.fasta\.gz";
		 if(-f "$genomes\/$n1"){}
		 elsif(-f "$genomes\/$n2"){}
		 elsif(-f "$genomes\/$n3"){}
		 elsif(-f "$genomes\/$n4"){}
		 elsif(-f "$genomes\/$refspecies"){}
		 else{print "No genome file for reference species in $genomes (option -f)\n"; exit 1;}
}
elsif($options{a}){}
else{print "No references species given! (option -f)\n"; exit 1;}


##similarity parameter
if ($options{s}){$seqsim = $options{s};}
else{$seqsim=0.9;}
if ($options{t}){$strucsim = $options{t};}
else{$strucsim=0.9;}
if($seqsim == -1 && $strucsim == -1){
    print "Similarity thresholds for both, sequence and structure is -1! At least one has to be > 0! (Use parameter -s for sequence and -t for structure to specify the value) \n"; exit 1;
}
if($seqsim==-1){$pseqsim = -1;}
elsif ($options{u}){$pseqsim = $options{u};}
else{$pseqsim=0.7;}
if($strucsim==-1){$pstrucsim = -1;}
elsif($options{v}){$pstrucsim = $options{v};}
else{$pstrucsim=0.7;}
if($pseqsim == -1 && $pstrucsim == -1){
    #then, pseudogenes are NOT included, thus no range between seqsim and pseqsim or strucsim and pstrucsim
    $pstrucsim=$strucsim;
    $pseqsim = $seqsim;
}


##CM
my $cmoptstr="";
my $usingCM = 1; ##take care for the option if either cm is used or a gene list is given. this changes the result page later on
if($options{c}){$cmfile = $options{c};$cmoption = "-sg $cmfile"; $cmoptstr = "-c $cmfile";}
elsif($options{l}){$genefile = $options{l};$cmoption = "-og $genefile"; $cmoptstr = "-l $genefile"; $usingCM = 0;}
elsif($options{a}){}
else{print "Please give either a cm file (option -c) or a list of genes (option -l)! \n"; exit 1;}


##Percentage of the lowest scoring blocks based on the MAF scores to be removed (between 0 and 100)
if ($options{q}){$perc = "-q $options{q}";}
else{$perc = "";}


##infernal parameter
my $inclopt="";
if ($options{w}){$inclopt = "-incE $options{w}";}
if ($options{b}){$inclopt = "-incT $options{b}";}


##run again with genes folder already created
my $doitagainstr = "";
my $doitagain="";
if ($options{a}){
    $doitagain = "$options{a}";
    $genefile="$options{a}";
    $usingCM=0;
    $doitagainstr="-a $doitagain";
    if(-e $doitagain){}
    else{print "Option -a given but argument folder $doitagain doesn't exist! \n"; exit 1;}
}

##Pathes
##Python
if ($options{y}){$pythonpath = $options{y};}
else{print "No path to python3 given! (option -y)\n"; exit 1;}

##Perl
if ($options{e}){$perlpath = $options{e};}
else{print "No path to perl given! (option -e)\n"; exit 1;}

##R (not needed anymore)
if ($options{r}){$rpath = $options{r};}
else{print "No path to R given! (option -r)\n"; exit 1;}

#Infernal, only needed when option -c
if ($options{i}){$infernalpath = $options{i};}
elsif($options{l}){$infernalpath = "";}
else{print "No path to infernal given! (option -i)\n"; exit 1;}

if ($options{x}){$checkgraphs = 0;}
if ($options{z}){$createalns = 0;}

my $newickstr = "";
if($options{n} && $createalns == 1){$newicktree = $options{n};$newickstr = "-n $newicktree";}
elsif((! $options{n}) && $createalns == 1){print "Please give a newick tree as input file (option -n)\n"; exit 1;}
else{}
my $idstr = "";
if($options{d} && $createalns == 1){$ids = $options{d};$idstr = "-d $ids";}
elsif((! $options{d}) && $createalns == 1){print "No file given to translate newick IDs into MAF IDs (option -d). Please make sure that IDs fit!\n";}
else{}

##Program call
my $optstr = "$outpathstr $genomesstr $mafstr $refspeciesstr $newickstr $idstr -s $seqsim -t $strucsim $doitagainstr $cmoptstr $perc $inclopt -y $pythonpath -e $perlpath -r $rpath -i $infernalpath";
print "program called with: $optstr \n";


##DEBUG
##pipe all possible errors in errorfile
system("touch $outpath\/errorsPart2.txt");
my $err = "$outpath\/errorsPart2.txt";
system("touch $outpath\/errorsPart1.txt");
my $err0 = "$outpath\/errorsPart1.txt";
##other debug prints
system("touch $outpath\/debug.txt");
my $db = "$outpath\/debug.txt";




##start the program
my $start_string = strftime "%a %b %e %H:%M:%S %Y", localtime;


#write summary file with analysis of cographs and noncographs and the text below
my $cmd33 = "touch $outpath\/summary.txt 2>>$err";
my @out33 = readpipe("$cmd33");
my $sumfile = "$outpath\/summary.txt";
open(my $outs,">>",$sumfile);
print $outs "Program started on $start_string\n";
print $outs "Program call: $optstr \n";
print $outs "\n\n";
print $outs "======Output statistics======\n";
print $outs "\n";


##Cam's part (sort the genetic elements into genomic anchors based on the maf blocks)
my $genesfolder="";
if(! $options{a}){
    print "analysis of maf files started (this might take a while)..\n";
    if($perc eq ""){
	open(PROG,"$pythonpath\/python3 $scripts_cam\/main.py $cmoption $inclopt $genomes $mafs $infernalpath $outpath $dirname $refspecies 2>>$err0 |") or die "Couldn't start prog!";
	while(<PROG>){print "$_";}
    }
    else{
	open(PROG,"$pythonpath\/python3 $scripts_cam\/main.py $perc $cmoption $inclopt $genomes $mafs $infernalpath $outpath $dirname $refspecies 2>>$err0 |") or die "Couldn't start prog!";
	while(<PROG>){print "$_";}
    }
        print "Done!\n";
    $genesfolder = "$outpath\/genes";
}
else{
    $genesfolder = $doitagain;
}


##create summaries folder to collect summaries and lateron summarize them
my $cmdsummary = "mkdir $outpath\/summaries 2>>$err";
my @outsummary = readpipe("$cmdsummary");

my $summarypath = "$outpath\/summaries";

my $sumcollectcluster = "$summarypath\/Summary_collectCluster.txt";
my $sumgetnumbers = "$summarypath\/Summary_getNumbers.txt";

##Construct clusters
my $cmd1 = "mkdir $outpath\/clusters 2>>$err";
my $cmd2 = "ls $genesfolder\/*.bed \> $genesfolder\/specieslist 2>>$err";
my $cmd3 = "$perlpath\/perl $scripts_sarah\/collectCluster\.pl $genesfolder\/specieslist $genesfolder $outpath\/clusters $sumcollectcluster 2>>$err";
my $cmd4 = "ls $outpath\/clusters/*.clus > $outpath\/clusters/precluslist 2>>$err";
my $cmd5 = "$perlpath/perl $scripts_sarah\/getNumbers.pl $outpath\/clusters/precluslist $outpath\/clusters $sumgetnumbers 2>>$err";
print "construct clusters..";
my @out1 = readpipe("$cmd1");
my @out2 = readpipe("$cmd2");
my @out3 = readpipe("$cmd3");
my @out4 = readpipe("$cmd4");
my @out5 = readpipe("$cmd5");

append2file($outs,$sumcollectcluster);
append2file($outs,$sumgetnumbers);
print "Done!\n";

#sort clusters without specific coordinates (none cluster)
my $cmd6 = "mkdir $outpath\/clusters\/NoneCluster 2>>$err";
my $cmd7 = "mv $outpath\/clusters\/cluster-None\* $outpath\/clusters\/NoneCluster 2>>$err";
my $cmd7a = "touch $outpath\/clusters\/NoneCluster/nonecluslist 2>>$err";;
my $cmd7b = "ls $outpath\/clusters\/NoneCluster/*.clus > $outpath\/clusters\/NoneCluster/nonecluslist 2>>$err";
##report the number of elements for each species that are contained in none clusters

print "sort clusters..";
my @out6 = readpipe("$cmd6");
my @out7 = readpipe("$cmd7");
my @out7a = readpipe("$cmd7a");
my @out7b = readpipe("$cmd7b");
print "Done!\n";

#create bedfile about clusters (without none clusters)
my $sumallclusters = "$summarypath\/Summary_allClusters.txt";
my $outname = "allClusters.bed";
my $cmd8 = "ls $outpath\/clusters\/\*\.clus \> $outpath\/clusters\/clusList 2>>$err";
my $cmd9 = "$perlpath\/perl $scripts_sarah\/writeBED\.pl $outpath\/clusters\/clusList $outpath $outname $sumallclusters 2>>$err";
print "create BED file..";
my @out8 = readpipe("$cmd8");
my @out9 = readpipe("$cmd9");
append2file($outs,$sumallclusters);
print "Done!\n";

#sort list of clusters by start coordinate
my $cmd10 = "$perlpath\/perl -pi.bak -e 's/-/ /g' $outpath\/clusters/clusList 2>>$err";
my $cmd11 = "sort -k2 -n $outpath\/clusters/clusList > $outpath\/clusters/clusList_sorted 2>>$err";
my $cmd12 = "$perlpath\/perl -pi.bak -e 's/ /-/g' $outpath\/clusters/clusList_sorted 2>>$err";

#join clusters and create list about which clusters were joined
my $cmd13 = "$perlpath\/perl $scripts_sarah\/joinClusters.pl $outpath\/clusters/clusList_sorted $outpath\/clusters $outpath\/clusters/joinlist 2>>$err";
#create list or current joined clusters
my $cmd14 = "ls $outpath\/clusters/*.clus > $outpath\/clusters/cluslist_joined 2>>$err";

#create bedfile about joined clusters (without none clusters)
my $sumallclustersjoined = "$summarypath\/Summary_allClusters_joined.txt";
my $outname2 = "allClusters_joined.bed";
my $cmd9a = "$perlpath\/perl $scripts_sarah\/writeBED\.pl $outpath\/clusters\/cluslist_joined $outpath $outname2 $sumallclustersjoined 2>>$err";

#sort out singletons
my $cmd15 = "mkdir $outpath\/clusters/singletons 2>>$err";
my $cmd16 = "$perlpath\/perl $scripts_sarah\/sortCluster.pl $outpath\/clusters/cluslist_joined $outpath\/clusters/singletons 2>>$err";
my $cmd17 = "ls $outpath\/clusters/*.clus > $outpath\/clusters/cluslist_nosingles 2>>$err";

print "join clusters..";
my @out10 = readpipe("$cmd10");
my @out11 = readpipe("$cmd11");
my @out12 = readpipe("$cmd12");
my @out13 = readpipe("$cmd13");
my @out14 = readpipe("$cmd14");
my @out9a = readpipe("$cmd9a");
my @out15 = readpipe("$cmd15");
my @out16 = readpipe("$cmd16");
my @out17 = readpipe("$cmd17");
append2file($outs,$sumallclustersjoined);
print "Done!\n";

##this string contains $species\-$num_singletons\=$species\-...has to be separated like this in order to be able to hand it over as a parameter
##give this to create alignments in order to add the singleton count to the genetic events list (as insertions)

my $singletoncount;
if($out16[0] eq ""){$singletoncount = "=";}
else{$singletoncount = "$out16[0]";}
print "singletoncount: $singletoncount \n";


#create graphs
my $sumbuildedges = "$summarypath\/Summary_buildedges.txt";
my $cmd18 = "mkdir $outpath\/graphs 2>>$err";
my $cmd19 = "$perlpath\/perl $scripts_sarah\/buildEdgeList.pl $outpath\/clusters/cluslist_nosingles $outpath\/clusters $outpath\/graphs $altnwpath $seqsim $strucsim $sumbuildedges 2>>$err";
my $cmd20 = "ls $outpath\/graphs/*.edli > $outpath\/graphs/edlilist 2>>$err";
##no edge graphs are graphs with node from only one species, as all other graphs have a completely connected graph (except same species)
#my $cmd21 = "mkdir $outpath\/graphs/noEdgeGraphs 2>>$err";   
#my $cmd22 = "$perlpath\/perl $scripts_sarah\/sortEdli.pl $outpath\/graphs/edlilist $outpath\/graphs $outpath\/graphs/noEdgeGraphs 2>>$err";
#my $cmd22 = "ls $outpath\/graphs/*.edli > $outpath\/graphs/edlilist_edgegraphs 2>>$err";
print "create graphs (this might take a while)..";
my @out18 = readpipe("$cmd18");
my @out19 = readpipe("$cmd19");
my @out20 = readpipe("$cmd20");
append2file($outs,$sumbuildedges);
#my @out21 = readpipe("$cmd21");
#my @out22 = readpipe("$cmd22");
#my @out22 = readpipe("$cmd22");
print "Done!\n";


if($checkgraphs == 1){
    ##check graph structure
    my $sumcheckgraph = "$summarypath\/Summary_checkgraph.txt";
    my $cmd24 = "touch $outpath\/graphs/cographs 2>>$err";
    my $cmd23 = "touch $outpath\/graphs/list-noEdgeGraphs.txt 2>>$err";
    my $cmd23a = "touch $outpath\/graphs/list-EdgeGraphs.txt 2>>$err";
    my $cmd25 = "touch $outpath\/graphs/noncographs 2>>$err";
    my $cmd26 = "mkdir $outpath\/graphs/showGraphs 2>>$err";
    ##TODO set the pseqsim and pstruclim if we have a solution for pseudogenes
    my $cmd27 = "$perlpath\/perl $scripts_sarah\/checkGraph.pl $outpath\/graphs/edlilist $outpath\/graphs $outpath\/graphs/showGraphs $seqsim $strucsim -1 -1 $outpath\/graphs/cographs $outpath\/graphs/noncographs $outpath\/graphs/list-noEdgeGraphs.txt $outpath\/graphs/list-EdgeGraphs.txt $sumcheckgraph >>$db 2>>$err";
    print "analyse graphs..";
    my @out23 = readpipe("$cmd23");
    my @out23a = readpipe("$cmd23a");
    my @out24 = readpipe("$cmd24");
    my @out25 = readpipe("$cmd25");
    my @out26 = readpipe("$cmd26");
    my @out27 = readpipe("$cmd27");
    append2file($outs,$sumcheckgraph);
    print "Done!\n";
}
#not needed at the moment
#my $cmd31 = "ls $outpath\/graphs/showGraphs/*.gr > $outpath\/graphs/showGraphs/graphsToDraw 2>>$err";
#my @out31 = readpipe("$cmd31");

if($createalns == 1){
    #create duplication alignments for each graph, thus take care for the similarity thresholds
    my $sumcreatealn = "$summarypath\/Summary_createAlignments.txt";
    my $cmd28 = "mkdir $outpath\/graphs/alignments 2>>$err";
    my $cmd29 = "touch $outpath\/geneticEvents.txt 2>>$err";
    my $cmd29a = "mkdir $outpath\/graphs/GainLoss 2>>$err";
    #my $cmd30 = "$perlpath\/perl $scripts_sarah\/getDuplication.pl $outpath\/graphs/showGraphs/graphsToDraw $outpath\/graphs/showGraphs $outpath\/graphs/alignments $altnwpath $outpath\/geneticEvents.txt 2>>$err";
    ##TODO set the pseqsim and pstruclim if we have a solution for pseudogenes
    my $cmd30 = "$perlpath\/perl $scripts_sarah\/createAlignments.pl $outpath\/graphs/edlilist $outpath\/graphs/alignments $altnwpath $seqsim $strucsim -1 -1 $singletoncount $outpath\/graphs/GainLoss $outpath\/geneticEvents.txt $sumcreatealn 2>>$err";
    
    print "create duplication alignments..";
    print "cmd createalns: $cmd30 \n";
    my @out28 = readpipe("$cmd28");
    my @out29 = readpipe("$cmd29");
    my @out29a = readpipe("$cmd29a");
    my @out30 = readpipe("$cmd30");
    append2file($outs,$sumcreatealn);
    print "Done!\n";
}



##commented draw graphs as this should go in an extra script to only print the graphs that the user wants to print

#draw graphs
##packages required: qgraph, igraph

##pipe R output to another file
#system("touch $outpath\/drawGraphs.Rout");
#my $rout = "$outpath\/drawGraphs.Rout";


#my $numrealgraphs;
#my $realgraphspath = "$outpath\/graphs/showGraphs/graphsToDraw";
#if(-z $realgraphspath){ $numrealgraphs=0;}
#else{
#    print "draw graphs..";
#    open CF,"<$outpath\/graphs/showGraphs/graphsToDraw" or die "can't open $outpath\/graphs/showGraphs/graphsToDraw\n";
#    while(<CF>){
#	chomp;
#	my $line = $_;
#	my $cmd32 = "Rscript --vanilla $scripts_sarah\/drawGraphs.R $line 2>>$rout";
#	my @out32 = readpipe("$cmd32");
#    }
#
#    print "Done!\n";
#    my $cmd38b = "wc -l $outpath\/graphs/showGraphs/graphsToDraw 2>>$err";
#    my @out38b = readpipe("$cmd38b");

#    my @tmp38b = split " ", $out38b[0];
#    $numrealgraphs = $tmp38b[0]; 
    
#}

#print "write summary..";
#my $cmd39 = "wc -l $outpath\/clusters\/clusList 2>>$err";
#my $cmd35 = "wc -l $outpath\/clusters/cluslist_nosingles 2>>$err";
#my $cmd36 = "ls $outpath\/clusters/singletons/*.clus > $outpath\/clusters/singletons/list_singletons 2>>$err";
#my $cmd36a = "wc -l $outpath\/clusters/singletons/list_singletons 2>>$err";
##my $cmd37 = "ls $outpath\/graphs/noEdgeGraphs/*.edli > $outpath\/graphs/noEdgeGraphs/list_noedgegraphs 2>>$err";
##my $cmd37a = "wc -l $outpath\/graphs/noEdgeGraphs/list_noedgegraphs 2>>$err";
#my $cmd38 = "wc -l $outpath\/graphs/edlilist 2>>$err";
#my $cmd38a = "wc -l $outpath\/clusters\/NoneCluster/nonecluslist 2>>$err";
#my $cmd39a = "wc -l $outpath\/graphs/list-noEdgeGraphs.txt 2>>$err";

#my @out39 = readpipe("$cmd39");
#my @out39a = readpipe("$cmd39a");
#my @out35 = readpipe("$cmd35");
#my @out36 = readpipe("$cmd36");
#my @out36a = readpipe("$cmd36a");
##my @out37 = readpipe("$cmd37");
##my @out37a = readpipe("$cmd37a");
#my @out38 = readpipe("$cmd38");
#my @out38a = readpipe("$cmd38a");

#my @tmp39 = split " ", $out39[0];
#my $numClus = $tmp39[0]; 

#my @tmp35 = split " ", $out35[0];
#my $numJoinClus = $tmp35[0]; 

#my @tmp36 = split " ", $out36a[0];
#my $numSingles = $tmp36[0]; 

#my @tmp39a = split " ", $out39a[0];
#my $numnoEdgeGr = $tmp39a[0]; 

#my @tmp38 = split " ", $out38[0];
#my $numGraphs = $tmp38[0]; 

#my @tmp38a = split " ", $out38a[0];
#my $numnoneclus = $tmp38a[0]; 

#if($genefile eq ""){$genefile = ".";}
#my $cmd34 = "$perlpath\/perl $scripts_sarah\/doSummary.pl $outpath\/summary.txt $outpath\/graphs/cographs $outpath\/graphs/noncographs $genesfolder\/specieslist $numClus $numJoinClus $numSingles $numnoEdgeGr $numGraphs $numnoneclus $numrealgraphs $usingCM $genefile $outpath 2>>$err";
#print "doSummary: $cmd34 \n";
#my @out34 = readpipe("$cmd34");
#print "Done!\n";


print $outs "===END===\n";
my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
print $outs "Program finished on $now_string\n";
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
#print $outs "Further analysis can be done ...TODO \n";
print $outs "\n";
print $outs "If you use this data, please cite:\n  $citation \n";
print $outs "If you have any further questions please write to:\n  $contact \n";

close($outs);


#print summary and where which result is located
print "\nEND!\n";
print "Summary and further explanations can be found here: $outpath\/summary.txt \n\n";
print "If you use this data, please cite:\n  $citation \n";
print "If you have any further questions please write to:\n  $contact \n";


##TODO further analysis???
##TODO version citation help page....


sub append2file{

    my @inp=@_;
    my $filehand = $inp[0];
    my $file = $inp[1];
    #OPEN FILE B.txt for READING (CHECK FOR FAILURES)
    open ( FL, "<", $file ) 
	or die "Could not open file $file \n";
    while ( my $line = <FL> ) {
	print $filehand $line;
    }
}


sub createHelp{
    my @inp=@_;
    my $toolname = $inp[0];
    my $version = $inp[1];
    my $citation = $inp[2];
    my $contact = $inp[3];

    
my $outstr = "";
    
    $outstr = "$outstr \nThis is the help page of $toolname.
Use this is a tool to detect homologous genetic elements in different species 
based on a model or a lit of the genetic elements and a multiple sequence 
alignment of the target species. You can choose if secondary structure 
information should be included in the analysis or not. \n";
    $outstr = "$outstr \n";
    $outstr = "$outstr Program use:
perl main.pl -o OutputFolder -g GenomesFolder -m MafFolder -f refspecies -n newicktree -d IDtranslation
{-c CM | -l genelist} -y python -e perl -r R -i infernal
[-s seqSimilarity] [-t strucSimilarity] [-u pseudoseqSimilarity] [-v pseudostrucSimilarity] 
[-q threshold] [-a genesFolder] [-x] [-z]  [-w incE] [-b incT] \n\n";
    $outstr = "$outstr Parameter explanations:\n";
    $outstr = "$outstr -o Outputfolder:
please specify where the Output files should be created. 
The folder has to exist. \n";
    $outstr = "$outstr -g GenomesFolder:
path to folder where genomes are located. The names of the genome files must 
be the same as the species' name in the MSA. The genome files should have 
fasta format and can be gzipped. Not required if option -a is taken. \n";
    $outstr = "$outstr -m MafFolder:
path to the folder where multiple sequence alignments in maf format are 
located. The names of the species in the maf files have to be the same as the 
names of the genome files. Format should be: .maf.gz, .maf.Z, .maf.bz2.  Not required if option -a is taken. \n";
    $outstr = "$outstr -f refspecies:
the maf files were created by using one species as the reference species. 
Please provide the name of the reference species as it can be found in the 
multiple sequence alignment and as a genome file. Not required if option -a is taken. \n";
    $outstr = "$outstr -n newicktree:
Specify a newick tree that represents the species used in the analysis. Numbers of occuring genetic
events will be written on the tree in an output file. This option is not needed when -z is activated \n";
    $outstr = "$outstr -d IDtranslation
As species names are different in newick trees and MAF alignments, we need a file showing both names to be able
to translate. If you are sure that the names in your newick tree are the same as in the MAF alignments, this
option is not needed. \n";
    $outstr = "$outstr -a pathToGenesFolder:
Give this options in case creating the clusters and graphs should be redone without analyzing the maf files
again. In this case, the option -a will specify the folder with .bed files for each species including information
about adjacent maf blocks (e.g. the genes folder from a former run could be stored and reused).\n";
    $outstr = "$outstr -s seqSimilarity and -t strucSimilarity:
when testing for homology, only sequences with a sequence similarity above 
seqSimilarity and a secondary structure similarity above strucSimilarity will 
be considered as homologous whereas the value is measured as percentage, thus 
1 means the sequences are the same. Default for both : 0.9; if for any of the 
two parameter -1 is given, this information will be skipped and the analysis 
will only be based on one of both, sequence or structure information. \n";
    $outstr = "$outstr -u pseudoseqSimilarity and -v pseudostrucSimilarity:
when testing for homology, only sequences with a sequence similarity above 
pseudoseqSimilarity and a secondary structure similarity above pseudostrucSimilarity will 
be considered as pseudogenes and thus counted as a pseudogenization event. 
Default for both : 0.7; if for any of the 
two parameter -1 is given, this information will be skipped and the analysis 
will only be based on one of both, sequence or structure information. 
If both are -1, pseudogenization is not considered in the analysis. 
If -t (strucsim) is set to -1, -v (pseudostrucsim) won't be considered either,
and analogous for the sequence similarities -s and -u. \n";
    $outstr = "$outstr -q threshold:
remove a percentage of the lowest scoring blocks based on the MAF scores, 
whereas the value is between 0 and 100 \n";
    $outstr = "$outstr -c CM or -l genelist:
either give a covariance model (cm file) for the gene of interest or give a 
list of genes (format: .bed, .bed.gz, .bed.Z, .bed.bz2). The gene list should 
contain one line per gene with tab separated elements as follows: chromosome     species     startCoord     endCoord     '+'or'-'     sequence     secondary_structure(optional). 
If either secondary structure or sequence are not given or should not be taken 
into account, please use the parameters -s and -t and set the corresponding 
one to -1. Not required if option -a is taken. \n";
    $outstr = "$outstr -w for inclusion E-value threshold or -b for an inclusion bit score threshold as parameter for the infernal run:
As described in the infernal manual (to be found here: eddylab.org/infernal), inclusion thresholds control which hits are considered to be significant.
In case you want to include more putative sequences into our pipeline as significant hits based on the given covariance model, please either increase the e-value threshold using parameter -w (default = 0.01) or decrease the bitscore threshold with parameter -b (this value is usually not set when running infernal with default parameters). Infernal will consider sequences significant if their e-value <= the e-value threshold or their bit score >= the bit score threshold. The README_output file will explain where to find the infernal output in order to see which sequences were included in the analysis.\n";
    $outstr = "$outstr if option -x is activated, graphs are NOT analysed for structural properties (cographs/noncographs).\n";
    $outstr = "$outstr if option -z is activated, duplication alignments are NOT created, thus genetic events are NOT counted.\n";
    $outstr = "$outstr -y python:
please specify the path to where python is installed. Version should be >= 3.0\n";
    $outstr = "$outstr -e perl:
please specify the path to where perl is installed. Version should be >= 5.0\n";
    $outstr = "$outstr -r R:
please specify the path to where R is installed. Version should be >= 3.2.0 \n";
    $outstr = "$outstr -i Infernal:
please specify the path to where Infernal is installed. Version should be >= 1.1.1 . Infernal is only needed when using option -c. With option -l Infernal is not used. \n";
    $outstr = "$outstr \n";
    $outstr = "$outstr Further parameter:\n";
    $outstr = "$outstr choose -k for contact information and citation, -v to see the 
version information and -h to display this help message.\n\n";


return $outstr;



}

__END__

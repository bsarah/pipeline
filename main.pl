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
#use Getopt::Std;
use Getopt::Long qw(GetOptions);
use POSIX qw(strftime);
use File::Basename;
use File::Find;

# declare the perl command line flags/options we want to allow
#my %options=();
#getopts("hVKo:g:m:f:s:t:c:l:q:w:b:a:y:e:r:i:xzu:v:n:d:", \%options);


##define variables
my $outpath;
my $genomes="";
my $mafs="";
my $refspecies="";
my $seqsim =0.9;
my $strucsim =0.9;
my $pseudoscore=-1;
#my $pstrucsim =-1;
my $cmfile="";
my $genefile="";
my $cmoption;
#my $infscore="";
my $perc="";
my $pathtocam; #path to output of cams program (genes folder)
my $pythonpath;
my $evalin=0.01;
my $bitvalin="";
my $perlpath;
my $rpath;
my $infernalpath;
#my $extravalue;
my $newicktree;
my $ids;
my $help;
my $cont;
my $vers;
my $cit;
my $skipg;
my $skipa;


my $dirname = dirname(__FILE__);


#define other information
my $toolname = "FindAFancyAbbrevation";
my $scripts_cam = "$dirname/scripts_cam";
my $scripts_sarah = "$dirname/scripts_sarah";
my $epope = "$dirname/epope";
my $version = "0.1";
my $citation = "Beerinformatics Leipzig";
my $contact = "bsarah at bioinf dot uni-leipzig dot de or anneh at bioinf.uni-leipzig dot de";
my $altnwpath = "$dirname/scripts_sarah";
##test if it work when altNW executable is in scripts folder. (altNW is only for 64 bit!)
#"/scr/gin/bsarah/software/alternativeNW/dist/build/altNW/";

my $helpstr = createHelp($toolname,$version,$citation,$contact);

GetOptions(
    'out|o=s' => \$outpath,
    'genomes=s' => \$genomes,
    'maf=s' => \$mafs,
    'ref=s' => \$refspecies,
    'seqsim=s' => \$seqsim,
    'strucsim=s' => \$strucsim,    
    'pseudo=s' => \$pseudoscore,
#    'pstrucsim=s' => \$pstrucsim,
    'cm=s' => \$cmfile,
    'genes=s' => \$genefile,
#    'infscore=s' => \$infscore,
    'filter=s' => \$perc,
    'again=s' => \$pathtocam,
    'incE=s' => \$evalin,
    'incT=s' => \$bitvalin,
    'newick=s' => \$newicktree,
    'id=s' => \$ids,
    'perl=s' => \$perlpath,
    'rpath=s' => \$rpath,
    'python=s' => \$pythonpath,
    'infernal=s' => \$infernalpath,
    'nograph' => \$skipg,
    'noaln' => \$skipa,
    'help|h' => \$help,
    'version|v' => \$vers,
    'citation' => \$cit,
    'contact' => \$cont
    ) or die "Usage: \n $helpstr \n";




##TODO: test if program can be called with -h -v together

##Help page
if ($help){print $helpstr;exit 0;}
##Contact
if($cont){print "If you have further questions, please contact $contact \n";exit 0;}
##Citation
if($cit){print "If you use this program, please cite: $citation \n"; exit 0;}
##Version
if($vers){print "Program version $version \n"; exit 0;}


###General parameters

my $checkgraphs =1;
my $createalns = 1;


##Outpath
my $outpathstr = "";
if ($outpath){$outpathstr = "-o $outpath";}
else{print "No output path given (option -o)!\n"; exit 1;}
##if output folder doesnt exist, create it!
if(-e $outpath){}
else{
    my $cmd42= "mkdir $outpath";
    my @out42 = readpipe("$cmd42");
}


##Refspecies
my $refspeciesstr = "";
if ($refspecies){$refspeciesstr="--ref $refspecies";
		 #check if file exists
		 if($genomes){
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
}
else{print "No references species given! (option --ref)\n"; exit 1;}


##Parameters for modes for skipping graph analysis or alignments
my $skipgstr = "";
my $skipastr = "";
if ($skipg){$checkgraphs = 0;$skipgstr = "--skipg";}
if ($skipa){$createalns = 0;$skipastr = "--skipa";}


##Newicktree and IDs
my $newickstr = "";
if($newicktree && $createalns == 1){$newickstr = "--newick $newicktree";}
elsif((! $newicktree) && $createalns == 1){print "Please give a newick tree as input file (option --newick)\n"; exit 1;}
else{}
my $idstr = "";
if($ids && $createalns == 1){$idstr = "--id $ids";}
elsif((! $ids) && $createalns == 1){print "No file given to translate newick IDs into MAF IDs (option --id). Please make sure that IDs fit!\n";}
else{}



##similarity parameter, default = 0.9, inactivated with -1
my $simstr = "--seqSim $seqsim --strucSim $strucsim --pseudo $pseudoscore";
if($seqsim == -1 && $strucsim == -1){
    print "Similarity thresholds for both, sequence and structure is -1! At least one has to be > 0! (Use parameter -s for sequence and -t for structure to specify the value) \n"; exit 1;
}
##pseudovalue, default = -1, thus inactivated


##Paths: perl, python, r
my $pystr = "";
my $pestr = "";
my $rstr = "";
##Python
if ($pythonpath){$pystr = "--python $pythonpath";}
else{print "No path to python3 given! (option --python)\n"; exit 1;}

##Perl
if ($perlpath){$pestr = "--perl $perlpath";}
else{print "No path to perl given! (option --perl)\n"; exit 1;}

##R (not needed at the moment)
if ($rpath){$rstr = "--rpath $rpath";}
else{print "No path to R given! (option --rpath)\n"; exit 1;}



###Parameters for infernal mode

##CM
my $cmoptstr="";
if($cmfile){$cmoption = "-sg $cmfile"; $cmoptstr = "--cm $cmfile";}


##Genomes
my $genomesstr="";
if ($genomes){ 
    $genomesstr = "--genomes genomes";
    $cmoption = "$cmoption $genomes";
    #check if folder is not empty
    if(-e $genomes){} else{print "Genomes folder is empty! (option -g)\n"; exit 1;}
}


##Maf
my $mafstr = "";
if ($mafs){$mafstr = "--maf $mafs";}


#Infernal path
my $infstr = "";
if ($infernalpath){$infstr = "--infernal $infernalpath";}


##infernal parameter
my $inclopt="";
my $incloptstr = "";
if ($evalin){$inclopt = "-incE $evalin"; $incloptstr = "$incloptstr\--incE $evalin";}
if ($bitvalin){$inclopt = "-incT $bitvalin";$incloptstr = "$incloptstr\--incT $bitvalin";}


##Percentage of the lowest scoring blocks based on the MAF scores to be removed (between 0 and 100)
my $filterstr = "";
if ($perc){$filterstr = "--filter $perc";}


if($cmfile){
    if($genomes && $mafs && $infernalpath){}
    else{print "Infernal mode but some required parameters are missing! See --help for more information! \n"; exit 1;}
}


###Genelist mode

if($genefile){$cmoption = "-og $genefile"; $cmoptstr = "--genes $genefile";}

##maf and filter should already be reported

if($genefile){
    if($mafs){}
    else{print "Genelist mode but --maf parameter is missing! See --help for more information! \n"; exit 1;}
}


###Repetition mode
my $doitagainstr = "";
if($pathtocam){
    $doitagainstr="--again $pathtocam";
    $genefile=$pathtocam;
    if(-e $pathtocam){}
    else{print "Option -a given but argument folder $pathtocam doesn't exist! \n"; exit 1;}
}




##Program call
my $optstr = "$outpathstr $cmoptstr $mafstr $refspeciesstr $newickstr $idstr $simstr $doitagainstr $filterstr $incloptstr $pystr $pestr $rstr $infstr";
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
print $outs "======Program parameters======\n";
print $outs "$optstr \n";
print $outs "\n\n";
print $outs "======Output statistics======\n";
print $outs "\n";


##Cam's part (sort the genetic elements into genomic anchors based on the maf blocks)
my $genesfolder="";
if(! $pathtocam){
    print "analysis of maf files started (this might take a while)..\n";
    if($perc eq ""){
	open(PROG,"$pythonpath\/python3 $scripts_cam\/main.py $cmoption $inclopt $mafs $infernalpath $outpath $dirname $refspecies 2>>$err0 |") or die "Couldn't start program!";
	while(<PROG>){print "$_";}
    }
    else{
	open(PROG,"$pythonpath\/python3 $scripts_cam\/main.py $perc $cmoption $inclopt $mafs $infernalpath $outpath $dirname $refspecies 2>>$err0 |") or die "Couldn't start program!";
	while(<PROG>){print "$_";}
    }
        print "Done!\n";
    $genesfolder = "$outpath\/genes";
}
else{
    $genesfolder = $pathtocam;
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
my $cmd3 = "$perlpath\/perl $scripts_sarah\/collectCluster\.pl $genesfolder\/specieslist $genesfolder $pseudoscore $outpath\/clusters $sumcollectcluster 2>>$err";
my $cmd4 = "ls $outpath\/clusters/*.clus > $outpath\/clusters/precluslist 2>>$err";
my $cmd5 = "$perlpath/perl $scripts_sarah\/getNumbers.pl $outpath\/clusters/precluslist $outpath\/clusters $sumgetnumbers 2>>$err";
print "construct clusters..";
my @out1 = readpipe("$cmd1");
my @out2 = readpipe("$cmd2");
my @out3 = readpipe("$cmd3");
my @out4 = readpipe("$cmd4");
my @out5 = readpipe("$cmd5");

##looks like: $species\-$num_elems\=$species\-...
my $totelemnumstr = $out3[0];
print "$totelemnumstr\n";

append2file($outs,$sumcollectcluster);
append2file($outs,$sumgetnumbers);
print "Done!\n";

#sort clusters without specific coordinates (none cluster)
my $cmd6 = "mkdir $outpath\/clusters\/NoneCluster 2>>$err";
my $cmd7 = "mv $outpath\/clusters\/cluster-None\* $outpath\/clusters\/NoneCluster 2>>$err";
my $cmd7a = "touch $outpath\/clusters\/NoneCluster/nonecluslist 2>>$err";;
my $cmd7b = "ls $outpath\/clusters\/NoneCluster/*.clus > $outpath\/clusters\/NoneCluster/nonecluslist 2>>$err";
my $cmd7c = "$perlpath\/perl $scripts_sarah\/countNones.pl $outpath\/clusters\/NoneCluster/nonecluslist";
##report the number of elements for each species that are contained in none clusters

print "sort clusters..";
my @out6 = readpipe("$cmd6");
my @out7 = readpipe("$cmd7");
my @out7a = readpipe("$cmd7a");
my @out7b = readpipe("$cmd7b");
my @out7c = readpipe("$cmd7c");
print "Done!\n";

my $nonestr = $out7c[0];

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
##sortCluster should output a summaryfile for singletons...
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

##add another summary file that includes the singleton and none cluster count
my $singlecounts = "$summarypath\/Singleton_Counts.txt";
my $cmdallsingles = "ls $outpath\/clusters/singletons\/*.clus > $outpath\/clusters/singletons\/singletonlist";
my $cmdsingles = "$perlpath\/perl $scripts_sarah\/countElems.pl $outpath\/clusters/singletons/singletonlist $summarypath\/Singleton_Counts.txt 2>>$err";
my $nonecounts = "$summarypath\/None_Counts.txt";
my $cmdnones = "$perlpath\/perl $scripts_sarah\/countElems.pl $outpath\/clusters\/NoneCluster/nonecluslist $summarypath\/None_Counts.txt 2>>$err";

readpipe("$cmdallsingles");
readpipe("$cmdsingles");
readpipe("$cmdnones");
append2file($outs,$singlecounts);
append2file($outs,$nonecounts);



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
    my $cmd281 = "mkdir $outpath\/data_iTOL 2>>$err";
    my $cmd291 = "touch $outpath\/matches.txt 2>>$err";
    my $cmd292 = "touch $outpath\/duplications.txt 2>>$err";
    my $cmd293 = "touch $outpath\/insertions.txt 2>>$err";
    my $cmd294 = "touch $outpath\/pseudogenes.txt 2>>$err";
    my $cmd29a = "mkdir $outpath\/graphs/GainLoss 2>>$err";
#    my $cmd29b = "touch $outpath\/graphs/showGraphs/graphsToDraw 2>>$err";
##getDuplication not needed anymore as it was replaced by create alignments
#    my $cmd30 = "$perlpath\/perl $scripts_sarah\/getDuplication.pl $outpath\/graphs/edlilist $outpath\/graphs/showGraphs $outpath\/graphs/alignments $altnwpath $outpath\/matches.txt $outpath\/duplications.txt $outpath\/insertions.txt $outpath\/pseudogenes.txt 2>>$err";
    my $cmd301 = "touch $outpath\/tree.out 2>>$err";
    my $cmd302 = "touch $outpath\/geneticEvents.txt 2>>$err";
    ##TODO set the pseqsim and pstruclim if we have a solution for pseudogenes
    my $cmd30b = "$perlpath\/perl $scripts_sarah\/createAlignments.pl $outpath\/graphs/edlilist $outpath\/graphs/alignments $altnwpath $seqsim $strucsim -1 $singletoncount $outpath\/graphs/GainLoss $outpath\/matches.txt $outpath\/duplications.txt $outpath\/insertions.txt $outpath\/pseudogenes.txt $sumcreatealn 2>>$err";
    my $cmd30a = "$perlpath\/perl $scripts_sarah\/countEvents.pl $newicktree $outpath\/matches.txt $outpath\/duplications.txt $outpath\/insertions.txt $outpath\/pseudogenes.txt $outpath\/tree.out $outpath\/geneticEvents.txt $totelemnumstr $nonestr $outpath\/data_iTOL 2>>$err";
    
    print "create duplication alignments..";
    my @out28 = readpipe("$cmd28");
    my @out281 = readpipe("$cmd281");
    my @out291 = readpipe("$cmd291");
    my @out292 = readpipe("$cmd292");
    my @out293 = readpipe("$cmd293");
    my @out294 = readpipe("$cmd294");
    my @out29a = readpipe("$cmd29a");
#    my @out29b = readpipe("$cmd29b");
#    my @out30 = readpipe("$cmd30");
    my @out301 = readpipe("$cmd301");
    my @out302 = readpipe("$cmd302");
    my @out30b = readpipe("$cmd30b");
    my @out30a = readpipe("$cmd30a");
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
    $outstr = "$outstr ##Program use:
         perl main.pl --out OutputFolder --ref refspecies --newick newicktree [--id IDtranslation]
                     [--seqsim seqSim] [--strucsim strucSim] 
                     [--pseudo pseudoScore] [--skipg] [--skipa]
                      --python pythonpath --perl perlpath --rpath Rpath 
{mode options 
(a)infernal mode:     --maf MafFolder --genomes GenomesFolder --cm Cmfile --infernal infernalpath 
                     [--incE e-value] [--incT bitscore] [--filter threshold],
(b)genelist mode:     --maf MafFolder --genes GeneList [--filter threshold],
(c)repetition mode:   --again pathToPreviousOutput,
} 

Please choose one of the modes and use its parameters!
  \n\n";
    $outstr = "$outstr \n";
    $outstr = "$outstr Input options:

--out|-o Outputfolder  specify the location of the output files. If folder doesn't exist,
                       it will be newly created.
--ref refSpecies       specify the reference species. The identifier must be the same as given in all
                       other input files, e.g. maf files, genelist, newicktree of IDtranslation file.
                       The specified reference species must be the same as in the maf alignments.
--newick newickTree    tree in newick format showing the phylogenetic relation between target species.
                       In case the species identifier in the tree do not fit to the maf alignments,
                       an ID translation file can be specified (--id).
                       Inner nodes should have a unique name.
--id IDtranslation     file including both, the species identifier in the tree and in the maf alignments
                       to be able to translate. Not needed if the names in the tree fit.
                       Tree and translation file can be downloaded from NCBI and should have NCBI's format.
--seqsim seqSim        optional, default=0.9, value between 0 and 1. Specifies the sequence similarity 
                       (percentage) that is used to detect orthologuous sequences (>= seqSim). 
                       if set to -1, only structure similarity is used.
                       If seqsim and strucsim are used, genetic elements are orthologuous if both thresholds
                       fit (>= seqSim AND >= strucSim).
--strucsim strucSim    optional, default=0.9, value between 0 and 1. Specifies the structure similarity
                       (percentage) that is used to detect orthologuous sequences. (>= strucSim)
                       if set to -1, only sequence similarity is used.
                       If seqsim and strucsim are used, genetic elements are orthologuous if both thresholds
                       fit (>= seqSim AND >= strucSim).
--pseudo pseudoScore   optional, default = -1 (inactivated). Based on the infernal score, all sequences
                       with infernal score < pseudoScore are considered pseudogenes, values above active genes.
                       After running the program for the first time with default parameters, the range of 
                       infernal scores will be included in the output. Using the repetition mode, the program
                       can be rerun specifying several parameters for the orthologuous elements detection. 
--skipg                If --skipg is given, graph analysis step will be inactivated, thus NO analysis of
                       structural properties. This is useful when repeating a run to save time. (optional)
--skipa                If --skipa is given, creation of duplication alignments is inactivated, thus NO counting of
                       genetic events. This is useful when repeating a run to save time. (optional)
--python pythonpath    Specifies the path to the installed python version (must be at least python 3.0)
--perl perlpath        Specifies the path to the installed perl version (must be at least perl 5.0)
--rpath Rpath          Specifies the path to the installed R version (must be at least R 3.2)
\n";
    $outstr = "$outstr ##The program modes:

a) Infernal mode: Enter a covariance model of your target elements, the genomes to be
                  scanned and the path to infernal. The pipeline will use the infernal 
                  output and maf alignments as input for further analysis. 
                  Input options (additional to the general ones):

--maf MafFolder           path to the folder where multiple sequence alignments in maf format are 
                          located. The names of the species in the maf files have to be the same as the 
                          names of the genome files. Format should be: .maf.gz, .maf.Z, .maf.bz2.
--genomes GenomesFolder   path to folder where genomes are located. The names of the genome files must 
                          be the same as the species' name in the MSA. The genome files should have 
                          fasta format and can be gzipped.
--cm CMfile               Give a covariance model (cm file) for the gene of interest created with cmbuild
                          a program of the infernal suite. Infernal manual: eddylab.org/infernal
--infernal infernalpath   path to the installed version of infernal (cmsearch) with version at least 1.1.1
                          Infernal manual: eddylab.org/infernal
--incE e-value            optional, inclusion e-value threshold for infernal run, default = 0.01. Specifies 
                          which infernal hits are considered to be significant. To catch more hits, increase
                          the inclusion e-value threshold. Infernal manual: eddylab.org/infernal
--incT bitscore           optional, inclusion bitscore threshold for infernal run, default = unset. Specifies 
                          which infernal hits are considered to be significant. To catch more hits, decrease
                          the inclusion bitscore threshold. Infernal manual: eddylab.org/infernal
--filter threshold        optional, remove a percentage of the lowest scoring blocks based on the MAF scores, 
                          whereas the value is between 0 and 100
\n
b) Genelist mode: Enter a list with target genes that will be sorted into the maf blocks
                  and used for further analysis. Input options (additional to the general ones):

--maf MafFolder           path to the folder where multiple sequence alignments in maf format are 
                          located. The names of the species in the maf files have to be the same as the 
                          names of the genome files. Format should be: .maf.gz, .maf.Z, .maf.bz2.
--genes GeneList          Give a list of genes (format: .bed, .bed.gz, .bed.Z, .bed.bz2). The gene list should 
                          contain one line per gene with tab separated elements as follows: 
                          chromosome, species, startCoord, endCoord, '+'or'-', sequence, 
                          secondary_structure, score.
                          If some columns are of no importance, please fill in some dummy values. Thus, there should
                          be 10 tab-separated columns in the input gene lists. 
                          If either secondary structure or sequence are not given or should not be taken 
                          into account, please use the parameters --seqSim and --strucSim and set the corresponding 
                          one to -1.
--filter threshold        optional, remove a percentage of the lowest scoring blocks based on the MAF scores, 
                          whereas the value is between 0 and 100
\n
(c) Repetition mode: After running the program once with default parameters, a second run can be started without 
                     analyzing the maf files again. Similarity parameters can be different than in the first run 
                     to change the output to the preferred direction.   

--again pathToPreviousOut Give this options in case the analysis should be redone without reading the maf
                          alignments again. This option can be used to adapt the similarity or scoring 
                          parameters after having seen a first output based on default parameters.
                          The option will specify the folder with .bed files for each species including information
                          about adjacent maf blocks (e.g. the genes folder from a former run could be stored and reused).
\n
";
    $outstr = "$outstr \n";
    $outstr = "$outstr Further parameter:
--help|-h                 print this help page
--version|-v              print version information
--contact                 print contact information
--citation                print citation information
\n";

    return $outstr;

}

__END__

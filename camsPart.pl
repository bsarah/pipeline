#!/usr/bin/perl -w

#this subprogram is for only running cameron's part, thus to sort the genes
#from the genes list or cm in between the maf blocks and output
#the temp, bed and genes folder in output

use Data::Dumper;
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use POSIX qw(strftime);
use File::Basename;
use File::Find;

my $outpath;
my $genomes="";
my $mafs="";
my $refspecies="";
#my $seqsim =0.9;
#my $strucsim =0.9;
#my $pseudoscore=-1;
#my $pstrucsim =-1;
my $cmfile="";
my $genefile="";
my $cmoption;
#my $infscore="";
my $perc="";
my $pathtocam; #path to output of cams program (genes folder)
my $pythonpath;
my $perlpath;
my $infernalpath;
my $evalin=0.01;
my $bitvalin="";


my $pwdcmd = "pwd";
my @outpwd = readpipe("$pwdcmd");

my $dirname = "/homes/biertank/bsarah/Documents/projects/trnaevo/pipeline";
print STDERR "dirname: $dirname \n";

my $toolname = "FindAFancyAbbrevation";
my $scripts_cam = "$dirname/scripts_cam";
my $scripts_sarah = "$dirname/scripts_sarah";



GetOptions(
    'out|o=s' => \$outpath,
    'genomes=s' => \$genomes,
    'maf=s' => \$mafs,
    'ref=s' => \$refspecies,
#    'seqsim=s' => \$seqsim,
#    'strucsim=s' => \$strucsim,    
#    'pseudo=s' => \$pseudoscore,
#    'pstrucsim=s' => \$pstrucsim,
    'cm=s' => \$cmfile,
    'genes=s' => \$genefile,
#    'infscore=s' => \$infscore,
    'filter=s' => \$perc,
    'again=s' => \$pathtocam,
    'incE=s' => \$evalin,
    'incT=s' => \$bitvalin,
#    'newick=s' => \$newicktree,
#    'id=s' => \$ids,
    'perl=s' => \$perlpath,
#    'rpath=s' => \$rpath,
    'python=s' => \$pythonpath,
    'infernal=s' => \$infernalpath,
#    'nograph' => \$skipg,
#    'noaln' => \$skipa,
#    'help|h' => \$help,
#    'version|v' => \$vers,
#    'citation' => \$cit,
#    'contact' => \$cont
    ) or die "Usage: \n blubb \n";

##Outpath
my $outpathstr = "";
if ($outpath){$outpathstr = "-o $outpath ";}
else{print "No output path given (option -o)!\n"; exit 1;}
##if output folder doesnt exist, create it!
if(-e $outpath){}
else{
    my $cmd42= "mkdir $outpath";
    my @out42 = readpipe("$cmd42");
}


##Refspecies
my $refspeciesstr = "";
if ($refspecies){$refspeciesstr="--ref $refspecies ";
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

my $pystr = "";
my $pestr = "";
##Python
if ($pythonpath){$pystr = "--python $pythonpath ";}
else{print "No path to python3 given! (option --python)\n"; exit 1;}

##Perl
if ($perlpath){$pestr = "--perl $perlpath ";}
else{print "No path to perl given! (option --perl)\n"; exit 1;}

my $mode = -1;

###Parameters for infernal mode

##CM
my $cmoptstr="";
if($cmfile){$cmoption = "-sg $cmfile"; $cmoptstr = "--cm $cmfile "; $mode = 0;}


##Genomes
my $genomesstr="";
if ($genomes){ 
    $genomesstr = "--genomes $genomes ";
    $cmoption = "$cmoption $genomes";
    #check if folder is not empty
    if(-e $genomes){} else{print "Genomes folder is empty! (option -g)\n"; exit 1;}
}

#infernal
my $infstr = "";
if ($infernalpath){
    $infstr = "--infernal $infernalpath ";
    $cmoption = "$cmoption $infernalpath";
}

##Maf
my $mafstr = "";
if ($mafs){$mafstr = "--maf $mafs ";}


##infernal parameter
my $inclopt="";
my $incloptstr = "";
if ($evalin && $cmfile){$inclopt = "-incE $evalin"; $incloptstr = "$incloptstr\--incE $evalin ";}
if ($bitvalin && $cmfile){$inclopt = "-incT $bitvalin";$incloptstr = "$incloptstr\--incT $bitvalin ";}


##Percentage of the lowest scoring blocks based on the MAF scores to be removed (between 0 and 100)
my $filterstr = "";
if ($perc){$filterstr = "--filter $perc ";}


if($cmfile){
    if($genomes && $mafs && $infernalpath){}
    else{print "Infernal mode but some required parameters are missing! See --help for more information! \n"; exit 1;}
}


###Genelist mode

if($genefile){$cmoption = "-og $genefile"; $cmoptstr = "--genes $genefile "; $mode = 1;}

##maf and filter should already be reported

if($genefile){
    if($mafs){}
    else{print "Genelist mode but --maf parameter is missing! See --help for more information! \n"; exit 1;}
}

my $err0 = "$outpath\/errorsPart1.txt";

if(! $pathtocam){
    print "analysis of maf files started (this might take a while)..\n";
    print STDERR "call to cam's prog:\n $pythonpath\/python3 $scripts_cam\/main.py $cmoption $inclopt $mafs $outpath $dirname $refspecies $pythonpath 2>>$err0 | \n";
    if($perc eq ""){
	open(PROG,"$pythonpath\/python3 $scripts_cam\/main.py $cmoption $inclopt $mafs $outpath $dirname $refspecies $pythonpath 2>>$err0 |") or die "Couldn't start program!";
	while(<PROG>){print "$_";}
    }
    else{
	open(PROG,"$pythonpath\/python3 $scripts_cam\/main.py $perc $cmoption $inclopt $mafs $outpath $dirname $refspecies $pythonpath 2>>$err0 |") or die "Couldn't start program!";
	while(<PROG>){print "$_";}
    }
    print "Done!\n";
#    $genesfolder = "$outpath\/genes";
}
else{
#    $genesfolder = $pathtocam;
}

#!/usr/bin/perl -w

#this should be a faster version of the pipeline, with almost no output files

#Input:
#path to camerons output
#path to output folder
#mode: if mode == 1, include types ad pseudogenes
#path to needleman wunsch
#similarity thresholds for edges
#newicktree
#path to temp folder of cameron

use Data::Dumper;
use strict;
use warnings;

my $inpath = shift; #here are the following folders: genes, bed, temp
my $outpath = shift;
my $scriptpath = shift;
my $mode = shift;
my $pathtonw = shift;
my $seqsim = shift;
my $strucsim = shift;
my $nwtree = shift;
#my $path2temp = shift;

my $specieslist = "$outpath\/fast_species";
my $gencmd = "ls  $inpath\/genes \> $specieslist";
my @outgncmd = readpipe("$gencmd");

my %blocks = ();
my %species = ();
my %pseuspecies = ();
my %nones = ();
my %pseudonones = ();
my %nonesNT = ();
my %pseudononesNT = ();

my $eventlist = "$outpath\/allClusters\.txt";
my $listcmd = "touch $eventlist";
readpipe("$listcmd");


open FA,"<$specieslist" or die "can't open $specieslist\n";

while(<FA>){
    chomp;
    my $precurfile = $_;
    my $curfile = "$inpath\/genes\/$precurfile";
    #print "CURFILE: $curfile\n";

    #get species name
    my @S1 = split '\/', $curfile;
    my @S2 = split '\.', $S1[-1];
    my $spec = $S2[0];
    #insert in %specis
    $species{$spec} = 0;
    $pseuspecies{$spec} = 0;
    #sort file and remove blanks
    my $curfile_nb = "$outpath\/$spec\_nb\.bed";
    my $curfile_sorted = "$outpath\/$spec\_sorted\.bed";
    my $cmdnoblanks = "awk \'NF\' $curfile \> $curfile_nb";
    my $cmdsort = "sort -n -k 6,7 $curfile_nb \> $curfile_sorted";

#    print "CMDNOBLANKS: $cmdnoblanks\n";
#    print "CMDSORT: $cmdsort\n";
    
    readpipe("$cmdnoblanks");
    readpipe("$cmdsort");
    #read file to blocks/nones
    open CF, "<$curfile_sorted" or die "can't open $curfile_sorted\n";
    while(<CF>){
	chomp;
	my $line = $_;
	if($line=~/^#/){next;}
	if($line eq ""){next;}
	my @F = split '\t', $line;
	if($mode==1 && ($F[-2] eq "t" || $F[-2] eq "T" ||$F[-2] eq "TRUE" || $F[-2] eq "True" || $F[-2] eq "true" || $F[-2] eq "1")){
	    $pseuspecies{$spec}++;
	}
	else{
	    $species{$spec}++;
	}
	my $none = "None";
	if($F[5] eq $none || $F[6] eq $none){
	    my $nonekey;
	    if($mode==1){
		my $type = $F[-3];
		$nonekey = "$spec\_$type";
		my $pseu = $F[-2];
		if($pseu eq "TRUE" || $pseu eq "True" || $pseu eq "true" || $pseu eq "1"){
		    if(exists($pseudonones{$nonekey})){
			$pseudonones{$nonekey}++;
		    }
		    else{
			$pseudonones{$nonekey}=1;
		    }
		    if(exists($pseudononesNT{$spec})){
			$pseudononesNT{$spec}++;
		    }
		    else{
			$pseudononesNT{$spec}=1;
		    }
		}
		else{
		    if(exists($nones{$nonekey})){
			$nones{$nonekey}++;
		    }
		    else{
			$nones{$nonekey}=1;
		    }
		    if(exists($nonesNT{$spec})){
			$nonesNT{$spec}++;
		    }
		    else{
			$nonesNT{$spec}=1;
		    }
		}
	    }
	    else{
		$nonekey = "$spec";
		if(exists($nones{$nonekey})){
		    $nones{$nonekey}++;
		}
		else{
		    $nones{$nonekey}=1;
		}
		if(exists($nonesNT{$spec})){
		    $nonesNT{$spec}++;
		}
		else{
		    $nonesNT{$spec}=1;
		}
	    }
	    next;
	}
	else{
	    #first, hash by left block variable
	    #then check the largest right block value and join until then in the hash
	    my $key;
	    if($F[5] < $F[6]){
		$key = "$F[5]\_$F[6]";
	    }
	    else{
		$key = "$F[6]\_$F[5]";
	    }
	    if(exists($blocks{$key})){
		my $tmp = $blocks{$key};
		$blocks{$key} = "$tmp\;$line";
	    }
	    else{
		$blocks{$key} = "$line";
	    }
	}
    }
}

#my $orignum = scalar (keys %blocks);
#print "ORIG CLUSNUM: $orignum \n";

##NO JOINING of blocks based on block number as they are only consecutive in the reference!

#join cluster by hash num
#my %joinedblocks = ();

#my $curblock = "";
#my $curstart = -5;
#my $curend = -1;
#foreach my $lb (sort { $a <=> $b} keys %blocks){
#    my $curb = $blocks{$lb};
#    my $rb = getMaxBlock($curb);
#    print "LB: $lb, RB: $rb \n";
#
#    if($lb <= $curend){
#	#join the cluster
#	$curblock="$curblock\;$curb";
#	if($curend < $rb){$curend = $rb;}
#    }
#    else{
#	#add the current cluster to joinedblocks and start a new one
#	if($curblock ne ""){
#	    my $key = "$curstart\_$curend";
#	    print "NEWKEY: $key \n";
#	    $joinedblocks{$key} = $curblock;
#	}
#	$curstart = $lb;
#	$curend = $rb;
#	$curblock = $curb;
#    }
#}




##print None hash to file
my $noneout = "$outpath\/nones.txt";
my $pseunoneout = "$outpath\/pseunones.txt";
my $nonestr = "";
my $pseudononestr = "";
foreach my $k (keys %nones){
    $nonestr = "$nonestr$k\t$nones{$k}\n";
}
foreach my $k (keys %pseudonones){
    $pseudononestr = "$pseudononestr$k\t$pseudonones{$k}\n";
}
open(my $outnone,">>",$noneout);
print $outnone $nonestr;
close $outnone;
open(my $outpnone, ">>", $pseunoneout);
print $outpnone $pseudononestr;
close $outpnone;



my $allnonestr = "";
foreach my $k (keys %nonesNT){
    $allnonestr = "$allnonestr$k\-$nonesNT{$k}\=";
}
$allnonestr = "$allnonestr\!";
foreach my $k (keys %pseudononesNT){
    $allnonestr = "$allnonestr$k\-$pseudononesNT{$k}\=";
}

my $allspecstr = "";
foreach my $t (keys %species){
    $allspecstr = "$allspecstr$t\-$species{$t}\=";
}
$allspecstr = "$allspecstr\!";
foreach my $t (keys %pseuspecies){
    $allspecstr = "$allspecstr$t\-$pseuspecies{$t}\=";
}

#print "totelemnum: $allspecstr \n";

my %dupevents =();
my %matevents =();
my %insevents =();
my %misevents =(); #missing data=no anchors
my %delevents = ();#deletions
my %psevents = (); #insertions of pseudogenes
my %psedels = ();
my %psemis = ();
my %pseins = ();




my %singletons = (); #hash species -> number
my %pseusingles = (); #include types
my %singletonsNT = (); #hash species -> number
my %pseusinglesNT = (); #include types

##go on with cluster hash

##sort keys of blocks and join adjacent ones?
#my $clusnum = scalar (keys %blocks);
#print "NUM CLUSTER: $clusnum \n";
open(my $outg, ">>",$eventlist);
my $cluscount = 0;
my $sumelems = 0;
foreach my $k (keys %blocks){
    my @A = split '_', $k;
    my $leftanchor = $A[0];
    my $rightanchor = $A[1];
    my @B = split ';', $blocks{$k};
    my $clussize = scalar @B;
    #write to allClustersList
    print $outg ">clus$cluscount $clussize\n";
    $cluscount++;
    $sumelems+=$clussize;
    for(my $cb = 0;$cb < $clussize;$cb++){
	print $outg "$B[$cb]\n";
    }
    if(scalar @B == 1){
	##sort single element cluster into hash
	my @F = split '\t', $B[0];
	my @C = split '_', $F[1];
	my $curspec = $C[0];
	my $singlekey;
	if($mode==1){
	    my $type = $F[-3];
	    $singlekey = "$curspec\_$type";
	    my $pseu = $F[-2];
	    if($pseu eq "T" || $pseu eq "TRUE" || $pseu eq "True" || $pseu eq "true" || $pseu eq "1"){
		if(exists($pseusingles{$singlekey})){
		    $pseusingles{$singlekey}++;
		}
		else{
		    $pseusingles{$singlekey}=1;
		}
		if(exists($pseusinglesNT{$curspec})){
		    $pseusinglesNT{$curspec}++;
		}
		else{
		    $pseusinglesNT{$curspec}=1;
		}
	    }
	    else{
		if(exists($singletons{$singlekey})){
		    $singletons{$singlekey}++;
		}
		else{
		    $singletons{$singlekey} = 1;
		}
		if(exists($singletonsNT{$curspec})){
		    $singletonsNT{$curspec}++;
		}
		else{
		    $singletonsNT{$curspec} = 1;
		}
		
	    }
	}
	else{
	    $singlekey = "$curspec";
	
	    if(exists($singletons{$singlekey})){
		$singletons{$singlekey}++;
	    }
	    else{
		$singletons{$singlekey} = 1;
	    }
	    if(exists($singletonsNT{$curspec})){
		$singletonsNT{$curspec}++;
	    }
	    else{
		$singletonsNT{$curspec} = 1;
	    }
	}
	next;
    }
    #build graph
    my $graphcmd = "perl $scriptpath\/buildEdgeList_fast.pl \"$blocks{$k}\" $mode $pathtonw $strucsim $seqsim";
#    print "GRAPHCMD: $graphcmd \n";
    my @outgraph = readpipe("$graphcmd");
    my $len = scalar @outgraph;
#    print "num edges: $len\n";
    my $graphstr = "";
    for(my $g=0;$g<scalar @outgraph;$g++){
	chomp($outgraph[$g]);
	$graphstr = "$graphstr\=$outgraph[$g]";
    }
#    print "GRAF: $graphstr\n";
    #check graph for cograph or not and edit slightly
    my $cglist = "$outpath\/list\_cographs\.txt";
    my $ncglist = "$outpath\/list\_noncographs\.txt";
    my $cgcmd = "touch $cglist";
    my $ncgcmd = "touch $ncglist";
    readpipe("$cgcmd");
    readpipe("$ncgcmd");
    my $checkcmd = "perl $scriptpath\/checkGraph_fast.pl \"$graphstr\" $k $seqsim $strucsim $mode $cglist $ncglist";
    my @newoutgraph = readpipe("$checkcmd");
    my $newgraphstr = $newoutgraph[0];
    
    #sort the non edges graphs

    #create alignment
    my $alncmd = "perl $scriptpath\/createAlignments_fast.pl \"$newgraphstr\" $outpath $pathtonw $seqsim $strucsim $mode 0 $nwtree $inpath\/temp $leftanchor $rightanchor";
#    print "ALNCMD: $alncmd \n";
    my @outaln = readpipe("$alncmd"); #this array contains: dup mat insertion pseins pseudomatch deletion missinganchor missinanchorpseu deletionpseu
    my $alen = scalar @outaln;
#    print "OUTALN: $alen\n";
    my $alnput = join('_',@outaln);
#    print "ALNOUT: $alnput \n";
    #the elements in each pos are separated by = and key and value are separated by -
    #starting with =!
    #if there are no counts, the string only consists of =, thats why, we start with 1 in the loops

    my @dups = split '=', $outaln[0];
    for(my $d=1;$d<scalar @dups;$d++){
	my @dt = split '-', $dups[$d];
	if(exists($dupevents{$dt[0]})){$dupevents{$dt[0]}+=$dt[1];}
	else{$dupevents{$dt[0]} = $dt[1];}
    }
    my @mats = split '=', $outaln[1];
    for(my $a=1;$a<scalar @mats;$a++){
	my @at = split '-', $mats[$a];
	if(exists($matevents{$at[0]})){$matevents{$at[0]}+=$at[1];}
	else{$matevents{$at[0]} = $at[1];}	
    }
    my @ins = split '=', $outaln[2];
    for(my $n=1;$n<scalar @ins;$n++){
	my @nt = split '-', $ins[$n];
	if(exists($insevents{$nt[0]})){$insevents{$nt[0]}+=$nt[1];}
	else{$insevents{$nt[0]} = $nt[1];}	
    }
    my @pins = split '=', $outaln[3];
    for(my $pn=1;$pn<scalar @pins;$pn++){
	my @pnt = split '-', $pins[$pn];
	if(exists($pseins{$pnt[0]})){$pseins{$pnt[0]}+=$pnt[1];}
	else{$pseins{$pnt[0]} = $pnt[1];}	
    }
    my @pmats = split '=', $outaln[4];
    for(my $pa=1;$pa<scalar @pmats;$pa++){
	my @pat = split '-', $pmats[$pa];
	if(exists($psevents{$pat[0]})){$psevents{$pat[0]}+=$pat[1];}
	else{$psevents{$pat[0]} = $pat[1];}	
    }
    my @dels = split '=', $outaln[5];
    for(my $dd=1;$dd<scalar @dels;$dd++){
	my @ddt = split '-', $dels[$dd];
	if(exists($delevents{$ddt[0]})){$delevents{$ddt[0]}+=$ddt[1];}
	else{$delevents{$ddt[0]} = $ddt[1];}
    }
    my @mis = split '=', $outaln[6];
    for(my $m=1;$m<scalar @mis;$m++){
	my @ms = split '-', $mis[$m];
	if(exists($misevents{$ms[0]})){$misevents{$ms[0]}+=$ms[1];}
	else{$misevents{$ms[0]} = $ms[1];}
    }
    my @pmis = split '=', $outaln[7];
    for(my $pm=1;$pm<scalar @pmis;$pm++){
	my @pms = split '-', $pmis[$pm];
	if(exists($psemis{$pms[0]})){$psemis{$pms[0]}+=$pms[1];}
	else{$psemis{$pms[0]} = $pms[1];}
    }
    my @pdels = split '=', $outaln[8];
    for(my $pdd=1;$pdd<scalar @pdels;$pdd++){
	my @pddt = split '-', $pdels[$pdd];
	if(exists($psedels{$pddt[0]})){$psedels{$pddt[0]}+=$pddt[1];}
	else{$psedels{$pddt[0]} = $pddt[1];}
    }
    
}
close $outg;
##print singletons to file and create allsinglestr
my $singleout = "$outpath\/singletons.txt";
my $pseusingleout = "$outpath\/pseusingletons.txt";
my $singlestr = "";
my $psinglestr = "";
foreach my $s (keys %singletons){
    $singlestr = "$singlestr$s\t$singletons{$s}\n";
}
foreach my $s (keys %pseusingles){
    $psinglestr = "$psinglestr$s\t$pseusingles{$s}\n";
}
open(my $outs, ">>", $singleout);
open(my $outt, ">>", $pseusingleout);
print $outs $singlestr;
print $outt $psinglestr;
close $outs;
close $outt;
my $allsinglestr = "";
foreach my $s (keys %singletonsNT){
    $allsinglestr = "$allsinglestr$s\-$singletonsNT{$s}\=";
}
$allsinglestr = "$allsinglestr\!";
foreach my $s (keys %pseusinglesNT){
    $allsinglestr = "$allsinglestr$s\-$pseusinglesNT{$s}\=";
}



#get output files which are the input to countEvents
my $matchout = "$outpath\/matches\.txt";
my $duplout = "$outpath\/duplications\.txt";
my $insout = "$outpath\/insertions\.txt";
my $pseout = "$outpath\/pseudomatches\.txt";
my $psemisout = "$outpath\/pseudomissing\.txt";
my $psedelout = "$outpath\/pseudodeletions\.txt";
my $pseinsout = "$outpath\/pseudoinsertions\.txt";
my $delout = "$outpath\/deletions\.txt";
my $misout = "$outpath\/missinganchors\.txt";


open(my $outd,">>",$duplout);

foreach my $du (sort keys %dupevents) {
    my @devs = split ',', $du;
    @devs = sort @devs;
    my $dustr = join(',',@devs);
    print $outd "$dustr\t$dupevents{$du}\n";
}
close $outd;

open(my $outm,">>",$matchout);
foreach my $ma (sort keys %matevents) {
    my @mevs = split ',', $ma;
    @mevs = sort @mevs;
    my $mastr = join(',',@mevs);
    print $outm "$mastr\t$matevents{$ma}\n";
}
close $outm;

open(my $outi,">>",$insout);
foreach my $in (sort keys %insevents) {
    my @ievs = split ',', $in;
    @ievs = sort @ievs;
    my $instr = join(',',@ievs);
    print $outi "$instr\t$insevents{$in}\n";
}
close $outi;


open(my $outpi,">>",$pseinsout);
foreach my $pin (sort keys %pseins) {
    my @pievs = split ',', $pin;
    @pievs = sort @pievs;
    my $pinstr = join(',',@pievs);
    print $outpi "$pinstr\t$pseins{$pin}\n";
}
close $outpi;


open(my $outp,">>",$pseout);
foreach my $mi (sort keys %psevents) {
    my @mevs = split ',', $mi;
    @mevs = sort @mevs;
    my $mistr = join(',',@mevs);
    print $outp "$mistr\t$psevents{$mi}\n";
}
close $outp;

open(my $oute,">>",$delout);
foreach my $dl (sort keys %delevents) {
    my @dls = split ',', $dl;
    @dls = sort @dls;
    my $dlstr = join(',',@dls);
    print $oute "$dlstr\t$delevents{$dl}\n";
}
close $oute;

open(my $outn,">>",$misout);
foreach my $ms (sort keys %misevents) {
    my @miss = split ',', $ms;
    @miss = sort @miss;
    my $msstr = join(',',@miss);
    print $outn "$msstr\t$misevents{$ms}\n";
}
close $outn;

open(my $outpm,">>",$psemisout);
foreach my $pm (sort keys %psemis) {
    my @pmiss = split ',', $pm;
    @pmiss = sort @pmiss;
    my $pmsstr = join(',',@pmiss);
    print $outpm "$pmsstr\t$psemis{$pm}\n";
}
close $outpm;

open(my $outpd,">>",$psedelout);
foreach my $pd (sort keys %psedels) {
    my @pdls = split ',', $pd;
    @pdls = sort @pdls;
    my $pdlstr = join(',',@pdls);
    print $outpd "$pdlstr\t$psedels{$pd}\n";
}
close $outpd;


#count events
##singletoncount = spec-num=spec-num...!spec-pnum=spec-pnum=...
##totelem and nonestr as singleton

my $err = "$outpath\/errors";
my $cmderr = "touch $err";
readpipe("$cmderr");

##create summarypath and itolout
my $cmditol = "mkdir $outpath\/data_iTOL";
readpipe("$cmditol");
my $cmdsum = "touch $outpath\/geneticEvents\.txt";
readpipe("$cmdsum");
my $cmdtree = "touch $outpath\/OutTree\.txt";
readpipe("$cmdtree");
my $countcmd = "perl $scriptpath\/countEvents.pl $nwtree $allsinglestr $matchout $duplout $insout $pseout $psemisout $psedelout $pseinsout $delout $misout $outpath\/OutTree\.txt $outpath\/geneticEvents\.txt $allspecstr $allnonestr $outpath\/data_iTOL 2>> $err";
#print "COUNTING: $countcmd\n";
readpipe("$countcmd");



print "Program run ended.\n";
print "Number of clusters: $cluscount\n";
my $avnum = $sumelems/$cluscount;
print "Average number of elements per cluster: $avnum\n";

sub getMaxBlock{
    my @inp = @_;
    my $cluster = $inp[0];
    my @c = split ';', $cluster;
    my $curmax = 0;
    for(my $i=0;$i<scalar @c;$i++){
	my @F = split '\t', $c[$i];
	if($F[5] > $F[6] && $F[5] > $curmax){$curmax = $F[5];}
	elsif($F[6] > $F[5] && $F[6] > $curmax){$curmax = $F[6];}
	else{}
    }
    return $curmax;
}

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
my %leftblocks = ();
my %species = ();
my %pseuspecies = ();
my %nones = ();
my %pseudonones = ();
my %nonesNT = ();
my %pseudononesNT = ();

my $eventlist = "$outpath\/allClusters\.txt";
my $listcmd = "touch $eventlist";
readpipe("$listcmd");
my $remlist = "$outpath\/remoldings\.txt";
my $inremlist = "$outpath\/inremoldings\.txt";
my $remcmd = "touch $remlist";
my $inremcmd = "touch $inremlist";
readpipe("$remcmd");
readpipe("$inremcmd");

open(my $outsr,">>", $remlist);
open(my $outsi,">>", $inremlist);


my %allTypes = (); #hash with spec_type -> num
my %pallTypes = (); #hash with spec_type -> num

my %block2spec = ();
my @allspec = ();

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
    push @allspec, $spec;
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

	if($mode==1)
	{
	    my $tkey = "$spec\_$F[-3]";
	    if($F[-2] eq "t" || $F[-2] eq "T" ||$F[-2] eq "TRUE" || $F[-2] eq "True" || $F[-2] eq "true" || $F[-2] eq "1"){
		$pseuspecies{$spec}++;
		##add to allTypes hash
		if(exists($pallTypes{$tkey})){$pallTypes{$tkey}++;}else{$pallTypes{$tkey}=1;}
	    }
	    else{
		$species{$spec}++;
		##add to allTypes hash
		if(exists($allTypes{$tkey})){$allTypes{$tkey}++;}else{$allTypes{$tkey}=1;}
	    }
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
	    #grep the block and species from temp and store blocks that appear in element context
	    my $bl1 = "$spec\_$F[5]";
	    my $bl2 = "$spec\_$F[6]";
#	    print "BLOCK1,2: $bl1,$bl2\n";
	    my $grepcmdleft = "zcat $inpath\/temp\/$spec\_temp\_sorted\.bed\.gz \| grep -w \"$bl1\" ";
	    my $grepcmdright = "zcat $inpath\/temp\/$spec\_temp\_sorted\.bed\.gz \| grep -w \"$bl2\" ";
	    my @outleft = readpipe("$grepcmdleft");
	    my $lstr = join(';',@outleft);
#	    print "LEFTGREP $lstr \n";
	    if(scalar @outleft == 0){print STDERR "anchor does not exist: $bl1\n";}
	    else{
		my @Fleft = split '\t', $outleft[0];
		my $newstr1 = "$bl1\=$Fleft[0]\=$Fleft[2]\=$Fleft[3]\=$Fleft[4]";#bl1 chr start len strand
		#check if the block already exists for THIS species!
		if(exists($block2spec{$F[5]}) && (index($block2spec{$F[5]},$bl1) == -1)){$block2spec{$F[5]} = "$block2spec{$F[5]}\=$newstr1";}
		else{$block2spec{$F[5]} = $newstr1;}
	    }
	    my @outright = readpipe("$grepcmdright");
	    my $rstr = join(';',@outright);
#	    print "RIGHTGREP $rstr \n";
	    if(scalar @outright == 0){print STDERR "anchor does not exist: $bl2\n";}
	    else{
		my @Fright = split '\t', $outright[0];
		my $newstr2 = "$bl2\=$Fright[0]\=$Fright[2]\=$Fright[3]\=$Fright[4]";
		if(exists($block2spec{$F[6]}) && (index($block2spec{$F[6]},$bl2)== -1)){$block2spec{$F[6]} = "$block2spec{$F[6]}\=$newstr2";}
		else{$block2spec{$F[6]} = $newstr2;}
	    }
	    


	    #first, hash by left block variable
	    #then check the largest right block value and join until then in the hash
	    my $key;
	    my $leftkey;
	    if($F[5] < $F[6]){
		$key = "$F[5]\_$F[6]";
		$leftkey = "$F[5]";
	    }
	    else{
		$key = "$F[6]\_$F[5]";
		$leftkey = "$F[6]";
	    }
	    if(exists($blocks{$key})){
		my $tmp = $blocks{$key};
		$blocks{$key} = "$tmp\;$line";
	    }
	    else{
		$blocks{$key} = "$line";
	    }
	    if(exists($leftblocks{$leftkey})){
		my $ltmp = $leftblocks{$leftkey};
		$leftblocks{$leftkey} = "$ltmp\;$line";
	    }
	    else{
		$leftblocks{$leftkey} = "$line";
	    }
	}
    }
}

my $orignum = scalar (keys %blocks);
print "ORIG CLUSNUM: $orignum \n";


my $numspec = scalar @allspec;
print "Number of species: $numspec \n";

my $blocktable = "$outpath\/blocktable\.tab";
my $nostr = "\=\t\=\t\=\t\=\t\=";
open(my $outbt, ">>", $blocktable);
foreach my $an (sort { $a <=> $b} keys %block2spec){
    #print species in the order as they are in allspec
    #if a species doesnt exist, put =
    #in this way we get a nice table to print
    my @tmp = split '=', $block2spec{$an};
    my $tabline = "$an";
    my $j=0;
    for(my $i=0;$i<scalar @allspec;$i++){
	if($j*5 >= scalar @tmp){
	    $tabline = "$tabline\t$nostr";
	    next;
	}
	my $tocheck = $tmp[(5*$j)+0];
	my $curs = "$allspec[$i]\_$an";
	if($tocheck eq $curs){
	    $tabline = "$tabline\t$tmp[(5*$j)+0]\t$tmp[(5*$j)+1]\t$tmp[(5*$j)+2]\t$tmp[(5*$j)+3]\t$tmp[(5*$j)+4]";
	    $j++;
	}
	else{
	    $tabline = "$tabline\t$nostr";
	}
    }
    print $outbt "$tabline\n";
}
close $outbt;

my %spec2file = (); #species and its sorted blockfile
for(my $s=0;$s<scalar @allspec;$s++){
    my $col1 = 5*$s+3;
    my $col2 = 5*$s+4;
    my $sortfile = "$outpath\/$allspec[$s]\_blocktable\.tab";
    my $sortcmd = "sort -k$col1,$col1 -nk$col2,$col2 $blocktable > $sortfile";
    readpipe("$sortcmd");
    $spec2file{$s} = $sortfile;
}


##NO JOINING of blocks based on block number as they are only consecutive in the reference!
##but the leftblocks are a first step and then, we check if all species are adjacent in the blocktable
#sorted for each species
my $blockfile = "$outpath\/Leftblocks\.txt";
open(my $outb,">>",$blockfile);

#join cluster by hash num
my %joinedblocks = ();

my $curblock = "";
my $curstart = -5;#a
my $curend = -1;#a
foreach my $lb (sort { $a <=> $b} keys %leftblocks){
    my $curb = $leftblocks{$lb};
    my $rb = getRightBlock($curb); #this is the right block, -1 if there are several rightblocks, then it doesnt work 

    print $outb "$lb\t$rb\t$curb\n";
    
    if($curstart == -5){
	$curblock = $curb;
	$curstart = $lb;
	$curend = $rb;
	next;
    }
    #as keys are sorted, if the blocks are adjacent and joinable, lb should be the same as curend
    #adjacent doesn't have to be directly adjacent but all elements in both clusters should be adjacent in all species
    my $isgood = 1;
    #change checkBlocks based on the created blocktable
    #grep for blocknum in each sorted file (for each species) and check if the blocknums are adjacent in all of them
    if($rb>0 && $curend > 0){
	foreach my $sk (keys %spec2file){
	    my $tmpgrep = "grep -A 1 -w $lb $spec2file{$sk}";
	    my @outtmpgrep = readpipe("$tmpgrep"); #output should contain two lines
	    if(scalar @outtmpgrep >= 2 && $outtmpgrep[1] =~ /^$rb/){ #the adjacent block is rb
		$isgood = 1;
	    }
	    else{
		$isgood = 0;
		last;
	    }
	    
	}
    }
    if($rb>0 && $isgood == 1 && $curend > 0){
	#should be joinable
	print STDERR "join a,b,c: $curstart, $curend, $rb\n";
	$curend = $rb;
	$curblock = "$curblock\;$curb";
    }
    else{
	#add the old block to joined clusters
	if($curend == -1){
	    #check if the elements are still neighbors
	    #if not separately add each element in curblock
	    my %rb2clus = ();
	    my @F = split ';', $curblock;
	    for(my $f=0;$f<scalar @F;$f++){
		my @G = split '\t', $F[$f];
		my $rkey;
		if($G[5] < $G[6]){$rkey = $G[6]}
		else{$rkey = $G[5]}
		if(exists($rb2clus{$rkey})){$rb2clus{$rkey}="$rb2clus{$rkey}\;$F[$f]";}		  
		else{$rb2clus{$rkey} = $F[$f];}
	    }
	    print STDERR "curend=-1, lb: $curstart, rb: $rb\n";
	    print STDERR Dumper(\%rb2clus);
	    print STDERR "\n";
	    #curstart = $lb
	    my $curclus = "";
	    my $curend = -1;
	    my $neighborsum = 0; #if this is equals scalar @rks, then we don't have neighbors
	    my @rks = sort { $a <=> $b} keys %rb2clus;
	    for(my $r=0;$r<scalar @rks;$r++){
		my $noneighbor = 0;
		for(my $k=$r+1;$k<scalar @rks;$k++){
		    foreach my $sk (keys %spec2file){
			my $tmpgrep = "grep -A 1 -w $rks[$r] $spec2file{$sk}";
			my @outtmpgrep = readpipe("$tmpgrep"); #output should contain two lines
			if(scalar @outtmpgrep >= 2 && $outtmpgrep[1] =~ /^$rks[$k]/){}
			else{
			    $neighborsum += 1;
			    last;
			}
		    }
		}
	    }
	    
	    if($neighborsum < scalar @rks){ ##check again, but we first check the prints
		print STDERR "join overlapping cluster\n"
		$joinedblocks{$rks[-1]} = join(';',values %rb2clus);
	    }
	    else{
		for(my $e=0;$e<scalar @rks;$e++){
		    my $newkey = "$lb\_$rks[$e]";
		    if(exists($joinedblocks{$newkey})){
			$joinedblocks{$newkey}="$joinedblocks{$newkey}\;$rb2clus{$rks[$e]}";
		    }
		    else{$joinedblocks{$newkey}="$rb2clus{$rks[$e]}";}
		}
	    }
	}
	else{
	    #all elements in curblock are in one cluster
	    my $newkey2 = "$curstart\_$curend";
	    if(exists($joinedblocks{$newkey2})){
		$joinedblocks{$newkey2}="$joinedblocks{$newkey2}\;$curblock";
	    }
	    else{$joinedblocks{$newkey2}=$curblock;}
	}
	#set the new block
	$curblock = $curb;
	$curstart = $lb;
	$curend = $rb;	
    }
    
}

my $newclusnum = scalar (keys %joinedblocks);
print "JOINED CLUSNUM: $newclusnum \n";



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

#print allTypes to file
my $typesout = "$outpath\/allTypes\.txt";
my $ptypesout = "$outpath\/allPseudoTypes\.txt";
open(my $tyout, ">>", $typesout);
foreach my $ty (keys %allTypes){
    print $tyout "$ty\t$allTypes{$ty}\n";
}
close $tyout;
open(my $ptyout, ">>", $ptypesout);
foreach my $pty (keys %pallTypes){
    print $ptyout "$pty\t$pallTypes{$pty}\n";
}
close $ptyout;


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
foreach my $k (keys %joinedblocks){
    my @A = split '_', $k;
    my $leftanchor = $A[0];
    my $rightanchor = $A[1];
    my @B = split ';', $joinedblocks{$k};
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
    #have a tmpfile for blocks string
    my $tmpfile0 = "$outpath\/tmpdata0";
    open(my $outtmp0, ">>", $tmpfile0);
    print $outtmp0 $joinedblocks{$k};
    close $outtmp0;
    my $tmpfile1 = "$outpath\/tmpdata1";
    #build graph, tmpfile1 will contain the output
    my $graphcmd = "perl $scriptpath\/buildEdgeList_fast.pl $tmpfile0 $mode $pathtonw $strucsim $seqsim $tmpfile1";
#    print "GRAPHCMD: $graphcmd \n";
    readpipe("$graphcmd");
    

    #check graph for cograph or not and edit slightly
    my $cglist = "$outpath\/list\_cographs\.txt";
    my $ncglist = "$outpath\/list\_noncographs\.txt";
    my $cgcmd = "touch $cglist";
    my $ncgcmd = "touch $ncglist";
    readpipe("$cgcmd");
    readpipe("$ncgcmd");

    my $checkcmd = "perl $scriptpath\/checkGraph_fast.pl $tmpfile1  $k $seqsim $strucsim $mode $cglist $ncglist";

    ##argument list is getting too long...create a temporary file

    my @newoutgraph = readpipe("$checkcmd");
    my $newgraphstr = $newoutgraph[0];
    
    #sort the non edges graphs

    #create alignment
    my $alncmd = "perl $scriptpath\/createAlignments_fast.pl $tmpfile1 $outpath $pathtonw $seqsim $strucsim $mode 0 $nwtree $inpath\/temp $leftanchor $rightanchor";

#    print "ALNCMD: $alncmd \n";
    my @outaln = readpipe("$alncmd"); #this array contains: dup mat insertion pseins pseudomatch deletion missinganchor missinanchorpseu deletionpseu
    my $rmcmd = "rm $tmpfile1 $tmpfile0";
    readpipe("$rmcmd");
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
    my @currems = split '=', $outaln[9];
    for(my $sr =0;$sr<scalar @currems;$sr++){
	print $outsr "$currems[$sr]\n";
    }
    my @curinrems = split '=', $outaln[10];
    for(my $ri =0;$ri<scalar @curinrems;$ri++){
	print $outsi "$curinrems[$ri]\n";
    }
    
    
}
close $outsi;
close $outsr;
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

    sub getRightBlock{
	#get THE right block. if it is several different ones, return -1 as the cluster seems to be too diverse
    my @inp = @_;
    my $cluster = $inp[0];
    my @c = split ';', $cluster;
    my $curmin = -1;
    for(my $i=0;$i<scalar @c;$i++){
	my @F = split '\t', $c[$i];
	my $rightanchor;
	if($F[5] < $F[6]){$rightanchor = $F[6];}
	else{$rightanchor = $F[5];}
	if($curmin==-1){$curmin = $rightanchor;}
	elsif($curmin != $rightanchor){return -1;}
	else{}
    }
    return $curmin;
}


sub checkBlocks{
    my @inp = @_;
    ##check for intersecting species
    my %spec1 = ();
    my @block1 = split ';', $inp[0];
    for(my $i=0;$i<scalar @block1;$i++){
	my @F = split '\ţ', $block1[$i];
	my @S = split '_', $F[1];
	my $spec = $S[0];
	$spec1{$spec}=1;
    }
    my @block2 = split ';', $inp[1];
    for(my $j=0;$j<scalar @block2;$j++){
	my @G = split '\t', $block2[$j];
	my @T = split '_', $G[1];
	my $spec2 = $T[0];
	if(exists($spec1{$spec2})){return 1;}
    }
    return 0;
}

#!/usr/bin/perl -w

## perl createAlignments.pl graphfile outpath pathtonw secsim strsim mode numdifftypes outfolder match dupl ins pseudo pseumis pseudels pseuins dels missing newicktree pathoTemp summary remodlingsout inremoldingsout

##program will produce sequences of letters, where the same letter means similar sequences, depending on the threshold.
##thus one letter for each connected component.
##then sort the elements by coordinate for each species
##get the sequence of one letters and put it in altnw for each pair of species.
##give own path for summary file because you dont want to have it in the same folder as all the alignments (many files).
## secsim and strucsim give the similarity thresholds for the sequence and structure, if one of them (and only one!) is -1, do not pay attention to this one, take only the other one into account
## singletoncount is the numbers of elements for each species that were included in the singleton clusters that were not part of the graph analysis and further. Now, they will be included in the geneticEvents file.


use Data::Dumper;
use strict;
use warnings;
use List::Util qw(any none first min max);
use POSIX qw(strftime);


##    my $cmd30 = "$perlpath\/perl $scripts_sarah\/createAlignments.pl $outpath\/graphs/edlilist $outpath\/graphs/alignments $altnwpath $seqsim $strucsim -1 -1 $singletoncount $outpath\/graphs/GainLoss $outpath\/geneticEvents.txt $sumcreatealn 2>>$err";

my $file = shift;
my $outpath=shift;
my $pathtonw = shift;
my $seqlim = shift;
my $strlim = shift;
#my $nonecount = shift; #I don't need those counts here
#my $pseudocount = shift;
#my $singletoncount = shift;
my $mode = shift;
my $numdifftypes = shift;
#my $outfolder = shift; ##write files for ePoPE #not used at the moment!
#my $matchout = shift;
#my $duplout = shift;
#my $insout = shift;
#my $pseout = shift;
#my $psemisout = shift;
#my $psedelout = shift;
#my $pseinsout = shift;
#my $delout = shift;
#my $misout = shift;
my $nwtree = shift;
#my $path2Temp = shift; #temp folder which contains gipped lists about anchors for each species
#my $summary = shift; #usual summary after running the script
#my $sumrems = shift;
#my $suminrems = shift;
my $leftanchor = shift;
my $rightanchor = shift;
my $del2check = shift;
my $pseudel2check = shift;
my $toprint = shift;
my $alnfile = shift;

my $aout;
if($toprint == 1){
    open($aout,">>",$alnfile);
}

#print STDERR "createALN input tree: $nwtree\n";

#create a hash with strings that show spec1_spec2,spec1,spec3,spec5 as key whereas
#spec1 is the current species and the others are in the same cluster.
#in this way, genetic events occuring in just this combination of species can
#be counted

my %dupevents =();
my %matevents =();
my %insevents =();
my %misevents =(); #missing data=no anchors
my %delevents = ();#deletions
my %psevents = (); #insertions of pseudogenes
my %psedels = ();
my %psemis = ();
my %pseins = ();


my @remoldings = ();
# The pairs of elements (separated with ':') are defined 
#as orthologs based on the similarity score but have distinct types according to the input";

my @inremoldings = ();
#The pairs of elements are defined as orthologs as they have the same types but the similarity is below the orthology threshold.
#this is only reported IFF there are more than just one type!

my %elemcount = (); #count how many elements of which species occur during the creation of alignments

open(my $outd2c,">>", $del2check);
open(my $outp2c,">>", $pseudel2check);

##filter which is the biggest alignment created
#my $maxlen = 0;
#my $maxlenname="";


#curfile are the .edli files, such that the user can run the analysis again and change the thresholds
#if(-z $file){
#    print $outs "no sequences found that fulfill the similarity thresholds, thus, only graphs without edges and no alignments. Lower the similarity thresholds for sequences and/or secondary structure to see alignments (options -s, -t). \n";
#    exit 0;
#}


##Counting: build graphs based on duplication alignments
##get connected components and use them as counting "entity"
##find neighbouring species in the tree and try to add them
##walk along the tree (top to bottom) and add insertions as far as possible, add losses later on if some species do not exist
##every node got a name, write numbers into a hash

my $alncount = 0;

##analyse files that are output for ePoPE, thus count them and now how many CCs is the most in which aln.
my $totALNnum = 0; ##num of alignments that are checked for ePoPE
my $totCCnum = 0; ##num of CCs for ePoPE in total
my $maxCC = 0;
my $maxCCaln = "";
my $maxCCnum = 0;
my $maxCCnumaln = "";


    

#    print $outgr ">$almostname\n";
    $alncount++;
    ##one letter codes needed
    ##are all the chars ok? ##no - or ~ as they appear in the alignment
    my @letters = ('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','0','1','2','3','4','5','6','7','8','9','!','@','#','$','%','^','&','*','(',')','=','+','>','<','?','[',']','|'); 
    my $letcount = 0;
    my $letnum = scalar @letters;
   
    my %node2letter = ();
    my %species = ();
    my %start2node = ();
    my %pgenes = (); ##pseudogenes
    my %spec2pseudo = ();   
    my @vertices = ();
    my @arcs = ();


#    print "file: $curfile \n";
    open FA,"<$file" or die "can't open $file\n";

    while(<FA>){
	chomp;

	my $line = $_;
	if($line eq ""){next;}
#	print STDERR "LINE: $line \n";
	my @F = split ' ', $line;
	my $n1 = $F[0]; #node 1
	my $n2 = $F[1]; #node 2
	my $seqsim = $F[2];
	my $strsim = $F[3];

	my @G = split '_', $n1;
	my @G2 = split '_', $n2;
	my $isp1 = 0; #check if genes are pseudo or not
	my $isp2 = 0;

	if($mode == 0){
	    #nodes look like: chr_spec_id_start_end_strand_pseudo
	    my $spec = $G[(scalar @G) - 6];
	    my $spec2 = $G2[(scalar @G2) - 6];
	    $species{$spec} = "";
	    $species{$spec2} = "";
	    $spec2pseudo{$spec} = "";
	    $spec2pseudo{$spec2} = "";
	    ##if we find a pseudogenes, add it to the pseudogene vector and remove it from the analysis
	    my $p1 = $G[(scalar @G) - 1];	    
	    if($p1 eq "P"){if(exists($pgenes{$n1})){}else{$pgenes{$n1}=1;} $isp1 = 1;}
#	    else{
		my $startvec = $G[(scalar @G) - 4] + (0.0001 * $G[(scalar @G) - 5]);
		if(exists $start2node{$startvec}){}
		else{$start2node{$startvec} = $n1;}
#	    }
	    my $p2 = $G2[(scalar @G2) - 1];
	    if($p2 eq "P"){if(exists($pgenes{$n2})){}else{$pgenes{$n2}=1;}$isp2 = 1;}
#	    else{	
		my $startvec2 = $G2[(scalar @G2) - 4] + (0.0001 * $G2[(scalar @G2) - 5]);
		if(exists $start2node{$startvec2}){}
		else{$start2node{$startvec2} = $n2;}
#	    }
	}
	else{
	    #nodes look like: chr_spec_id_start_end_strand_type_pseudo
	    #	    chr2A_gorGor3_276_9333523_9333524_-_Undet_TRUE
#	    chr6_GL000252v2_alt_hg38_1276_10814083_10814084_-_Ala_FALSA
	    my $spec = $G[(scalar @G) - 7];
	    $species{$spec} = "";
	    my $spec2 = $G2[(scalar @G2) - 7];
	    $species{$spec2} = "";
	    $spec2pseudo{$spec} = "";
	    $spec2pseudo{$spec2} = "";
	    ##if we find a pseudogenes, add it to the pseudogene vector and remove it from the analysis
	    my $p1 = $G[(scalar @G) - 1];
	    my $p2 = $G2[(scalar @G2) - 1];
	    if($p1 eq "T" || $p1 eq "t" || $p1 eq "True" || $p1 eq "true"
	       || $p1 eq "1" || $p1 eq "TRUE"){if(exists($pgenes{$n1})){}else{$pgenes{$n1}=1;}$isp1 = 1;}
#	    else{
		my $startvec = $G[(scalar @G) - 5] + (0.0001 * $G[(scalar @G) - 6]);
		if(exists $start2node{$startvec}){}
		else{$start2node{$startvec} = $n1;}
#	    }
	    if($p2 eq "T" || $p2 eq "t" || $p2 eq "True" || $p2 eq "true"
	       || $p2 eq "1" || $p2 eq "TRUE"){if(exists($pgenes{$n2})){}else{$pgenes{$n2}=1;} $isp2 = 1;}
#	    else{	    
		my $startvec2 = $G2[(scalar @G2) - 5] + (0.0001 * $G2[(scalar @G2) - 6]);
		if(exists $start2node{$startvec2}){}
		else{$start2node{$startvec2} = $n2;}
#	    }
	    ##check types:
	    my $t1 = $G[(scalar @G) - 2];
	    my $t2 = $G2[(scalar @G2) - 2];
	    if($numdifftypes > 0){
		if($numdifftypes > 1 && $seqsim >= $seqlim && $strsim >= $strlim && $t1 ne $t2 && $isp1 == 0 && $isp2 == 0){
		    ##sequence similarity is high but types seem to be different
		    my $remstr = "";
		    if($n1 le $n2){$remstr = "$n1\:$n2";}
		    else{$remstr = "$n2\:$n1";}
		    if(grep( /^$remstr$/, @remoldings)){}
		    elsif($remstr ne ""){push @remoldings, $remstr;}
		    else{}
		}
		if($numdifftypes > 1 && $t1 eq $t2 && $seqsim < $seqlim && $isp1 == 0 && $isp2 == 0){
		    ##equal types but low sequence similarity
		    my $inremstr = "";
		    if($n1 le $n2){$inremstr = "$n1\:$n2";}
		    else{$inremstr = "$n2\:$n1";}
		    if(grep( /^$inremstr$/, @inremoldings)){}
		    elsif($inremstr ne ""){push @inremoldings, $inremstr;}
		    else{}
		    
		}
	    }
	}

#	if($isp1 == 1 && $isp2 == 1){next;}
#	elsif($isp1 == 1 && $isp2 == 0){
#	    if(exists $node2letter{$n2}){}
#	    else{
#		$node2letter{$n2} = $letters[$letcount];
#		$letcount++;
#	    }
#	    next;
#	}
#	elsif($isp1 == 0 && $isp2 == 1){
#	    if (exists $node2letter{$n1}) {}
#	    else{
#		$node2letter{$n1} = $letters[$letcount];
#		$letcount++;
#	    }
#	    next;
#	}
#	else{}
#	##the following is the else case, as we do 'next' above, this is ok
#	##both genes are not pseudo

	if($letcount>=scalar @letters){$letcount = 0;}
	
	if($seqsim >= $seqlim && $strsim >= $strlim){
	    ##similarity fits, nodes get the same letter
	    if (exists $node2letter{$n1}) {
		if(exists $node2letter{$n2}){
		    if($node2letter{$n1} eq $node2letter{$n2}){}
		    else{$node2letter{$n2} = $node2letter{$n1};} ##should fit as we iterate over all nodes.
		}
		else{
		    $node2letter{$n2} = $node2letter{$n1};
		}
	    } 
	    else{
		if(exists $node2letter{$n2}){
		    $node2letter{$n1} = $node2letter{$n2};
		}
		else{
		    $node2letter{$n1} = $letters[$letcount];
		    $node2letter{$n2} = $letters[$letcount];
		    if($letcount > $letnum){$letcount=0;}
		    else{$letcount++;}
		}
	    }
	}
	else{
	    ##similarity does not fit, nodes get different letters.
	    if (exists $node2letter{$n1}) {}
	    else{$node2letter{$n1} = $letters[$letcount];
		 if($letcount > $letnum){$letcount=0;}
		 else{$letcount++;}
	    }
	    if (exists $node2letter{$n2}) {}
	    else{$node2letter{$n2} = $letters[$letcount];
		 if($letcount > $letnum){$letcount=0;}
		 else{$letcount++;}
	    }
	}

    }


    
    
#    print $outp "pgenes\n";
#    my $pgenestr = join(",",@pgenes);
#    print $outp "$pgenestr \n";

    my $specstr = "";
    $specstr = join(',',(keys %species));


    my $null = "0";
    my $eins = "1";


##do a special case if there is only one species!
if(scalar (keys %species) == 1){
    #start2node contains all current nodes
    #pgenes contains all pseudogenes
    #now add the pseudogenes to pseudinsertions and the remaining genes to insertions and return
    my $numpgenes = scalar (keys %pgenes);
    my $numgenes = scalar (keys %start2node);
    my $numrealgenes = $numgenes - $numpgenes;
    my @SP = keys %species;
    if(exists($insevents{$SP[0]})){$insevents{$SP[0]}+=$numrealgenes;}
    else{$insevents{$SP[0]}=$numrealgenes;}
    if(exists($pseins{$SP[0]})){$pseins{$SP[0]}+=$numpgenes;}
    else{$pseins{$SP[0]}=$numpgenes;}
}
else{



    ##for every species, sort nodes in the correct order and create the sequences to be aligned
    foreach my $k1 (sort { $a <=> $b } keys(%start2node) ) {
	my $curnode = $start2node{$k1};
#	print "curnode: $curnode \n";
	my @H = split "_", $curnode;
	my $spec = "";
	my $state = "";
	if($mode == 0){
	    $spec = $H[(scalar @H) - 6];
	    $state = $H[(scalar @H) - 1];
	}
	else{
	    $spec = $H[(scalar @H) - 7];
	    $state = $H[(scalar @H) - 1];
	}
#	print "state: $state \n";
	my $pstr = $spec2pseudo{$spec};
#	print "pstr: $pstr \n";
	##write numbers in the same order as chars in sequence, 0 for normal, 1 for pseudo
	if($mode == 0){ ##cm mode or genelist mode
	    if($state eq "N"){
		#	    print "blubb \n";
		$spec2pseudo{$spec} = "$pstr$null";
	    }
	    else{
		$spec2pseudo{$spec} = "$pstr$eins";
	    }
	}
	else{
	    if($state eq "T" || $state eq "t" || $state eq "True" || $state eq "true"
	       || $state eq "1" || $state eq "TRUE"){
		$spec2pseudo{$spec} = "$pstr$eins";
	    }
	    else{$spec2pseudo{$spec} = "$pstr$null";}
	}
	my $sstr = $species{$spec};
#	print "sstr: $sstr \n";
	$species{$spec} = "$sstr$node2letter{$curnode}";
    }
    ##add extra element (x) to spec2pseudo for being sure that it is not getting out of range
##do we really need that? it makes more problems than it helps...
#    my $minone = "x";
#    foreach my $psk (keys %spec2pseudo){
#	my $permstr = $spec2pseudo{$psk};
#	$spec2pseudo{$psk}="$permstr$minone";
#    }
    

    ##create graphs out of the duplciation alignments, edges are where the columns fit
    ##nodes are: letter_spec_id!!!
    ##edges are alphabetically sorted e1 < e2 (unique edges), written with e1 e2

    
    
    
    my $lseq1sum = 0;
    my $jump = 0;
    ##get alignments for each species against each other
    foreach my $k (keys %species){
	my $lseq1 = $species{$k};
	my $lseq1len = length($lseq1);
	$lseq1sum += $lseq1len;
	if(exists($elemcount{$k})){$elemcount{$k}+=$lseq1len;}
	else{$elemcount{$k}=$lseq1len;}
	##skip if it is the same species as there cannot be an edge
	##problem, if the cluster consists of only nodes of the same species!
	##thus, add those elements directly to the insertions
	if(scalar (keys %species) == 1){
	    my $psseq = $spec2pseudo{$k};
	    my $zcount = $psseq =~ tr/0/0/; #normal genes
	    my $ecount = $psseq =~ tr/1/1/; #pseudogenes
	    if(exists $insevents{$k}){$insevents{$k} += $zcount;}
	    else{$insevents{$k} = $zcount;}
	    if(exists($psevents{$k})){$psevents{$k}+=$ecount;}
	    else{$psevents{$k}=$ecount;}
	    $jump=1;
	    last;
	}
	if($toprint == 1){print $aout "\@$k\t$species{$k}\n";}
#	print "\@$k\t$species{$k}\n";
	my $match=0; #m
	my $pseu=0; #s
	my $dupl=0; #d
	my $ins=0; #i
	my $del=0; #l
	foreach my $k2 (keys %species){
	    if($k eq $k2) {next;}
	    my $lseq2 = $species{$k2};
	    my $lseq2len = length($lseq2);
	    my $cmd1 = "$pathtonw/altNW 1 1 \"$lseq1\" \"$lseq2\"";
#	    print "altnwcmd: $cmd1 \n";
	    my @out1 = readpipe("$cmd1");
	    if($toprint == 1){print $aout ">$k2\t";}
	    if($toprint == 1){print $aout "@out1";}
#	    print "aln:\n";
#	    print "@out1";
#	    print "\n";
	    chomp(@out1);
#	    if(length $out1[1] > $maxlen){
#		$maxlen = length $out1[1];
#		$maxlenname = $newname;
#	    }
	    ##count duplication etc for each pos with ~ or - and write it down
	    ##such that we take the max for event (all letters together)
	    my $m=0;
	    my $l=0;
	    my $s=0;
	    my $i=0;
#	    my $p =0;
	    my @ref = split '',$out1[1];
	    my @rseq = split '', $lseq1;
	    my @oth = split '',$out1[2];
	    my @oseq = split '', $lseq2;
	    if(scalar @ref != scalar @oth){print STDERR "alignment does not fit!\n";}
	    my $rcount =0;
	    my $ocount=0;
	    for(my $r=0;$r<scalar @ref;$r++){
		if($ref[$r] eq '-' || $ref[$r] eq '~'){}
		else{$rcount++;}
		if($oth[$r] eq '-' || $oth[$r] eq '~'){}
		else{$ocount++;}
	    }
	    #compare sequence lengths with seq len in alignment
	    if($lseq1len != $rcount){print STDERR "lengths don't fit! sequence:$lseq1len, alnseq: $rcount\n";}
	    if($lseq2len != $ocount){print STDERR "lengths don't fit! sequence:$lseq2len, alnseq: $ocount\n";}

#	    print "spec2pseudo k,k2: $spec2pseudo{$k},$spec2pseudo{$k2} \n";
	    my @sp12pseudo = split '', $spec2pseudo{$k};  ##this is the sequence for species k for pseudogene or not
	    my @sp22pseudo = split '', $spec2pseudo{$k2};	    
	    my $tild = "~";
	    my $mins = "-";
	    my $rz = 0;
	    my $oz = 0;
	    #should have the same length as it is an alignment
	    #build a graph between the chars of the aligned sequences
	    for(my $z=0;$z < scalar @ref;$z++){
		#create nodes and edges
		my $curc = $ref[$z];
		my $v1;
		my $v2;
		my $ed;
		#do not declare pe1 and p2 here, use the array index directly when needed?
		my $pe1 = $sp12pseudo[$rz];
		my $pe2 = $sp22pseudo[$oz];
		
#		print "pe1: $pe1, pe2: $pe2 \n";
		if($curc eq $oth[$z])
		{
		    ##add $z behind the nodes names to make them unique
		    $v1 = "$curc\_$k\_$pe1\_$rz";
		    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		    if($v1 lt $v2){$ed = "$v1 $v2";}
		    else{$ed = "$v2 $v1";}
		    if(none {$_ eq $ed} @arcs){push @arcs, $ed;}
		    $rz++;
		    $oz++;

		}
		elsif($curc eq $tild){
		    ##add an edge from the last letter in the current to the current letter in the other

#		    my $zz = $rz;
		    my $pe1p = $sp12pseudo[$rz-1];
#		    print STDERR "Duplication of ref, $k: $out1[1], rz: $rz, oth, $k2: $out1[2], oz: $oz, z: $z\n ";
		    my $smrz = $rz-1;
		    $v1 = "$rseq[$smrz]\_$k\_$pe1p\_$smrz"; ##v1 should already exist!
		    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1; print STDERR "v1 should exist! $v1 \n";}
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		    if($v1 lt $v2){$ed = "$v1 $v2";}
		    else{$ed = "$v2 $v1";}
		    if(none {$_ eq $ed} @arcs){push @arcs, $ed;}
		    $oz++;

#		    while($zz >= 0){
#			if($ref[$zz] ne $tild){
#			    last;
#			}
#			$zz--;
#		    }
		}
		elsif($curc eq $mins){

		    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		    $oz++;
		}
		elsif($oth[$z] eq $mins){

		    $v1 = "$curc\_$k\_$pe1\_$rz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		    $rz++;
		}
		elsif($oth[$z] eq $tild){

		    #		    my $zz2 = $oz;
		    my $pe2p = $sp22pseudo[$oz-1];
#		    print STDERR "Duplication of oth, $k2: $out1[2]; oz: $oz, ref, $k: $out1[1], rz: $rz, z: $z\n ";
		    $v1 = "$curc\_$k\_$pe1\_$rz";
		    my $smoz = $oz-1;
		    $v2 = "$oseq[$smoz]\_$k2\_$pe2p\_$smoz"; ##v2 should already exist!
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2; print STDERR "v2 should exist! $v2 \n";}
		    if($v1 lt $v2){$ed = "$v1 $v2"; }
		    else{$ed = "$v2 $v1";}
		    if(none {$_ eq $ed} @arcs){push @arcs, $ed;}
		    $rz++;

#		    while($zz2 >= 0){
#			if($oth[$zz2] ne $tild){
#			    my $pe2p = $sp22pseudo[$zz2];
#			    last;
#			}
#			$zz2--;
#		    }
		}
		else{
		    $v1 = "$curc\_$k\_$pe1\_$rz";
		    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		    $rz++;
		    $oz++;

		}
	    }
	}

    }
    ##if jump==1, no graph could be built
    if($jump==1){
	print STDERR "No graph could be build..\n";
    }

    
    ##graph for this cluster is built.
    ##get connected components
    my $vertices = join(',',@vertices);
    my $curvertnum = scalar @vertices;
    if($lseq1sum != $curvertnum){print STDERR "Sequence sum does not fit vertex num: seq: $lseq1sum, vertex: $curvertnum, vertices: $vertices\n";}
    my $arcs = join(',',@arcs);
 #   print "vertices: $vertices \n";
    #   print "arcs: $arcs \n";


    my @CCs = connectedComponents($vertices,$arcs);
    
    ##count events for each CC, build ePoPE input for each CC
    my $perccnum = scalar @CCs; #add all nodes of each cc for this graph to see if it fits the curvertnum
    #   print "ccnum: $ccnum \n";
    my $ccnodecount = 0;
    my $cccount = 0;
    for(my $ii = 0; $ii < scalar @CCs; $ii++){
#	print "cc: $CCs[$ii] \n";
	my %spe2count = ();
	my %pse2count = ();
	my @verts = split ' ', $CCs[$ii];
	$ccnodecount += scalar @verts;
#	my $outpo;
#	if(scalar @verts >= 2){
#	    my $outname = "$newname\_$cccount\.ali";
#	    my $path2 = "$outfolder\/$outname";
#	    open($outpo,">>",$path2);
#	    if($cccount == 0){$totALNnum++;}
#	    $cccount++;
#	    $totCCnum++;
#	    if($cccount > $maxCCnum){
#		$maxCCnum = $cccount;
#		$maxCCnumaln = $newname;
#	    }
#	    if(scalar @verts > $maxCC){
#		$maxCC = scalar @verts;
#		$maxCCaln = $outname;
#	    }
#	}	
	##what about singleton clusters?
	##In order to make the program faster, 
	##only multigene clusters should be used with ePoPE, 
	##everything else can be added by hand.
	##in the counting file, everything from the graphs is included
#	print "fill pse2count \n";
	for(my $jj=0;$jj< scalar @verts; $jj++){
#	    print "$verts[$jj] \n";
	    my @SP = split '_', $verts[$jj];
	    my $spe = $SP[1];
	    my $pseudi = $SP[2];
#	    if(scalar @verts >= 2){	
#		my $word = "sequence";
#		my $outst = "$spe $word\n";
#		print $outpo $outst;
	    #	    }


	    if($pseudi eq "1"){
		if(exists $pse2count{$spe}){
		    $pse2count{$spe}++;
		}
		else{
		    $pse2count{$spe}=1;
		}
	    }
	    else{
		if(exists $spe2count{$spe}){
		    $spe2count{$spe}++;
		}
		else{
		    $spe2count{$spe}=1;
		}
	    }
	    
	}


#	my @kp = keys %pse2count;
#	my $kplen = scalar @kp;
#	print "size of pse2count: $kplen \n";
#	if(scalar @kp >0){
#	    print $outp "hash pse2count \n";
#	    print $outp Dumper(\%pse2count);
#	    print $outp "\n";
#	}
#	print "count now!\n";
	
	##counting: only matches and dup (and pseudogenes?), losses later when adding to the tree
	##add elements from the singleton clusters that were sorted
	##and specify the number of elements per species for the None cluster
	my $spstr = join(',',sort (keys %spe2count));
	my $spnum = scalar (keys %spe2count);	
#	print "spstr: $spstr \n";
	my @vals = values %spe2count;
	my $vnum = scalar @vals;
	my $mini = min @vals;
#	print "hash spe2count \n";
#	print Dumper(\%spe2count);
#	print "\n";
	#	print "num vals: $vnum\n";
	#	my @smallvals = splice(@vals,1);
	if(scalar @vals == 0){}
	elsif(scalar @vals == 1){##singleton/insertion
	    if(exists $insevents{$spstr}){$insevents{$spstr} += $vals[0];}
	    else{$insevents{$spstr} = $vals[0];}
	}
	elsif(none {$_ != $vals[0]} @vals)#check if all values are equal
	{
	    if(exists $matevents{$spstr}){$matevents{$spstr} += $vals[0];}
	    else{$matevents{$spstr} = $vals[0];}
	}
	else{##count duplications

#	    print "min vals $mini\n";
	    if(exists $matevents{$spstr}){$matevents{$spstr} += $mini;}
	    else{$matevents{$spstr} = $mini;}
	    foreach my $kk (keys %spe2count){
		my $diff = $spe2count{$kk} - $mini; 
#		print "diff: $diff \n";
		if($diff > 0){
#		    print "spe2count k and v: $kk, $spe2count{$kk}\n";
		    if(exists $dupevents{$kk}){$dupevents{$kk} += $diff;}
		    else{$dupevents{$kk} = $diff;}
#		    print "duplication events: $kk,$dupevents{$kk}\n";
		}
	    }
	}
	my @pseuvals = values %pse2count;
	my $psnum = scalar @pseuvals;
	my $psestr = "";
	my $pmin = min @pseuvals;
	if(scalar @pseuvals > 0){
	    $psestr = join(',',sort (keys %pse2count));
	    ##do the same distinguishing for the pseudogenes as for matching
	    if(scalar @pseuvals == 1){#singleton
		if(exists $pseins{$psestr}){$pseins{$psestr} += $pseuvals[0];}
		else{$pseins{$psestr} = $pseuvals[0];}
	    }
	    elsif(none {$_ != $pseuvals[0]} @pseuvals)#all the same
	    {
		if(exists $psevents{$psestr}){$psevents{$psestr} += $pseuvals[0];}
		else{$psevents{$psestr} = $pseuvals[0];}
	    }
	    else{#add the minimum common value for the complete psestr and the remainings as singletons, as the duplications have been counted already
		$pmin = min @pseuvals;
		if(exists $psevents{$psestr}){$psevents{$psestr} += $pmin;}
		else{$psevents{$psestr} = $pmin;}
		foreach my $ppp (keys %pse2count){
		    my $pdiff = $pse2count{$ppp} - $pmin;
		    if($pdiff > 0){
			if(exists $psevents{$ppp}){$psevents{$ppp} += $pdiff; }
			else{$psevents{$ppp} = $pdiff;}	    
		    }
		}
	    

	    }
	}
	#do check for missing data or deletions here
	
	##if cluster >= 2 nodes && num_species >= 2
	##check if all species below the LCA appear in this cluster
	##if a species does not appear, check if
	##the missing species have the anchoring blocks
	##if yes: write the combination of species as a deletion(like matches)
	##if no: ignore this element in the species (add to missing data list)
	
	#do this for every CC as we are looking at genetic events of homologs
	#steps:
	#collect set of species
	#missing species = find LCA and missing species below (sub)
	#check if anchors exist
	#thus here: counting of deletion events and missing data
	if($spnum > 1){
	    #species: comma-separated in $spstr
	    my @missingspecs = getMissingSpecs($spstr,$nwtree);
	    if(scalar @missingspecs > 0){
		my $missstr = join(',',@missingspecs);
		my $d2cstr = "$leftanchor\_$rightanchor\t$missstr\t$mini\n";
		print $outd2c $d2cstr;

	    }
	}

	#do the same for pseudogenes in order to detect missing data and deletions for pseudogenes
	#extra files for pseudo: -singletons, ins, del, ex, mis


	if($psnum > 0){
	    	    #species: comma-separated in $psestr
	    my @pmissingspecs = getMissingSpecs($psestr,$nwtree);
	    if(scalar @pmissingspecs > 0){
		my $pmissstr = join(',',@pmissingspecs);
		my $p2cstr = "$leftanchor\_$rightanchor\t$pmissstr\t$pmin\n";
		print $outp2c $p2cstr;
	    }
	}
	

	
    }#done check for every CC
    #check if ccnodecount == curvertnum
    if($ccnodecount != $curvertnum){
	print STDERR "numbers don't fit! previous nodenum: $curvertnum, now: $ccnodecount\n";
    }

    ##pseudogenes here!
#    foreach my $k (keys %pgenes){
#	my @Ktmp = split '_', $k;
#	my $kspec;
#	if($mode == 0){
#	    $kspec = $Ktmp[-6];
#	}
#	else{
#	    $kspec=$Ktmp[-7];
#	}
#	if(exists($psevents{$kspec})){$psevents{$kspec}++;}
#	else{$psevents{$kspec}=1;}
 #   }
    


#my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;



##shift the singletoncount to countEvents.pl
##include insertion events from the singleton clusters that were not included in the graph analysis, as they were sorted out before
#my @SC = split "=", $singletoncount;
#for(my $s = 0; $s < scalar @SC; $s++){
#    if($SC[$s] eq ""){next;}
#    my @stmp = split "-", $SC[$s];
#    my $tmpnum = $insevents{$stmp[0]};
#    if(exists $insevents{$stmp[0]}){$insevents{$stmp[0]} = $tmpnum + $stmp[1];}
#    else{$insevents{$stmp[0]} = $stmp[1];}
#}



#my @PC = split "=", $pseudocount;
#for(my $pc = 0;$pc < scalar @PC; $pc++){
#    if($PC[$pc] eq ""){next;}
#    my @ptmp = split "-", $PC[$pc];
#    if(exists($psevents{$ptmp[0]})){$psevents{$ptmp[0]}+=$ptmp[1];}
#    else{$psevents{$ptmp[0]}=$ptmp[1];}
#}

}#end else if more than one species

##sort the entries in each entry for the hashes alphabetically

my $outstring = "";
#instead of writing to files, we put all information in a string and return it
#there is not too much information as its only the counts for one graph

foreach my $du (sort keys %dupevents) {
    my @devs = split ',', $du;
    @devs = sort @devs;
    my $dustr = join(',',@devs);
    $outstring = "$outstring\=$dustr\-$dupevents{$du}";
}
$outstring = "$outstring\n";

foreach my $ma (sort keys %matevents) {
    my @mevs = split ',', $ma;
    @mevs = sort @mevs;
    my $mastr = join(',',@mevs);
    $outstring = "$outstring\=$mastr\-$matevents{$ma}";
}
$outstring = "$outstring\n";


foreach my $in (sort keys %insevents) {
    my @ievs = split ',', $in;
    @ievs = sort @ievs;
    my $instr = join(',',@ievs);
    $outstring = "$outstring\=$instr\-$insevents{$in}";
}
$outstring = "$outstring\n";


foreach my $pin (sort keys %pseins) {
    my @pievs = split ',', $pin;
    @pievs = sort @pievs;
    my $pinstr = join(',',@pievs);
    $outstring = "$outstring\=$pinstr\-$pseins{$pin}";
}
$outstring = "$outstring\n";


foreach my $mi (sort keys %psevents) {
    my @mevs = split ',', $mi;
    @mevs = sort @mevs;
    my $mistr = join(',',@mevs);
    $outstring = "$outstring\=$mistr\-$psevents{$mi}";
}
$outstring = "$outstring\n";

foreach my $dl (sort keys %delevents) {
    my @dls = split ',', $dl;
    @dls = sort @dls;
    my $dlstr = join(',',@dls);
    $outstring = "$outstring\=$dlstr\-$delevents{$dl}";
}
$outstring = "$outstring\n";

foreach my $ms (sort keys %misevents) {
    my @miss = split ',', $ms;
    @miss = sort @miss;
    my $msstr = join(',',@miss);
    $outstring = "$outstring\=$msstr\-$misevents{$ms}";
}
$outstring = "$outstring\n";


foreach my $pm (sort keys %psemis) {
    my @pmiss = split ',', $pm;
    @pmiss = sort @pmiss;
    my $pmsstr = join(',',@pmiss);
    $outstring = "$outstring\=$pmsstr\-$psemis{$pm}";
}
$outstring = "$outstring\n";

foreach my $pd (sort keys %psedels) {
    my @pdls = split ',', $pd;
    @pdls = sort @pdls;
    my $pdlstr = join(',',@pdls);
    $outstring = "$outstring\=$pdlstr\-$psedels{$pd}";
}
$outstring = "$outstring\n";



if(scalar @remoldings > 0){
#    print $outsr "The following pairs of elements (separated with ':') are defined 
#as orthologs based on the similarity score but have distinct types according to the input:\n";
    for(my $i=0;$i< scalar @remoldings; $i++){
	$outstring = "$outstring$remoldings[$i]\=";
    }
}
$outstring="$outstring\n";

if(scalar @inremoldings > 0){
#    print $outsi "The following pairs of elements (separated with ':') are defined 
#as orthologs as they have the same types but they sequence similarity is below the given threshold:\n";
    for(my $i=0;$i< scalar @inremoldings; $i++){
	$outstring = "$outstring$inremoldings[$i]\=";
    }
}

$outstring="$outstring\n";

print $outstring;

sub connectedComponents{

    my @inp = @_;
    my @uniqedges = split ',', $inp[1];
    my @permnodes1 = split ',', $inp[0];
    my @permnodes2=();

    for(my $i=0;$i<scalar @uniqedges;$i++){
	my @E = split ' ', $uniqedges[$i];
	my $n1 = $E[0];
	my $n2 = $E[1];
#	print "n1,n2: $n1, $n2\n";
	my $id1 =-1;
	my $id2 =-1;
#	print "curedge: $uniqedges[$i]\n";
	for(my $j=0;$j<scalar @permnodes1;$j++)
	{
#	    print "curnode: $permnodes1[$j]\n";
	    
	    if(index($permnodes1[$j],$n1)!= -1) #remember current position of node
	    {
		$id1 = $j;
	    }
	    if(index($permnodes1[$j],$n2)!= -1) #remember current position of node
	    {
		$id2 = $j;
	    }
	}
	#	print "id1,id2: $id1, $id2\n";
	for(my $k=0;$k<scalar @permnodes1;$k++)
	{
	    if($k != $id1 && $k != $id2)#node is not in this edge
	    {
		push @permnodes2, $permnodes1[$k];
	    }
	}
	if($id1 != $id2){#node occur in this edge and are in one cc
	    my $str = "$permnodes1[$id1] $permnodes1[$id2]";
	    push @permnodes2, $str;
	}
	else{#node do not belong to the same cc
	    push @permnodes2, $permnodes1[$id1];
	}

	@permnodes1 = ();
	@permnodes1 = @permnodes2;
	@permnodes2 = ();
    }
    
    return @permnodes1;
}



#find LCA
#get diffset of leafs_under_LCA - spstr
sub getMissingSpecs{
    my @inp = @_;
    my $spstr = $inp[0];
    my $treefile = $inp[1];

    open TF,"<$treefile" or die "can't open $treefile\n";
    
    my $tree="";
    
    while(<TF>){
	chomp;
	$tree = $_;
	last;
    }


    
    my @output = ();
    my @species = split ',', $spstr;
    if(scalar @species <= 1){return @output;}


    if($tree eq ""){print STDERR "createAlignments: tree format doesn't fit!\n"; exit 1;}
    
    #print $outs "tree: $tree \n";
    
    ##split the tree into an array of its elements
    my @T = (); ##tree with its components
    my @N = (); ##at each position there is a number showing opening brackets - closing brackets before this position, except for ( ) , ; then -2
    my @L = (); #leaves
    my @Lids = ();
    my @tr = split '', $tree;
    my $brackets = 0;
    my $tmp = "";
    for(my $i = 0; $i < scalar @tr; $i++)
    {
	if($tr[$i] eq ')' || $tr[$i] eq '(' || $tr[$i] eq ',' || $tr[$i] eq ';'){
	    if($tmp ne ""){
		push @T, $tmp; 
		push @N, $brackets;
		#	    print "leaves? $T[(scalar @T) -2] \n";
		if($T[(scalar @T) -2] ne ")"){ #leaves
		    push @L, $tmp;
		    push @Lids, (scalar @T) -1;
		}
		$tmp="";
	    }
	    push @T, $tr[$i];
	    push @N, -2;
	    if($tr[$i] eq '('){$brackets++;}
	    if($tr[$i] eq ')'){$brackets--;}
	}
	else{
	    $tmp = "$tmp$tr[$i]";
	}
    }

    my @leafbutnospec = (); #insert ids of leaves in T: Lids
    for(my $l=0;$l< scalar @L;$l++){
	my $isnotspec = 1;
	for(my $s = 0;$s < scalar @species;$s++){
	    if($species[$s] eq $L[$l]){
		$isnotspec = 0;
		last;
	    }
	}
	if($isnotspec == 1){
	    push @leafbutnospec, $Lids[$l];
	}
    }
    
    my $treestr = join('=',@T);
    my $leafbuitnotspecstr = join('=',@leafbutnospec);
    #print STDERR "createALN findLCA leafbutnotspec: $leafbuitnotspecstr \n";
    #print STDERR "createALN findLCA: $spstr \n";
    
    #print STDERR "createALN findLCA: $treestr \n";
    my $lca = findLCA($treestr,$spstr);
    #print STDERR "createALNs lca: $lca \n";
    #lca is the id in T
    #check for all elements in leafbutnospec if lca is on their path to the root
    #if yes, put in output
    for(my $z=0;$z<scalar @leafbutnospec;$z++){
	my $num = $N[$leafbutnospec[$z]];
	for(my $t=$leafbutnospec[$z];$t<scalar @T;$t++){
	    if($N[$t] != $num){next;}
	    if($T[$lca] eq $T[$t]){
		push @output, $T[$leafbutnospec[$z]];
		last;
	    }
	    $num--;
	    if($num < 0){last;}
	}
    }

    my $outputstr = join('=',@output);
    #print STDERR "createAlns missingspecs $outputstr\n";
    return @output;

}





sub findLCA{
    my @inp= @_;
    #tree and leaf string is separated by =
    my @T = split '=', $inp[0];
    my @L = split ',', $inp[1]; #names of the species
    my @Ltmp = split ',', $inp[1];
    
    my @output = ();
    my $rootid = -1;
    
    my @N = ();
    my @allleaves = ();
    my $brackets = 0;
    my $maxbracket = 0;
    for(my $i = 0; $i < scalar @T; $i++)
    {
	if($T[$i] eq ')' || $T[$i] eq '(' || $T[$i] eq ',' || $T[$i] eq ';'){
	    
	    push @N, -2;
	    if($T[$i] eq '('){$brackets++;}
	    if($T[$i] eq ')'){$brackets--;}
	    if($brackets > $maxbracket){$maxbracket = $brackets;}
	}
	else{
	    push @N, $brackets;
	    if($brackets == 0){$rootid = $i;}
	    if($i>0){
		if($T[$i-1] eq ')'){
		    ##this is an inner node
		}
		else{
		    push @allleaves, $T[$i];
		}
	    }	
	}
    }

    if($rootid == -1){return @output;}
    
    #print join(" ",@T);
    #print "\n";
    #print join(" ",@N);
    #print "\n";
    #print join(" ",@allleaves);
    #print "\n";
    my @Lids = (); #ids of the leaves in T and N
    for(my $l=0;$l<scalar @L;$l++){
	for(my $t=0;$t<scalar @T;$t++){
	    if($L[$l] eq $T[$t]){
		push @Lids, $t;
	    }
	}	
    }

    my @Ltmpids = @Lids;
    my @L2tmp = ();
    my @L2tmpids = ();
    
    #######################################
    my $pstr1="";
    my $pstr2="";

    ##try to find curlca by always looking at two species at the same time
    ##then, get set of nodes from the pathes from the current leaves to the root and intersect them
    ##if there are several nodes in the intersection, take the one with the highest N
    my $curlca = $Lids[0]; #start with species at $S[0]
    for(my $k=1; $k < scalar @L; $k++){
	my @path1 = ();
	my @path2 = ();
	push @path1, $curlca;
	push @path2, $Lids[$k];
	my $num1 = $N[$curlca];
	for(my $t4=$curlca;$t4<scalar @T;$t4++){
	    #got up in the tree and check the bracketnums
	    if($N[$t4] != $num1){next;}
	    $num1--; #bracketnum
	    if($num1 < 0){push @path1, $rootid;last;}
	    push @path1, $t4;
	}

	my $num2 = $N[$Lids[$k]];
	for(my $t5=$Lids[$k];$t5<scalar @T;$t5++){
	    if($N[$t5] != $num2){next;}
	    $num2--;
	    if($num2 < 0){push @path2, $rootid;last;}
	    push @path2, $t5;
	}

	if(scalar @path1 == 0){
#	    print $outs "path1 zero: $F[0]\n";
	}
	else{$pstr1 = join("=",@path1);}

	if(scalar @path2 == 0){
#	    print $outs "path2 zero: $F[0]\n";
	}

	if(scalar @path1 == 0 || scalar @path2 == 0){
	    next;
	}
	else{$pstr2 = join("=",@path2);}

	


	#do the intersection of both
	##as we go from leaves to root, we take the first element that we have in common
	my @pdiff = ();
	for(my $a1 = 0; $a1 < scalar @path1;$a1++){
	    for(my $b1 = 0; $b1 < scalar @path2; $b1++){
		if($path1[$a1] == $path2[$b1]){
		    push @pdiff, $path2[$b1];
#		    last; #this can be deleted later on in order to get all overlapping nodes
		}
	    }

	}

	my $pdiffstr = join('=',@pdiff);

	#print STDERR "createALNs path1: $pstr1\n";	
	#print STDERR "createALNs path2: $pstr2\n";
	#print STDERR "createALNs pdiff: $pdiffstr\n";

	if(scalar @pdiff == 0){
	    print STDERR "son mist pathes: $pstr1; $pstr2\n";
	    #	    print $outs "son mist pathes: $pstr1; $pstr2\n";
	}
	elsif(scalar @pdiff == 1){
	    $curlca = $pdiff[0];
	}
	else{
	    #choose the curlca with the highest N
	    my $curplca = $N[$pdiff[0]];
	    my $curpid = $pdiff[0];
	    for(my $pd=1;$pd<scalar @pdiff;$pd++){
		if($N[$pdiff[$pd]] >= $curplca){
		    $curplca = $N[$pdiff[$pd]];
		    $curpid = $pdiff[$pd];
		}
	    }
	    $curlca = $curpid;
	}
    }

    if($N[$curlca] < 0){
	print STDERR "createAlignments: strange curlca\n";
    }



    return $curlca;

    


}

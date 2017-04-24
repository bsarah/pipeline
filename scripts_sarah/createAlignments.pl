#!/usr/bin/perl -w

## perl createAlignments.pl edlilist outpath pathtonw secsim strsim pseudoscore singletoncount mode numdifftypes outfolder match dupl ins pseudo summary remodlingsout inremoldingsout

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
my $pseudosc = shift;
my $singletoncount = shift;
my $mode = shift;
my $numdifftypes = shift;
my $outfolder = shift; ##write files for ePoPE
my $matchout = shift;
my $duplout = shift;
my $insout = shift;
my $pseout = shift;
my $summary = shift; #usual summary after running the script
my $sumrems = shift;
my $suminrems = shift;


#create a hash with strings that show spec1_spec2,spec1,spec3,spec5 as key whereas
#spec1 is the current species and the others are in the same cluster.
#in this way, genetic events occuring in just this combination of species can
#be counted

my %dupevents =();
my %matevents =();
my %insevents =();
##are mismatches needed? use pseudogenes INSTEAD of mistmatches
#my %misevents =();
my %psevents = ();
my @remoldings = ();
# The pairs of elements (separated with ':') are defined 
#as orthologs based on the similarity score but have distinct types according to the input";

my @inremoldings = ();
#The pairs of elements are defined as orthologs as they have the same types but the similarity is below the orthology threshold.
#this is only reported IFF there are more than just one type!

my %elemcount = (); #count how many elements of which species occur during the creation of alignments


open(my $outm,">>",$matchout);
open(my $outd,">>",$duplout);
open(my $outi,">>",$insout);
open(my $outp,">>",$pseout);
open(my $outs,">>", $summary);

open(my $outsr,">>", $sumrems);
open(my $outsi,">>", $suminrems);


##filter which is the biggest alignment created
my $maxlen = 0;
my $maxlenname="";


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


open FA,"<$file" or die "can't open $file\n";
while(<FA>){
    chomp;
    my $curfile = $_;
    my @D = split '\.', $curfile;
    my $prename = $D[(scalar @D)-2];
    my @E = split '\/', $prename;
    my $almostname = $E[(scalar @E)-1];
    my $newname="$almostname\.aln";
    my $path = "$outpath\/$newname";
    open(my $outgr,">>",$path);

    print $outgr ">$almostname\n";
    $alncount++;
    ##one letter codes needed
    ##are all the chars ok?
    my @letters = ('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','0','1','2','3','4','5','6','7','8','9','!','@','#','$','%','^','&','*','(',')','=','+','-','~'); 
    my $letcount = 0;
   
    my %node2letter = ();
    my %species = ();
    my %start2node = ();
    my @pgenes = (); ##pseudogenes
    my %spec2pseudo = ();   
    my @vertices = ();
    my @arcs = ();


#    print "file: $curfile \n";
    open CF,"<$curfile" or die "can't open $curfile\n";
    while(<CF>){
	chomp;
	my $line = $_;
#	print "line: $line \n";
	my @F = split ' ', $line;
	my $n1 = $F[0];
	my $n2 = $F[1];
	my $seqsim = $F[2];
	my $strsim = $F[3];

	my @G = split '_', $n1;
	my @G2 = split '_', $n2;
	if($mode == 0){
	    #nodes look like: chr_spec_id_start_end_strand_pseudo
	    my $spec = $G[(scalar @G) - 6];
	    my $spec2 = $G2[(scalar @G2) - 6];
	    $species{$spec} = "";
	    $species{$spec2} = "";
	    $spec2pseudo{$spec} = "";
	    my $startvec = $G[(scalar @G) - 4] + (0.0001 * $G[(scalar @G) - 5]);
	    if(exists $start2node{$startvec}){}
	    else{$start2node{$startvec} = $n1;}
	    my $p1 = $G[(scalar @G) - 1];
	    
	    my $p2 = $G2[(scalar @G2) - 1];
	    if($p1 eq "P"){push @pgenes, $n1;}
	    if($p2 eq "P"){push @pgenes, $n2;}
	}
	else{
	    #nodes look like: chr_spec_id_start_end_strand_type_pseudo
	    my $spec = $G[(scalar @G) - 7];
	    $species{$spec} = "";
	    my $spec2 = $G2[(scalar @G2) - 7];
	    $species{$spec2} = "";
	    $spec2pseudo{$spec} = "";
	    my $startvec = $G[(scalar @G) - 5] + (0.0001 * $G[(scalar @G) - 6]);
	    if(exists $start2node{$startvec}){}
	    else{$start2node{$startvec} = $n1;}
	    my $p1 = $G[(scalar @G) - 1];
	    my $p2 = $G2[(scalar @G2) - 1];
	    if($p1 eq "T" || $p1 eq "t" || $p1 eq "True" || $p1 eq "true"
	       || $p1 eq "1" || $p1 eq "TRUE"){push @pgenes, $n1;}
	    if($p2 eq "T" || $p2 eq "t" || $p2 eq "True" || $p2 eq "true"
	       || $p2 eq "1" || $p2 eq "TRUE"){push @pgenes, $n2;}

	    ##check types:
	    my $t1 = $G[(scalar @G) - 2];
	    my $t2 = $G2[(scalar @G2) - 2];
	    if($numdifftypes > 0){
		if($numdifftypes > 1 && $seqsim >= $seqlim && $strsim >= $strlim && $t1 ne $t2){
		    ##sequence similarity is high but types seem to be different
		    my $remstr = "";
		    if($n1 le $n2){$remstr = "$n1\:$n2";}
		    else{$remstr = "$n2\:$n1";}
		    if(grep( /^$remstr$/, @remoldings)){}
		    else{push @remoldings, $remstr;}
		}
		if($numdifftypes > 1 && $t1 eq $t2 && $seqsim < $seqlim){
		    ##equal types but low sequence similarity
		    my $inremstr = "";
		    if($n1 le $n2){$inremstr = "$n1\:$n2";}
		    else{$inremstr = "$n2\:$n1";}
		    if(grep( /^$inremstr$/, @inremoldings)){}
		    else{push @inremoldings, $inremstr;}
		    
		}
	    }
	}
	    
	
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
		    $letcount++;
		}
	    }
	}
	else{
	    ##similarity does not fit, nodes get different letters.
	    if (exists $node2letter{$n1}) {}
	    else{$node2letter{$n1} = $letters[$letcount];$letcount++;}
	    if (exists $node2letter{$n2}) {}
	    else{$node2letter{$n2} = $letters[$letcount];$letcount++;}
	}

    }


    
    
#    print $outp "pgenes\n";
#    my $pgenestr = join(",",@pgenes);
#    print $outp "$pgenestr \n";

    my $specstr = "";
    $specstr = join(',',(keys %species));


    my $null = "0";
    my $eins = "1";

    
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
    ##add extra 0 to spec2pseudo for being sure that it is not getting out of range
    foreach my $psk (keys %spec2pseudo){
	my $permstr = $spec2pseudo{$psk};
	$spec2pseudo{$psk}="$permstr$null";
    }
    

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
	    if(exists $insevents{$k}){$insevents{$k} += $lseq1len;}
	    else{$insevents{$k} = $lseq1len;}
	    $jump=1;
	    last;
	}
	print $outgr "\@$k\t$species{$k}\n";
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
	    print $outgr ">$k2\t";
	    print $outgr "@out1";
#	    print "aln:\n";
#	    print "@out1";
#	    print "\n";
	    chomp(@out1);
	    if(length $out1[1] > $maxlen){
		$maxlen = length $out1[1];
		$maxlenname = $newname;
	    }
	    ##count duplication etc for each pos with ~ or - and write it down
	    ##such that we take the max for event (all letters together)
	    my $m=0;
	    my $l=0;
	    my $s=0;
	    my $i=0;
#	    my $p =0;
	    my @ref = split '',$out1[1];
	    my @oth = split '',$out1[2];
	    if(scalar @ref != scalar @oth){print $outs "alignment does not fit! file: $curfile \n";}
	    my $rcount =0;
	    my $ocount=0;
	    for(my $r=0;$r<scalar @ref;$r++){
		if($ref[$r] eq '-' || $ref[$r] eq '~'){}
		else{$rcount++;}
		if($oth[$r] eq '-' || $oth[$r] eq '~'){}
		else{$ocount++;}
	    }
	    #compare sequence lengths with seq len in alignment
	    if($lseq1len != $rcount){print $outs "lengths don't fit! sequence:$lseq1len, alnseq: $rcount, file: $curfile \n";}
	    if($lseq2len != $ocount){print $outs "lengths don't fit! sequence:$lseq2len, alnseq: $ocount, file: $curfile \n";}

#	    print "spec2pseudo k,k2: $spec2pseudo{$k},$spec2pseudo{$k2} \n";
	    my @sp12pseudo = split '', $spec2pseudo{$k};
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
		    $rz++;
		    $oz++;
		    $v1 = "$curc\_$k\_$pe1\_$rz";
		    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		    if($v1 lt $v2){$ed = "$v1 $v2";}
		    else{$ed = "$v2 $v1";}
		    if(none {$_ eq $ed} @arcs){push @arcs, $ed;}
		}
		elsif($curc eq $tild){
		    ##add an edge from the last letter in the current to the current letter in the other
		    $oz++;
		    my $zz = $z-1;
		    while($zz >= 0){
			if($ref[$zz] ne $tild){
			    $v1 = "$ref[$zz]\_$k\_$pe1\_$rz"; ##v1 should already exist!
			    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
			    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
			    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
			    if($v1 lt $v2){$ed = "$v1 $v2";}
			    else{$ed = "$v2 $v1";}
			    if(none {$_ eq $ed} @arcs){push @arcs, $ed;}
			    last;
			}
			$zz--;
		    }
		}
		elsif($curc eq $mins){
		    $oz++;
		    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		}
		elsif($oth[$z] eq $mins){
		    $rz++;
		    $v1 = "$curc\_$k\_$pe1\_$rz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		}
		elsif($oth[$z] eq $tild){
		    $rz++;
		    my $zz2 = $z-1;
		    while($zz2 >= 0){
			if($oth[$zz2] ne $tild){
			    $v1 = "$curc\_$k\_$pe1\_$rz"; 
			    $v2 = "$oth[$zz2]\_$k2\_$pe2\_$oz"; ##v2 should already exist!
			    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
			    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
			    if($v1 lt $v2){$ed = "$v1 $v2";}
			    else{$ed = "$v2 $v1";}
			    if(none {$_ eq $ed} @arcs){push @arcs, $ed;}
			    last;
			}
			$zz2--;
		    }
		}
		else{
		    $rz++;
		    $oz++;
		    $v1 = "$curc\_$k\_$pe1\_$rz";
		    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		}
	    }
	}

    }
    ##if jump==1, no graph could be built
    if($jump==1){
	next;
    }

    
    ##graph for this cluster is built.
    ##get connected components
    my $vertices = join(',',@vertices);
    my $curvertnum = scalar @vertices;
    if($lseq1sum != $curvertnum){print STDERR "Sequence sum does not fit vertex num: seq: $lseq1sum, vertex: $curvertnum, vertices: $vertices, file: $curfile \n";}
    my $arcs = join(',',@arcs);
 #   print "vertices: $vertices \n";
    #   print "arcs: $arcs \n";

    ##check if all species appear in this cluster
    ##if a species does not appear, check if
    ##it is a singleton ok
    ##the missing species have the anchoring blocks
    ##if yes: write the combination of species as a deletion(like matches)
    ##if no: ignore this element in the species (add to missing data list)
    
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
	my $outpo;
	if(scalar @verts >= 2){
	    my $outname = "$newname\_$cccount\.ali";
	    my $path2 = "$outfolder\/$outname";
	    open($outpo,">>",$path2);
	    if($cccount == 0){$totALNnum++;}
	    $cccount++;
	    $totCCnum++;
	    if($cccount > $maxCCnum){
		$maxCCnum = $cccount;
		$maxCCnumaln = $newname;
	    }
	    if(scalar @verts > $maxCC){
		$maxCC = scalar @verts;
		$maxCCaln = $outname;
	    }
	}	
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
	    if(scalar @verts >= 2){	
		my $word = "sequence";
		my $outst = "$spe $word\n";
		print $outpo $outst;
	    }
	    if(exists $spe2count{$spe}){
		$spe2count{$spe}++;
	    }
	    else{
		$spe2count{$spe}=1;
	    }
#	    if($pseudi eq "0"){print "0\n";}
	    if($pseudi eq "1"){
		if(exists $pse2count{$spe}){
		    $pse2count{$spe}++;
		}
		else{
		    $pse2count{$spe}=1;
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
#	print "spstr: $spstr \n";
	my @vals = values %spe2count;
	my $vnum = scalar @vals;
#	print "hash spe2count \n";
#	print Dumper(\%spe2count);
#	print "\n";
	#	print "num vals: $vnum\n";
#	my @smallvals = splice(@vals,1);
	if(scalar @vals == 1){##singleton/insertion
	    if(exists $insevents{$spstr}){$insevents{$spstr} += $vals[0];}
	    else{$insevents{$spstr} = $vals[0];}
	}
	elsif(none {$_ != $vals[0]} @vals)#check if all values are equal
	{
	    if(exists $matevents{$spstr}){$matevents{$spstr} += $vals[0];}
	    else{$matevents{$spstr} = $vals[0];}
	}
	else{##count duplications
	    my $mini = min @vals;
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
	if(scalar @pseuvals > 0){
	    my $psestr = join(',',sort (keys %pse2count));
	    ##do the same distinguishing for the pseudogenes as for matching
	    if(scalar @pseuvals == 1){#singleton
		if(exists $psevents{$psestr}){$psevents{$psestr} += $pseuvals[0];}
		else{$psevents{$psestr} = $pseuvals[0];}
	    }
	    elsif(none {$_ != $pseuvals[0]} @pseuvals)
	    {
		if(exists $psevents{$psestr}){$psevents{$psestr} += $pseuvals[0];}
		else{$psevents{$psestr} = $pseuvals[0];}
	    }
	    else{#add the minimum common value for the complete psestr and the remainings as singletons, as the duplications have been counted already
		my $pmin = min @pseuvals;
		if(exists $psevents{$psestr}){$psevents{$psestr} += $pmin;}
		else{$psevents{$psestr} = $pmin;}
		foreach my $ppp (keys %pse2count){
		    my $pdiff = $pse2count{$ppp} - $pmin;
		    if($pdiff > 0){
			if(exists $psevents{$psestr}){$psevents{$psestr} += $pdiff;}
			else{$psevents{$psestr} = $pdiff;}	    
		    }
		}
	    

	    }
	}	
    }#done check for every CC
    #check if ccnodecount == curvertnum
    if($ccnodecount != $curvertnum){
	print $outs "numbers don't fit! previous nodenum: $curvertnum, now: $ccnodecount\n";
    }
    
    ##how to get those events into the tree?
    
}

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



##sort the entries in each entry for the hashes alphabetically

foreach my $du (sort keys %dupevents) {
    my @devs = split ',', $du;
    @devs = sort @devs;
    my $dustr = join(',',@devs);
    print $outd "$dustr\t$dupevents{$du}\n";
}


foreach my $ma (sort keys %matevents) {
    my @mevs = split ',', $ma;
    @mevs = sort @mevs;
    my $mastr = join(',',@mevs);
    print $outm "$mastr\t$matevents{$ma}\n";
}

foreach my $in (sort keys %insevents) {
    my @ievs = split ',', $in;
    @ievs = sort @ievs;
    my $instr = join(',',@ievs);
    print $outi "$instr\t$insevents{$in}\n";
}


foreach my $mi (sort keys %psevents) {
    my @mevs = split ',', $mi;
    @mevs = sort @mevs;
    my $mistr = join(',',@mevs);
    print $outp "$mistr\t$psevents{$mi}\n";
}





#my $totALNnum = 0; ##num of alignments that are checked for ePoPE
#my $totCCnum = 0; ##num of CCs for ePoPE in total
#my $maxCC = 0;
#my $maxCCaln = "";

#print $outs "===============Element counts\===============\n";
#for my $kel (keys %elemcount){
#    print $outs "$kel $elemcount{$kel}\n";
#}

print $outs "===============Duplication alignments\===============\n";
print $outs "Duplication alignments and genetic events information: \n";
print $outs "Number of clusters: $alncount. \n";
print $outs "The longest alignment has length $maxlen and is in file $maxlenname . \n";
print $outs "\n";
print $outs "For the Gain/Loss analysis, $totALNnum alignments are used that
are subdivided in total in $totCCnum connected components.
The alignment with the maximal number of connected components contains $maxCCnum connected
components and is called $maxCCnumaln. 
The connected component with most nodes contains $maxCC nodes and can be found in $maxCCaln. \n";
print $outs "The Gain/Loss analysis for multi gene cluster is done based on the files
that can be found in $outfolder \n";
print $outs "\n";
print $outs "The summary file contains information about how many genetic events were 
counted in a certain combination of species. This can be used to draw a 
phylogenetic tree with genetic events at its nodes. \n";
print $outs "The files containing the duplications alignments (.aln) for 
each cluster and pairs of species can be found here: 
$outpath \n";
print $outs "\n";
print $outs "Format of alignment files (.aln): At first the cluster and 
its species are defined. Then, connected nodes get mapped to the same 
one-letter-code (needed for the alignment). Afterwards for each species, its 
genetic elements are sorted by coordinate and depicted by the letter code. The
 letter code is aligned to the ones of each other species in the cluster. '~' 
stand for duplications, '-' for insertions or deletions in the alignment. \n";
print $outs "\n";

if(scalar @remoldings > 0){
    print $outsr "The following pairs of elements (separated with ':') are defined 
as orthologs based on the similarity score but have distinct types according to the input:\n";
    for(my $i=0;$i< scalar @remoldings; $i++){
	print $outsr "$remoldings[$i] \n";
    }
}

if(scalar @inremoldings > 0){
    print $outsi "The following pairs of elements (separated with ':') are defined 
as orthologs as they have the same types but they sequence similarity is below the given threshold:\n";
    for(my $i=0;$i< scalar @inremoldings; $i++){
	print $outsi "$inremoldings[$i] \n";
    }
}



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

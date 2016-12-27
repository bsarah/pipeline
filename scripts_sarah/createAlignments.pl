#!/usr/bin/perl -w

## perl createAlignments.pl edlilist outpath pathtonw secsim strsim psecsim pstrsim singletoncount outfolder geneticEvents summary

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
my $pseqlim = shift;
my $pstrlim = shift;
my $singletoncount = shift;
my $outfolder = shift; ##write files for ePoPE
my $sumary = shift; #geneticEvents
my $summary = shift; #usual summary after running the script



##set the lowest similarity score as threshold
my $gseqlim;
my $gstrlim;

if($pseqlim >= 0 && $pseqlim <= $seqlim){
    $gseqlim = $pseqlim;}
else{$gseqlim = $seqlim;}

if($pstrlim >= 0 && $pstrlim <= $strlim){
    $gstrlim = $pstrlim;}
else{$gstrlim = $strlim;}

#create a hash with strings that show spec1_spec2,spec1,spec3,spec5 as key whereas
#spec1 is the current species and the others are in the same cluster.
#in this way, genetic events occuring in just this combination of species can
#be counted

my %dupevents =();
my %matevents =();
my %delevents =();
my %insevents =();
##are mismatches needed? use pseudogenes INSTEAD of mistmatches
my %misevents =();
my %psevents = ();

open(my $outs,">>", $sumary);

open(my $outm,">>", $summary);

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
	my $spec = $G[(scalar @G) - 5];
	$species{$spec} = "";
	my $startvec = $G[(scalar @G) - 3] + (0.0001 * $G[(scalar @G) - 4]);
	if(exists $start2node{$startvec}){}
	else{$start2node{$startvec} = $n1;}
	
	
	if($seqsim >= $gseqlim && $strsim >= $gstrlim){
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
    

    my $specstr = "";
    $specstr = join(',',(keys %species));


    
    
    ##for every species, sort nodes in the correct order and create the sequences to be aligned
    foreach my $k1 (sort { $a <=> $b } keys(%start2node) ) {
	my $curnode = $start2node{$k1};
	my @H = split "_", $curnode;
	my $spec = $H[(scalar @H) - 5];
	my $sstr = $species{$spec};
	$species{$spec} = "$sstr$node2letter{$curnode}";
    }


    ##create graphs out of the duplciation alignments, edges are where the columns fit
    ##nodes are: letter_spec_id!!!
    ##edges are alphabetically sorted e1 < e2 (unique edges), written with e1 e2
    

    ##get alignments for each species against each other
    foreach my $k (keys %species){
	print $outgr "\@$k\t$species{$k}\n";
#	print "\@$k\t$species{$k}\n";
	my $match=0; #m
	my $mismatch=0; #s
	my $dupl=0; #d
	my $ins=0; #i
	my $del=0; #l
	foreach my $k2 (keys %species){
	    if($k eq $k2) {next;}
	    my $lseq1 = $species{$k};
	    my $lseq2 = $species{$k2};
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
	    my $d=0;
	    my @ref = split '',$out1[1];
	    my @oth = split '',$out1[2];
	    my $tild = "~";
	    my $mins = "-";
	    my $rz = 0;
	    my $oz = 0;
	    #should have the same length as it is an alignment
	    for(my $z=0;$z < scalar @ref;$z++){
		#create nodes and edges
		my $curc = $ref[$z];
		my $v1;
		my $v2;
		my $ed;
		if($curc eq $oth[$z])
		{
		    ##add $z behin the nodes names to make them unique
		    $rz++;
		    $oz++;
		    $v1 = "$curc\_$k\_$rz";
		    $v2 = "$oth[$z]\_$k2\_$oz";
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
			    $v1 = "$ref[$zz]\_$k\_$rz"; ##v1 should already exist!
			    $v2 = "$oth[$z]\_$k2\_$oz";
			    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
			    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
			    if($v1 lt $v2){$ed = "$v1 $v2";}
			    else{$ed = "$v2 $v1";}
			    if(none {$_ eq $ed} @arcs){push @arcs, $ed;}
			}
			$zz--;
		    }
		}
		elsif($curc eq $mins){
		    $oz++;
		    $v2 = "$oth[$z]\_$k2\_$oz";
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		}
		elsif($oth[$z] eq $mins){
		    $rz++;
		    $v1 = "$curc\_$k\_$rz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		}
		elsif($oth[$z] eq $tild){
		    $rz++;
		    my $zz2 = $z-1;
		    while($zz2 >= 0){
			if($oth[$zz2] ne $tild){
			    $v1 = "$curc\_$k\_$rz"; 
			    $v2 = "$oth[$zz2]\_$k2\_$oz"; ##v2 should already exist!
			    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
			    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
			    if($v1 lt $v2){$ed = "$v1 $v2";}
			    else{$ed = "$v2 $v1";}
			    if(none {$_ eq $ed} @arcs){push @arcs, $ed;}
			}
			$zz2--;
		    }
		}
		else{
		    $rz++;
		    $oz++;
		    $v1 = "$curc\_$k\_$rz";
		    $v2 = "$oth[$z]\_$k2\_$oz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		}
	    }
	}

    }
    
    ##graph for this cluster is built.
    ##get connected components
    my $vertices = join(',',@vertices);
    my $arcs = join(',',@arcs);
 #   print "vertices: $vertices \n";
 #   print "arcs: $arcs \n";
    my @CCs = connectedComponents($vertices,$arcs);

    my %spe2count = ();
    ##count events for each CC, build ePoPE input for each CC
    my $ccnum = scalar @CCs;
 #   print "ccnum: $ccnum \n";
    my $cccount = 0;
    for(my $ii = 0; $ii < scalar @CCs; $ii++){
#	print "cc: $CCs[$ii] \n";
	my @verts = split ' ', $CCs[$ii];
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
	for(my $jj=0;$jj< scalar @verts; $jj++){
	    my @SP = split '_', $verts[$jj];
	    my $spe = $SP[1];
	    if(scalar @verts >= 2){	
		my $word = "sequence";
		my $outst = "$spe\t$word\n";
		print $outpo $outst;
	    }
	    if(exists $spe2count{$spe}){
		$spe2count{$spe}++;
	    }
	    else{
		$spe2count{$spe}=1;
	    }
	}

	

#	print "count now!\n";
	
	##counting: only matches and dup (and pseudogenes?), losses later when adding to the tree
	##add elements fromt the singleton clusters that were sorted
	##and specify the number of elements per species for the None cluster
	my $spstr = join(',',sort (keys %spe2count));
#	print "spstr: $spstr \n";
	my @vals = values %spe2count;
	my $vnum = scalar @vals;
#	print "hash spe2count \n";
#	print Dumper(\%spe2count);
#	print "\n";
#	print "num vals: $vnum\n";
	if(scalar @vals == 1){##singleton/insertion
	    if(exists $insevents{$spstr}){$insevents{$spstr} += $vals[0];}
	    else{$insevents{$spstr} = $vals[0];}
	}
	elsif(none {$_ != $vals[0]} @vals)
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
    }
    
    ##how to get those events into the tree?
    
}

my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;

#print "now: $now_string \n";

print $outs "File created on $now_string\n";
print $outs "The following output shows the number of events that occur in the specified combination of species relative to the species written first.\n";
print $outs "Thus, summary for each event is a tab separated table with: reference_species species_in_cluster number_of_event\n";
print $outs "Events are Duplication, Matches, Insertion, Deletion, Mismatches.\n";
print $outs "The corresponding alignment files can be found in $outpath \n";
print $outs "\n\n";


##include insertion events from the singleton clusters that were not included in the graph analysis, as they were sorted out before
my @SC = split "=", $singletoncount;
for(my $s = 0; $s < scalar @SC; $s++){
    if($SC[$s] eq ""){next;}
    my @stmp = split "-", $SC[$s];
    if(exists $insevents{$stmp[0]}){$insevents{$stmp[0]} += $stmp[1];}
    else{$insevents{$stmp[0]} = $stmp[1];}
}



##sort the entries in each entry for the hashes alphabetically

print $outs "EVENT: Duplication\n";
foreach my $du (sort keys %dupevents) {
    my @devs = split ',', $du;
    @devs = sort @devs;
    my $dustr = join(',',@devs);
#    my $realdupnum = $dupevents{$du}/2;
    print $outs "$dustr\t$dupevents{$du}\n";
}
print $outs "\n\n";

print $outs "EVENT: Matches\n";
foreach my $ma (sort keys %matevents) {
    my @mevs = split ',', $ma;
    @mevs = sort @mevs;
    my $mastr = join(',',@mevs);
    print $outs "$mastr\t$matevents{$ma}\n";
}
print $outs "\n\n";

print $outs "EVENT: Insertion\n";
foreach my $in (sort keys %insevents) {
    my @ievs = split ',', $in;
    @ievs = sort @ievs;
    my $instr = join(',',@ievs);
    print $outs "$instr\t$insevents{$in}\n";
}
print $outs "\n\n";

print $outs "EVENT: Deletion\n";
foreach my $de (sort keys %delevents) {
    my @deevs = split ',', $de;
    @deevs = sort @deevs;
    my $destr = join(',',@deevs);
    print $outs "$destr\t$delevents{$de}\n";
}
print $outs "\n\n";

##exchange mismatches by pseudogenization?

print $outs "EVENT: Mismatches\n";
foreach my $mi (sort keys %misevents) {
    my @mevs = split ',', $mi;
    @mevs = sort @mevs;
    my $mistr = join(',',@mevs);
    print $outs "$mistr\t$misevents{$mi}\n";
}


print $outs "\n\n";



#my $totALNnum = 0; ##num of alignments that are checked for ePoPE
#my $totCCnum = 0; ##num of CCs for ePoPE in total
#my $maxCC = 0;
#my $maxCCaln = "";
	

print $outm "===============Duplication alignments\===============\n";
print $outm "Duplication alignments and genetic events information: \n";
print $outm "Number of clusters: $alncount. \n";
print $outm "The longest alignment has length $maxlen and is in file $maxlenname . \n";
print $outm "\n";
print $outm "The results for the analysis of genetic events are written to
 $sumary .
Here, all events are counted, including the insertion events due to the singleton
clusters that were not included in the graph analysis, as they only contain one node.\n";
print $outm "For the Gain/Loss analysis, $totALNnum alignments are used that
are subdivided in total in $totCCnum connected components.
The alignment with the maximal number of connected components contains $maxCCnum connected
components and is called $maxCCnumaln. 
The connected component with most nodes contains $maxCC nodes and can be found in $maxCCaln. \n";
print $outm "The Gain/Loss analysis for multi gene cluster is done based on the files
that can be found in $outfolder \n";
print $outm "\n";
print $outm "The summary file contains information about how many genetic events were 
counted in a certain combination of species. This can be used to draw a 
phylogenetic tree with genetic events at its nodes. \n";
print $outm "The files containing the duplications alignments (.aln) for 
each cluster and pairs of species can be found here: 
$outpath \n";
print $outm "\n";
print $outm "Format of alignment files (.aln): At first the cluster and 
its species are defined. Then, connected nodes get mapped to the same 
one-letter-code (needed for the alignment). Afterwards for each species, its 
genetic elements are sorted by coordinate and depicted by the letter code. The
 letter code is aligned to the ones of each other species in the cluster. '~' 
stand for duplications, '-' for insertions or deletions in the alignment. \n";
print $outm "\n";




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
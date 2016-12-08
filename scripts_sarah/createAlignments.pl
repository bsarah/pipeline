#!/usr/bin/perl -w

## perl createAlignments.pl edlilist outpath pathtonw secsim strsim geneticEvents summary

##program will produce sequences of letters, where the same letter means similar sequences, depending on the threshold.
##thus one letter for each connected component.
##then sort the elements by coordinate for each species
##get the sequence of one letters and put it in altnw for each pair of species.
##give own path for summary file because you dont want to have it in the same folder as all the alignments (many files).
## secsim and strucsim give the similarity thresholds for the sequence and structure, if one of them (and only one!) is -1, do not pay attention to this one, take only the other one into account



use Data::Dumper;
use strict;
use warnings;
use List::Util qw(any none first);
use POSIX qw(strftime);


my $file = shift;
my $outpath=shift;
my $pathtonw = shift;
my $seqlim = shift;
my $strlim = shift;
my $sumary = shift; #geneticEvents
my $summary = shift; #usual summary after running the script


#create a hash with strings that show spec1_spec2,spec1,spec3,spec5 as key whereas
#spec1 is the current species and the others are in the same cluster.
#in this way, genetic events occuring in just this combination of species can
#be counted

my %dupevents =();
my %matevents =();
my %delevents =();
my %insevents =();
my %misevents =();

open(my $outs,">>$sumary");

open(my $outm,">>$summary");



#curfile are the .edli files, such that the user can run the analysis again and change the thresholds
if(-z $file){
    print $outs "no sequences found that fulfill the similarity thresholds, thus, only graphs without edges and no alignments. Lower the similarity thresholds for sequences and/or secondary structure to see alignments (options -s, -t). \n";
    exit 0;
}

open FA,"<$file" or die "can't open $file\n";
while(<FA>){
    chomp;
    my $curfile = $_;
    my @D = split '\.', $curfile;
    my $prename = $D[(scalar @D)-2];
    my @E = split '\/', $prename;
    my $almostname = $E[(scalar @E)-1];
    my $newname="$almostname\.aln";
    open(my $outgr,">>$outpath\/$newname");

    print $outgr ">$almostname\n";

    ##one letter codes needed
    ##are all the chars ok?
    my @letters = ('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','0','1','2','3','4','5','6','7','8','9','!','@','#','$','%','^','&','*','(',')','=','+','_','-','~'); 
    my $letcount = 0;
   
    my %node2letter = ();
    my %species = ();
    my %start2node = ();
    
    open CF,"<$curfile" or die "can't open $curfile\n";
    while(<CF>){
	chomp;
	my $line = $_;
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
    

    my $specstr = "";
    $specstr = join(',',(keys %species));


    ##for every species, sort nodes in the correct order and create the sequences to be aligned
    foreach my $k1 (sort { $a <=> $b } keys(%start2node) ) {
	my $curnode = $start2node{$k1};
	my @H = split "_", $curnode;
	my $spec = $H[1];
	#my $sstr = $species{$spec};
	$species{$spec} = "$node2letter{$curnode}";
    }
    ##get alignments for each species against each other
    foreach my $k (keys %species){
	print $outgr "\@$k\t$species{$k}\n";
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
	    my @out1 = readpipe("$cmd1");
	    print $outgr ">$k2\t";
	    print $outgr "@out1";
#	    print "aln:\n";
#	    print "@out1";
#	    print "\n";
	    chomp(@out1);
	    ##count duplication etc for each pos with ~ or - and write it down
	    ##such that we take the max for event (all letters together)
	    my $m =0;
	    my $l = 0;
	    my $s=0;
	    my $i=0;
	    my $d=0;
	    my @ref = split '',$out1[1];
	    my @oth = split '',$out1[2];
	    my $tild = "~";
	    my $mins = "-";
	    #should have the same length as it is an alignment
	    for(my $z=0;$z < scalar @ref;$z++){
		my $curc = $ref[$z];
#		print "cur char: $curc \n";
		if($curc eq $oth[$z]){$m++;}
		elsif($curc eq $tild){}
		elsif($curc eq $mins){$l++;}
		elsif($oth[$z] eq $mins){$i++;}
		elsif($oth[$z] eq $tild){$d++;}
		else{$s++;}
	    }
	    if($m > $match){$match = $m;}
	    if($s > $mismatch){$mismatch = $s;}
	    if($d > $dupl){$dupl = $d;}
	    if($i > $ins){$ins = $i;}
	    if($l > $del){$del = $l;}
	}
	my $ms = "M";
	my $ss = "S";
	my $ds = "D";
	my $is = "I";
	my $ls = "L";
	print $outgr "#$k\t$ms$match$ds$dupl$is$ins$ls$del$ss$mismatch\t$specstr\n";

	#write summary hashes
	my $sumstr = "$k\t$specstr";
	#duplication
	if($dupl>0){
	    if(exists $dupevents{$sumstr})
	    {
		my $ntmp = $dupevents{$sumstr};
		$dupevents{$sumstr} = $ntmp + $dupl;
	    }
	    else{$dupevents{$sumstr} = $dupl;}
	}
	if($match>0){
	    #matches
	    if(exists $matevents{$sumstr})
	    {
		my $ntmp = $matevents{$sumstr};
		$matevents{$sumstr} = $ntmp + $match;
	    }
	    else{$matevents{$sumstr} = $match;}

	}
	if($ins>0){
	    #insertion
	    if(exists $insevents{$sumstr})
	    {
		my $ntmp = $insevents{$sumstr};
		$insevents{$sumstr} = $ntmp + $ins;
	    }
	    else{$insevents{$sumstr} = $ins;}
	}
	if($del>0){
	    #deletion
	    if(exists $delevents{$sumstr})
	    {
		my $ntmp = $delevents{$sumstr};
		$delevents{$sumstr} = $ntmp + $del;
	    }
	    else{$delevents{$sumstr} = $del;}
	}
	if($mismatch>0){
	    #mismatch
	    if(exists $misevents{$sumstr})
	    {
		my $ntmp = $misevents{$sumstr};
		$misevents{$sumstr} = $ntmp + $mismatch;
	    }
	    else{$misevents{$sumstr} = $mismatch;}
	}
    }
}

my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;

print $outs "File created on $now_string\n";
print $outs "The following output shows the number of events that occur in the specified combination of species relative to the species written first.\n";
print $outs "Thus, summary for each event is a tab separated table with: reference_species species_in_cluster number_of_event\n";
print $outs "Events are Duplication, Matches, Insertion, Deletion, Mismatches.\n";
print $outs "The corresponding alignment files can be found in $outpath \n";
print $outs "\n\n";


##sort the entries in each entry for the hashes alphabetically

print $outs "EVENT: Duplication\n";
foreach my $du (sort keys %dupevents) {
    my @devs = split ',', $dupevents{$du};
    @devs = sort @devs;
    my $dustr = join(',',@devs);
    print $outs "$du\t$dustr\n";
}
print $outs "\n\n";

print $outs "EVENT: Matches\n";
foreach my $ma (sort keys %matevents) {
    my @mevs = split ',', $matevents{$ma};
    @mevs = sort @mevs;
    my $mastr = join(',',@mevs);

    print $outs "$ma\t$mastr\n";
}
print $outs "\n\n";

print $outs "EVENT: Insertion\n";
foreach my $in (sort keys %insevents) {
    my @ievs = split ',', $insevents{$in};
    @ievs = sort @ievs;
    my $instr = join(',',@ievs);
    print $outs "$in\t$instr\n";
}
print $outs "\n\n";

print $outs "EVENT: Deletion\n";
foreach my $de (sort keys %delevents) {
    my @deevs = split ',', $delevents{$de};
    @deevs = sort @deevs;
    my $destr = join(',',@deevs);
    print $outs "$de\t$destr\n";
}
print $outs "\n\n";

print $outs "EVENT: Mismatches\n";
foreach my $mi (sort keys %misevents) {
    my @mevs = split ',', $misevents{$mi};
    @mevs = sort @mevs;
    my $mistr = join(',',@mevs);
    print $outs "$mi\t$mistr\n";
}


print $outs "\n\n";

	

print $outm "===============Duplication alignments\===============\n";
print $outm "Duplication alignments and genetic events information: \n";
print $outm "The results for the analysis of genetic events are written to
 $sumary .
The summary file contains information about how many genetic events were 
counted in a certain combination of species. This can be used to draw a 
phylogenetic tree with genetic events at its nodes. \n";
print $outm "The files containing the duplications alignments (.aln) for 
each cluster and pairs of species can be found here: 
$outpath \n";
print $outm "Format of alignment files (.aln): At first the cluster and 
its species are defined. Then, connected nodes get mapped to the same 
one-letter-code (needed for the alignment). Afterwards for each species, its 
genetic elements are sorted by coordinate and depicted by the letter code. The
 letter code is aligned to the ones of each other species in the cluster. '~' 
stand for duplications, '-' for insertions or deletions in the alignment. \n";
print $outm "\n";



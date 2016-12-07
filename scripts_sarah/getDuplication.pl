#!/usr/bin/perl -w

## perl getDuplication edlilist inpath outpath pathtonw summaryfile

##program will produce sequences of letters, where the same letter means similar sequences, depending on the threshold.
##thus one letter for each connected component.
##then sort the elements by coordinate for each species
##get the sequence of one letters and put it in altnw for each pair of species.
##give own path for summary file because you dont want to have it in the same folder as all the alignments (many files).
## secsim and strucsim give the similarity thresholds for the sequence and structure, if one of them (and only one!) is -1, do not pay attention to this one, take only the other one into account


## Only use graphs with edges as the alignments are only done if the thresholds fit!!!
##BUT: this doesn't include single elements in the clusters which count as insertion/deletion!

#thus take the weighted edge lists as input for graphs and use the thresholds

use Data::Dumper;
use strict;
use warnings;
use List::Util qw(any none first);
#use Time::localtime;
use POSIX qw(strftime);

my $file = shift;
my $inpath = shift;
my $outpath=shift;
#my $secsim = shift;
#my $strucsim = shift;
my $pathtonw = shift;
my $sumary = shift;

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

#curfile are the .gr files, thus only showing the possible edges
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
   
    my @edges = ();
    my @nodes= ();
    my @uniqedges =();
    my @species=();

    open CF,"<$curfile" or die "can't open $curfile\n";
    while(<CF>){
	chomp;
	my $line = $_;
	my @F = split ' ', $line;
	my $n1 = $F[0];
	my $n2 = $F[1];
#	print "nodes: $n1,$n2\n";
	my @F1 = split '_', $n1;
#	my @F2 = split '_', $n2;
	my $spec1 = $F1[(scalar @F1)-5];
#	my $spec2 = $F2[(scalar @F2)-5];
	if(none {$_ eq $spec1} @species){
#	    print "pushspec\n";
	    push @species, $spec1;
	}
#all possible edges exist, thus each species will be spec1 at some point
#	if(none {$_ eq $spec2} @species){
#	    push @species, $spec2;
#	}


	my $ed = "$n1 $n2";
	push @edges, $ed;
	
	my $ed0;
	if($n1 lt $n2){$ed0 = "$n1 $n2";}
	else{$ed0 = "$n2 $n1";}
	if(none {$_ eq $ed0} @uniqedges){
	    push @uniqedges, $ed0;
	}
	
	if(none {$_ eq $n1} @nodes){
#	    print "pushnode\n";
	    push @nodes,$n1;
	}
	
    }

#    print "curfile:\n";
#    print "$curfile\n";
    
#    print "len nodes:\n";
#    print scalar @nodes;
#    print "\n";
    
    my $specstr = join(',',@species);

#    print "len species:\n";
#    print scalar @species;
#    print "\n";
	
    my @ccs = GetCCs(join(' ',@nodes),@uniqedges);


#    print "len nodes2:\n";
#    print scalar @nodes;
#    print "\n";

#    print "len species2:\n";
#    print scalar @species;
#    print "\n";

#    print "ccs:\n";
#    print join(" ",@ccs);
#    print "\n";
    
#    print "len ccs:\n";
#    print scalar @ccs;
#    print "\n";

    
    my @letters = ('A','B','C','D','E','F','G','H'.'I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z');
    
    my %n2l=(); #hash to translate each node to a letter

    #assign letter to each connected component
    for(my $c=0;$c< scalar @ccs;$c++){
	my @N = split " ", $ccs[$c];
	for(my $n=0;$n< scalar @N;$n++){
	    $n2l{$N[$n]}=$letters[$c];
	}
	print $outgr "$letters[$c]\t$ccs[$c]\n";
    }

#    print "n2l hash\n";
#    print Dumper(\%n2l);
#    print "\n";
    
    my %sp2nod = ();
    for(my $p=0;$p<scalar @species;$p++){

	$sp2nod{$species[$p]}="";
    }

    #write string of nodes for each species
    for(my $q=0;$q < scalar @nodes;$q++){
	my @T= split "_", $nodes[$q];
	my $spec = $T[(scalar @T)-5];
	my $spstr = $sp2nod{$spec};
	$spstr = "$spstr $nodes[$q]";
	$sp2nod{$spec}=$spstr;
    }

    ##get letter seq for each species
    my %sp2seq=();
    foreach my $key (keys %sp2nod){
	my $sstr = $sp2nod{$key};
	my @A= split " ", $sstr;
	my $letseq = "";
	my %tmp = ();
	for(my $a=0;$a< scalar @A;$a++){
	    my @B = split '_', $A[$a];
	    my $start = $B[(scalar @B)-3];
	    $tmp{$start} = $A[$a];
	}

	
#	print "tmp hash\n";
#	print Dumper(\%tmp);
#	print "\n";

	
	foreach my $key2 (sort { $a <=> $b } keys(%tmp) ) {
	    my $let = $tmp{$key2};
	    my $bet = $n2l{$let};
	    $letseq="$letseq$bet";
	}
	$sp2seq{$key}=$letseq;
    }

    ##get alignments for each species against each other
    foreach my $k (keys %sp2seq){
	print $outgr "\@$k\t$sp2seq{$k}\n";
	my $match=0; #m
	my $mismatch=0; #s
	my $dupl=0; #d
	my $ins=0; #i
	my $del=0; #l
	foreach my $k2 (keys %sp2seq){
	    if($k eq $k2) {next;}
	    my $lseq1 = $sp2seq{$k};
	    my $lseq2 = $sp2seq{$k2};
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


#write summary file

#my $dt = DateTime->today;
#my $time = $dt->date;
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

    
sub GetCCs{ #arguments is nodestring, uniqedgelist, input is a connected graph
    
    #first: get complement

    my @edges = @_;
    my @nodes = split ' ', $edges[0];
    $edges[0]="";
    
    my @uniqedges=(); #unique edges from complement
 
#    print "getComplement \n";

    for(my $a=0;$a<scalar @nodes;$a++){
	for(my $b=0;$b<scalar @nodes;$b++){
	    if($a == $b){next;}
	    my $tmpstr;
	    if($nodes[$a] lt $nodes[$b]){
		$tmpstr = "$nodes[$a] $nodes[$b]";	       
	    }
	    else{
		$tmpstr = "$nodes[$b] $nodes[$a]";
	    }
	    if(none {$tmpstr eq $_ } @edges){
		if(none {$tmpstr eq $_ } @uniqedges){
		    push @uniqedges, $tmpstr;
		}
	    }
	}

    }
#    my $ednum = scalar @uniqedges;
#    print "complement edge num: $ednum \n";
#    my $edtmp = join(",",@uniqedges);
#    print "complement edges: $edtmp \n";

    
    #check if we can find more than one cc
    
    my @permnodes1 = @nodes;
    my @permnodes2=();

#    my $nodstr = join(',',@permnodes1);
#    print "nodes: $nodstr\n";
    
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
	    if($k != $id1 && $k != $id2)
	    {
		push @permnodes2, $permnodes1[$k];
	    }
	}
	if($id1 != $id2){
	    my $str = "$permnodes1[$id1] $permnodes1[$id2]";
	    push @permnodes2, $str;
	}
	else{
	    push @permnodes2, $permnodes1[$id1];
	}
	@permnodes1=();
	@permnodes1 = @permnodes2;
	@permnodes2=();

#	my $permnodestr = join(',',@permnodes1);
#	print "curr conn comp: $permnodestr\n";
    }

    return @permnodes1;

}

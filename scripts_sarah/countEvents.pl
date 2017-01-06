#!/usr/bin/perl -w

##coutEvents.pl newickTree(with correct identifiers) matches dupl insertions pseudo treeout summary totelelemstr nonestr iTOLfolder
##get temporary files containing the genetic events and the corresponding numbers
##species,species... number
##output: new geneticEvents file including deletions, new newick tree including numbers at nodes


use Data::Dumper;
use strict;
use warnings;
use List::Util qw(any none first);
#use Getopt::Std;
use POSIX qw(strftime);
#use File::Basename;
#use File::Find;



my $treefile = shift;
my $matches = shift;
my $dupl = shift;
my $ins = shift;
my $pseudo = shift;
my $treeout = shift;
my $summary = shift;
my $totelemstr = shift; ##for later, printing iTOL files
my $nonestr = shift; ##for later, printing iTOL files
my $iTOLout = shift; ##for later, printing iTOL files

open(my $outt,">>",$treeout);
open(my $outs,">>",$summary);

my %plusnodes = (); ##nodes of the tree with numbers (+)
my %minusnodes = (); ##nodes of the tree with numbers (-)

my %duplications = ();
my %insertions = ();
my %pseudos = ();

open TF,"<$treefile" or die "can't open $treefile\n";

my $tree="";

while(<TF>){
    chomp;
    $tree = $_;
    last;
}
    
if($tree eq ""){print "tree format doesn't fit!\n"; exit 1;}

##split the tree into an array of its elements
my @T = (); ##tree with its components
my @N = (); ##at each position there is a number showing opening brackets - closing brackets before this position, except for ( ) , ; then -2
my @L = (); #leaves
my @tr = split '', $tree;
my $brackets = 0;
my $tmp = "";
for(my $i = 0; $i < scalar @tr; $i++)
{
    if($tr[$i] eq ')' || $tr[$i] eq '(' || $tr[$i] eq ',' || $tr[$i] eq ';'){
	if($tmp ne ""){
	    push @T, $tmp; 
	    push @N, $brackets;
	    $plusnodes{$tmp}=0;
	    $minusnodes{$tmp}=0;
	    $duplications{$tmp} = 0;
	    $insertions{$tmp} = 0;
	    $pseudos{$tmp} = 0;
#	    print "leaves? $T[(scalar @T) -2] \n";
	    if($T[(scalar @T) -2] ne ")"){ #leaves
		push @L, $tmp;
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


#print join(" ",@T);
#print "\n";
#print join(" ",@N);
#print "\n";
#print join(" ",@L);
#print "\n";



##MATCHES
##analyse match nodes first as they include several species


open FA,"<$matches" or die "can't open $matches\n";

while(<FA>){
    chomp;
    my $line = $_;
    my @F = split '\t', $line; ##should have at least 2! entries
    if(scalar @F < 2){print "misformatted input line, will be skipped! \n"; next;}
    
    my @S = split ',', $F[0]; ##species
    my $num = $F[1];
    my @nums = (); ##get bracket numbers for each species
    my @ids = (); ##indices in the tree array

    
    
 #   print "species: $F[0] \n";
    my $maxbracket = 0;
    
    ##Find lowest common ancestor (lca) for the species    
    ##sort numbers, indices, species
    for(my $j = 0;$j < scalar @S; $j++){
	my $idx;
	for(my $jj = 0; $jj < scalar @T; $jj++){
	    if($T[$jj] eq $S[$j]){
		$idx = $jj;
		last;
	    }
	}
#	print "index of $S[$j] : $idx \n";
	if($N[$idx] > $maxbracket){$maxbracket = $N[$idx];}
	push @ids, $idx;
	push @nums, $N[$idx];
    }
    
#    print join(" ",@ids);
#    print "\n";
#    print join(" ",@nums);
#    print "\n";
    
    ##find lca
    my $curlca = $ids[0]; #start with species at $S[0]
    for(my $k=1; $k < scalar @S; $k++){
#	print "curlca: $curlca, nums: $N[$curlca] \n";
	my $startidx;
	my $startnum;
	my $otheridx;
	my $othernum;
#	print "k: $k, id: $ids[$k], num: $N[$ids[$k]] ";
	if($N[$curlca] >= $N[$ids[$k]]){
	    $startidx = $curlca;
	    $startnum = $N[$curlca];
	    $otheridx = $ids[$k];
	    $othernum = $N[$ids[$k]]
	}
	else{
	    $startidx = $ids[$k];
	    $startnum = $N[$ids[$k]];
	    $otheridx = $curlca;
	    $othernum = $N[$curlca]
	}
#	print "before in between startidx: $startidx, startnum: $startnum, otheridx: $otheridx, othernum: $othernum \n";
	if($startnum == 0 && $othernum == 0)##both are the root
	{
	    $curlca = $startidx;
	    next;
	}
	my $pnum = $startnum;
	my $pidx = $startidx;
#	print "in between pidx: $pidx, otheridx: $otheridx, pnum: $pnum, othernum: $othernum \n";
	while($pnum != $startnum-1){
	    if($pnum == 0){last;}
	    $pidx++;
	    $pnum = $N[$pidx];
#	    last;
	}
#	print "in between1a pidx: $pidx, pnum: $pnum, otheridx: $otheridx, othernum: $othernum, curlca: $curlca \n";
    
	if($pidx == $otheridx){$curlca = $pidx;}
	else{
	    my $go = 1;
	    while($go){
#		last;
		my $minidx; ##smaller in the tree but higher number
		my $minnum;
		my $maxidx;
		my $maxnum;
		if($pnum >= $othernum){
		    $minnum = $pnum;
		    $minidx = $pidx;
		    $maxnum = $othernum;
		    $maxidx = $otheridx;
		}
		else{
		    $maxnum = $pnum;
		    $maxidx = $pidx;
		    $minnum = $othernum;
		    $minidx = $otheridx;
		}
#		print "min idx num: $minidx, $minnum; max idx num: $maxidx, $maxnum \n ";
		if($minnum == 0 && $maxnum==0){
		    $curlca = $minidx;
		    $go = 0;
		    next;
		}
		$pidx = $minidx;
		$pnum = $N[$pidx];#==$minnum
#		print "in between pidx: $pidx, pnum: $pnum \n";
		while($pnum != $minnum-1 && $pidx < (scalar @N) -2){ 
		    $pidx++;
		    $pnum = $N[$pidx];
#		    print "in between2 pidx: $pidx, pnum: $pnum \n";
#		    last;
		}
		if($pidx == $maxidx){$curlca = $pidx; $go = 0;}
		
		if($pidx == (scalar @T) -2){$curlca = $pidx; $go=0;}##make sure that there is an end
		$otheridx = $maxidx;
		$othernum = $maxnum;
	    }
	}
    }


    ##curlca of all species for this event found
    ##add corresponding number to this event and do deletions for all species not involved but child of lca
    $plusnodes{$T[$curlca]} = $num;
    ##find all nodes under curlca and check if it is leaves
    my @leaves = (); 
    for(my $m = $curlca -1; $m >=0; $m--){
	if($N[$m] == $N[$curlca]){last;}
	if($N[$m] > $N[$curlca]){push @leaves, $m;}
    }

#    print "leaves: \n";
#    print join(" ", @leaves);
#    print "\n";

    my @lili = ();
    for(my $p = 0; $p < scalar @leaves; $p++){
	push @lili, $T[$leaves[$p]];
    }

    
    my %species = map{$_=>1} @S;
    my %lili = map{$_=>1} @lili;
    my %L = map{$_=>1} @L;

    ##intersection of leaves and L
    my @realL = grep( $L{$_}, @lili);

    ##intersection of species and realL
    my @intersect = grep( $species{$_}, @realL);
    my %intersect = map{$_=>1} @intersect;
    
    ##diff between realL and intersect
    ##should give you only the leaves that are real leaves and not contained in species!
    my @diff = grep(!defined $intersect{$_}, @realL);
    
    
#    print "diff: \n";
#    print join (" ", @diff);
#     print "\n";
    ##todo: summarize the deleted leaves such that the number appears at common ancestors if all children are deleted
    ##easiest way: take all leaves that have a higher brackets count, find their lca and write number there. all other leaves, use as singleton deletions.
    my @leaf = (); ##take leaves that have a higher bracket count
    my @leafidx = ();
    my @singledels = ();


    for(my $d = 0;$d < scalar @diff; $d++){
	my $idx2;
	for(my $dd = 0; $dd < scalar @T; $dd++){
	    if($T[$dd] eq $diff[$d]){
		$idx2 = $dd;
		last;
	    }
	}
	if($N[$idx2] > $maxbracket)
	{push @leaf, $diff[$d];
	 push @leafidx, $idx2;}
	else{push @singledels, $diff[$d];
	     $minusnodes{$T[$idx2]}= $num;
	}
    }
    
#    print "leaf: \n";
#    print join (" ", @leaf);
#    print "\n";

    if(scalar @leaf > 0){
    ##find lca of leaves in leaf
        ##find lca
    my $llca = $leafidx[0]; #start with species at $S[0]
    for(my $k2=1; $k2 < scalar @leaf; $k2++){
	my $startidx2;
	my $startnum2;
	my $otheridx2;
	my $othernum2;
	if($N[$llca] >= $N[$k2]){
	    $startidx2 = $llca;
	    $startnum2 = $N[$llca];
	    $otheridx2 = $leafidx[$k2];
	    $othernum2 = $N[$leafidx[$k2]]
	}
	else{	    
	    $startidx2 = $leafidx[$k2];
	    $startnum2 = $N[$leafidx[$k2]];
	    $otheridx2 = $llca;
	    $othernum2 = $N[$llca]
	}
	if($startnum2 == 0 && $othernum2 == 0)##both are the root
	{
	    $llca = $startidx2;
	    next;
	}
	my $pnum2 = $startnum2;
	my $pidx2 = $startidx2;
#	print "in between pidx2: $pidx2, otheridx2: $otheridx2, pnum2: $pnum2 \n";
	while($pnum2 != $startnum2-1){
	    $pidx2++;
	    $pnum2 = $N[$pidx2];
	}
#	print "in between pidx2: $pidx2, otheridx2: $otheridx2, llca: $llca \n";
	if($pidx2 == $otheridx2){$llca = $pidx2;}
	else{
	    my $go3 = 1;
	    while($go3){
#		last;
		my $minidx2; ##smaller in the tree but higher number
		my $minnum2;
		my $maxidx2;
		my $maxnum2;
		if($pnum2 >= $othernum2){
		    $minnum2 = $pnum2;
		    $minidx2 = $pidx2;
		    $maxnum2 = $othernum2;
		    $maxidx2 = $otheridx2;
		}
		else{
		    $maxnum2 = $pnum2;
		    $maxidx2 = $pidx2;
		    $minnum2 = $othernum2;
		    $minidx2 = $otheridx2;
		}
#		print "min2 idx2 num2: $minidx2, $minnum2; max idx2 num2: $maxidx2, $maxnum2 \n ";
		if($minnum2 == 0 && $maxnum2==0){
		    $llca = $minidx2;
		    $go3 = 0;
		    next;
		}
		$pidx2 = $minidx2;
		$pnum2 = $N[$pidx2];#==$minnum
#		print "in between pidx2: $pidx2, pnum2: $pnum2 \n";
		while($pnum2 != $minnum2-1 && $pidx2 < scalar @N){ 
#		    last;
		    $pidx2++;
		    $pnum2 = $N[$pidx2];
#		    print "in between2 pidx2: $pidx2, pnum2: $pnum2 \n";
		}
		if($pidx2 == $maxidx2){$llca = $pidx2; $go3 = 0;}
		$otheridx2 = $maxidx2;
		$othernum2 = $maxnum2;
	    }
	}
    }
    ###################
    ##llca for all is now defined, = $T[$llca]
    $minusnodes{$T[$llca]} = $num;
#    print "LLCA: $T[$llca] \n";
    }
    
}


#print Dumper(\%plusnodes);
#print "\n";
#print Dumper(\%minusnodes);
#print "\n";

#other events:
#only one species per entry, thus easy
#might be changed later




##DUPLICATIONS
open FD,"<$dupl" or die "can't open $dupl\n";

while(<FD>){
    chomp;
    my $dline = $_;
    my @D = split '\t', $dline; ##should have at least 2! entries
    if(scalar @D < 2){print "misformatted input line in duplications, will be skipped! \n"; next;}
    $duplications{$D[0]} += $D[1];
}

    
##INSERTIONS
open FI,"<$ins" or die "can't open $ins\n";

while(<FI>){
    chomp;
    my $iline = $_;
    my @I = split '\t', $iline; ##should have at least 2! entries
    if(scalar @I < 2){print "misformatted input line in insertions, will be skipped! \n"; next;}
    $insertions{$I[0]} += $I[1];
}


##PSEUDOGENES
open FP,"<$pseudo" or die "can't open $pseudo\n";

while(<FP>){
    chomp;
    my $pline = $_;
    my @P = split '\t', $pline; ##should have at least 2! entries
    if(scalar @P < 2){print "misformatted input line in pseudogenes, will be skipped! \n"; next;}
    $pseudos{$P[0]} += $P[1];
}



##enter all the numbers in the tree
##each node has a tuple [a+i,-b,dx,py]
##with a = matches, i = insertions, b = deletions, 
##x = duplications whereas d stays as letter for indication, 
##y as pseudogenes whereas p stays as letter for indication

my @newT = ();
for(my $t=0; $t < scalar @N; $t++){
    my $vert = $T[$t];
    if($N[$t] >= 0){ ##not a -2, thus a node in the tree
	my $allplus = $plusnodes{$vert} + $insertions{$vert};
	my $newname = "$T[$t]\[$allplus\|-$minusnodes{$vert}\|d$duplications{$vert}\|p$pseudos{$vert}\]";
	push @newT, $newname;
    }
    else{
	push @newT, $T[$t];
    }
}

my $tstr = join ("", @newT);
#print "final tree:\n";
#print "$tstr \n";

print $outt "Newick tree with numbers specifying event counts at its nodes.
Each node name contains the following numbers in event specification:
[a|-b|dx|py] where:
a are insertions,
b are deletions, 
x are duplications (indicated by d) and 
y pseudogenes (indicated by p).\n

Tree:\n";
print $outt "$tstr\n";
close $outt;


#write summary file

#my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;

#print $outs "File created on $now_string\n";
print $outs "The following output shows the number of events that occur at the specified node in the tree.\n";
print $outs "Thus, summary for each event is a tab separated table with: tree_node number_of_events\n";
print $outs "Events are Insertion, Deletion, Duplication, Pseudogenization if chosen at the beginning.\n";
print $outs "\n\n";

##sort the entries in each entry for the hashes alphabetically

print $outs "EVENT: Insertion\n";
foreach my $in (sort keys %insertions) {
    my $sum = $insertions{$in} + $plusnodes{$in};
    print $outs "$in\t$sum\n";
}
print $outs "\n\n";


print $outs "EVENT: Deletion\n";
foreach my $de (sort keys %minusnodes) {
    print $outs "$de\t$minusnodes{$de}\n";
}
print $outs "\n\n";


print $outs "EVENT: Duplication\n";
foreach my $du (sort keys %duplications) {
    print $outs "$du\t$duplications{$du}\n";
}
print $outs "\n\n";


print $outs "EVENT: Pseudogenization\n";
foreach my $mi (sort keys %pseudos) {
    print $outs "$mi\t$pseudos{$mi}\n";
}
print $outs "\n\n";


####print output files for iTOL input
my $treeout2 = "$iTOLout\/tree.txt";
my $namesout = "$iTOLout\/tree_labels.txt";
my $totout = "$iTOLout\/tree_labels_tot.txt";
my $insout = "$iTOLout\/tree_labels_ins.txt";
my $delout = "$iTOLout\/tree_labels_del.txt";
my $dupout = "$iTOLout\/tree_labels_dup.txt";
my $pseout=  "$iTOLout\/tree_labels_pse.txt";
my $nonout=  "$iTOLout\/tree_labels_non.txt";

open(my $outt2,">>",$treeout2);
open(my $outn,">>",$namesout);
open(my $outx,">>",$totout);
open(my $outi,">>",$insout);
open(my $outd,">>",$delout);
open(my $outu,">>",$dupout);
open(my $outp,">>",$pseout);
open(my $outo,">>",$nonout);


print $outt2 $tree;

my $header_names = 
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL, Internal labels
COLOR,#000000\n";
    
my $header_tot = 
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL, Sum of elements & Legend
COLOR,#000000\n";

my $header_ins =
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL,Insertions
COLOR,#00ff00\n";

my $header_del =
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL,Deletions
COLOR,#ff0000\n";

my $header_dup =
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL,Duplications
COLOR,#0000ff\n";

my $header_pse =
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL,Pseudogenizations
COLOR,#ff00ff\n";

my $header_non =
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL,Excluded
COLOR,#00ffff\n";

my $inbetweentext =
"#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#
 
#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#
    
#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap. Used only for text labels which are displayed on the outside
MARGIN,1.0
    
#applies to external text labels only; if set, text labels associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
SHOW_INTERNAL,0
   
#Rotate all labels by the specified angle
LABEL_ROTATION,0
    
#applies to external text labels only; If set to 1, labels will be displayed in arcs aligned to the tree (in circular mode) or vertically (in normal mode). All rotation parameters (global or individual) will be ignored.
ALIGN_TO_TREE,0
    
#font size factor; For external text labels, default font size will be slightly less than the available space between leaves, but you can set a multiplication factor here to increase/decrease it (values from 0 to 1 will decrease it, values above 1 will increase it)
SIZE_FACTOR,1
    
#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the \"DATA\" keyword              #
#=================================================================#
#the following fields are possible for each node:
#ID,label,position,color,style,size_factor,rotation
    
#position defines the position of the pie chart on the tree:
#  -1 = external label
#  a number between 0 and 1 = internal label positioned at the specified value along the node branch (for example, position 0 is exactly at the start of node branch, position 0.5 is in the middle, and position 1 is at the end)
#style can be 'normal',''bold','italic' or 'bold-italic'
#size factor will be multiplied with the standard font size\n";

my $legend =
"
LEGEND_TITLE,Element Counts
LEGEND_SHAPES,1,1,1,1,1,1
LEGEND_COLORS,#000000,#00ff00,#ff0000,#0000ff,#ff00ff,#00ffff
LEGEND_LABELS,Total,Insertions,Deletions,Duplications,Pseudogenizations,Excluded
";



print $outn "$header_names$inbetweentext\n";
print $outn "\n\nDATA\n";

print $outx "$header_tot$inbetweentext\n";
print $outx "$legend\n";
print $outx "\n\nDATA\n";


print $outi "$header_ins$inbetweentext";
print $outi "\n\nDATA\n";

print $outd "$header_del$inbetweentext";
print $outd "\n\nDATA\n";

print $outu "$header_dup$inbetweentext";
print $outu "\n\nDATA\n";

print $outp "$header_pse$inbetweentext";
print $outp "\n\nDATA\n";

print $outo "$header_non$inbetweentext";
print $outo "\n\nDATA\n";

my %totnumbers = ();
my @GT = split '=', $totelemstr;
for(my $gt=0;$gt < scalar @GT; $gt++){
    my @HT = split '-', $GT[$gt];
    $totnumbers{$HT[0]} = $HT[1];
}

my %nonenums = ();
my @GG = split '=', $nonestr;
for(my $gg=0;$gg<scalar @GG;$gg++){
    my @HG = split '-', $GG[$gg];
    $nonenums{$HG[0]} = $HG[1];
}



for(my $tt=0;$tt < scalar @T;$tt++){
    if($N[$tt] > 0){##we do not include the artificial root 
	my $n1 = $T[$tt];
	my $isleaf = 0;
	for(my $ll=0;$ll < scalar @L;$ll++){
	    if($L[$ll] eq $T[$tt]){
		$isleaf = 1;
		last;
	    }
	}
	if($isleaf){
	    ##leaf node, different line printing
	    my $nline = "$totnumbers{$T[$tt]}";
	    print $outx "$T[$tt],$nline,-1,#000000,normal,1,0\n";
	    my $sumi = $insertions{$T[$tt]} + $plusnodes{$T[$tt]};
	    print $outi "$T[$tt],$sumi,-1,#00ff00,normal,1,0\n";
	    print $outd "$T[$tt],$minusnodes{$T[$tt]},-1,#ff0000,normal,1,0\n";
	    print $outu "$T[$tt],$duplications{$T[$tt]},-1,#0000ff,normal,1,0\n";
	    print $outp "$T[$tt],$pseudos{$T[$tt]},-1,#ff00ff,normal,1,0\n";
	    my $nonum=0; 
	    if(exists $nonenums{$T[$tt]}){$nonum = $nonenums{$T[$tt]};}
	    print $outo "$T[$tt],$nonum,-1,#00ffff,normal,1,0\n";		
	}
	else{
	    ##inner nodes, different line printing
	    print $outn "$T[$tt],$T[$tt],0,#000000,normal,1,0\n";
	    my $sumi2 = $insertions{$T[$tt]} + $plusnodes{$T[$tt]};
	    print $outi "$T[$tt],$sumi2,0.6,#00ff00,normal,1,0\n";
	    print $outd "$T[$tt],$minusnodes{$T[$tt]},0.7,#ff0000,normal,1,0\n";
	    print $outu "$T[$tt],$duplications{$T[$tt]},0.8,#0000ff,normal,1,0\n";
	    print $outp "$T[$tt],$pseudos{$T[$tt]},0.9,#ff00ff,normal,1,0\n";		
	}
    }
}


##Find father of a node n in a newick tree:
##give each node in the newick string a positive number
##the number of opening gaps - closing gaps before the node
##the father of the current nodes is the next node in the newick
##string with a number one smaller than the one of the current node

##find lca: if both have the same number:
##simultaneously find the father and compare them
##if they have different numbers: start with the one with
##the higher number and check if the detected father is the node
##with the lower number. if the numbers become the same and the father is
##not found, go to the first case

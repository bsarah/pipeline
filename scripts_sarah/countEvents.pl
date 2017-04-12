#!/usr/bin/perl -w

##coutEvents.pl newickTree(with correct identifiers) singletoncount matches dupl insertions pseudo treeout summary totelelemstr nonestr totpseudostr iTOLfolder
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
my $singletoncount = shift;
my $matches = shift;
my $dupl = shift;
my $ins = shift;
my $pseudo = shift;
my $treeout = shift;
my $summary = shift;
my $totelemstr = shift; ##for later, printing iTOL files
my $nonestr = shift; ##for later, printing iTOL files
my $totpseudostr = shift; ##for later, printing iTOL files, if = then nothing, this is the string for singleton pseudogenes, the others are listed in the input file pseudo.txt
my $iTOLout = shift; ##for later, printing iTOL files

open(my $outt,">>",$treeout);
open(my $outs,">>",$summary);

my %plusnodes = (); ##nodes of the tree with numbers (+)
my %minusnodes = (); ##nodes of the tree with numbers (-)

my %duplications = ();
my %insertions = ();
my %pseudos = ();
my %singletons = ();

open TF,"<$treefile" or die "can't open $treefile\n";

my $tree="";

while(<TF>){
    chomp;
    $tree = $_;
    last;
}
    
if($tree eq ""){print "tree format doesn't fit!\n"; exit 1;}

print $outs "tree: $tree \n";

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


#define root ids
my $rootid=-1;
for(my $nn = (scalar @N) -1;$nn>=0;$nn--){
    if($N[$nn]==0){$rootid =$nn;last;}
}

if($rootid == -1){
    print "no root in the tree!\n"; exit;
}

print join(" ",@T);
print "\n";
print join(" ",@N);
print "\n";
print join(" ",@L);
print "\n";


##pseudogenes

open PS,"<$pseudo" or die "can't open $pseudo\n";

while(<PS>){
    chomp;
    my $psline = $_;
    my @PS = split '\t', $psline;
    if(scalar @PS < 2){print "misformatted input line, will be skipped! \n"; next;}
    my @SS = split ',', $PS[0];
    my $psnum = $PS[1];
    if(scalar @SS == 1){
	$pseudos{$SS[0]} += $psnum;
    }
    else{
	my $Tstr = join('=', @T);
	my $Sstr = join('=',@SS);
#	print $outs "findParentOfAll($Tstr,$Sstr)\n";
	my @pslcas = findParentOfAll($Tstr,$Sstr);
	for(my $i=0;$i<scalar @pslcas;$i++){
	    $pseudos{$pslcas[$i]} += $psnum;
	}
    }
}




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
    my $pstr1="";
    my $pstr2="";

    ##try to find curlca by always looking at two species at the same time
    ##then, get set of nodes from the pathes from the current leaves to the root and intersect them
    ##if there are several nodes in the intersection, take the one with the highest N
    my $curlca = $ids[0]; #start with species at $S[0]
    for(my $k=1; $k < scalar @S; $k++){
	my @path1 = ();
	my @path2 = ();
	push @path1, $curlca;
	push @path2, $ids[$k];
	my $num1 = $N[$curlca];
	for(my $t4=$curlca;$t4<scalar @T;$t4++){
	    if($N[$t4] != $num1){next;}
	    $num1--;
	    if($num1 < 0){push @path1, $rootid;last;}
	    print "t4: $t4 \n";
	    push @path1, $t4;
	}

	my $num2 = $N[$ids[$k]];
	for(my $t5=$ids[$k];$t5<scalar @T;$t5++){
	    if($N[$t5] != $num2){next;}
	    $num2--;
	    if($num2 < 0){push @path2, $rootid;last;}
	    push @path2, $t5;
	}

	if(scalar @path1 == 0){
	    print $outs "path1 zero: $F[0]\n";
	}
	else{$pstr1 = join("=",@path1);}

	if(scalar @path2 == 0){
	    print $outs "path2 zero: $F[0]\n";
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
		    last; #this can be deleted later on in order to get all overlapping nodes
		}
	    }

	}

	if(scalar @pdiff == 0){
	    print $outs "son mist pathes: $pstr1; $pstr2\n";
	}
	elsif(scalar @pdiff == 1){
	    $curlca = $pdiff[0];
	}
	else{
	    #choose the curlca with the highest N
	    my $curplca = $N[$pdiff[0]];
	    my $curpid = 0;
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
	print $outs "strange curlca $T[$curlca], species: $F[0], pathes: $pstr1; $pstr2\n";
    }

    
=for comment    
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
	else{#if N[ids[k]] < N[curlca] then idsk is above curlca
	    $curlca = $ids[$k];
	    $startidx = $ids[$k];
	    $startnum = $N[$ids[$k]];
	    $otheridx = $curlca;
	    $othernum = $N[$curlca]
	}
#	print "before in between startidx: $startidx, startnum: $startnum, otheridx: $otheridx, othernum: $othernum \n";
	if($startnum == 0 && $othernum == 0)##both are the root
	{
	    $curlca = $startidx;
	    last;
	}	
    

        #try to find the common parent of both current nodes in order to get the curlca
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
=cut

    ##curlca of all species for this event found
    ##add corresponding number to this event and do deletions for all species not involved but child of lca
    $plusnodes{$T[$curlca]} += $num;
    ##find all nodes under curlca and check if it is leaves
    my @leaves = (); 
    for(my $m = $curlca -1; $m >=0; $m--){
	if($N[$m] == $N[$curlca]){last;}
	if($N[$m] > $N[$curlca]){push @leaves, $m;}
    }

#    print "leaves: \n";
#    print join(" ", @leaves);
#    print "\n";

    my @lili = (); #names of leaves
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
    ##thus, those leaves need to be deleted
    my @diff = grep(!defined $intersect{$_}, @realL);
    if(scalar @diff == 0 && scalar @intersect == 0 && scalar @realL == 0){
	print $outs "no elements with $T[$curlca]; pathes: $pstr1; $pstr2\n";
    }

    if(scalar @diff >= scalar @L){
	print $outs "all leaves deleted?\n";
    }
    
    if(scalar @diff == 1){
	$minusnodes{$diff[0]} += $num;
    }
    elsif(scalar @diff == 0){next;}
    else{
	my $Treestr = join('=', @T);
	my $Leafstr = join('=',@diff);
	my @leafsToDel = findParentOfAll($Treestr,$Leafstr);
	if(scalar @leafsToDel == 0){
	    print $outs "no parent of all: $Leafstr\n";
	}
	for(my $ltd=0;$ltd < scalar @leafsToDel;$ltd++){
	    $minusnodes{$leafsToDel[$ltd]} += $num;
	}
    }

=for comment
   
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
	while($pnum2 != $startnum2-1 && $pidx2++ < scalar @N){
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

=cut


#    print "LLCA: $T[$llca] \n";
}


#print Dumper(\%plusnodes);
#print "\n";
#print Dumper(\%minusnodes);
#print "\n";

#other events:
#only one species per entry, thus easy
#might be changed later


#singletoncount
my @SC = split "=", $singletoncount;
for(my $s = 0; $s < scalar @SC; $s++){
    if($SC[$s] eq ""){next;}
    my @stmp = split "-", $SC[$s];
    $singletons{$stmp[0]} = $stmp[1];
}


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


##not needed anymore, see above
##PSEUDOGENES
#open FP,"<$pseudo" or die "can't open $pseudo\n";
#
#while(<FP>){
#    chomp;
#    my $pline = $_;
#    my @P = split '\t', $pline; ##should have at least 2! entries
#    if(scalar @P < 2){print "misformatted input line in pseudogenes, will be skipped! \n"; next;}
#    $pseudos{$P[0]} += $P[1];
#}



##enter all the numbers in the tree
##each node has a tuple [a+i,-b,dx,py]
##with a = matches, i = insertions, b = deletions, 
##x = duplications whereas d stays as letter for indication, 
##y as pseudogenes whereas p stays as letter for indication



my %totnumbers = ();
my @GT = split '=', $totelemstr;
for(my $gt=0;$gt < scalar @GT; $gt++){
    if($GT[$gt] eq ""){next;}
    my @HT = split '-', $GT[$gt];
    $totnumbers{$HT[0]} = $HT[1];
}

my %nonenums = ();
my @GG = split '=', $nonestr;
for(my $gg=0;$gg<scalar @GG;$gg++){
    if($GG[$gg] eq ""){next;}
    my @HG = split '-', $GG[$gg];
    $nonenums{$HG[0]} = $HG[1];
}

my @GP = split '=', $totpseudostr;
for(my $gp=0;$gp<scalar @GP;$gp++){
    if($GP[$gp] eq ""){next;}
    my @HP = split '-', $GP[$gp];
    $pseudos{$HP[0]} += $HP[1];
}







my @newT = ();
for(my $t=0; $t < scalar @N; $t++){
    my $vert = $T[$t];
    if($N[$t] >= 0){ ##not a -2, thus a node in the tree
	##check that the values we wanna get from the hash really exist
	if(exists($plusnodes{$vert})){}else{$plusnodes{$vert}=0;}
	if(exists($insertions{$vert})){}else{$insertions{$vert}=0;}
	my $allplus = $plusnodes{$vert} + $insertions{$vert};
	if(exists($minusnodes{$vert})){}else{$minusnodes{$vert}=0;}
	if(exists($totnumbers{$vert})){}else{$totnumbers{$vert} = 0;}
	if(exists($duplications{$vert})){}else{$duplications{$vert}=0;}
	if(exists($singletons{$vert})){}else{$singletons{$vert}=0;}
	if(exists($pseudos{$vert})){}else{$pseudos{$vert}=0;}
	if(exists($nonenums{$vert})){}else{$nonenums{$vert}=0;}
	my $newname = "$T[$t]\[t$totnumbers{$vert}|i$allplus\|l$minusnodes{$vert}\|d$duplications{$vert}\|s$singletons{$vert}\|p$pseudos{$vert}|n$nonenums{$vert}\]";
	push @newT, $newname;
    }
    else{
	push @newT, $T[$t];
    }
}

my %numdiffs = ();

###############Test numbers, thus sum up all the numbers for each species and check the difference to the total number of elems
for(my $ll=0;$ll<scalar @L;$ll++){
    #go through all leaves
    #add up all the numbers from nodes that have a N-val < N[$ll] and are not leaves, each number appears once!!!
    #find ID of leave in T
    my $curid;
    for(my $t0=0;$t0<scalar @T;$t0++){
	if($L[$ll] eq $T[$t0]){
	    $curid = $t0;
	    last;
	}
    }
    my $allelems = 0;
    my $elemsum = 0;
    my $curval = $N[$curid];
    #go through Tnew and N and add up all the numbers
    for(my $t1 = $curid;$t1<scalar @newT;$t1++){
	if($N[$t1] != $curval){next;}
	$curval--;
	my $curel = $newT[$t1];
	my @split1 = split '\[', $curel;
	my @split1b = split '\]', $split1[1];
	my @split2 = split '\|', $split1b[0];
	for(my $s2 = 0;$s2 < scalar @split2; $s2++){
	    #letters are: t,i,-l,d,n,s
	    my $curnum = $split2[$s2];
	    if($curnum =~ /^t/){$allelems += substr($curnum,1);}#total
	    if($curnum =~ /^i/){$elemsum += substr($curnum,1);}#ins
	    if($curnum =~ /^l/){$elemsum -= substr($curnum,1);}#loss
	    if($curnum =~ /^d/){$elemsum += substr($curnum,1);}#dupl
	    if($curnum =~ /^s/){$elemsum += substr($curnum,1);}#single	    
	    if($curnum =~ /^p/){$elemsum += substr($curnum,1);}#pseudo
	    if($curnum =~ /^n/){$elemsum += substr($curnum,1);}#none
	}
    }
    my $diff = $allelems-$elemsum;
    $numdiffs{$L[$ll]} = $diff;
}

#add the diff numbers to the tree?

my $tstr = join ("", @newT);
#print "final tree:\n";
#print "$tstr \n";

print $outt "Newick tree with numbers specifying event counts at its nodes.
Each node name contains the following numbers in event specification:
[ta|ib|lc|dx|se|py|nz] where:
a are total number of elements,
b are insertions,
c are deletions, thus losses, 
x are duplications,
e are singletons, 
y pseudogenes and
z are the elements that were excluded from the analysis.\n

Tree:\n";
print $outt "$tstr\n";

#print diff numbers first here:
foreach my $dif (keys %numdiffs){
    print $outt "$dif $numdiffs{$dif}\n";
}


close $outt;


#write summary file

#my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;

#print $outs "File created on $now_string\n";
print $outs "The following output shows the number of events that occur at the specified node in the tree.\n";
print $outs "Thus, summary for each event is a tab separated table with: tree_node number_of_events\n";
print $outs "Events are Insertion, Deletion, Duplication, Pseudogenization if chosen at the beginning.\n";
print $outs "\n\n";

##sort the entries in each entry for the hashes alphabetically


my %allins = ();
foreach my $in1 (sort keys %insertions) {
    if(exists($allins{$in1})){$allins{$in1} += $insertions{$in1};}
    else{$allins{$in1} = $insertions{$in1};}
}
foreach my $in2 (sort keys %plusnodes) {
    if(exists($allins{$in2})){$allins{$in2} += $plusnodes{$in2};}
    else{$allins{$in2} = $plusnodes{$in2};}
}

print $outs "EVENT: Insertion\n";
foreach my $in (sort keys %allins) {
    print $outs "$in\t$allins{$in}\n";
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


print $outs "EVENT: Singletons\n";
foreach my $si (sort keys %singletons) {
    print $outs "$si\t$singletons{$si}\n";
}
print $outs "\n\n";


print $outs "EVENT: No Anchors\n";
foreach my $ni (sort keys %nonenums) {
    print $outs "$ni\t$nonenums{$ni}\n";
}
print $outs "\n\n";


print $outs "EVENT: Others\n";
foreach my $oi (sort keys %numdiffs) {
    print $outs "$oi\t$numdiffs{$oi}\n";
}
print $outs "\n\n";


##TODO add extra files for singletons and numdiffs?

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



for(my $tt=0;$tt < scalar @T;$tt++){
    if($N[$tt] >= 0){##we do include the artificial root to debug
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
	    print $outd "$T[$tt],$minusnodes{$T[$tt]},0.8,#ff0000,normal,1,0\n";
	    #duplications do not occur at inner nodes
	    #print $outu "$T[$tt],$duplications{$T[$tt]},0.8,#0000ff,normal,1,0\n";
	    print $outp "$T[$tt],$pseudos{$T[$tt]},0.95,#ff00ff,normal,1,0\n";		
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


##given a list of leaves and a tree, find the lca of the leaves s.t. there is NO child of the parent that is NOT in the given leaf list
##This function will be used to determine the deletions and the pseudogenes
##the return value is again a list of nodes in the tree
sub findParentOfAll{
    my @inp= @_;
    #tree and leaf string is separated by =
    my @T = split '=', $inp[0];
    my @L = split '=', $inp[1];
    my @Ltmp = split '=', $inp[1];
    
    my @output = ();
    
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
    
    #print join(" ",@T);
    #print "\n";
    #print join(" ",@N);
    #print "\n";
    #print join(" ",@allleaves);
    #print "\n";
    
    
    while(scalar @Ltmp > 0){
	
	for(my $j = 0; $j<scalar @Ltmp; $j++){
	    my $id = 0; ##index of leaf in T and N
	    for(my $l=0;$l < scalar @T;$l++){
		if($T[$l] eq $Ltmp[$j]){
		    $id = $l;
		    last;
		}
	    }
	    
	    my $nullLtmp = $Ltmp[0];
	    ##find brothers in tree
	    my @bros = (); ##ids of brothers in T, N
	    my $curnum = $N[$id];
	    push @bros, $id;
	    my $down=$id-1;
	    while($down >=0){
		if($N[$down] == -2){$down--;next;}
		if($N[$down] < $curnum){last;}#break
		elsif($N[$down] > $curnum){}#do nothing
		else{#equality, this is a bro!
		    push @bros, $down;
		}
		$down--;
	    }
	    my $up = $id+1;
	    while($up < scalar @N){
		if($N[$up]==-2){$up++;next;}
		if($N[$up] < $curnum){last;}#break
		elsif($N[$up] > $curnum){}#do nothing
		else{#equality, this is a bro!
		    push @bros, $up;
		}
		$up++;
	    }
	    
	    ##check, if all bros are in L
	    ##if yes, get their father and remove them from L
	    ##if not, check if they are leaves
	    ##if yes, break, no parent of all
	    ##if no, find the children
	    my $allL = 1;
	    for(my $b = 0; $b<scalar @bros; $b++){
		if(none {$_ eq $T[$bros[$b]]} @L){
		    $allL = 0;
		    last;
		}
	    }
	    
	    if($allL){
		##find parent and delete bros from input Ltmp, add parent to L and Ltmp
		my $parent = 0;
		my $bronum = $N[$bros[0]];
		my $up = $bros[0];
		while($up < scalar @N){
		    if($N[$up] == $bronum -1 ){
			$parent = $up;
			last;
		    }#break
		    $up++;
		}
		#delete entries in Ltmp
		for(my $b = 0; $b<scalar @bros; $b++){
		    my $index = -1;
		    for(my $ii=0;$ii < scalar @Ltmp; $ii++){
			if($Ltmp[$ii] eq $T[$bros[$b]]){$index = $ii;last;}
		    }
		    if($index>=0){
			splice(@Ltmp,$index,1);
		    }
		}
		push @Ltmp, $T[$parent];
		push @L, $T[$parent];
		
	    }
	    else{
		##check if they are leaves
		my @noleaves = ();
		for(my $b = 0; $b<scalar @bros; $b++){
		    if(none {$_ eq $T[$bros[$b]]} @allleaves){
			push @noleaves, $bros[$b];
		    }
		    #check if leaves are in L, if they are not in L, delete them from Ltmp
		    else{
			my $index1 = -1;
			my $index2 = -1;
			for(my $i1 = 0; $i1 < scalar @L;$i1++){
			    if($L[$i1] eq $T[$bros[$b]]){$index1 = $i1;last;}
			} 
			if($index1 < 0){
			    for(my $i2 = 0; $i2 < scalar @Ltmp;$i2++){
				if($Ltmp[$i2] eq $T[$bros[$b]]){$index2 = $i2;last;}
			    }
			    if($index2 >0){
				splice(@Ltmp,$index2,1);
			    }
			}
			
		    }
		}
		
		
		if(scalar @noleaves == 0){##all in L are leaves
		    print "case noleaves\n";
		    ##there is no parent of all, put the bros that were in input L into output list and delete from input list Ltmp
		    for(my $b = 0; $b<scalar @bros; $b++){
			my $index1 = -1;
			my $index2 = -1;
			for(my $i1 = 0; $i1 < scalar @L;$i1++){
			    if($L[$i1] eq $T[$bros[$b]]){$index1 = $i1;last;}
			} 
			if($index1 >= 0){
			    for(my $i2 = 0; $i2 < scalar @Ltmp;$i2++){
				if($Ltmp[$i2] eq $T[$bros[$b]]){$index2 = $i2;last;}
			    }
			    if($index2 >= 0){
				splice(@Ltmp,$index2,1);
			    }
			    if(none {$_ eq $T[$bros[$b]]} @output){
				push @output, $T[$bros[$b]];
			    }
			}
			else{
			    my $index3 = -1;
			    for(my $i3 = 0; $i3 < scalar @Ltmp;$i3++){
				if($Ltmp[$i3] eq $T[$bros[$b]]){$index3 = $i3;last;}
			    }
			    if($index3 >= 0){
				splice(@Ltmp,$index3,1);
			    }
			    
			    
			}
		    }
		}
		else{
		    ##find the children of the ones which are not leaves
		    #if the children are all leaves and are all in L, push their parent in to output and delete from Ltmp, if the children are not leaves, push them into Ltmp
		    for(my $k = 0; $k < scalar @noleaves; $k++){
			for(my $llk=0;$llk<scalar @Ltmp; $llk++){
			    if($Ltmp[$llk] eq $noleaves[$k]){
				splice(@Ltmp,$llk,1);
				last;
			    }
			}
			
			
			my $curnum2 = $N[$noleaves[$k]];
			my $down=$noleaves[$k]-1;
			my $notall = 0;
			my @leafchildrentoout = ();
			while($down >=0){
			    if($N[$down] == -2){$down--;next;}
			    if($N[$down] <= $curnum2){last;}#break
			    elsif($N[$down] = $curnum2 + 1){
				if(none {$_ eq $T[$down]} @allleaves){#check if children are leaves
				    if(none {$_ eq $T[$down]} @L){
					push @Ltmp, $T[$down];
				    }
				}
				elsif(none {$_ eq $T[$down]} @L){#check if leafchild is not in L
				    $notall = 1;
				}
				else{#leaf child is in L, thats good
				    push @leafchildrentoout, $T[$down];
				}
			    }#that's the child!
			    else{ ##bigger than +1, nothing
			    }
			    $down--;
			}
			
			
			if($notall == 0){#push current notleaf to output, delete from Ltmp
			    if(none {$_ eq $T[$noleaves[$k]]} @L){}
			    else{
				if(none {$_ eq $T[$noleaves[$k]]} @output){
				    push @output, $T[$noleaves[$k]];
				}
				my $index15 = -1;
				for(my $i15 = 0; $i15 < scalar @Ltmp; $i15++){
				    if($Ltmp[$i15] eq $T[$noleaves[$k]]){$index15 = $i15;last;}
				}
				if($index15 >= 0){
				    splice(@Ltmp,$index15,1);
				}
			    }
			}
			else{#push all children in @leafchildrentoout in output and delete from Ltmp
			    if(scalar @leafchildrentoout == 0){next;}
			    for(my $lcto = 0;$lcto < scalar @leafchildrentoout;$lcto++){
				if(none {$_ eq $leafchildrentoout[$lcto]} @output){
				    push @output, $leafchildrentoout[$lcto];
				}
				my $index16 = -1;
				for(my $i16 = 0; $i16 < scalar @Ltmp; $i16++){
				    if($Ltmp[$i16] eq $leafchildrentoout[$lcto]){$index16 = $i16;last;}
				}
				if($index16 >= 0){
				    splice(@Ltmp,$index16,1);
				}
				
				
			    }
			}
			
		    }
		    
		}
		##what happens if all runs through and curnum is still in Ltmp?
		
		if(scalar @Ltmp > 0 && $nullLtmp eq $Ltmp[0]){##put Ltmp0 in output and delete it from Ltmp
		    if(none {$_ eq $Ltmp[0]} @output){
			if(none {$_ eq $Ltmp[0]} @L){}
			else{
			    push @output, $Ltmp[0];
			}
		    }
		    splice(@Ltmp,0,1);
		}
	    }
	    
	}
    }
    
    
    
    
    
    return @output;
    
} 

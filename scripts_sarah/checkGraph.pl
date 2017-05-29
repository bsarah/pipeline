#!/usr/bin/perl -w


## perl checkGraph.pl .edlilist inpath outpath secsim strucsim mode cographlist noncographlist noEdgeGraohs_list edgegraphs summary

##outpath is for files that will be visualized with R (no weights), extra folder
## .edli is the edgelist file with node1 node2 seqval strucval (separated with space)
## secsim is the lower bound for the similarity of sequences, e.g. 90%, thus 0.9
## strucsim is the lower bound for the similarity of secondary structure, e.g. 90%, thus 0.9
## if secsim or strucsim (just one of both) is -1, do not take for this value!
## edge is only drawn if both values are higher than lower bound
##output R compatible edge list to draw graph (thus node1 node2 AND node2 node1)
##cographlist and noncographlist are tables with information to fill in for each graph, depending if it is a cograph or not. tables look like:
## clustername(filename) nodenum edgenum speclist duplications

#handle pseudogenes as if they were usual homologos genes and include them in the graphs to test

use Data::Dumper;
use strict;
use warnings;
use List::Util qw(any none first);
use File::Basename;


my $file = shift;
my $inpath = shift;
my $outpath=shift;
my $secsim = shift;
my $strucsim = shift;
my $mode = shift; ##0 is cmmode, 1 is genelist
#my $psecsim = shift;
#my $pstrucsim = shift;
my $cglist = shift;
my $ncglist = shift;
my $noedges = shift;
my $edgelist = shift;
my $summary = shift;

open(my $outcg, ">>",$cglist);
open(my $outn,">>",$ncglist);

open(my $outo, ">>",$noedges);
open(my $oute, ">>",$edgelist);

open(my $outs, ">>",$summary);

my $seqlim = $secsim;
my $struclim = $strucsim;

#if($psecsim >= 0 && $psecsim <= $secsim){
#    $seqlim = $psecsim;
#}
#else{$seqlim = $secsim;}
#if($pstrucsim >= 0 && $pstrucsim <= $strucsim){
#    $struclim = $pstrucsim;
#}
#else{$struclim = $strucsim;}

my @changedGraphs = ();

my $numcg = 0;
my $numcliques = 0;
my $numng = 0;
my $sumnodescg = 0;
my $sumnodesng = 0;
my $sumedgescg = 0;
my $sumedgesng = 0;
my $sumdenscg = 0;
my $sumdensng = 0;

my $maxncg = 0;
my $maxecg = 0;
my $minncg = 3000000;
my $minecg = 3000000;


my $maxnng = 0;
my $maxeng = 0;
my $minnng = 3000000;
my $mineng = 3000000;

open FA,"<$file" or die "can't open $file\n";
while(<FA>){
    chomp;
    my $curfile = $_;

    my $name;
    my $path;
    my $suffix;
    ($name,$path,$suffix) = fileparse($curfile);
    

    my @D = split '\.', $name;
    my $almostname = $D[0];
    my $newname="$almostname.gr";

#    print "checkgraph args: almostname: $almostname, secsim: $secsim, strucsim: $strucsim \n";
    
#    print "filename: $newname \n";

    my @edges = ();
    my @nodes= ();
    my @uniqedges =();
    my %species = ();
    my %seqscore = ();
    my %strscore= ();

    open CF,"<$curfile" or die "can't open $curfile\n";
    while(<CF>){
	chomp;
	my $line = $_;
	my @F = split ' ', $line;
	my $n1 = $F[0];
	my $n2 = $F[1];
	my $s1 = $F[2];
	my $s2 = $F[3];
	#sort out entries with same species
	my @F1 = split '_', $n1;
	my @F2 = split '_', $n2;
	my $spec1;
	my $spec2;
	if($mode == 0){
	    $spec1 = $F1[(scalar @F1)-5];
	    $spec2 = $F2[(scalar @F2)-5];
	}
	else{
	    $spec1 = $F1[(scalar @F1)-7];
	    $spec2 = $F2[(scalar @F2)-7];
	}
	if($spec1 eq $spec2){next;}
	if($seqlim == -1){
	    if($s2 >= $struclim){
		my $ed = "$n1 $n2";
		push @edges, $ed;
		
		my $ed0;
		if($n1 lt $n2){$ed0 = "$n1 $n2";}
		else{$ed0 = "$n2 $n1";}
		if(none {$_ eq $ed0} @uniqedges){
		    push @uniqedges, $ed0;
		}
		if(exists($strscore{$ed0})){}
		else{$strscore{$ed0}=$s2;}
	    }
	}
	elsif($struclim == -1){
	    if($s1 >= $seqlim){
		my $ed = "$n1 $n2";
		push @edges, $ed;
		
		my $ed0;
		if($n1 lt $n2){$ed0 = "$n1 $n2";}
		else{$ed0 = "$n2 $n1";}
		if(none {$_ eq $ed0} @uniqedges){
		    push @uniqedges, $ed0;
		}
		if(exists($seqscore{$ed0})){}
		else{$seqscore{$ed0}=$s1;}
	    }
	}
	else{
	    if($s1 >= $seqlim && $s2 >= $struclim){
		my $ed = "$n1 $n2";
		push @edges, $ed;
		
		my $ed0;
		if($n1 lt $n2){$ed0 = "$n1 $n2";}
		else{$ed0 = "$n2 $n1";}
		if(none {$_ eq $ed0} @uniqedges){
		    push @uniqedges, $ed0;
		}
		if(exists($seqscore{$ed0})){}
		else{$seqscore{$ed0}=$s1;}
		if(exists($strscore{$ed0})){}
		else{$strscore{$ed0}=$s2;}
	    }
	}
	
	if(none {$_ eq $n1} @nodes){
	    push @nodes,$n1;
	    #count how many different species
	    if(exists($species{$spec1})){$species{$spec1}++;}
	    else{$species{$spec1}=1;}
	}
	
    }
    
    ##do not write graph files at the moment, as they are not visualized
    ##write graph file to be able to visualize with R
#    if(scalar @edges > 0){
#	print "outpath: $outpath";
#	open(my $outgr,">>$outpath\/$newname");
#	for(my $e = 0;$e < scalar @edges;$e++){
#	    print $outgr "$edges[$e]\n";
#	}
#	close($outgr);
#    }
    
    my $numn = scalar @nodes;
    my $nume = scalar @uniqedges;
    my $numsp = scalar (keys %species);
    
    if($numn == 0 || $nume == 0 || $numsp <= 1){
	print $outo "$almostname \n";
	next;
    }
    else{
	print $oute "$almostname \n";
    }

    
    #check if it is a clique first
    #this only works if we only have one species
    #for several species, check the cliquenum edgenum +
    #edges for a complete graph between samespeciesedges
    my $sumextraeds = 0;
    foreach my $s (keys %species){
	my $nodenum = $species{$s};
	if($nodenum > 1){
	    my $addnodes = ($nodenum * ($nodenum-1))/2;
	    $sumextraeds += $addnodes;
	}
    }

    my $numenew = $nume + $sumextraeds;

    #calculate density with numenew
    my $density = 0;
#    if($nume > 0){
	$density = $numn/$numenew;
    #    }

    #density in a clique: 2/(n-1)


    
    my $cliquenum=(($numn-1)*$numn)/2;
    if($nume == $cliquenum || $numenew == $cliquenum){
	$numcg++;
	$numcliques++;
	$sumnodescg = $sumnodescg + $numn;
	$sumedgescg = $sumedgescg + $numenew;
	$sumdenscg = $sumdenscg + $density;
	if($numn > $maxncg){$maxncg = $numn;}
	if($numn < $minncg){$minncg = $numn;}
	if($numenew > $maxecg){$maxecg = $numenew;}
	if($numenew < $minecg){$minecg = $numenew;}
	print $outcg "$almostname\t$numn\t$nume\t$numenew\t$density\n";
	next;
    }
    
    my @noCGCC = ();
    my @CGCC = ();
    my @singlenodes = ();
    ##get CCs for this graph and then check if this is a cograph
    my $vertices = join(',',@nodes);
    my $arcs = join(',',@uniqedges);
    my @CCs = connectedComponents($vertices,$arcs);
    #get induced subgraph and check if this is a cograph
    for(my $c=0;$c<scalar @CCs;$c++){
	my @V = split ' ', $CCs[$c];
	if(scalar @V == 1){
	    push @singlenodes, $V[0];
	}
	my @edges = ();
	for(my $v = 0;$v<(scalar @V)-1;$v++){
	    for(my $w = $v+1;$w<scalar @V;$w++){
		my $n1 = $V[$v];
		my $n2 = $V[$w];
		my $ed;
		if($n1 lt $n2){$ed = "$n1 $n2";}
		else{$ed = "$n2 $n1";}
		if(exists($strscore{$ed}) || exists($seqscore{$ed})){
		    push @edges, $ed;
		}
	    }
	}
	my $cg = IsCograph(join(' ',@V),@edges);
	if($cg == 1){push @CGCC, $CCs[$c];}
	else{push @noCGCC, $CCs[$c];}

    }
    if(scalar @noCGCC == 0){
	#is a cograph
	print $outcg "$almostname\t$numn\t$nume\t$numenew\t$density\n";
	$numcg++;
	$sumnodescg = $sumnodescg + $numn;
	$sumedgescg = $sumedgescg + $numenew;
	$sumdenscg = $sumdenscg + $density;
	if($numn > $maxncg){$maxncg = $numn;}
	if($numn < $minncg){$minncg = $numn;}
	if($numenew > $maxecg){$maxecg = $numenew;}
	if($numenew < $minecg){$minecg = $numenew;}
    }
    else{
	#no cograph
	$numng++;
	$sumnodesng = $sumnodesng + $numn;
	$sumedgesng = $sumedgesng + $numenew;
	$sumdensng = $sumdensng + $density;
	if($numn > $maxnng){$maxnng = $numn;}
	if($numn < $minnng){$minnng = $numn;}
	if($numenew > $maxeng){$maxeng = $numenew;}
	if($numenew < $mineng){$mineng = $numenew;}
	my $ediff = $numn*($numn-1)/2 - $numenew;
	print $outn "$almostname\t$numn\t$nume\t$numenew\t$density\t$ediff\n";
	my @allnewedges = ();
	#edit only CCs that are noncographs
	for(my $z=0;$z<scalar @noCGCC;$z++){
	    #get induced subgraph for CCs
	    my @N = split ' ', $noCGCC[$z];
	    my @ccedges = getInducedSubgraph(join('=',@N),join('=',@uniqedges));
	    #check if this is almost a clique, if yes, do not delete any edges as they automatically will be added later
	    my %tmpspnum = ();
	    for(my $y=0;$y<scalar @N;$y++){
		my @NN = split '_', $N[$y];
		my $sp;
		if($mode ==0){$sp = $NN[-5];}
		else{$sp = $NN[-7];}
		if(exists($tmpspnum{$sp})){$tmpspnum{$sp}++;}
		else{$tmpspnum{$sp}=1;}
	    }
	    my $sumextra = 0;
	    foreach my $sc (keys %tmpspnum){
		my $nnu = $tmpspnum{$sc};
		if($nnu > 1){
		    my $addn = ($nnu * ($nnu-1))/2;
		    $sumextra += $addn;
		}
	    }
	    my $newcceds = scalar @uniqedges + $sumextra;
#	    my $cliquenum=(($numn-1)*$numn)/2;
	    #cliquenum is max possible edgenum
	    if($newcceds/$cliquenum > 0.5){
		push @allnewedges, @uniqedges;
		next;
	    }
	    
	    my $strkeys = join('=',keys %strscore);
	    my $strvals = join('=',values %strscore);
	    my $seqkeys = join('=',keys %seqscore);
	    my $seqvals = join('=',values %seqscore);
	    my (@newnodes,@newedges) = cographEditing(join('=',@N),join('=',@ccedges),$strkeys,$strvals,$seqkeys,$seqvals);
	    push @allnewedges, @newedges;
	}
	for(my $y=0;$y< scalar @CGCC;$y++){
	    my @M = split ' ', $CGCC[$y];
	    my @tmpedges = getInducedSubgraph(join('=',@M),join('=',@uniqedges));
	    push @allnewedges, @tmpedges;
	}
	#now print the edges in allnewedges in a file
	my $newcurfile = "$curfile\_old";
	my $mvcmd = "mv $curfile $newcurfile";
	#print STDERR "mvcommand: $mvcmd\n";
	readpipe("$mvcmd");
	push @changedGraphs, $curfile;
	open(my $outc, ">>",$curfile);
	#take care that nodes don't get lost, if they are not contained in an edge!
	for(my $n1 = 0; $n1 < scalar @nodes-1;$n1++){
	    for(my $n2=$n1+1;$n2 < scalar @nodes;$n2++){
		my $tmped;
		if($nodes[$n1] lt $nodes[$n2]){$tmped = "$nodes[$n1] $nodes[$n2]";}
		else{$tmped = "$nodes[$n2] $nodes[$n1]";}
		if(none {$tmped eq $_ } @allnewedges){
		    print $outc "$nodes[$n1] $nodes[$n2] 0 0\n";
		    print $outc "$nodes[$n2] $nodes[$n1] 0 0\n";
		}
		else{
		    print $outc "$nodes[$n1] $nodes[$n2] 1 1\n";
		    print $outc "$nodes[$n2] $nodes[$n1] 1 1\n";
		}
	    }
	}
    }
}

    
    #print "nodenum: $numn, edgenum: $nume\n";
    #my $edtmp = join(",",@uniqedges);
    #print "edges: $edtmp \n";
#    if($density == 0 || $density == 2){
#	print $outcg "$almostname\t$numn\t$nume\t$density\n";
#	$numcg++;
	#$numcliques++;
#	$sumnodescg = $sumnodescg + $numn;
#	$sumedgescg = $sumedgescg + $nume;
#	$sumdenscg = $sumdenscg + $density;
#	if($numn > $maxncg){$maxncg = $numn;}
#	if($numn < $minncg){$minncg = $numn;}
#	if($nume > $maxecg){$maxecg = $nume;}
#	if($nume < $minecg){$minecg = $nume;}
#    }
#    else{
=for
	my $cg = IsCograph(join(' ',@nodes),@uniqedges);
	if($cg ==1){
	    print $outcg "$almostname\t$numn\t$nume\t$numenew\t$density\n";
	    $numcg++;
	    $sumnodescg = $sumnodescg + $numn;
	    $sumedgescg = $sumedgescg + $numenew;
	    $sumdenscg = $sumdenscg + $density;
	    if($numn > $maxncg){$maxncg = $numn;}
	    if($numn < $minncg){$minncg = $numn;}
	    if($numenew > $maxecg){$maxecg = $numenew;}
	    if($numenew < $minecg){$minecg = $numenew;}
	}
	else{

	    $numng++;
	    $sumnodesng = $sumnodesng + $numn;
	    $sumedgesng = $sumedgesng + $numenew;
	    $sumdensng = $sumdensng + $density;
	    if($numn > $maxnng){$maxnng = $numn;}
	    if($numn < $minnng){$minnng = $numn;}
	    if($numenew > $maxeng){$maxeng = $numenew;}
	    if($numenew < $mineng){$mineng = $numenew;}
	    my $ediff = $numn*($numn-1)/2 - $numenew;
	    print $outn "$almostname\t$numn\t$nume\t$numenew\t$density\t$ediff\n";
	    #implement a cograph editing with weighted edges
	    #thus hand over the file name
	    #read the file
	    #build graph and its connected components
	    #check for each connected component
	    #delete edge with smallest score (or random?)
	    #check again for connected components and if cograph
	    #iterate until cograph
	    #write the new graph and replace the old file
	    my $strkeys = join('=',keys %strscore);
	    my $strvals = join('=',values %strscore);
	    my $seqkeys = join('=',keys %seqscore);
	    my $seqvals = join('=',values %seqscore);
	    my (@newnodes,@newedges) = cographEditing(join('=',@nodes),join('=',@uniqedges),$strkeys,$strvals,$seqkeys,$seqvals);
	    my $newnodesstr = join(';',@newnodes);
	    my $newedgesstr = join(';',@newedges);
	    print STDERR "changed graph: nodes $newnodesstr\n edges:$newedgesstr\n";
	    #delete the old file and write this into a new one,
	    #thus write a complete graph with score 0 for nonedges and 1 for edges
	    my $newcurfile = "$curfile\_old";
	    my $mvcmd = "mv $curfile $newcurfile";
	    #print STDERR "mvcommand: $mvcmd\n";
	    readpipe("$mvcmd");
	    push @changedGraphs, $curfile;
	    open(my $outc, ">>",$curfile);
	    for(my $i=0;$i<(scalar @newnodes)-1; $i++){
		for(my $j=$i+1;$j<scalar @newnodes;$j++){
		    my $ed0;
		    my $ed1;
		    if($newnodes[$i] lt $newnodes[$j]){
			$ed0 = "$newnodes[$i] $newnodes[$j]";
			$ed1 = "$newnodes[$j] $newnodes[$i]";
		    }
		    else{
			$ed0 = "$newnodes[$j] $newnodes[$i]";
			$ed1 = "$newnodes[$i] $newnodes[$j]";
		    }
		    my $found = 0;
		    for(my $k=0;$k<scalar @newedges;$k++){
			if($newedges[$k] eq $ed0){
			    print $outc "$ed0 1 1\n";
			    print $outc "$ed1 1 1\n";
			    $found = 1;
			    last;
			}
		    }
		    
		    if($found == 0){
			print $outc "$ed0 0 0\n";
			print $outc "$ed1 0 0\n";
		    }
		}
	    }
	    
	}
    #    }

}    

=cut

my $avnumncg = 0;
my $avnumecg = 0;
my $avdenscg = 0;

if($numcg > 0){
    $avnumncg = sprintf("%.2f",$sumnodescg/$numcg);
    $avnumecg = sprintf("%.2f",$sumedgescg/$numcg);
    $avdenscg = sprintf("%.2f",$sumdenscg/$numcg);
}
my $avnumnng = 0;
my $avnumeng = 0;
my $avdensng = 0;

if($numng > 0){
    $avnumnng = sprintf("%.2f",$sumnodesng/$numng);
    $avnumeng = sprintf("%.2f",$sumedgesng/$numng);
    $avdensng = sprintf("%.2f",$sumdensng/$numng);
}

print $outs "===============Graph analysis\===============\n";
print $outs "Number of cographs: $numcg \n";
print $outs "Number of cliques (included in cographs): $numcliques \n";
print $outs "Number of noncographs: $numng \n";
print $outs "Cographs: 
  average node number: $avnumncg;
  average edge number: $avnumecg; 
  Max\/min number of nodes: $maxncg\/$minncg; 
  Max\/min number of edges: $maxecg\/$minecg; 
  average density:  $avdenscg \n";
print $outs "Noncographs: 
  average node number: $avnumnng; 
  average edge number: $avnumeng; 
  Max\/min number of nodes: $maxnng\/$minnng; 
  Max\/min number of edges: $maxeng\/$mineng; 
  average density:  $avdensng \n";
print $outs "\n";
if(scalar @changedGraphs > 0){
    print $outs "The following graphs were edited in order to obtain a cograph structure:\n";
    my $changedstr = join("\n",@changedGraphs);
    print $outs "$changedstr\n";
    print $outs "\n";
}


sub getInducedSubgraph{
    my @inp = @_;
    my @N = split '=', $inp[0]; #nodes of CC
    my @uniqedges = split '=', $inp[1]; #all edges of the graph
    
    my @ccedges = ();
    #get induced subgraph for CCs
    for(my $m=0;$m<scalar @N -1;$m++){
	for(my $n=$m+1;$n<scalar @N;$n++){
	    my $eddy;
	    if($N[$m] lt $N[$n]){$eddy = "$N[$m] $N[$n]";}
	    else{$eddy = "$N[$n] $N[$m]";}
	    if(none {$eddy eq $_ } @uniqedges){next;}		    
	    else{push @ccedges, $eddy;}
	}
    }

    return @ccedges;

}




#edit the noncographs
sub cographEditing{#($outpath, $struclim, $seqlim, $curfile);
    my @input = @_;
    my @nodes = split '=', $input[0];
    my @uniqedges = split '=', $input[1];
    my @strkeys = split '=', $input[2];
    my @strvals = split '=', $input[3];
    my @seqkeys = split '=', $input[4];
    my @seqvals = split '=', $input[5];

#    my $outpath = $input[0];
#    my $struclim = $input[1];
#    my $seqlim = $input[2];
#    my $curfile = $input[3];

#    my @newgraph = ();
    my @newnodes = ();
    my @newedges = ();

#    my $grapheditlist = "$outpath\/list\_edited\_graphs";
#    open(my $oute, ">>", $grapheditlist);


    my $vertices = join(',',@nodes);
#    print STDERR "vertices: $vertices\n";
    my $arcs = join(',',@uniqedges);

    my @CCs = connectedComponents($vertices,$arcs);
    
       	#CCnodes = connectedComponents(nodes,uniqedges)
	#separated with ','; edges n1 n2 with n1 < n2
	#CCnodes also one string per CC, separated with ','
    for(my $i=0;$i<scalar @CCs; $i++){
	my @V = split ' ', $CCs[$i];
	if(scalar @V == 1){
	    push @newnodes, $V[0];
	    next;
	}
	#get induced subgraph with vertices in V
	my %curseq = ();
	my %curstr = ();
	my @edges = ();
	for(my $v = 0;$v<(scalar @V)-1;$v++){
	    for(my $w = $v+1;$w<scalar @V;$w++){
		my $n1 = $V[$v];
		my $n2 = $V[$w];
		my $ed;
		if($n1 lt $n2){$ed = "$n1 $n2";}
		else{$ed = "$n2 $n1";}
		for(my $sc = 0;$sc < scalar @strkeys;$sc++){
		    if($strkeys[$sc] eq $ed){
			$curstr{$ed}=$strvals[$sc];
		    }
		}
		for(my $so = 0;$so < scalar @seqkeys;$so++){
		    if($seqkeys[$so] eq $ed){
			$curstr{$ed}=$seqvals[$so];
		    }
		}
		#print STDERR "edge: $ed\n";
		#edges are all correct
		push @edges, $ed;
	    }
	}
	my $cg = IsCograph(join(' ',@V),@edges);
	if($cg==1){
	    for(my $j = 0;$j<scalar @V;$j++){
		push @newnodes, $V[$j];
	    }
	    for(my $k=0;$k<scalar @edges;$k++){
		push @newedges, $edges[$k];
	    }
	    next;
	}
	else{
	    my @curstrkeys = keys %curstr;
	    my @curstrvals = values %curstr;
	    my @curseqkeys = keys %curseq;
	    my @curseqvals = values %curseq;
	    #delete edge with lowest val
	    #if two or more have the same val, the first one is deleted
	    my $ed2del;
	    my $lowval = 3; #something bigger than 1
	    if(scalar @curseqkeys == 0){
		#check in str
		for(my $st = 0;$st<scalar @curstrvals;$st++){
		    if($curstrvals[$st] < $lowval){
			$lowval = $curstrvals[$st];
			$ed2del = $curstrkeys[$st];
		    }
		}
	    }
	    else{
		for(my $se = 0;$se<scalar @curseqvals;$se++){
		    if($curseqvals[$se] < $lowval){
			$lowval = $curseqvals[$se];
			$ed2del = $curseqkeys[$se];
		    }
		}
	    }
	    print STDERR "edge to delete: $ed2del\n";
	    #delete ed2del
	    my @nextedges = ();
	    my @nextstrvals = ();
	    my @nextstrkeys = ();
	    my @nextseqvals = ();
	    my @nextseqkeys = ();
	    for(my $u = 0; $u<scalar @edges;$u++){
		if($edges[$u] eq $ed2del){next;}
		push @nextedges, $edges[$u];
	    }
	    for(my $sr=0;$sr<scalar @curstrkeys;$sr++){
		if($curstrkeys[$sr] eq $ed2del){next;}
		push @nextstrkeys, $curstrkeys[$sr];
		push @nextstrvals, $curstrvals[$sr];
	    }
	    for(my $sq = 0;$sq<scalar @curseqkeys;$sq++){
		if($curseqkeys[$sq] eq $ed2del){next;}
		push @nextseqkeys, $curseqkeys[$sq];
		push @nextseqvals, $curseqvals[$sq];
	    }
	    #do cographediting again
	    my $vstr = join('=',@V);
	    my $edstr = join('=',@nextedges);
	    print STDERR "cograph editing: vertex:$vstr\n edges:$edstr\n";
	    
	    my (@tmpnodes,@tmpedges) = cographEditing(join('=',@V),join('=',@nextedges),join('=',@nextstrkeys),join('=',@nextstrvals),join('=',@nextseqkeys),join('=',@nextseqvals));
	    #add to return arrays
	    for(my $j = 0;$j<scalar @tmpnodes;$j++){
		push @newnodes, $tmpnodes[$j];
	    }
	    for(my $k=0;$k<scalar @tmpedges;$k++){
		push @newedges, $tmpedges[$k];
	    }
		
	}

    }
    return (@newnodes,@newedges);

}
##check if it is a cograph
    
sub IsCograph{ #arguments is nodestring, uniqedgelist, input is a connected graph
    
    #first: get complement

    my @edges = @_;
    my @nodes = split ' ', $edges[0];
    $edges[0]="";
    if(scalar @edges <= 3){return 1;}
    
    my @uniqedges=(); #unique edges from complement
    if(scalar @nodes <4){return 1;}

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
   
    
    if(scalar @permnodes1 == 1){return 0;} ##no cograph
    elsif(scalar @permnodes1 == scalar @nodes){return 1;} ##cograph
    else{ ##recurse
	for(my $h=0;$h< scalar @permnodes1;$h++)
	{
	    my @nextedges=();
	    ##each entry is a string with nodes, get edges for it
	    my @K = split " ", $permnodes1[$h];
	    if(scalar @K <= 1){next;}
	    for(my $r=0;$r<(scalar @K -1);$r++){
		for(my $s=1;$s<scalar @K;$s++){
		    my $tmpedge;
		    if($K[$r] lt $K[$s])
		    {
			$tmpedge = "$K[$r] $K[$s]";
		    }
		    else{
			$tmpedge = "$K[$s] $K[$r]";
		    }
		    if(any{$tmpedge eq $_ } @uniqedges)
		    {
			push @nextedges, $tmpedge;
		    }
		}
	    }
	
	    my $blubb = IsCograph($permnodes1[$h], @nextedges);
	    if($blubb == 0){return 0;}
	}
	##if the loop runs through, all returned 1
	return 1;
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

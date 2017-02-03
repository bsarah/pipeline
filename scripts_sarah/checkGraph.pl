#!/usr/bin/perl -w


## perl checkGraph.pl .edlilist inpath outpath secsim strucsim psecsim pstrucsim cographlist noncographlist noEdgeGraohs_list edgegraphs summary

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
my $psecsim = shift;
my $pstrucsim = shift;
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

my $seqlim;
my $struclim;

if($psecsim >= 0 && $psecsim <= $secsim){
    $seqlim = $psecsim;
}
else{$seqlim = $secsim;}
if($pstrucsim >= 0 && $pstrucsim <= $strucsim){
    $struclim = $pstrucsim;
}
else{$struclim = $strucsim;}



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
my $minncg = 3000;
my $minecg = 3000;


my $maxnng = 0;
my $maxeng = 0;
my $minnng = 3000;
my $mineng = 3000;

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
	my $spec1 = $F1[(scalar @F1)-5];
	my $spec2 = $F2[(scalar @F2)-5];
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
	    }
	}
	
	if(none {$_ eq $n1} @nodes){
	    push @nodes,$n1;
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

    if($numn == 0 || $nume == 0){
	print $outo "$almostname \n";
	next;
    }
    else{
	print $oute "$almostname \n";
    }

    #max possible density is 2!
    my $density = 0;
    if($nume > 0){
	$density = $numn/$nume;
    }
    
    #check if it is a clique first
    my $cliquenum=(($numn-1)*$numn)/2;
    if($nume == $cliquenum){
	$numcg++;
	$numcliques++;
	$sumnodescg = $sumnodescg + $numn;
	$sumedgescg = $sumedgescg + $nume;
	$sumdenscg = $sumdenscg + $density;
	if($numn > $maxncg){$maxncg = $numn;}
	if($numn < $minncg){$minncg = $numn;}
	if($nume > $maxecg){$maxecg = $nume;}
	if($nume < $minecg){$minecg = $nume;}
	print $outcg "$almostname\t$numn\t$nume\t2.0\n";
	next;
    }
    

    #print "nodenum: $numn, edgenum: $nume\n";
    #my $edtmp = join(",",@uniqedges);
    #print "edges: $edtmp \n";
    if($density == 0 || $density == 2){
	print $outcg "$almostname\t$numn\t$nume\t$density\n";
	$numcg++;
	#$numcliques++;
	$sumnodescg = $sumnodescg + $numn;
	$sumedgescg = $sumedgescg + $nume;
	$sumdenscg = $sumdenscg + $density;
	if($numn > $maxncg){$maxncg = $numn;}
	if($numn < $minncg){$minncg = $numn;}
	if($nume > $maxecg){$maxecg = $nume;}
	if($nume < $minecg){$minecg = $nume;}
    }
    else{
	my $cg = IsCograph(join(' ',@nodes),@uniqedges);
	if($cg ==1){
	    print $outcg "$almostname\t$numn\t$nume\t$density\n";
	    $numcg++;
	    $sumnodescg = $sumnodescg + $numn;
	    $sumedgescg = $sumedgescg + $nume;
	    $sumdenscg = $sumdenscg + $density;
	    if($numn > $maxncg){$maxncg = $numn;}
	    if($numn < $minncg){$minncg = $numn;}
	    if($nume > $maxecg){$maxecg = $nume;}
	    if($nume < $minecg){$minecg = $nume;}
	}
	else{
	    print $outn "$almostname\t$numn\t$nume\t$density\n";
	    $numng++;
	    $sumnodesng = $sumnodesng + $numn;
	    $sumedgesng = $sumedgesng + $nume;
	    $sumdensng = $sumdensng + $density;
	    if($numn > $maxnng){$maxnng = $numn;}
	    if($numn < $minnng){$minnng = $numn;}
	    if($nume > $maxeng){$maxeng = $nume;}
	    if($nume < $mineng){$mineng = $nume;}

	}
    }
}    

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

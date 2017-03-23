#!/usr/bin/perl -w
## call: buildEdgeList filelist pathtoinfiles mode outpath pathtoaltnw seqsim strucsim pseudoscore summary

##if seqsim or strucsim is -1 (only one of them) then just check the other one for similarity as input to altnw
##output: weighted graph, thus node1 node2 seq_diff struc_diff

##TODO: add parameter which excludes the sequence comparison and only checks equality between types for orthology measure


use Data::Dumper;
use strict;
use warnings;

my $filename = shift;
my $inpath=shift;
my $mode = shift; #1=genelist, 0=cm
my $outpath=shift;
my $pathtonw=shift;
my $seqsim = shift;
my $strucsim = shift;
my $pseudoscore = shift;
my $summary = shift;

open FA,"<$filename" or die "can't open $filename\n";

open(my $outs, ">>$summary");

#print "file: $filename \n";

my $sumgraphs = 0;
my $sumrealedges = 0;
my $sumrealedgegraphs = 0;
my $totsumsecsim = 0;
my $totsumstrsim = 0;
my $sumalledges = 0;
my $sumspecies = 0;
my $sumnodes = 0;



while(<FA>){
    chomp;
    my $curfile=$_;

#    print "curfile: $curfile \n";
    
    open CF,"<$curfile" or die "can't open $curfile\n";

    my @cursplit = split '\.', $curfile;
    my $len = (scalar @cursplit) -1;
    my $curname=$cursplit[$len-1];
    my @newsplit = split '\/', $curname;
    my $newname = $newsplit[scalar @newsplit -1];
    my $outname="$outpath\/$newname.edli";


#    print "outfile: $outname \n";

    my @nodes=();
    my @seqs=();
    my @strucs = ();
    my @specs=();
    
    
    while(<CF>){
	my $line = $_;
	$line=~s/ /_/g;
	chomp $line;
	my @F = split '\t', $line; #chr spec_geneID start end strand blockleft blockright secstruc seq_lowercase seq_orig score

	my $nodename = "";
	my $seq = "";
	my $struc = "";
	my $spec = "";
	if($mode == 0){
	    my $chr = $F[0];
	    my $prespec = $F[1];
	    my @G = split "_", $prespec;
	    $spec = $G[0];
	    my $id = $G[1];
	    my $start = $F[2];
	    my $end = $F[3];
	    my $strand = $F[4];
	    my $leftblock = $F[5];
	    my $rightblock = $F[6];
	    $struc = $F[7];
	    $seq = $F[8];
	    my $prescore = $F[10]; ##seems to look like score_, why is there a '_'?
	    my @SC = split '_', $prescore;
	    my $score = $SC[0];
	    my $suffix = "N";
	    if($pseudoscore >= 0) {
		##check if element is pseudogenized
		if($score < $pseudoscore){
		    $suffix = "P";
		}
	    }
	    $nodename = "$chr\_$prespec\_$start\_$end\_$strand\_$suffix";
	}
	else{
    	    my $chr = $F[0];
	    my $prespec = $F[1];
	    my @G = split "_", $prespec;
	    $spec = $G[0];
	    my $id = $G[1];
	    my $start = $F[2];
	    my $end = $F[3];
	    my $strand = $F[4];
	    my $leftblock = $F[5];
	    my $rightblock = $F[6];
	    $struc = $F[7];
	    $seq = $F[8];
	    my $type = $F[9]; 
	    my $pseudogene = $F[10];
	    my $comment = $F[11];
	    $nodename = "$chr\_$prespec\_$start\_$end\_$strand\_$type\_$pseudogene";
	}



	
	push @nodes,$nodename;
	push @seqs,$seq;
	push @strucs,$struc;
	push @specs, $spec;
    }
    
    
    open(my $outf, ">>$outname");

    if((scalar @nodes)==1) {
	my $res = "$nodes[0]";
	print $outf "$res\n";
    }

    ##divide all those numbers by 2 as every edge appears twice in the edgelist!
    my $edgecount = 0; #number of edges fulfilling the similarity thresholds
    my $sumsecsim = 0;
    my $sumstrsim = 0;
    my $totedges = 0;

    for(my $i=0 ; $i < (scalar @nodes) ; $i++){
	my $res = "$nodes[$i]";
	for(my $j=0 ; $j < scalar @nodes ; $j++){
	    ##create complete weighted graphs
	    if($nodes[$i] eq $nodes[$j]){next;}
	    my $cmd1 = "$pathtonw/altNW 0 1 \"$seqs[$i]\" \"$seqs[$j]\"";
	    my $cmd2 = "$pathtonw/altNW 0 1 \"$strucs[$i]\" \"$strucs[$j]\"";

	    my @out1 = ();
	    my @out2 = ();
	    
	    if($seqsim == -1){
		@out2 = readpipe("$cmd2");
		@out1 = ("","","");
	    }
	    elsif($strucsim == -1){
		@out1 = readpipe("$cmd1");
		@out2 = ("","","");
	    }
	    else{
		@out1 = readpipe("$cmd1");
		@out2 = readpipe("$cmd2");
	    }	
	    chomp(@out1);
	    chomp(@out2);
	    my $len = length($out1[1]);
	    my $diff1 = 0;
	    if($len > 0){
		$diff1 = $out1[0]/$len;
	    }
	    my $diff_rounded = sprintf("%.2f", $diff1);
	    my $len2 = length($out2[1]);
	    my $diff2 = 0;
	    if($len2 >0){
		$diff2 = $out2[0]/$len2;
	    }
	    my $diff2_rounded = sprintf("%.2f", $diff2);
	    my $res2 = "$res $nodes[$j] $diff_rounded $diff2_rounded\n";
	    print $outf "$res2";
	    if($diff_rounded >= $seqsim && $diff2_rounded >= $strucsim){
		$edgecount++;
	    }
	    $totedges++;
	    $sumsecsim = $sumsecsim+$diff_rounded;
	    $sumstrsim = $sumstrsim + $diff2_rounded;
	}	
    }

    $sumgraphs++;
    if($edgecount > 0){$sumrealedgegraphs++;}
    $sumrealedges = $sumrealedges + $edgecount;
    $sumalledges = $sumalledges + $totedges;
    my $avsecsim = sprintf("%.2f",$sumsecsim/$sumalledges);
    my $avstrsim = sprintf("%.2f",$sumstrsim/$sumalledges);
    $totsumsecsim = $totsumsecsim + $avsecsim;
    $totsumstrsim = $totsumstrsim + $avstrsim;
    $sumspecies = $sumspecies + scalar @specs;
    $sumnodes = $sumnodes + scalar @nodes;
    
}


my $avedgenum = sprintf("%.2f",($sumalledges/2)/$sumgraphs);
my $avnodenum = sprintf("%.2f",$sumnodes/$sumgraphs);
my $avspecnum = sprintf("%.2f",$sumspecies/$sumgraphs);

my $avtotsecsim = sprintf("%.2f",($totsumsecsim/2)/$sumgraphs);
my $avtotstrsim = sprintf("%.2f",($totsumstrsim/2)/$sumgraphs);

my $avrealedgenum = sprintf("%.2f",($sumrealedges/2)/$sumrealedgegraphs);

my $percrealgraphs = sprintf("%.2f",$sumrealedgegraphs/$sumgraphs);

print $outs "===============Graph information\===============\n";
print $outs "Number of graphs: $sumgraphs\n";
print $outs "Average number of nodes per graph: $avnodenum\n";
print $outs "Average number of edges per graph: $avedgenum\n";
print $outs "Average number of species per graph: $avspecnum\n";

print $outs "Average sequence similarity score of average score for each graph: $avtotsecsim \n";
print $outs "Average structure similarity score of average score for each graph: $avtotstrsim \n";

print $outs "Number of graphs with edges fulfilling the thresholds (real edges): $sumrealedgegraphs\n";
print $outs "Percentage of graphs with real edges: $percrealgraphs \n";
print $outs "Average number of real edges in graphs with real edges: $avrealedgenum \n";
print $outs "In the graphs without edges, no edges were added because of 
either not enough nodes (thus, including singletons) or all nodes
 belonging to the same species or the similarity thresholds were not fulfiled.
 Use input parameters -s and/or -t to change the similarity thresholds. \n";
print $outs "Graph files (.edli) are located in $outpath
which contains a file listing graphs without any edges. \n";
print $outs "Graph files are edgelists whereas an edge (a,b) appears as (a,b) 
and (b,a). As graphs are weighted (similarity threshold) the graph files are 
tab separated files containing: node1, node2, sequences similarity, structure 
similarity.\n";
print $outs "If the number of graphs with edges is 0, no sequences fulfilled 
the similarity thresholds for sequence and/or structure. In this way, no graphs
 were drawn or alignments were built. If you want to get graphs with edges, 
similarity thresholds can be lowered using parameters -s and/or -t. The 
 values currently set are $seqsim and $strucsim, the default values are 0.9 meaning 
a similarity of 90%.\n";
print $outs "\n";

__END__
my $sumgraphs = 0;
my $sumrealedges = 0;
my $sumrealedgegraphs = 0;
my $totsumsecsim = 0;
my $totsumstrsim = 0;
my $sumalledges = 0;
my $sumspecies = 0;
my $sumnodes = 0;

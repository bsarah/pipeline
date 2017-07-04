#!/usr/bin/perl -w
##get one cluster list and build the edge list out of it and print out (return)

## call: buildEdgeList filelist pathtoinfiles mode outpath pathtoaltnw seqsim strucsim pseudoscore summary

##if seqsim or strucsim is -1 (only one of them) then just check the other one for similarity as input to altnw
##output: weighted graph, thus node1 node2 seq_diff struc_diff

##TODO: add parameter which excludes the sequence comparison and only checks equality between types for orthology measure


use Data::Dumper;
use strict;
use warnings;

my $instring = shift;
my $mode = shift; #1=genelist, 0=cm
my $pathtonw=shift;
my $strucsim = shift;
my $seqsim = shift;


my @Elements = split ';', $instring;

my @nodes=();
my @seqs=();
my @strucs = ();
my @specs=();
    
for(my $i=0;$i<scalar @Elements;$i++){

	my $line = $Elements[$i];
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
#	    if($pseudoscore >= 0) {
#		##check if element is pseudogenized
#		if($score < $pseudoscore){
#		    $suffix = "P";
#		}
#	    }
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
    
    
#    open(my $outf, ">>$outname");

    if((scalar @nodes)==1) {
	my $res = "$nodes[0]";
	print "$res\n";
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
	    print "$res2";
	}	
}

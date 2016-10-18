#!/usr/bin/perl -w
## call: buildEdgeList filelist pathtoinfiles outpath pathtoaltnw seqsim strucsim

##if seqsim or strucsim is -1 (only one of them) then just check the other one for similarity as input to altnw
##output: weighted graph, thus node1 node2 seq_diff struc_diff


use Data::Dumper;
use strict;
use warnings;

my $filename = shift;
my $inpath=shift;
my $outpath=shift;
my $pathtonw=shift;
my $seqsim = shift;
my $strucsim = shift;

open FA,"<$filename" or die "can't open $filename\n";

#print "file: $filename \n";

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
	my @F = split '\t', $line; #chr spec_geneID start end strand blockleft blockright secstruc seq_lowercase seq_orig
	my $chr = $F[0];
	my $prespec = $F[1];
	my @G = split "_", $prespec;
	my $spec = $G[0];
	my $id = $G[1];
	my $start = $F[2];
	my $end = $F[3];
	my $strand = $F[4];
	my $leftblock = $F[5];
	my $rightblock = $F[6];
	my $struc = $F[7];
	my $seq = $F[8];
	my $nodename = "$chr\_$prespec\_$start\_$end\_$strand";
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
    
    for(my $i=0 ; $i < (scalar @nodes) ; $i++){
	my $res = "$nodes[$i]";
	for(my $j=0 ; $j < scalar @nodes ; $j++){
	    #include edges of the same species for the duplication detection, remove them later in the checkGraph
	    #if($specs[$i] eq $specs[$j]){next;}
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
	}
	
    }
    
}

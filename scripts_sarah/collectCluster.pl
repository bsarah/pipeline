 
#!/usr/bin/perl -w
## call: collectCluster filelist pathtoinfiles pathtooutfiles
##filelist = consists of a list of files called e.g. example_genes.bed
## input file format:
## Chromosome (tab) Species_IDnum (tab) startPosition (tab) endposition
## (tab) '+' or '-' for the strand (tab) 5'blocknum (tab) 3'blocknum
## (tab) secondary structure(including (,),,,<,>,.,_,:) (tab) sequence(including *,[,],numbers,-,upper case, lower case)


## further input parameter: path to where the clus files should be written?

##ignore empty lines and lines starting with #
##when reading sequence, delete all except letters, put everything lower case, store ints on information

## output:
## Chromosome (tab) Species_IDnum (tab) startPosition (tab) endposition
## (tab) '+' or '-' for the strand (tab) 5'blocknum (tab) 3'blocknum
## (tab) secondary structure(including (,),,,<,>,.,_,:) (tab) sequence(only lower case) (tab) (intron1 start,intron1 end);(intron2 start, intron2 end)...
## without intron: last col  = (0,0);


use Data::Dumper;
use strict;
use warnings;

my @clusters = ();

my $filename = shift;
my $inpath = shift;
my $outpath = shift;
my $outf;
my $outf2;
#my $cluscount=0;


open FA,"<$filename" or die "can't open $filename\n";

while(<FA>){
    chomp;
    my $curfile = $_;
    my @cursplit = split '\.', $curfile;
    my $curname=$cursplit[0];

#    print "current file: ",$curfile,"\n";
#    print join(", ",@clusters);
#    print "\n";
    open CF,"<$curfile" or die "can't open $curfile\n";
    while(<CF>){
	my $line = $_;
	chomp $line;
	##skip empty or starting with #
	if($line =~ /^\s*#/) {next;}
	if($line =~ /^$/) {next;}
#	print "line: $line\n";
	my @F = split '\t', $line;
	my $arrlen = scalar @F; ##should be 9
#	print "arrlen = $arrlen \n";
#	print join(", ",@F); #check if indices are ok!

	my $clusstart;
	my $clusend;
	#	if($F[$arrlen-1] eq "NA"){next;}
	## clusstart is always the smaller number
	if($F[$arrlen-4] eq "None")
	{
	    $clusstart = $F[$arrlen-4];
	    $clusend = $F[$arrlen-3];
	}
	elsif($F[$arrlen-3] eq "None"){
	    $clusstart = $F[$arrlen-3];
	    $clusend = $F[$arrlen-4];
	}
	elsif($F[$arrlen-4] < $F[$arrlen-3]){
	    $clusstart = $F[$arrlen-4];
	    $clusend = $F[$arrlen-3];
	}
	else{
	    $clusstart = $F[$arrlen-3];
	    $clusend = $F[$arrlen-4];
	}

	#standardize secondary structure
	my $struc = $F[$arrlen-2];
	my $a = "<";
	my $b="(";
	my $c=">";
	my $d = ")";
	my $p = ".";
	my $k=",";
	my $dp = ":";
	my $us="_";
	my $til = "~";
	my $mini = "-";
	$struc=~s/$a/$b/g;
	$struc =~ s/$c/$d/g;
	$struc =~ s/$k/$p/g;
	$struc =~ s/$dp/$p/g;
	$struc =~ s/$us/$p/g;
	$struc =~ s/$til/$p/g;
	$struc =~ s/$mini/$p/g;
#	print "struc: $struc\n";
	$F[$arrlen-2]=$struc;
	
	# work on sequence, find intron, delete - and turn to lower case
	my $preseq = $F[$arrlen-1];
	my @introns = ();
	my @intronpos = ();
	my $pos1 = index($preseq,"*");
	my $seq;
	if($pos1<0){
	    my $noint = "(0,0);"; 
	    push @introns, $noint;
	    $seq=$preseq;
	}
	else{
	    push @intronpos, $pos1;
	    my $offset = $pos1;
	    my $position=0;
	    my $first = $pos1;
	    my $second = -1;
	    while ( $position >= 0 )
	    {
		$position = index($preseq, "*", $offset+1);
		last if ( $position < 0 );		
#		print "position= $position \n";
		if($first==-1){$first = $position; push @intronpos, $first;}
		else{
		    $second = $position+1;
		    push @intronpos, $second;
		    my $intron1 = "($first,$second);";
		    push @introns, $intron1;
		    $first = -1;
		    $second = -1;
		}

		$offset=$position;
	    }
	    push @intronpos, length($preseq);
#	    print @intronpos,"\n";
	    my @substring = ();
	    my $curpos =0;
	    for(my $i=0;$i<scalar @intronpos;){
#		print "curpos=$curpos,other=$intronpos[$i]\n";
		my $tempstr = substr $preseq, $curpos, ($intronpos[$i]-$curpos);
#		print "tempstr: $tempstr \n";
		push @substring, $tempstr;
		$curpos = $intronpos[$i+1];
		$i=$i+2;

	    }
	    $seq = join("",@substring);
	}
#	print "seq: $seq\n";

	my $minu = "-";
	$seq =~ s/$minu//g;
	my $seq2 = lc $seq;
#	print "seq2: $seq2\n";

#	my $intronvec = join ("",@introns);
	$F[$arrlen-1]=$seq2;
	push @F, $preseq;

	my $outline = join("\t",@F);
#	print "outline: $outline\n";
	
	#print "startend: $clusstart,$clusend \n";
	my $outname = "cluster-$clusstart-$clusend.clus";
#	print "outname: $outname \n";
	if ( grep( /^$outname$/, @clusters ) ) {
	    open($outf2, ">>$outpath/$outname");
	    print $outf2 "$outline \n"; 
	    close $outf2;
	}
	else {
	    push @clusters,$outname;
	    open($outf, ">$outpath/$outname");
	    print $outf "$outline \n";
	    close $outf;
	}
    }
}
##print join(", ",@clusters);
##print "\n";

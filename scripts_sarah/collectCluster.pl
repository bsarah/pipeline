 
#!/usr/bin/perl -w
## call: collectCluster filelist pathtoinfiles pseudo pathtooutfiles summary
##filelist = consists of a list of files called e.g. example_genes.bed
## input file format:
## Chromosome (tab) Species_IDnum (tab) startPosition (tab) endposition
## (tab) '+' or '-' for the strand (tab) 5'blocknum (tab) 3'blocknum
## (tab) secondary structure(including (,),,,<,>,.,_,:) (tab) sequence(including *,[,],numbers,-,upper case, lower case) (tab) Infernal score


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
my $pseudo = shift; #if -1, deactivated, if this is the infernal SCORE, all values that are SMALLER than $pseudo count for pseudo genes, if this is an evalue, it is the other way round.
my $outpath = shift;
my $summary = shift;
my $outf;
my $outf2;
my $scoresum = 0;
my $minscore = 300;
my $maxscore = 0;
my $numseqs = 0;
#my $cluscount=0;

my %species = ();
my %pseudos = ();
open(my $outs, ">>$summary");


open FA,"<$filename" or die "can't open $filename\n";

while(<FA>){
    chomp;
    my $curfile = $_;
    my @cursplit = split '\.', $curfile;
    my $curname=$cursplit[0]; #species
    
    my $elementnum = 0;
    open CF,"<$curfile" or die "can't open $curfile\n";
    while(<CF>){
	my $line = $_;
	chomp $line;
	##skip empty or starting with #
	if($line =~ /^\s*#/) {next;}
	if($line =~ /^$/) {next;}
	$elementnum++;
	my @F = split '\t', $line;
	my $arrlen = scalar @F; ##should be 10
	my $len = scalar @F;
	if($len != 10){
	    print "incorrect format of bed file $curfile! Line has $len entries and will be skipped!\n";
	    next;
	}
	
	my $clusstart;
	my $clusend;

	## clusstart is always the smaller number
	if($F[$arrlen-5] eq "None")
	{
	    $clusstart = $F[$arrlen-5];
	    $clusend = $F[$arrlen-4];
	}
	elsif($F[$arrlen-4] eq "None"){
	    $clusstart = $F[$arrlen-4];
	    $clusend = $F[$arrlen-5];
	}
	elsif($F[$arrlen-5] < $F[$arrlen-4]){
	    $clusstart = $F[$arrlen-5];
	    $clusend = $F[$arrlen-4];
	}
	else{
	    $clusstart = $F[$arrlen-4];
	    $clusend = $F[$arrlen-5];
	}

	#standardize secondary structure
	my $struc = $F[$arrlen-3];
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
	$F[$arrlen-3]=$struc;
	
	# work on sequence, find intron, delete - and turn to lower case
	my $preseq = $F[$arrlen-2];
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
	    my @substring = ();
	    my $curpos =0;
	    for(my $i=0;$i<scalar @intronpos;){
		my $tempstr = substr $preseq, $curpos, ($intronpos[$i]-$curpos);
		push @substring, $tempstr;
		$curpos = $intronpos[$i+1];
		$i=$i+2;

	    }
	    $seq = join("",@substring);
	}

	my $minu = "-";
	$seq =~ s/$minu//g;
	my $seq2 = lc $seq;

	$F[$arrlen-2]=$seq2;
	my $score = $F[$arrlen-1];
	$F[$arrlen-1] = $preseq;
	push @F, $score; ##now F got 11 entries

	$numseqs++;
	if($pseudo >= 0){
	    $scoresum = $scoresum + $score;
	    if($score > $maxscore){$maxscore = $score;}
	    if($score < $minscore){$minscore = $score;}
	    ##check for pseudogenes here and write into hash that is returned and included later again
	    ##most simple way of counting pseudogenes (extend lateron into the graphs structure)
	    if($score < $pseudo){
		if(exists $pseudos{$curname}){$pseudos{$curname}++;}
		else{$pseudos{$curname}=1;}
	    }

	}
	
	my $outline = join("\t",@F);
	my $outname = "cluster-$clusstart-$clusend.clus";
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
    $species{$curname}=$elementnum;
}

my $spstr = "";
my $numspec = 0;
foreach my $key (sort (keys(%species))) {
    my @L = split '\/', $key;
    my @M = split '\.', $L[(scalar @L) -1];
    $spstr = "$spstr$M[0] $species{$key}\n";
    $numspec++;
}


my $psestr = "";
foreach my $spi (sort (keys(%pseudos))) {
    my @L2 = split '\/', $spi;
    my @M2 = split '\.', $L2[(scalar @L2) -1];
    $psestr = "$psestr$M2[0] $pseudos{$spi}\n";
}
my $nopseudos = 0;
if($psestr eq ""){$psestr = "No pseudogenes detected.\n";$nopseudos = 1;}

if($minscore == 300){$minscore = 0;}

my $avscore = sprintf("%.2f",$scoresum/$numseqs);

print $outs "===============Species information\===============\n";
print $outs "Number of Species: $numspec\n 
Species Number_of_genetic_elements
$spstr \n";
print $outs "===============Infernal information\===============\n";
print $outs "Total number of sequences detected with infernal: $numseqs
Maximal infernal score: $maxscore
Minimal infernal score: $minscore
Average infernal score: $avscore\n";
print $outs "\n";
print $outs "===============Pseudogenes\===============\n";
print $outs "Species Number_of_pseudogenes
$psestr \n";
print $outs "\n";




my $astr="\n";
my $bstr = "=";
my $cstr = "-";
my $dstr = " ";

$spstr=~s/$astr/$bstr/g;
$spstr=~s/$dstr/$cstr/g;

if($nopseudos){$psestr = "";}
else{
    $psestr=~s/$astr/$bstr/g;
    $psestr=~s/$dstr/$cstr/g;
}

my $outstr = "$spstr\!$psestr";

print $outstr;

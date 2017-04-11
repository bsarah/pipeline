 
#!/usr/bin/perl -w
## call: collectCluster filelist mode pathtoinfiles pseudo pathtooutfiles summary
##filelist = consists of a list of files called e.g. example_genes.bed

##mode input: 1 = genelist, 0=cm

## input file format for cm mode:
## Chromosome (tab) Species_IDnum (tab) startPosition (tab) endposition
## (tab) '+' or '-' for the strand (tab) 5'blocknum (tab) 3'blocknum
## (tab) secondary structure(including (,),,,<,>,.,_,:) (tab) sequence(including *,[,],numbers,-,upper case, lower case) (tab) Infernal score

##input file format for genellist mode:
##Chromosome (tab) Species_IDnum (tab) startPosition (tab) endposition
## (tab) '+' or '-' for the strand (tab) 5'blocknum (tab) 3'blocknum
## (tab) sequence (tab) structure (tab) type (tab) pseudogene (tab) comment


## further input parameter: path to where the clus files should be written?

##ignore empty lines and lines starting with #
##when reading sequence, delete all except letters, put everything lower case, store ints on information

## output:
## Chromosome (tab) Species_IDnum (tab) startPosition (tab) endposition
## (tab) '+' or '-' for the strand (tab) 5'blocknum (tab) 3'blocknum
## (tab) secondary structure(including (,),,,<,>,.,_,:) (tab) sequence(only lower case) (tab) (intron1 start,intron1 end);(intron2 start, intron2 end)...
## without intron: last col  = (0,0);

#output for cm additional column: sequence_orig score
#output for genelist format got additional columns: type pseudogene comment


use Data::Dumper;
use strict;
use warnings;

my @clusters = ();

my $filename = shift;
my $mode = shift; #1=genelist, 0=cm
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
my %types = ();
my %difftypes = ();
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
	my $arrlen = scalar @F; ##should be 10 in cm mode, 12 in genelist mode
	my $len = scalar @F;
	if($len != 10 && $mode == 0){
	    print "incorrect format of bed file $curfile! Line has $len entries and will be skipped!\n";
	    next;
	}
	if($len != 12 && $mode == 1){
	    print "incorrect format of bed file $curfile! Line has $len entries and will be skipped!\n";
	    next;
	}
	
	my $clusstart;
	my $clusend;

	if($mode == 0){
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
	}
	else{#mode==1
	    ## clusstart is always the smaller number
	    if($F[$arrlen-7] eq "None")
	    {
		$clusstart = $F[$arrlen-7];
		$clusend = $F[$arrlen-6];
	    }
	    elsif($F[$arrlen-6] eq "None"){
		$clusstart = $F[$arrlen-6];
		$clusend = $F[$arrlen-7];
	    }
	    elsif($F[$arrlen-7] < $F[$arrlen-6]){
		$clusstart = $F[$arrlen-7];
		$clusend = $F[$arrlen-6];
	    }
	    else{
		$clusstart = $F[$arrlen-6];
		$clusend = $F[$arrlen-7];
	    }
	}
	#standardize secondary structure
	my $struc;
	if($mode==0){$struc = $F[$arrlen-3];}
	else{$struc = $F[$arrlen-4];}
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
	if($mode == 0){
	    $F[$arrlen-3]=$struc;
	}
	else{$F[$arrlen-4]=$struc;}
	    
	# work on sequence, find intron, delete - and turn to lower case

	my $preseq;
	if($mode == 0){$preseq = $F[$arrlen-2];}
	else{$preseq = $F[$arrlen-5];}
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

	if($mode ==0){$F[$arrlen-2]=$seq2;}
	else{$F[$arrlen-5]=$seq2;}

	if($mode == 0){
	    my $score = $F[$arrlen-1];
	    $F[$arrlen-1] = $preseq;
	    push @F, $score; ##now F got 11 entries
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

	}
	else{
	    my $pretype = $F[$arrlen-3];
	    if(exists $difftypes{$pretype}){$difftypes{$pretype}++;}
	    else{$difftypes{$pretype}=1;}
	    my @L3 = split '\/', $curname;
	    my @M3 = split '\.', $L3[(scalar @L3) -1];
	    my $type = "$M3[0]\_$pretype";
	    if(exists $types{$type}){$types{$type}++;}
	    else{$types{$type}=1;}

	    my $pseulabel = $F[$arrlen-2];
	    if($pseulabel eq "T" || $pseulabel eq "TRUE" ||  $pseulabel eq "true" || $pseulabel eq "t" || $pseulabel eq "1")
	    {
		if(exists $pseudos{$curname}){$pseudos{$curname}++;}
		else{$pseudos{$curname}=1;}
	    }

	    
	}
	
	$numseqs++;
	
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
    my $tmpnum = 0;
#    if(scalar @M >= 2){$tmpnum = -2;}
    $spstr = "$spstr$M[$tmpnum] $species{$key}\n";
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

my $typstr = "";
foreach my $ty (sort (keys(%types))) {
    $typstr = "$typstr$ty $types{$ty}\n";
}



if($minscore == 300){$minscore = 0;}



print $outs "===============Species information\===============\n";
print $outs "Number of Species: $numspec\n 
Species Number_of_genetic_elements
$spstr \n";
if($mode==0){
    my $avscore = sprintf("%.2f",$scoresum/$numseqs);
    print $outs "===============Infernal information\===============\n";
    print $outs "Total number of sequences detected with infernal: $numseqs
		Maximal infernal score: $maxscore
		Minimal infernal score: $minscore
		Average infernal score: $avscore\n";
    print $outs "\n";
}
else{
    #report types?
    print $outs "===============Element types\===============\n";
    print $outs "Species_Type Number_of_type\n";
    print $outs	"$typstr \n";
    print $outs "\n";
}
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

#include typeinformation in order to hand it over
my $numdifftypes = scalar (keys %difftypes);

$outstr = "$numdifftypes\!$outstr";

print $outstr;

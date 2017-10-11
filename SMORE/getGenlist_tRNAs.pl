#usage:
#perl getGenlist_tRNAs.pl /scr/k61san/trnaevo/tRNAscan/Leipzig_primates/<species>.sec <species> <species>.bed
#transform tRNAscan output to genelist format
#remove identical copies

#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $sec = shift;  	#<species>.sec
my $species = shift;	#<species>
my $output = shift;	#<species>.bed
my @array = ();
my $index = -1;

open IN, "< $sec" or die "can t open $sec\n";
open OUT,"> $output" or die "can t open $output\n";

while(<IN>){
	if($_=~m/(\S*).trna\d*\s\((\d*)-(\d*)\)/){
		$index++;
		$array[$index][0] = $1; #chromosome
		my $start = $2;
		my $end = $3;
		if ($start le $end){
			$array[$index][1] = $start; #start
			$array[$index][2] = $end; #start
			$array[$index][3] = "+"; #strand
		}
		else{
			$array[$index][2] = $start; #start
                        $array[$index][1] = $end; #start
                        $array[$index][3] = "-"; #strand
		}
		$array[$index][6] = "FALSE";  #pseudogene
	}
	elsif($_=~m/Type:\s(\S+)\s*Anticodon:\s(\S+).+\(\d*-\d*\)/){
		$array[$index][4] = $1; #tRNA type
		$array[$index][5] = $2; #anticodon
	}
	elsif($_=~m/pseudo/){
		$array[$index][6] = "TRUE";  #pseudogene
	}
	elsif($_=~m/Seq:\s(.+)/){
		$array[$index][7] = $1;  #sequence
		#print $1."\n";
        }
	elsif($_=~m/Str:\s(.+)/){
		$array[$index][8] = $1; #secondary structure
	}
}

for(my $i = 0; $i <= $index; $i++){
	#with anticodon and type
	#print OUT "$array[$i][0]\t$array[$i][1]\t$array[$i][2]\t$species\t$array[$i][3]\t$array[$i][4]_$array[$i][5]\t$array[$i][6]\t$array[$i][7]\t$array[$i][8]\tNA\n";
	#without anticodon, with type
	print OUT "$array[$i][0]\t$array[$i][1]\t$array[$i][2]\t$species\t$array[$i][3]\t$array[$i][4]\t$array[$i][6]\t$array[$i][8]\t$array[$i][7]\tNA\n";
	#type is just tRNA
	#print OUT "$array[$i][0]\t$array[$i][1]\t$array[$i][2]\t$species\t$array[$i][3]\ttRNA\t$array[$i][6]\t$array[$i][7]\t$array[$i][8]\tNA\n";

}




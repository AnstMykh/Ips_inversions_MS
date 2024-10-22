#!/usr/bin/perl
#The script to genotype inversions based on NGSadmix

use warnings;
use strict;

my $inds = $ARGV[0]; #list_of_individuals
my $admix = $ARGV[1]; #admix results 
my $out = $ARGV[2]; #output with classification

open(INDS, "< $inds") or die "file $inds not opened", $!, "\n";
open(ADMIX, "< $admix") or die "file $admix not opened", $!, "\n";
open(OUT, "> $out") or die "file $out not opened", $!, "\n";


my @inds=();
my @admix=();

while (<INDS>) {
	
	chomp $_; 
	#print $_;
	push @inds, $_;

}

while (<ADMIX>) {
	
	chomp $_; 
	#print $_;
	push @admix, $_;

}



for (my $i=0; $i<=$#inds; $i++){
	#print $admix[$i];
	my @splitted = split(/ /, $admix[$i]);
	my $ancestry = $splitted[0];
	#print $split[0];
	#print $ancestry;
	
	if ($ancestry >= 0.7){
			#print OUT $ancestry." it gonna be AA\n";
			print OUT $inds[$i]."\tAA\n";
	}
	if ($ancestry <= 0.3){
		#print OUT $ancestry." it gonna be BB\n";
		print OUT $inds[$i]."\tBB\n";
	}
	if ($ancestry >= 0.3 && $ancestry <= 0.7){
	#print OUT $ancestry." it gonna be AB\n";
	print OUT $inds[$i]."\tAB\n";
	}
	
	
}




close INDS;
close ADMIX;
close OUT;


	

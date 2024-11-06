#!/usr/bin/perl
#description
#calculate rho in 100kb windows

use warnings;
use strict;
use List::Util qw(sum);

my $file = $ARGV[0];
my $out = $ARGV[1];
my $win_size = $ARGV[2]; 
my $win_slide = $ARGV[3];##window size in kb!!

open(RHO, "< $file") or die "file $file not opened", $!, "\n";
open(OUT, "> $out") or die "file $out not opened", $!, "\n";


my @rho =();

while (<RHO>) {
	chomp $_;
	push @rho, $_;
}

my @last_snp_line= split(/\t/, $rho[$#rho]);
my $last_snp = $last_snp_line[0];
#my $contig = $last_snp_line[0];

my $number_of_windows=$last_snp/$win_slide;
#print "$number_of_windows";


my $win_start=0;
my $win_end=$win_size;


for (my $j=0; $j<=$number_of_windows; $j++){
	
	my @window=();
	
	for (my $s=0; $s<=$#rho; $s++){
	
		my @line = split(/\t/, $rho[$s]);
		if($line[0] >=$win_start && $line[0] <=$win_end ){
			push @window, $line[1];
		}
	}
	my $snp=$#window;
	
	if(@window ne 0){
	my $mean_rho = sum(@window)/@window; # correct based on the all sites passing filters
	print OUT $win_start."\t".$win_end."\t".$mean_rho."\t".$snp."\n";
	$win_start=$win_start+$win_slide;
	$win_end=$win_end+$win_slide;
	}
	else{
	print OUT $win_start."\t".$win_end."\t0.00\to\n";
	$win_start=$win_start+$win_slide;
	$win_end=$win_end+$win_slide;
	}
	

}




close(RHO);
close(OUT);

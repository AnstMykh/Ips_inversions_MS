#!/usr/bin/perl
#

use warnings;
use strict;

my $eigenvec = $ARGV[0]; #eigenvec plink file
my $out = $ARGV[1];
my $pc3 = $ARGV[2];#do you want to genotype based on PCA3?

open(PCA, "< $eigenvec") or die "file $eigenvec not opened", $!, "\n";
open(OUT, "> $out") or die "file $out not opened", $!, "\n";


my @pca=();
my @file=split(/_/, $eigenvec);
my $contig=$file[0];
my $start=$file[1];
print $start;
print $contig;


while (<PCA>) {
	
	chomp $_; 
	#print $_;
	push @pca, $_;

}

if ($pc3 eq "no"){
	if ($contig eq "IpsContig22" and $start==677434){

		for (my $i=1; $i<=$#pca; $i++){
			#print $admix[$i];
			my @splitted = split(/ /, $pca[$i]);
			my $ind = $splitted[0];
			my $pca1= $splitted[2];
			my $pca2= $splitted[3];
			my $pca3= $splitted[4];
			
			if($pca2 =~ /^-/){
				my @splitted_pca2=split(/-/, $pca2);
				my $new_pca2=$splitted_pca2[1];
				if($pca1 =~ /^-/){
					my @splitted_pca1=split(/-/, $pca1);
					my $new_pca1=$splitted_pca1[1];
					print $new_pca1."\n";
	
				
					if ($new_pca1 > 0.1 && $new_pca2 > 0){

							print OUT $ind."\tAA\n";
					}
					
					if ($new_pca1 < 0.05 && $new_pca1 > 0 && $new_pca2 > 0){
						print OUT $ind."\tAB\n";
					}
				}
				
				if ($pca1 > 0 && $new_pca2 > 0){

					print OUT $ind."\tBB\n";
				}
				


			}
			
			if($pca2 >0 && $pca1 >0){

				if ($pca1 > 0 and $pca2 > 0  and $pca2 < 0.1){
					print OUT $ind."\tBC\n";
				}
				
				if ($pca1 > 0 and $pca2 > 0.15){
					print OUT $ind."\tCC\n";
				}
	
				
			}
			if ($pca2 >0 && $pca1 =~ /^-/){
				my @splitted_pca1=split(/-/, $pca1);
				my $new_pca1=$splitted_pca1[1];
				
				if ($new_pca1 < 0.05 and $new_pca1 > 0  and $pca2 > 0){
					print OUT $ind."\tAC\n";
				}	
				
			}
		}
	}

}

if ($pc3 eq "yes"){
	if ($contig eq "IpsContig22" and $start==677434){

		for (my $i=1; $i<=$#pca; $i++){
			#print $admix[$i];
			my @splitted = split(/ /, $pca[$i]);
			my $ind = $splitted[0];
			my $pca1= $splitted[2];
			my $pca2= $splitted[3];
			my $pca3= $splitted[4];
			
			if($pca3 >=0){
				print OUT $ind."\tAA\n";
			}
			
			if($pca3 =~ /^-/ && $pca1 =~ /^-/ ) {
				my @splitted_pca3=split(/-/, $pca3);
				my $new_pca3=$splitted_pca3[1];
				my @splitted_pca1=split(/-/, $pca1);
				my $new_pca1=$splitted_pca1[1];
				
				if ($new_pca3 <= 0.1 && $new_pca1 >=0 ){
				
						print OUT $ind."\tAA\n";
				}
				if ($new_pca3 > 0.15 && $new_pca3 < 0.25 && $new_pca1 > 0 ){
				
					print OUT $ind."\tAB\n";
				}
				

			
			}
			
			if ($pca3 =~ /^-/ && $pca1 >0){
				my @splitted_pca3=split(/-/, $pca3);
				my $new_pca3=$splitted_pca3[1];
				#print $new_pca3."\n";
				
				if ($new_pca3 > 0.25 ){
					print OUT $ind."\tBB\n";
				}
				if ($new_pca3 > 0.05 && $new_pca3 < 0.25 && $pca1 > 0 ){
				#print "TUUUUU";
					print OUT $ind."\tAB\n";
				}
				if ($new_pca3 < 0.05 ){
					print OUT $ind."\tAA\n";
				}				
				
			}
			

			
		}
	}

}

if ($contig eq "IpsContig22" and $start==17434){

	for (my $i=1; $i<=$#pca; $i++){
		my @splitted = split(/ /, $pca[$i]);
		my $ind = $splitted[0];
		my $pca1= $splitted[2];
		#print $pca1."\n";
		my $pca2= $splitted[3];
		#print $pca2."\n";
		my $pca3= $splitted[4];

		if($pca2 =~ /^-/){
			#print $admix[$i];
			my @splitted_pca2=split(/-/, $pca2);
			my $new_pca2=$splitted_pca2[1];
			
			if ($new_pca2 < 0.1 ){
				#print $pca2."> -0.1\n";
					print OUT $ind."\tAA\n";
			}
			if ($new_pca2 < 0.27 && $new_pca2 > 0.1 ){
			print $new_pca2."> -0.27 and < -0.1\n";
				print OUT $ind."\tAB\n";
			}
			
			if ($new_pca2 > 0.27 ){
				print $new_pca2."< -0.27\n";
				print OUT $ind."\tBB\n";
			}
	
		}
		
		
		
		else {

			
			if ($pca2 >= 0 ){
				#print $pca2."> -0.1\n";
					print OUT $ind."\tAA\n";
			}
	
		}
	}
}

if ($contig eq "IpsContig22" and $start==417434){
	


	for (my $i=1; $i<=$#pca; $i++){
		my @splitted = split(/ /, $pca[$i]);
		my $ind = $splitted[0];
		my $pca1= $splitted[2];
		#print $pca1."\n";
		my $pca2= $splitted[3];
		#print $pca2."\n";
		my $pca3= $splitted[4];

		if($pca2 =~ /^-/){
			#print $admix[$i];
			my @splitted_pca2=split(/-/, $pca2);
			my $new_pca2=$splitted_pca2[1];
			
			if ($new_pca2 < 0.1 ){
				#print $pca2."> -0.1\n";
					print OUT $ind."\tAA\n";
			}
			if ($new_pca2 < 0.27 && $new_pca2 > 0.1 ){
			print $new_pca2."> -0.27 and < -0.1\n";
				print OUT $ind."\tAB\n";
			}
			
			if ($new_pca2 > 0.27 ){
				print $new_pca2."< -0.27\n";
				print OUT $ind."\tBB\n";
			}
	
		}
		
		
		
		else {

			
			if ($pca2 >= 0 ){
				#print $pca2."> -0.1\n";
					print OUT $ind."\tAA\n";
			}
	
		}
	}
}


if ($contig eq "IpsContig22" and $start==1917434){

	for (my $i=1; $i<=$#pca; $i++){
		#print $admix[$i];
		my @splitted = split(/ /, $pca[$i]);
		my $ind = $splitted[0];
		my $pca1= $splitted[2];
		my $pca2= $splitted[3];
		my $pca3= $splitted[4];
		
		if($pca1 > 0){
			#print " $ind I am pca1 positive\n";
			print OUT $ind."\tAA\n";
			
		}
		
		if($pca1 =~ /^-/ and $pca2 =~ /^-/){
			
			print OUT $ind."\tAA\n";
			
		}
		
		if($pca1 =~ /^-/ and $pca2 >=0){
			#print " $ind I am pca1 negative and pca2 positive\n";
			if ($pca2 >=0 && $pca2 < 0.05 ){
			#print " $ind I am pca1 positive\n";
					print OUT $ind."\tAA\n";
			}
			if ($pca2 > 0.05 and $pca2 < 0.25 ){
			print "AAA";
				print OUT $ind."\tAB\n";
			}
			
			if ($pca2 > 0.25 ){
				print OUT $ind."\tBB\n";
			}
		}
		
	}
}

if ($contig eq "IpsContig14" and $start==2077682){

	for (my $i=1; $i<=$#pca; $i++){
		#print $admix[$i];
		my @splitted = split(/ /, $pca[$i]);
		my $ind = $splitted[0];
		my $pca1= $splitted[2];
		my $pca2= $splitted[3];
		my $pca3= $splitted[4];
		

		
		if($pca2 >= 0 && $pca1 >=0){
			#print " $ind I am pca1 positive\n";
			print OUT $ind."\tAA\n";
			
		}
		
		if($pca2 >= 0 && $pca1 =~ /^-/){
			#print " $ind I am pca1 positive\n";
			print OUT $ind."\tAA\n";
			
		}
		
		if($pca1 >0 and $pca2 =~ /^-/){
			my @splitted_pca2=split(/-/, $pca2);
			my $new_pca2=$splitted_pca2[1];
			
			if ($new_pca2 >0.3){
				print OUT $ind."\tBB\n";
			}
			
			if ($new_pca2 > 0.1 and $new_pca2 < 0.25 ){
			#print "AAA";
				print OUT $ind."\tAB\n";
			}
			
			if ($new_pca2 < 0.1 ){
				print OUT $ind."\tAA\n";
			}
			
			
			
		}
		
		
		if($pca1 =~ /^-/ and $pca2 =~ /^-/){
			#print " $ind I am pca1 negative and pca2 positive\n";
			my @splitted_pca2=split(/-/, $pca2);
			my $new_pca2=$splitted_pca2[1];
			my @splitted_pca1=split(/-/, $pca1);
			my $new_pca1=$splitted_pca1[1];
			
			if ($new_pca2 > 0.1 and $new_pca2 < 0.25 ){
			#print "AAA";
				print OUT $ind."\tAB\n";
			}
			
			if ($new_pca2 < 0.1 ){
				print OUT $ind."\tAA\n";
			}
		}
		
	}
}



if ($contig eq "IpsContig14" and $start==107682){

	for (my $i=1; $i<=$#pca; $i++){
		#print $admix[$i];
		my @splitted = split(/ /, $pca[$i]);
		my $ind = $splitted[0];
		my $pca1= $splitted[2];
		my $pca2= $splitted[3];
		my $pca3= $splitted[4];
		

		
		if($pca2 >= 0 && $pca1 >=0){
			if ($pca2 >0.3){
				print OUT $ind."\tBB\n";
			}
			
			if ($pca2 > 0.1 and $pca2 < 0.25 ){
			#print "AAA";
				print OUT $ind."\tAB\n";
			}
			
		}
		
		if($pca1 > 0 && $pca2 =~ /^-/){
			#print " $ind I am pca1 positive\n";
			print OUT $ind."\tAA\n";
			
		}
		
		if($pca2 >=0 and $pca1 =~ /^-/){
			my @splitted_pca1=split(/-/, $pca1);
			my $new_pca1=$splitted_pca1[1];
			
			if ($new_pca1 <0.05){
				print OUT $ind."\tAB\n";
			}
			
			if ($new_pca1 > 0.05 ){
			#print "AAA";
				print OUT $ind."\tAA\n";
			}
			
		}
		
		
		if($pca1 =~ /^-/ and $pca2 =~ /^-/){

			print OUT $ind."\tAA\n";

		}
		
	}
}


if ($contig eq "IpsContig16"){

	for (my $i=1; $i<=$#pca; $i++){
		#print $admix[$i];
		my @splitted = split(/ /, $pca[$i]);
		my $ind = $splitted[0];
		my $pca1= $splitted[2];
		my $pca2= $splitted[3];
		my $pca3= $splitted[4];
		
		if($pca2 < 0.1){
			#print " $ind I am pca1 positive\n";
			print OUT $ind."\tAA\n";
			
		}
		
		if($pca2 >0.1 and $pca2 <0.4){
			
			print OUT $ind."\tAB\n";
			
		}
		
		if($pca2 >0.4){
			print OUT $ind."\tBB\n";
		}
		
	}
}


if ($contig eq "IpsContig23"){
	


	for (my $i=1; $i<=$#pca; $i++){
		my @splitted = split(/ /, $pca[$i]);
		my $ind = $splitted[0];
		my $pca1= $splitted[2];
		#print $pca1."\n";
		my $pca2= $splitted[3];
		#print $pca2."\n";
		my $pca3= $splitted[4];

		if($pca2 =~ /^-/){
			#print $admix[$i];
			my @splitted_pca2=split(/-/, $pca2);
			my $new_pca2=$splitted_pca2[1];
			
			if ($new_pca2 < 0.1 ){
				#print $pca2."> -0.1\n";
					print OUT $ind."\tAA\n";
			}
			if ($new_pca2 > 0.1 && $new_pca2 < 0.4 ){
			
				print OUT $ind."\tAB\n";
			}
			
			if ($new_pca2 > 0.4 ){
				
				print OUT $ind."\tBB\n";
			}
	
		}
		
		
		
		else {

			
			if ($pca2 >= 0 ){
				#print $pca2."> -0.1\n";
					print OUT $ind."\tAA\n";
			}
	
		}
	}
}

#south
if ($contig eq "IpsContig2"){

		for (my $i=1; $i<=$#pca; $i++){
			#print $admix[$i];
			my @splitted = split(/ /, $pca[$i]);
			my $ind = $splitted[0];
			my $pca1= $splitted[2];
			my $pca2= $splitted[3];
			my $pca3= $splitted[4];
			
			if($pca1 =~ /^-/){
				my @splitted_pca2=split(/-/, $pca2);
				my $new_pca2=$splitted_pca2[1];
				print OUT $ind."\tAA\n";

			}
			
			if ($pca1 <0.1 && $pca1 >0){

				print OUT $ind."\tAA\n";
				
			}
			
			
			if($pca2 >0 && $pca1 >0.2){

				print OUT $ind."\tBB\n";	
			}
			
			if ($pca1 >0.1 && $pca1 <0.2){

				print OUT $ind."\tAB\n";
				
			}
		}
}





close PCA;
close OUT;

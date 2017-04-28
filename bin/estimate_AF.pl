#!/usr/bin/perl
use strict;
use warnings;

open F, $ARGV[0] or die;
while(<F>){
	chomp;
	my @col = split("\t", $_);
	my $het = 0;
	my $hom = 0;
	foreach(@col[2..$#col]){
		if($_ =~ /A/){
			$het++;
		}
		elsif($_ =~ /a/){
			$hom++;
		}
	}
	my $total = $het + $hom;
	my $het_f = $het/$total;
	my $hom_f = $hom/$total;
	if($hom_f > 0.8 || $het_f > 0.8){
		
		next;
	}
	else{
		print "$_\n";
	}
}
close F;

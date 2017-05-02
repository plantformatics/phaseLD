#!/usr/bin/perl
use strict;
use warnings;

my %hash;
open F, $ARGV[0] or die;
while(<F>){
	chomp;
	my @col = split("\t", $_);
	if($col[2] < 0.1){
		$hash{$col[0]} = 1;
	}
}
close F;

open G, $ARGV[1] or die;
while(<G>){
	chomp;
	my @col = split("\t", $_);
	my $id = join("_", @col[0..1]);
	if(exists $hash{$id}){
		print "$_\n";
	}
}
close G;

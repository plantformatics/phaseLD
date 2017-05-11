#!/usr/bin/perl
use strict;
use warnings;

open F, $ARGV[0] or die;
while(<F>){
	chomp;
	if($_ =~ /#/){
		next;
	}
	else{
			my @col = split("\t", $_);
			print "$col[0]\t$col[1]";
			foreach(@col[9..$#col]){
				if($_ =~ /0\/0/){
					print "\ta";
				}
				elsif($_ =~ /0\/1/){
					print "\tA";
				}
				elsif($_ =~ /1\/1/){
					print "\tA";
				}
				else{
					print "\t-";
				}
			}
			print "\n";
	}
}
close F;

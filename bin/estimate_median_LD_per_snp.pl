#!/usr/bin/perl
use strict;
use warnings;
use Sort::Naturally;
use Statistics::Distributions;

my %hash;
open F, $ARGV[0] or die;
my @file = <F>;
close F;

my $win = 25;
my $step = 5;

for (my $i = 0; $i < @file - $win; $i+=$step){
	for (my $x = $i; $x < $i + $win - 1; $x++){
		for (my $t = $x + 1; $t < $i + $win; $t++){
			chomp($file[$x]);
			chomp($file[$t]);
			my @col1 = split("\t", $file[$x]);
			my @col2 = split("\t", $file[$t]);
			my $id1 = join("_",@col1[0..1]);
			my $id2 = join("_",@col2[0..1]);
			my %temp;
			$temp{3}{AA} = 0;
			$temp{3}{Aa} = 0;
			$temp{3}{aA} = 0;
			$temp{3}{aa} = 0;
			$temp{1}{A} = 0;
			$temp{1}{a} = 0;
			$temp{2}{A} = 0;
			$temp{2}{a} = 0;
			my $total = 0;
			my $skip = 0;
			for (my $z = 2; $z < @col1; $z++){
				if($col1[$z] eq '-' || $col2[$z] eq '-'){
					$skip++;
					next;
				}
				else{
					$total++;
					my $hap = join("_",$col1[$z],$col2[$z]);
					$temp{1}{$col1[$z]}++;
					$temp{2}{$col2[$z]}++;
					$temp{3}{$hap}++; 
				}
			}
			print STDERR "$id1=>$id2\tskip=$skip\ttotal=$total\n";
			if($total < 10){
				push(@{$hash{$id1}},'NA');
				push(@{$hash{$id2}},'NA');
				next;
			}
			else{
				my @sort = sort {$temp{3}{$b} <=> $temp{3}{$a}} keys $temp{3};
				my @haps = split("_",$sort[0]);
				my $pA = $temp{1}{$haps[0]}/$total;
				my $pB = $temp{2}{$haps[1]}/$total;
				if($pA == 0 || $pB == 0 || $pA == 1 || $pB == 1){
					push(@{$hash{$id1}},'NA');
					push(@{$hash{$id2}},'NA');
					next;
				}
				my $D = ($temp{3}{$sort[0]}/$total) - (($temp{1}{$haps[0]}/$total)*($temp{2}{$haps[1]}/$total));
				my $denom = $pA*(1-$pA)*$pB*(1-$pB);
				if($denom == 0){
					push(@{$hash{$id1}},'NA');
					push(@{$hash{$id2}},'NA');
					print STDERR "Some how these got through... pA=$pA pB=$pB\n";
				}
				my $r2 = ($D**2)/$denom;
				my $chi2 = $r2*$total;
				my $p_val = Statistics::Distributions::chisqrprob(1, $chi2);
				push(@{$hash{$id1}},$p_val);
				push(@{$hash{$id2}},$p_val);
			}
		}
	}
}

use List::Util qw(sum);

my @keys = nsort keys %hash;
for (my $x = 0; $x < @keys; $x++){
	my @p_vals = @{$hash{$keys[$x]}};
	print "$keys[$x]";
	my @sub;
	foreach(@p_vals){
		if($_ ne 'NA'){
			push(@sub, $_);
		}
	}
	my $ave = 0;
	my $med = 0;
	if(@sub > 0){
		$ave = mean(@sub);
		$med = median(@sub);
	}
	print "\t$ave\t$med";
#	foreach(@p_vals){
#		print "\t$_";
#	}
	print "\n";
}


sub mean{
	return sum(@_)/@_;
}

sub median {
	my @sorted = sort { $a <=> $b } @_;
  	($sorted[$#sorted/2 + 0.1] + $sorted[$#sorted/2 + 0.6])/2;
}

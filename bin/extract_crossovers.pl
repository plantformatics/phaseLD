#!/usr/bin/perl
use strict;
use warnings;
use Sort::Naturally;
use PDL;
use PDL::Stats;

die "$0 [file.bayes] [file.raw_hap]\n" unless @ARGV == 2;
my %hash;
open F, $ARGV[0] or die;
print STDERR "Loading window-based haplotypes...\n";
my $rank = 0;
while(<F>){
	chomp;
	$rank++;
	my @col = split("\t",$_);
	for (my $i = 1; $i < @col; $i++){     # caicc commented
		my $id = join("-",$rank,$col[0]);   # caicc commented
		my $adj = $i - 1;                   # caicc commented
		my @feat = split(":",$col[$i]);     # caicc commented
		my @best = split("=",$feat[0]);     # caicc commented
		$hash{$id}{$adj} = $best[0];        # caicc commented
	}
}
close F;

print STDERR "Loading SNV-based haplotypes...\n";
my %haps;
open G, $ARGV[1] or die;
while(<G>){
	chomp;
	my @col = split("\t",$_);
	my $id = join("_",@col[0..1]);
	for (my $i = 5; $i < @col; $i++){
		my $adj = $i - 5;
		$haps{$id}{$adj} = $col[$i];
	}
}
close G;

my $log_file = 'co_log.txt';
open (my $log, '>', $log_file) or die;

my @wins = nsort keys %hash;
my @snps = nsort keys %haps;
print STDERR "Begin crossover detection...\n";
for (my $i = 0; $i < @wins - 1; $i++){
	my $j = $i + 1; 
	
	# verbose
	print STDERR " - running logistic regression for bins $wins[$i] - $wins[$j]\n";

	# specify windows
	my @win1a = split("-",$wins[$i]);
	my @win2a = split("-",$wins[$j]);
	my @win1 = split(":",$win1a[1]);
	my @win2 = split(":",$win2a[1]);
	my @chr_pos_win1_1 = split("_",$win1[0]);
	my @chr_pos_win1_2 = split("_",$win1[1]);
	my @chr_pos_win2_1 = split("_",$win2[0]);
	my @chr_pos_win2_2 = split("_",$win2[1]);
	my $start1 = $chr_pos_win1_1[1];
	my $end1 = $chr_pos_win1_2[1];
	my $start2 = $chr_pos_win2_1[1];
	my $end2 = $chr_pos_win2_2[1];
	my $ave1 = int(($start1 + $end1)/2);
	my $ave2 = int(($start2 + $end2)/2);
	my $distance = $ave2 - $ave1;
	if($distance < 0){
		next;
	}
	elsif($distance > 2000000 && ($start2 > $end1)){
		next;
	}
	elsif($start2 < $end1 || $distance < 2000000){
		my @inds = nsort keys %{$hash{$wins[$i]}};
		for (my $n = 0; $n < @inds; $n++){
			my $hap1 = $hash{$wins[$i]}{$inds[$n]}; # $call1
			my $hap2 = $hash{$wins[$j]}{$inds[$n]}; # $call2
			my $break = 0;
			my %temp;
			my %types;
			if($hap1 ne $hap2){
				for (my $t = 0; $t < @snps; $t++){
					my @ids = split("_",$snps[$t]);
					if($ids[1] >= $start1 && $ids[1] <= $end2){
						if($haps{$snps[$t]}{$inds[$n]} eq '-'){
							next;
						}
						else{
							$temp{$snps[$t]} = $haps{$snps[$t]}{$inds[$n]};
							$types{$temp{$snps[$t]}} = 1;
						}
					}
					elsif($ids[1] > $end2){
						last;
					}
				}
				my @uniqtypes = keys %types;
				if(@uniqtypes > 1){
					$break++;
				}
			}
			my $little_break;
			if($break > 0){
				my @sites = nsort keys %temp;
				my @calls;
				my @dummy;
				my $cos = 0;
				my @putative;
				my @perf1;
				my @perf2;
				for (my $x = 0; $x < @sites; $x++){
					my @coord = split("_",$sites[$x]);
					my $adj = $x + 1;
					push(@dummy, $adj);
					push(@calls, $temp{$sites[$x]});
					my $call1 = $temp{$sites[$x]};
					if($x < @sites - 1){
						my $nxt = $x + 1;
						my $call2 = $temp{$sites[$nxt]};
						if($call1 ne $call2){
							$cos++;
							push(@putative, $sites[$x], $sites[$nxt]);
							push(@perf1, $call1);
							push(@perf2, $call2);
						}	
					}
				}
				if($cos == 1){
					my @sites1 = split("_", $putative[0]);
					my @sites2 = split("_",$putative[1]);
					my $len = $sites2[1] - $sites1[1];
					print "$sites1[0]\t$sites1[1]\t$sites2[1]\t$perf1[0]=>$perf2[0]\t1\t$inds[$n]\t$len\n";
				}
				elsif($cos > 1){
					my $x = pdl @dummy;
					my $y = pdl @calls;
					my %results = $y->logistic($x);
					my @interpret = keys %results;
					my $prob = $results{y_pred};
					$prob =~ s/\[//g;
					$prob =~ s/\]//g;
					chomp($prob);
					my @probs = split/\s+/,$prob;
					my $p_num = @probs;
					my $m_num = @dummy;
					my %test_diffs;
					my %test_lengths;
					for (my $t = 0; $t < @sites - 1; $t++){
						for (my $j = $t + 1; $j < @sites; $j++){
							my $h1 = $temp{$sites[$t]};
							my $h2 = $temp{$sites[$j]};
							if($h1 ne $hap1 || $h2 ne $hap2){
								next;
							}
							elsif(!$probs[$t] || !$probs[$j]){
								next;
							}
							my $dif = abs($probs[$t] - $probs[$j]);
							my $track = join("_",$t,$j);
							my @ss1 = split("_",$sites[$t]);
							my @ss2 = split("_",$sites[$j]);
							my $len = abs($ss2[1]-$ss1[1]);
							$test_lengths{$track} = $len;
							$test_diffs{$track} = $dif;
						}
					}
					my @keepers;
					my $finalInt;
					my @intervals = sort {$test_lengths{$a} <=> $test_lengths{$b}} keys %test_lengths;
					foreach(@intervals){
						if($test_diffs{$_} > 0.9){
							push(@keepers, $_);
						}
					}
					my @ids;
					if(@keepers > 0){
						@ids = split("_", $keepers[0]);
						$finalInt = $keepers[0];
					}else{
						my @breaks = sort {$test_diffs{$b} <=> $test_diffs{$a}} keys %test_diffs;
						@ids = split("_",$breaks[0]);
						$finalInt = $breaks[0];
					}
					my $one = $ids[0];
					my $two = $ids[1];
					my @first = split("_", $sites[$one]);
					my @secon = split("_", $sites[$two]);
					my $length = $secon[1] - $first[1];
					print "$first[0]\t$first[1]\t$secon[1]\t$hap1=>$hap2\t$test_diffs{$finalInt}\t$inds[$n]\t$length\n";
					for (my $z = 0; $z < @calls; $z++){
						print $log "$inds[$n]\t$sites[$z]\t$calls[$z]\t$probs[$z]\n";
					}
				}
			}
		}
	}
}	
close $log;
	

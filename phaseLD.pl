#!/usr/bin/perl
use strict;
use warnings;
use Sort::Naturally;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(sum);

######################################################
##                   Arguments                      ##
######################################################

my $gen;
my $out = 'phased.out';
my $help;
my $win = 20;
my $step = 1;
my $bwin = 50;
my $bstep = 5;

GetOptions(
        'in=s'          => \$gen,
        'out:s'         => \$out,
        'win:i'         => \$win,
        'step:i'        => \$step,
	'bwin:i'	=> \$bwin,
	'bstep:i'	=> \$bstep,
        'help!'         => \$help
) or pod2usage(-verbose => 1) && exit;

pod2usage(-verbose => 1) && exit if defined $help;

if(!$gen){
        pod2usage(-verbose => 1) && exit;
}

####################################################
##		  Phase Gentoypes 		  ##
####################################################

open F, $gen or die;
my @file = <F>;
close F;
open (my $output, '>', $out), or die;

my $logs = 'log.txt';
open (my $log, '>', $logs), or die;

my %hash;
print STDERR "Calculating LD...\n";
for (my $i = 0; $i < @file; $i+=$step){
	my $extend = 0;
	my $end;
	if(@file > $i + $win){
		$end = $i + $win;
	}
	else{
		$end = @file;
	}
	for (my $j = $i+1; $j < $end; $j++){
		chomp($file[$i]);
		chomp($file[$j]);
		my @pos1 = split("\t",$file[$i]);
		my @pos2 = split("\t",$file[$j]);
		my $id1 = join("_",@pos1[0..1]);
		my $id2 = join("_",@pos2[0..1]);
		my @results = calc_ld(\@pos1,\@pos2);
		my $linkage = $results[0];
		if($linkage eq 'NA'){
			print $log "$id1\t$id2\tNA\n";
			next;
		}
		else{
			my $r2 = $results[1];
			if($r2 == 0){
				next;
			}
			else{
				$extend++;
				my %haps = %$linkage;
				my @keys = sort {$haps{$b} <=> $haps{$a}} keys %haps;
				my $att = join(":",$id1,$id2,$r2);
				print $log "$id1\t$id2\t$r2";
				for (my $t = 0; $t < @keys; $t++){
					$hash{$att}{$keys[$t]} = $haps{$keys[$t]};
					print $log "\t$keys[$t]=$haps{$keys[$t]}";
				}
				print $log "\n";
			}
		}
	}
}
close $log;
## analyze LD information for phasing

my %phased;
my %probs;
my @keys = nsort keys %hash;
my $check = 0;
print STDERR "Begin phasing algorithm...\n";
for (my $j = 0; $j < @keys; $j++){
	my @coords = split(":",$keys[$j]);
	my @haps = sort {$hash{$keys[$j]}{$b} <=> $hash{$keys[$j]}{$a}} keys $hash{$keys[$j]};
	my $coup = 0;
	my $rep = 0;
	foreach(@haps){
		if($_ eq 'AA' || $_ eq 'aa'){
			$coup = $coup + $hash{$keys[$j]}{$_};
		}
		elsif($_ eq 'Aa' || $_ eq 'aA'){
			$rep = $rep + $hash{$keys[$j]}{$_};
		}
	}
	my $total = $coup + $rep;
	my $c_freq = $coup / $total;
	my $r_freq = $rep / $total;
	if($c_freq < 0.6 && $r_freq < 0.6){
		next;
	}
	else{
		$check++;
		if($check == 1){
			if($c_freq > $r_freq){
				push(@{$phased{$coords[0]}{0}},'A');
				push(@{$phased{$coords[0]}{1}},'a');
				push(@{$phased{$coords[1]}{0}},'A');
				push(@{$phased{$coords[1]}{1}},'a');
				push(@{$probs{$coords[0]}},$c_freq); 
				push(@{$probs{$coords[1]}},$c_freq);
			}
			else{
				push(@{$phased{$coords[0]}{0}},'A');
				push(@{$phased{$coords[0]}{1}},'a');
				push(@{$phased{$coords[1]}{0}},'a');
				push(@{$phased{$coords[1]}{1}},'A');
				push(@{$probs{$coords[0]}},$r_freq);
				push(@{$probs{$coords[1]}},$r_freq);
			}
		}
		else{
			if(exists $phased{$coords[0]}){
				my $cnt_A = 0;
				my $cnt_a = 0;
				if(@{$phased{$coords[0]}{0}} > 1){
					my @array = @{$phased{$coords[0]}{0}};
					foreach(@array){
						if($_ eq 'A'){
							$cnt_A++;
						}
						else{
							$cnt_a++;
						}
					}
				}
				else{
					if($phased{$coords[0]}{0} eq 'A'){
						$cnt_A++;
					}
					else{
						$cnt_a++;
					}
				}
				my $break;
				if($cnt_A > $cnt_a){
					if($c_freq > $r_freq){
						push(@{$phased{$coords[1]}{0}},'A');
						push(@{$phased{$coords[1]}{1}},'a');
						push(@{$probs{$coords[1]}},$c_freq);
						push(@{$probs{$coords[0]}},$c_freq);
					}
					else{
						push(@{$phased{$coords[1]}{0}},'a');
						push(@{$phased{$coords[1]}{1}},'A');
						push(@{$probs{$coords[1]}},$r_freq);
						push(@{$probs{$coords[0]}},$r_freq);
					}
				}
				else{
					if($c_freq > $r_freq){
						push(@{$phased{$coords[1]}{0}},'a');
						push(@{$phased{$coords[1]}{1}},'A');
						push(@{$probs{$coords[0]}},$c_freq);
						push(@{$probs{$coords[1]}},$c_freq);
					}
					else{
						push(@{$phased{$coords[1]}{0}},'A');
						push(@{$phased{$coords[1]}{1}},'a');
						push(@{$probs{$coords[0]}},$r_freq);
						push(@{$probs{$coords[0]}},$r_freq);
					}
				}
			}
			else{
				next;
			}
		}
	}
}

## create new hash with haplotype assignments
my %calls;
my %ave_prob;
my @sites = nsort keys %probs;
for (my $t = 0; $t < @file; $t++){
	chomp($file[$t]);
	my @cols = split("\t",$file[$t]);
	my $id = join("_",@cols[0..1]);
	if(!exists $probs{$id}){
		next;
	}
	elsif(!exists $phased{$id}){
		next;
	}
	my @prob = @{$probs{$id}};
	my $ave = sprintf("%.3f",mean(@prob));
	$ave_prob{$id} = $ave;
	my $count_A = 0;
	my $count_a = 0;
	foreach(@{$phased{$id}{0}}){
		if($_ eq 'A'){
			$count_A++;
		}
		else{
			$count_a++;
		}
	}
	if($count_A > 0 && $count_a > 0){
		next;
	}
	my $hap0 = $count_A + $count_a;
	if($hap0 < 10){
		next;
	}
	my $count_A1 = 0;
	my $count_a1 = 0;
	foreach(@{$phased{$id}{1}}){
		if($_ eq 'A'){
			$count_A1++;
		}
		else{
			$count_a1++;
		}
	}
	if($count_A1 > 0 && $count_a1 > 0){
		next;
	}
	elsif($count_A1 + $count_a1 < 10){
		next;
	}
	my $phase_0;
	my $phase_1;
	if($count_A > $count_a){
		$phase_0 = 'A';
		$phase_1 = 'a';
	}
	else{
		$phase_0 = 'a';
		$phase_1 = 'A';
	}
	my $it = 0;
	foreach(@cols[2..$#cols]){
		$it++;
		if($_ eq $phase_0){
			$calls{$id}{$it}=0;		
		}
		elsif($_ eq $phase_1){
			$calls{$id}{$it}=1;
		}
		else{
			$calls{$id}{$it}='-';
		}
	}
}

## call window phases based on bayesian probabilities
print STDERR "Begin Bayes window based haplotype calling...\n";
my @snps = nsort keys %calls;
for (my $i = 0; $i < @snps; $i++){
	my $end;
	if($i + $win > @snps){
		$end = @snps;
	}
	else{
		$end = $i + $win;
	}
	my %temp;
	my %temp_probs;
	for (my $j = $i; $j < $end; $j++){
		my @inds = nsort keys $calls{$snps[$i]};
		$temp_probs{$snps[$j]} = $ave_prob{$snps[$j]};
		for (my $z = 0; $z<@inds;$z++){
			$temp{$inds[$z]}{$snps[$j]} = $calls{$snps[$i]}{$inds[$z]};
		}
	}
	my $bayes_calls = bayes(\%temp,\%temp_probs);
	my %bcalls = %$bayes_calls;
	my @keys = nsort keys %bcalls;
	for (my $t = 0; $t < @keys; $t++){
		print $output "$keys[$t]";
		my @inds = nsort keys $bcalls{$keys[$t]};
		for (my $z = 0; $z < @inds; $z++){
			print $output "\t0=$bcalls{$keys[$t]}{$inds[$z]}{0}:1=$bcalls{$keys[$t]}{$inds[$z]}{1}";
		}
		print $output "\n";
	}
}

close $output;
####################################################
## 		    Subroutines			  ##
####################################################

sub mean{
	sum(@_)/@_;
}

sub bayes {
	my ($ref,$ref1) = @_;
	my %hash = %$ref;
	my %hash_p = %$ref1;
	my @indivs = nsort keys %hash;
	my %bcalls;
	for (my $i = 0; $i < @indivs; $i++){
		my $hap0 = 0;
		my $hap1 = 0;
		my $missing = 0;
		my @snps = nsort keys $hash{$indivs[$i]};
		my $window_snps = join(":",$snps[0],$snps[$#snps]);
		my $ave_call = 0;
		for (my $x = 0; $x < @snps; $x++){
			my $hap = $hash{$indivs[$i]}{$snps[$x]};
			$ave_call = $ave_call + $hash_p{$snps[$x]};
			if($hap eq '0'){
				$hap0++;
			}
			elsif($hap eq '1'){
				$hap1++;
			}
			else{
				$missing++;
			}
		}
		my $n_call = $ave_call/@snps;
		my $total = $hap0 + $hap1;
		my $prob_0_r = binomial($total,$hap1)*((1-$n_call)**$hap1)*($n_call**$hap0)*0.5;
		my $prob_1_r = binomial($total,$hap0)*($n_call**$hap1)*((1-$n_call)**$hap0)*0.5;
		my $total_prob = $prob_0_r + $prob_1_r;
		my $prob_0 = $prob_0_r/$total_prob;
		my $prob_1 = $prob_1_r/$total_prob;
		$bcalls{$window_snps}{$indivs[$i]}{0} = $prob_0;
		$bcalls{$window_snps}{$indivs[$i]}{1} = $prob_1;
	}
	return(\%bcalls);
}

sub binomial {
	my ($n,$k)=@_;my$r=1;$r*=$n/($n-$k),$n--while$n>$k;$r;
};

sub calc_ld{
	my ($loc1,$loc2) = @_;
	my @snp1 = @$loc1;
	my @snp2 = @$loc2;
	my %temp;
	my %out;
	for (my $t = 2; $t < @snp1; $t++){
		if($snp1[$t] ne '-' && $snp2[$t] ne '-'){
			my $hap = join("",$snp1[$t],$snp2[$t]);
			$temp{3}{$hap}++;
			$temp{1}{$snp1[$t]}++;
			$temp{2}{$snp2[$t]}++;
		}
	}
	my @keysA = sort {$temp{1}{$b} <=> $temp{1}{$a}} keys $temp{1};
	my @keysB = sort {$temp{2}{$b} <=> $temp{2}{$a}} keys $temp{2};
	my @haps  = sort {$temp{3}{$b} <=> $temp{3}{$a}} keys $temp{3};
	if(@keysA != 2 || @keysB != 2){
		return('NA');
	}
	my $tA = $temp{1}{$keysA[0]} + $temp{1}{$keysA[1]};
	my $tB = $temp{2}{$keysB[0]} + $temp{2}{$keysB[1]};
	my $r2 = 0;
	if($tA == 0 || $tB == 0){
		return('NA');
	} 
	else{
		for (my $x = 0; $x < @haps; $x++){
			my $haplotype = $haps[$x];
			my @alleles = split("",$haplotype);
			my $fpA = $temp{1}{$alleles[0]} / $tA;
			my $fpB = $temp{2}{$alleles[1]} / $tB;
			my $fpAB = $temp{3}{$haplotype} / $tA;
			if($fpA == 1){
				$fpA = 0.99;
			}
			my $break;
			if($fpB == 1){
				$fpB = 0.99;
			}
			my $Dab = $fpAB - ($fpA*$fpB);
			my $denom = $fpA*(1-$fpA)*$fpB*(1-$fpB);
			$r2 = ($Dab**2)/$denom;
			$out{$haplotype} = $fpAB;	
		}
	}
	return(\%out,$r2);
}

####################################################
##                    Help Info                   ##  
####################################################

__END__

=head1 NAME

F1 Diploid Genotype Phasing, Alexandre Marand & Hainan Zhao

=head1 SYNOPSIS

phase_F1_LD.pl --in [OPTIONS]

=head1 OPTIONS

=over 8

=item B<--in>

Genotype file [String|REQUIRED]

=item B<--out>

Output file name [String|Default=phased.out]

=item B<--win>

Set minimum window size for estimating LD [Int|Default=20]

=item B<--step>

Set minimum step size for LD calculations in a sliding window, defaults to 1 [Int|Default=1]

=item B<--help>

Print a brief help message and exit

=item B<--bwin>

Set minimum window size for bayes haplotype calling [Int|Default=50]

=item B<--bstep>

Set minimum step size for bayes haplotype calling [Int|Default=5]

=back

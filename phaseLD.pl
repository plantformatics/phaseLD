#!/usr/bin/perl
use strict;
use warnings;
use Sort::Naturally;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(sum);
use Parallel::ForkManager;

######################################################
##                   Arguments                      ##
######################################################

my $gen;
my $out = 'phased';
my $help;
my $win = 20;
my $step = 1;
my $bwin = 50;
my $bstep = 5;
my $threads = 1;
my $filt;
my $r2_threshold = 0.25;

GetOptions(
        'in|i=s'        => \$gen,
        'out|o:s'       => \$out,
        'win|w:i'       => \$win,
        'step|s:i'      => \$step,
	'rpen|r:i'	=> \$r2_threshold,
	'bwin|n:i'	=> \$bwin,
	'bstep|t:i'	=> \$bstep,
	'threads|p:i'	=> \$threads,
	'filter|f:s'	=> \$filt,
        'help|h!'         => \$help
) or pod2usage(-verbose => 1) && exit;

pod2usage(-verbose => 1) && exit if defined $help;

if(!$gen){
        pod2usage(-verbose => 1) && exit;
}

####################################################
##           Multithreading Parameters            ##
####################################################

my %bayesian_calls;
my $pm = Parallel::ForkManager->new($threads);
## allow child process to return data structure (HASH reference) upon completion.
$pm->run_on_finish( sub{
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        if($data_structure_reference){
                my $ref_type = ref($data_structure_reference);
                my %ref_hash = %$data_structure_reference;
                my @keys = nsort keys %ref_hash;
                for (my $ii = 0; $ii < @keys; $ii++){
                        my @individual = nsort keys $ref_hash{$keys[$ii]};
                        for (my $xx = 0; $xx < @individual; $xx++){
                                my @haps = nsort keys $ref_hash{$keys[$ii]}{$individual[$xx]};
                                for (my $tt = 0; $tt < @haps; $tt++){
                                        $bayesian_calls{$keys[$ii]}{$individual[$xx]}{$haps[$tt]} = $ref_hash{$keys[$ii]}{$individual[$xx]}{$haps[$tt]};
                                }
                        }
               }
        }
        else{
                print STDERR "$pid did not send anything back\n";
        }
});

print STDERR "Running with $threads threads...\n";

####################################################
##		  Phase Gentoypes 		  ##
####################################################

open F, $gen or die;
my @file = <F>;
close F;

## remove bad markers if filter flag is one
my %remove;
if($filt){
	print STDERR "Removing markers from $filt before processing...\n";
	open G, $filt or die;
	while(<G>){
		chomp;
		my @col = split("\t",$_);
		$remove{$col[0]} = 1;
	}
	close G;
	my @file1;
	for (my $t = 0; $t < @file; $t++){
		my $line = $file[$t];
		chomp($file[$t]);
		my @col = split("\t",$file[$t]);
		my $id = join("_",@col[0..1]);
		if(exists $remove{$id}){
			next;
		}
		else{
			push(@file1,$line);
		}
	}
	@file = @file1;
}

## phase genotypes
my $out1 = $out . '.out';
my $logs = $out . '.log';
my $log_hap = $out . '.raw_haplotypes';
my $bad_log = $out . '.bad';
open (my $output, '>', $out1) or die;
open (my $log, '>', $logs) or die;
open (my $rawhap, '>', $log_hap) or die;
open (my $bad, '>', $bad_log) or die;

my %hash;
my %probs;
my $check = 0;
print STDERR "Calculating LD...\n";
for (my $i = 0; $i < @file; $i+=$step){
	my $end;
	if(@file > $i + $win){
		$end = $i + $win;
	}
	else{
		$end = @file;
	}
	my $current_marker_pass = 0;
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
			my %haps = %$linkage;
			my @keys = sort {$haps{$b} <=> $haps{$a}} keys %haps;
			print $log "$id1\t$id2\t$r2";
			my $coupling = 0;
			my $repulsion = 0;
			for (my $t = 0; $t < @keys; $t++){
				if($keys[$t] eq 'AA' || $keys[$t] eq 'aa'){
					$coupling = $coupling + $haps{$keys[$t]};
				}
				elsif($keys[$t] eq 'Aa' || $keys[$t] eq 'aA'){
					$repulsion = $repulsion + $haps{$keys[$t]};
				}
				print $log "\t$keys[$t]=$haps{$keys[$t]}";
			}
			print $log "\n";
			if($coupling < 0.7 && $repulsion < 0.7 || $r2 < $r2_threshold){
				$current_marker_pass++;
				if($current_marker_pass > $win/2){
					print $bad "$id1\tconsistently_poor_r2\n";
					last;
				}
				next;
			}
			elsif($r2 >= $r2_threshold){
				$check++;
			}
			my $break1;
			if($check == 1){
				if($coupling > $repulsion){
					push(@{$hash{$id1}{0}},'A');
					push(@{$hash{$id2}{0}},'A');
					push(@{$hash{$id1}{1}},'a');
					push(@{$hash{$id2}{1}},'a'); 
					push(@{$probs{$id1}},$coupling);
					push(@{$probs{$id2}},$coupling);
				}
				else{
					push(@{$hash{$id1}{0}},'A');
					push(@{$hash{$id2}{0}},'a');
					push(@{$hash{$id1}{1}},'a');
					push(@{$hash{$id2}{1}},'A');
					push(@{$probs{$id1}},$repulsion);
					push(@{$probs{$id2}},$repulsion);
				}
			}
			elsif($check > 0){
				if(exists $hash{$id1}){
					my $cnt_A = 0;
					my $cnt_a = 0;
					my @array = @{$hash{$id1}{0}};
					foreach(@array){
						if($_ eq 'A'){
							$cnt_A++;
						}
						else{
							$cnt_a++;
						}
					}
					if($coupling > $repulsion){
						if($cnt_A >= $cnt_a){
							push(@{$hash{$id2}{0}},'A');
							push(@{$hash{$id2}{1}},'a');
						}
						else{
							push(@{$hash{$id2}{0}},'a');
							push(@{$hash{$id2}{1}},'A');
						}
						push(@{$probs{$id1}},$coupling);
						push(@{$probs{$id2}},$coupling);
					}
					else{
						if($cnt_A >= $cnt_a){
							push(@{$hash{$id2}{0}},'a');
							push(@{$hash{$id2}{1}},'A');
						}
						else{
							push(@{$hash{$id2}{0}},'A');
							push(@{$hash{$id2}{1}},'a');
						}
						push(@{$probs{$id1}},$repulsion);
						push(@{$probs{$id2}},$repulsion);
					}		
				}
				elsif(!exists $hash{$id1} && exists $hash{$id2}){
					my $cnt_A = 0;
					my $cnt_a = 0;
					my @array = @{$hash{$id2}{0}};
					foreach(@array){
						if($_ eq 'A'){
							$cnt_A++;
						}
						else{
							$cnt_a++;
						}
					}
					if($coupling > $repulsion){
						if($cnt_A >= $cnt_a){
							push(@{$hash{$id1}{0}},'A');
							push(@{$hash{$id1}{1}},'a');
						}
						else{
							push(@{$hash{$id1}{0}},'a');
							push(@{$hash{$id1}{1}},'A');
						}
						push(@{$probs{$id1}},$coupling);
						push(@{$probs{$id2}},$coupling);
					}
					else{
						if($cnt_A >= $cnt_a){
							push(@{$hash{$id1}{0}},'a');
							push(@{$hash{$id1}{1}},'A');
						}
						else{
							push(@{$hash{$id1}{0}},'A');
							push(@{$hash{$id1}{1}},'a');
						}
						push(@{$probs{$id1}},$repulsion);
						push(@{$probs{$id2}},$repulsion);
					}
				}
			}
		}
	}
}
close $log;

################################################
## create new hash with haplotype assignments ##
################################################

print STDERR "Assigning haplotypes to SNP calls...\n";
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
	elsif(!exists $hash{$id}){
		next;
	}
	my @prob = @{$probs{$id}};
	my $ave = sprintf("%.3f",mean(@prob));
	if($ave < 0.7){
		print $bad "$id\tlow_hap_frequency\n";
		next;
	}
	$ave_prob{$id} = $ave;
	my $count_A = 0;
	my $count_a = 0;
	foreach(@{$hash{$id}{0}}){
		if($_ eq 'A'){
			$count_A++;
		}
		else{
			$count_a++;
		}
	}
	my $hap0 = $count_A + $count_a;
	my $freq_A = $count_A / $hap0;
	my $freq_a = $count_a / $hap0;
	if($freq_A < 0.8 && $freq_a < 0.8){
		print $bad "$id\tsite_freq_less_than_0.8\n";
		next;
	}
	elsif($hap0 < $win/(3*$step)){
		print $bad "$id\tlow_counts_corroboration\n";
		next unless $hap0 == 1;
	}
	my $count_A1 = 0;
	my $count_a1 = 0;
	foreach(@{$hash{$id}{1}}){
		if($_ eq 'A'){
			$count_A1++;
		}
		else{
			$count_a1++;
		}
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
	print $rawhap "$cols[0]\t$cols[1]\t$phase_0\t$phase_1\t$ave";
	foreach(@cols[2..$#cols]){
		$it++;
		if($_ eq $phase_0){
			$calls{$id}{$it}=0;		
			print $rawhap "\t1";
		}
		elsif($_ eq $phase_1){
			$calls{$id}{$it}=1;
			print $rawhap "\t2";
		}
		else{
			$calls{$id}{$it}='-';
			print $rawhap "\t3";
		}
	}
	print $rawhap "\n";
}
close $rawhap;
close $bad;

########################################################
## call window phases based on bayesian probabilities ##
########################################################

print STDERR "Begin Bayes window based haplotype calling...\n";
my @snps = nsort keys %calls;
for (my $i = 0; $i < @snps-$bwin+$bstep; $i+=$bstep){
	my $pid = $pm->start and next;
	my $end;
	if($i + $bwin > @snps){
		$end = @snps;
	}
	else{
		$end = $i + $bwin;
	}
	my $win_snp = join(":",$snps[$i],$snps[$end-1]);
	print STDERR "Current window => $win_snp, iteration $i in progress...\n";
	my %temp;
	my %temp_probs;
	for (my $j = $i; $j < $end; $j++){
		my @inds = nsort keys $calls{$snps[$j]};
		$temp_probs{$snps[$j]} = $ave_prob{$snps[$j]};
		for (my $z = 0; $z<@inds;$z++){
			$temp{$inds[$z]}{$snps[$j]} = $calls{$snps[$j]}{$inds[$z]};
		}
	}
	my $bayes_calls = bayes(\%temp,\%temp_probs);
	my %sub = %$bayes_calls;
	$pm->finish(0, \%sub);
}
$pm->wait_all_children;

my @w_calls = nsort keys %bayesian_calls;
for (my $z = 0; $z < @w_calls; $z++){
	print $output "$w_calls[$z]";
	my @inds = nsort keys $bayesian_calls{$w_calls[$z]};
	for (my $i = 0; $i < @inds; $i++){
		print $output "\t";
		my @hap = sort {$bayesian_calls{$w_calls[$z]}{$inds[$i]}{$b} <=> $bayesian_calls{$w_calls[$z]}{$inds[$i]}{$a}} keys $bayesian_calls{$w_calls[$z]}{$inds[$i]};
		for (my $j = 0; $j < @hap; $j++){
			print $output "$hap[$j]=$bayesian_calls{$w_calls[$z]}{$inds[$i]}{$hap[$j]}:";
		}
	}
	print $output "\n";
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
		my $prob_1_r = binomial($total,$hap1)*($n_call**$hap1)*((1-$n_call)**$hap0)*0.5;
		my $total_prob = $prob_0_r + $prob_1_r;
		my $prob_0 = $prob_0_r/$total_prob;
		my $prob_1 = $prob_1_r/$total_prob;
		$bcalls{$window_snps}{$indivs[$i]}{0} = sprintf("%.10f", $prob_0);
		$bcalls{$window_snps}{$indivs[$i]}{1} = sprintf("%.10f", $prob_1);
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

phaseLD.pl --in <foo.gen> [OPTIONS]

=head1 OPTIONS

=over 14

=item B<--in|-i>

Genotype file [String|REQUIRED]

=item B<--out|-o>

Prefix for output files [String|Default='phased']

=item B<--win|-w>

Set minimum window size for estimating LD [Int|Default=20]

=item B<--step|-s>

Set minimum step size for LD calculations in a sliding window, defaults to 1 [Int|Default=1]

=item B<--rpen|-r>

Set minimum r2 value needed to retain SNPs. This threshold helps to remove false positive SNP call. Defaults to 0.25. [Int|Default=0.25]

=item B<--bwin|-n>

Set minimum window size for bayes haplotype calling [Int|Default=50]

=item B<--bstep|-t>

Set minimum step size for bayes haplotype calling [Int|Default=5]

=item B<--threads|-p>

Set number of threads [Int|Default=1]

=item B<--filter|-f>

Optionally use file with suffix 'bad' generated by phaseLD to filter poor markers. 'prefix.bad' contains IDs of bad SNPs identified in the first pass. Useful to cleaning haplotype switch problems [String|Optional]

=item B<--help|-h>

Print a brief help message and exit

=back

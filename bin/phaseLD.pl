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
my $rpen = 0.25;
my $fast = '';
my $quick = '';

GetOptions(
        'in|i=s'        => \$gen,
        'out|o:s'       => \$out,
        'win|w:i'       => \$win,
        'step|s:i'      => \$step,
	'rpen|r:f'	=> \$rpen,
	'bwin|n:i'	=> \$bwin,
	'bstep|t:i'	=> \$bstep,
	'fast_mode|f'	=> \$fast,
	'quick_mode|q'	=> \$quick,
	'threads|p:i'	=> \$threads,
	'filter|f:s'	=> \$filt,
        'help|h!'       => \$help
) or pod2usage(-verbose => 99) && exit;

pod2usage(-verbose => 99) && exit if defined $help;

if(!$gen){
        pod2usage(-verbose => 99) && exit;
}

my $time1 = ts();

print STDERR "#########################\n";
print STDERR "## Selected Parameters ##\n";
print STDERR "#########################\n\n";
print STDERR "Time: $time1\n\n";
if(defined $gen){
	print STDERR "--in		= $gen\n";
}
if(defined $out){
	print STDERR "--out		= $out\n";
}
if(defined $win){
	print STDERR "--win 		= $win\n";
}
if(defined $step){
	print STDERR "--step 		= $step\n";
}
if(defined $rpen){
	print STDERR "--rpen 		= $rpen\n";
}
if(defined $bwin){
	print STDERR "--bwin		= $bwin\n";
}
if(defined $bstep){
	print STDERR "--bstep		= $bstep\n";
}
if(defined $threads){
	print STDERR "--threads	= $threads\n";
}
if(defined $filt){
	print STDERR "--filter	= $filt\n";
}
if($fast eq 1){
	print STDERR "--fast_mode	enabled\n";
}
if($quick eq 1){
	print STDERR "--quick_mode	enabled\n";
}
print STDERR "\n";

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


####################################################
##		  Phase Gentoypes 		  ##
####################################################

open F, $gen or die "Could not open file $gen\n";
my @file = <F>;
close F;

## remove bad markers if filter flag is on
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

## Output files
my $out1 = $out . '.bayes';
my $logs = $out . '.log';
my $log_hap = $out . '.raw_hap';
#my $bad_log = $out . '.bad';
open (my $output, '>', $out1) or die;
open (my $log, '>', $logs) or die;
open (my $rawhap, '>', $log_hap) or die;
#open (my $bad, '>', $bad_log) or die;

#####################
## Phase Genotypes ##
#####################

my %hash;
my %probs;
my $check = 0;
my $time2 = ts();
my %skipped;
print STDERR "$time2\t...Calculating LD...\n";
for (my $i = 0; $i < @file; $i+=$step){
	my $end;
	if(@file > $i + $win){
		$end = $i + $win;
	}
	else{
		$end = @file;
	}
	my $bingo = 0;
	my $its = 0;
	my @good_array = ();
	my $bounce = 0;
	for (my $j = $i+1; $j < $end; $j++){
		$its++;
		if($bounce > 0){
			last;
		}
		elsif($j == $end - 1){
			if($good_array[0]){
				$i = $good_array[0] - 1;
				last;
			}
			else{
				my $try_again = $end + 2;
				if(@file > $try_again){
					$end = $try_again;
					redo;
				}
				else{
					next;
				}
			}
		}
		chomp($file[$i]);
		chomp($file[$j]);
		my @pos1 = split("\t",$file[$i]);
		my @pos2 = split("\t",$file[$j]);
		my $id1 = join("_",@pos1[0..1]);
		my $id2 = join("_",@pos2[0..1]);
		my @results = calc_ld(\@pos1,\@pos2);
		my $linkage = $results[0];
		if($linkage eq 'NA'){
			print $log "$its...$id1\t$id2\tNA\n";
			$skipped{$id1}++;
			$skipped{$id2}++;
			next;
		}
		else{
			my %haps = %$linkage;
			my $r2 = $results[1];
			my $pr2 = sprintf("%.6f",$r2);
			my @keys = sort {$haps{$b} <=> $haps{$a}} keys %haps;
			print $log "A:$check..L:$its..B:$bingo\t\t$id1\t$id2\t$pr2\tT_$rpen";
			my $dist1 = 0;
			my $dist2 = 0;
			my $dist3 = 0;
			my $dist4 = 0;
                        foreach(@keys){
				my $val = sprintf("%.4f",$haps{$_});
                                print $log "\t$_=$val";
				if($_ eq 'AA' || $_ eq 'Aa'){
					$dist1 = $haps{$_} + $dist1;
				}
				if($_ eq 'AA' || $_ eq 'aA'){
					$dist2 = $haps{$_} + $dist2;
				}
				if($_ eq 'aa' || $_ eq 'Aa'){
					$dist3 = $dist3 + $haps{$_};
				}
				if($_ eq 'aa' || $_ eq 'aA'){
					$dist4 = $dist4 + $haps{$_};
				}
                        }
			if($r2 < $rpen){
				print $log "\tskipped\n";
				$skipped{$id1}++;
				$skipped{$id2}++;
				next;
			}
			elsif($haps{AA} > 0.75 || $haps{aa} > 0.75 || $haps{Aa} > 0.75 || $haps{aA} > 0.75){
				if($r2 < 0.95){
					print $log "\tskipped\n";
					$skipped{$id1}++;
					$skipped{$id2}++;
					next;
				}
			}
			elsif($dist1 > 0.8 || $dist2 > 0.8 || $dist3 > 0.8 || $dist4 > 0.8){
				print $log "\tskipped\n";
				$skipped{$id1}++;
				$skipped{$id2}++;
				next;
			}
			my $best = $keys[0];
			my $rev = $best;
			$rev =~ tr/Aa/aA/;
			my $dif = $haps{$best} - $haps{$rev};
			if($dif > 0.5){
				print $log "\tskipped\n";
				$skipped{$id1}++;
				$skipped{$id2}++;
				next;
			}
			elsif($haps{$rev} < 0.15){
				print $log "\tskipped\n";
				$skipped{$id1}++;
				$skipped{$id2}++;
				next;
			}
			elsif($haps{$rev} + $haps{$best} < 0.8){
				if($rev ne $keys[1]){
					print $log "\tskipped\n";
					$skipped{$id1}++;
					$skipped{$id2}++;
					next;
				}
				else{
					print $log "\tretained\n";
					$bingo++;
					push(@good_array,$j);
				}
			}
			else{
				print $log "\tretained\n";
				$bingo++;
				push(@good_array, $j);
			}
			my $coupling = 0;
			my $repulsion = 0;
			for (my $t = 0; $t < @keys; $t++){
				if($keys[$t] eq 'AA' || $keys[$t] eq 'aa'){
					$coupling = $coupling + $haps{$keys[$t]};
				}
				elsif($keys[$t] eq 'Aa' || $keys[$t] eq 'aA'){
					$repulsion = $repulsion + $haps{$keys[$t]};
				}
			}
			$check++;
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
				my $b;
				if($quick){
					if($fast){
						print STDERR "Error: --fast_mode and --quick_mode cannot be used together\n";
						exit;
					}
					else{
						if(($bingo > 10 && $win < 100) || ($bingo > $win/20 && $its < $win && $win > 100)){
							$i = $good_array[0] - 1;
							$bounce++;
							last;
						}
						elsif($its == $win && $bingo > 0){
							$i = $good_array[0] - 1;
							$bounce++;
							last;
						}
					}
				}
				if($fast){
					if($quick){
						print STDERR "Error: --fast_mode and --quick_mode cannot be used together\n";
						exit;
					}
					else{
						$i = $j - 1;
						$bounce++;
						last;
					}
				}
				if($bounce > 0){
					last;
				}
			}
		}
	}
}
close $log;

################################################
## create new hash with haplotype assignments ##
################################################

my $time3 = ts();
print STDERR "$time3\t...Assigning haplotypes to SNP calls...\n";
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
			print $rawhap "\t0";
		}
		elsif($_ eq $phase_1){
			$calls{$id}{$it}=1;
			print $rawhap "\t1";
		}
		else{
			$calls{$id}{$it}='-';
			print $rawhap "\t-";
		}
	}
	print $rawhap "\n";
}
close $rawhap;
#close $bad;

########################################################
## call window phases based on bayesian probabilities ##
########################################################

my $time4 = ts();
print STDERR "$time4\t...Begin Bayes window based haplotype calling...\n";
my @snps = nsort keys %calls;
for (my $i = 0; $i < @snps-$bwin+$bstep; $i+=$bstep){
	my $pid = $pm->start and next;
	my $end;
	my $time5 = ts();
	if($i + $bwin > @snps){
		$end = @snps;
	}
	else{
		$end = $i + $bwin;
	}
	my $win_snp = join(":",$snps[$i],$snps[$end-1]);
	print STDERR "$time5\t...Current window => $win_snp, iteration $i in progress...\n";
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
		if($total < 5){
			$bcalls{$window_snps}{$indivs[$i]}{0} = 0.5000000000;
	                $bcalls{$window_snps}{$indivs[$i]}{1} = 0.5000000000;
			next;
		}
		my $p0_bc = binomial($total,$hap1);
		my $p1_bc = binomial($total,$hap0);
		my $p0_top = $p0_bc * ((1-$n_call)**$hap1) * ($n_call**$hap0) * 0.5;
		my $p1_top = $p1_bc * ((1-$n_call)**$hap0) * ($n_call**$hap1) * 0.5;
		my $total_prob = $p0_top + $p1_top;
		my $prob_0 = $p0_top/$total_prob;
		my $prob_1 = $p1_top/$total_prob;
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
	my $total1 = 0;
	$temp{1}{A} = 0;
	$temp{1}{a} = 0;
	$temp{2}{A} = 0;
	$temp{2}{a} = 0;
	$temp{3}{AA} = 0;
	$temp{3}{aa} = 0;
	$temp{3}{Aa} = 0;
	$temp{3}{aA} = 0;
	for (my $t = 2; $t < @snp1; $t++){
		if($snp1[$t] ne '-' && $snp2[$t] ne '-'){
			$total1++;
			my $hap = join("",$snp1[$t],$snp2[$t]);
			$temp{3}{$hap}++;
			$temp{1}{$snp1[$t]}++;
			$temp{2}{$snp2[$t]}++;
		}
	}
	my @keysA = sort {$temp{1}{$b} <=> $temp{1}{$a}} keys $temp{1};
	my @keysB = sort {$temp{2}{$b} <=> $temp{2}{$a}} keys $temp{2};
	my @haps  = sort {$temp{3}{$b} <=> $temp{3}{$a}} keys $temp{3};
	my $r2 = 0;
	if($total1 == 0){
		return('NA');
	} 
	else{
		for (my $x = 0; $x < @haps; $x++){
			my $haplotype = $haps[$x];
			my @alleles = split("",$haplotype);
			my $fpA = $temp{1}{$alleles[0]} / $total1;
			my $fpa = 1 - $fpA;
			my $fpB = $temp{2}{$alleles[1]} / $total1;
			my $fpb = 1 - $fpB;
			my $fpAB = $temp{3}{$haplotype} / $total1;
			if($fpA == 0){
				$fpA = 1/$total1;
			}
			if($fpa == 0){
				$fpa = 1/$total1;
			}
			if($fpB == 0){
				$fpB = 1/$total1;
			}
			if($fpb == 0){
				$fpb = 1/$total1;
			}
			my $Dab = $fpAB - ($fpA*$fpB);
			my $denom = $fpA*$fpa*$fpB*$fpb;
			$r2 = ($Dab**2)/$denom;
			$out{$haplotype} = $fpAB;	
		}
	}
	return(\%out,$r2);
}

sub ts {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%02d|%02d|%04d-%02d:%02d:%02d",
                                   $mon+1,$mday,$year+1900,$hour,$min,$sec);
    return $nice_timestamp;
}

####################################################
##                    Help Info                   ##  
####################################################

__END__

=head1 NAME

phaseLD v.0.01, developed by Alexandre Marand & Hainan Zhao, 2017

=head1 SYNOPSIS

phaseLD.pl --in <file.gen> [OPTIONS]

=head1 DESCRIPTION

A complete list of parameter options. See 'https://github.com/plantformatics/phaseLD' for more information. 

=over 14

=back

=head2 Input/Output

=over 20

=item B<--in|-i>

Genotype file [String|REQUIRED]

=item B<--out|-o>

Prefix for output files [String|Default='phased']

=back

=head2 Performance 

=over 20

=item B<--win|-w>

Set minimum window size for estimating Pairwise LD (when --quick_mode or --fast_mode are False). If --quick_mode or --fast_mode are enabled, this distance becomes the search area for linked markers. Windows are dynamically extended if no markers within the specified window are in linkage with the current iteration. Markers which are not linked with neighboring markers are typically false positives. See prefix.log to identify markers with low linkage. [Int|Default=20]

=item B<--step|-s>

Set minimum step size for LD calculations in a sliding window, defaults to 1. [Int|Default=1]

=item B<--bwin|-n>

Set minimum window size for bayes haplotype calling. [Int|Default=50]

=item B<--bstep|-t>

Set minimum step size for bayes haplotype calling. [Int|Default=5]

=item B<--rpen|-r>

Set minimum r2 value needed to retain SNPs. This threshold helps to remove false positive SNP call. Defaults to 0.25. [Int|Default=0.25]

=back

=head2 Run Time Options

=over 20

=item B<--threads|-p>

Set number of threads. [Int|Default=1]

=item B<--quick_mode|-f>

Instead of calculating LD for every marker in the window, calculate LD for window_size/10. The next position is the first marker above the R2 threshold. This makes the computation 10-20X faster for the default window size. Time saved increases exponentially as window sizes get larger. [Boolean|Default=False]

=item B<--fast_mode|-q>

Rather than use pairwise LD calculations within a window, use LD chaining instead. Starting at marker X, the algorithm searches for the nearest marker, y, with R2 above the threshold. The next iteration starts at marker y and continues until all markers are exhausted. [Boolean|Default=False] 

=item B<--filter|-f>

Optionally use a file to mask markers. This file should contain the chromosome and position of the marker (example: chr01_1000) on a single line. Useful to cleaning haplotype switch problems or thinning data to specific markers. [String|Optional]

=back

=head2 Misc

=over 20

=item B<--help|-h>

Print a brief help message and exit.

=back

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

GetOptions(
        'in=s'          => \$gen,
        'out:s'         => \$out,
        'win:i'         => \$win,
        'step:i'        => \$step,
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
for (my $i = 0; $i < @file - $win; $i+=$step){
	my $extend = 0;
	for (my $j = $i+1; $j < $i + $win; $j++){
		chomp($file[$i]);
		chomp($file[$j]);
		my @pos1 = split("\t",$file[$i]);
		my @pos2 = split("\t",$file[$j]);
		my $id1 = join("_",@pos1[0..1]);
		my $id2 = join("_",@pos2[0..1]);
		my @results = calc_ld(\@pos1,\@pos2);
		my $linkage = $results[0];
		if($linkage eq 'NA'){
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
	if($extend < $win/2){
		$win = $win + $step;
		redo;
	}
	else{
		$win = $win;
	}
}
close $log;
## analyze LD information for phasing

my %phased;
my %probs;
my @keys = nsort keys %hash;
my $check = 0;
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
				$phased{$coords[0]}{0} = 'A';
				$phased{$coords[0]}{1} = 'a';
				$phased{$coords[1]}{0} = 'A';
				$phased{$coords[1]}{1} = 'a';
				push(@{$probs{$coords[0]}},$c_freq); 
				push(@{$probs{$coords[1]}},$c_freq);
			}
			else{
				$phased{$coords[0]}{0} = 'A';
				$phased{$coords[0]}{1} = 'a';
				$phased{$coords[1]}{0} = 'a';
				$phased{$coords[1]}{1} = 'A';
				push(@{$probs{$coords[0]}},$r_freq);
				push(@{$probs{$coords[1]}},$r_freq);
			}
		}
		else{
			if(exists $phased{$coords[0]}){
				if($phased{$coords[0]}{0} eq 'A'){
					if($c_freq > $r_freq){
						$phased{$coords[1]}{0} = 'A';
						$phased{$coords[1]}{1} = 'a';
						push(@{$probs{$coords[1]}},$c_freq);
						push(@{$probs{$coords[0]}},$c_freq);
					}
					else{
						$phased{$coords[1]}{0} = 'a';
						$phased{$coords[1]}{1} = 'A';
						push(@{$probs{$coords[1]}},$r_freq);
						push(@{$probs{$coords[0]}},$r_freq);
					}
				}
				else{
					if($c_freq > $r_freq){
						$phased{$coords[1]}{0} = 'a';
						$phased{$coords[1]}{1} = 'A';
						push(@{$probs{$coords[0]}},$c_freq);
						push(@{$probs{$coords[1]}},$c_freq);
					}
					else{
						$phased{$coords[1]}{0} = 'A';
						$phased{$coords[1]}{1} = 'a';
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
	my $ave = mean(@prob);
	print $output "$sites[$t]\t$phased{$id}{0}\t$phased{$id}{1}\t$ave";
	foreach(@cols[2..$#cols]){
		if($_ eq $phased{$id}{0}){
			print $output "\t0";
		}
		elsif($_ eq $phased{$id}{1}){
			print $output "\t1";
		}
		else{
			print $output "\t-";
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

sub assign_phase {
	
}

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

Set minimum step size for the sliding window, defaults to 1 [Int|Default=1]

=item B<--help>

Print a brief help message and exit

=back

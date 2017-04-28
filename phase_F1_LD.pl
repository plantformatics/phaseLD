#!/usr/bin/perl
use strict;
use warnings;
use Sort::Naturally;
use Getopt::Long;
use Pod::Usage;

######################################################
##                   Arguments                      ##
######################################################

my $gen;
my $out = 'phased.out';
my $help;
my $win = 20;
my $step = $win / 4;

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
			$extend++;
			my $r2 = $results[1];
			if($r2 == 0){
				next;
			}
			else{
				my %haps = %$linkage;
				my @keys = sort {$haps{$b} <=> $haps{$a}} keys %haps;
				my $att = join(":",$id1,$id2,$r2);
				for (my $t = 0; $t < @keys; $t++){
					$hash{$att}{$keys[$t]} = $haps{$keys[$t]};
				}
			}
		}
	}
	if($extend == 0){
		$win = $win + $step;
	}
	else{
		$win = $win;
	}
}
## analyze LD information for phasing

my %phased;
my %probs;
my @keys = nsort keys %hash;
for (my $j = 0; $j < @keys; $j++){
	my @coords = split(":",$keys[$j]);
	my @haps = sort {$hash{$keys[$j]}{$b} <=> $hash{$keys[$j]}{$a}} keys $hash{$keys[$j]};
	if($j == 0){
			
	}
	else{

	}
}


close $output;
####################################################
## 		    Subroutines			  ##
####################################################
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

Set minimum window size for estimating recombination frequencies [Int|Default=20]

=item B<--step>

Set minimum step size for the sliding window, defaults to a quarter of the window size [Int|Default=5]

=item B<--help>

Print a brief help message and exit

=back

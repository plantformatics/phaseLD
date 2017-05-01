# phaseLD
&copy; Alexandre Marand & Hanain Zhao, May 2017


A simple LD based method to quickly phase pseudo-oneway test cross populations. Users should use a two step approach, first identifying the most probable haplotypes, via Bayesian methods (within phaseLD.pl), followed by HMM refinement of called haplotypes with HMM_snp.pl. 

Please cite "" if you use this software.  

## Dependencies
phaseLD.pl and HMM_snp.pl require users to have the following perl modules available to @INC:

1) [Sort::Naturally](http://search.cpan.org/~bingos/Sort-Naturally-1.03/lib/Sort/Naturally.pm)
2) [Getopt::Long](http://perldoc.perl.org/Getopt/Long.html)
3) [Pod::Usage](http://search.cpan.org/~marekr/Pod-Usage-1.69/lib/Pod/Usage.pm)
4) [List::Util](http://search.cpan.org/~pevans/Scalar-List-Utils-1.47/lib/List/Util.pm)
5) [Parallel::ForkManager](search.cpan.org/~yanick/Parallel-ForkM…)

All modules are available from CPAN, and can be accessed with cpan/cpanm.

Click [here] for link to cpanm.

[here]: http://search.cpan.org/~miyagawa/Menlo-1.9004/script/cpanm-menlo

## Usage
```
Usage:
    phase_F1_LD.pl --in [OPTIONS]

Options:
    --in    Genotype file [String|REQUIRED]

    --out   Output file name [String|Default=phased.out]

    --win   Set minimum window size for estimating LD [Int|Default=20]

    --step  Set minimum step size for LD calculations in a sliding window,
            defaults to 1 [Int|Default=1]

    --help  Print a brief help message and exit

    --bwin  Set minimum window size for bayes haplotype calling
            [Int|Default=50]

    --bstep Set minimum step size for bayes haplotype calling
            [Int|Default=5]
```

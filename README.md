# phaseLD
&copy; Alexandre Marand & Hanain Zhao, 2017


A simple LD based method to quickly phase pseudo-oneway test cross populations. Users should use a two step approach, first identifying the most probable haplotypes, via Bayesian methods (within phaseLD.pl), followed by HMM refinement of called haplotypes with HMM_snp.pl. 

Please cite "" if you use this software.  

## Dependencies
phaseLD.pl and HMM_snp.pl require users to have the following perl modules available to @INC:

1) Sort::Naturally
2) Getopt::Long
3) Pod::Usage
4) List::Util
5) Parallel::ForkManager

## Usage
```
Usage:
    phaseLD.pl --in [OPTIONS]

Options:
    --in    Genotype file [String|REQUIRED]

    --out   Output file name [String|Default=phased.out]

    --win   Set minimum window size for estimating LD
            [Int|Default=20]

    --step  Set minimum step size for the sliding window, 
	    defaults to 1 [Int|Default=1]

    --help  Print a brief help message and exit
```

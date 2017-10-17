# phaseLD
Alexandre Marand, 2017

A simple LD based method to quickly phase pseudo-oneway test cross populations by identifying window haplotypes via Bayesian inference. 

Please cite the following article if you use this software. Click [here] for a link to the paper.  

[here]: https://genomebiology.biomedcentral.com

Marand et al. (2017). "Meiotic crossovers are associated with open chromatin and enriched with *Stowaway* transposons in potato." *Genome Biology*


## Dependencies
phaseLD.pl requires users to have the following perl modules available to @INC:

1) [Sort::Naturally](http://search.cpan.org/~bingos/Sort-Naturally-1.03/lib/Sort/Naturally.pm)
2) [Pod::Usage](http://search.cpan.org/~marekr/Pod-Usage-1.69/lib/Pod/Usage.pm)
3) [List::Util](http://search.cpan.org/~pevans/Scalar-List-Utils-1.47/lib/List/Util.pm)
4) [Parallel::ForkManager](search.cpan.org/~yanick/Parallel-ForkM…)

All modules are available from CPAN, and can be accessed with cpan/cpanm.

Click [here] for a link to cpanm.

[here]: http://search.cpan.org/~miyagawa/Menlo-1.9004/script/cpanm-menlo

## Installation
After installing dependencies, simply add ```PATH=$PATH:/path/to/phaseLD/bin``` to your PATH in .bashrc or .zshrc.

## Preparing Data
phaseLD expects pseudo-oneway test cross populations, segregation ratios of 1:1 of homozygous to heterozygous genotypes. ```vcf2gen.pl```, which is included in the directory ```bin```, will convert homozygous reference genotypes (0/0) to 'a' calls, and heterozygous genotypes (0/1) to 'A'. All homozygous alternate allele calls (1/1) are converted to heterozygous genotypes, 'A'. Users can adjust the vcf2gen script to accomodate markers which segregate 1:1 for homozygous alternate alleles (1/1) and heterozygous genotypes (0/1). 

Convert VCF file to gen file format.
```
perl vcf2gen.pl file.vcf > file.gen
```
Example of gen file format:
```
chr01	1000	A	-	-	a	A	a	a	A
chr01	1100	a	-	A	A	a	A	-	-
chr01	1200	a	a	A	A	a	-	A	a
```
Sample data can be found in the directory ```examples```

Example of filter file format:
```
chr01_15879     
chr01_104598    
chr01_158589    
chr01_158620    
```

## Output
phaseLD will output 3 files, with the suffix '.bayes', '.log', '.raw_hap'. The '.bayes' file contains the probabilities of both haplotypes for each individual in a given window. The '.log' file contains all the information from the LD calculations. The '.raw_hap' file contains the non-window haplotype, with the normalized probability of the correct phase assignment. 

## Usage
```
Name:
    phaseLD v.0.01, developed by Alexandre Marand, 2017

Usage:
    phaseLD.pl --in <file.gen> [OPTIONS]

Description:
    A complete list of parameter options. See 
    'https://github.com/plantformatics/phaseLD' for more information.

  Input/Output:
    --in|-i             Genotype file [String|REQUIRED]

    --out|-o            Prefix for output files [String|Default='phased']

  Performance:
    --win|-w            Set minimum window size for estimating Pairwise LD
                        (when --quick_mode or --fast_mode are False). If
                        --quick_mode or --fast_mode are enabled, this
                        distance becomes the search area for linked markers.
                        Windows are dynamically extended if no markers
                        within the specified window are in linkage with the
                        current iteration. Markers which are not linked with
                        neighboring markers are typically false positives.
                        See prefix.log to identify markers with low linkage.
                        [Int|Default=100]

    --step|-s           Set minimum step size for LD calculations in a
                        sliding window, defaults to 1. [Int|Default=1]

    --bwin|-n           Set minimum window size for bayes haplotype calling.
                        [Int|Default=50]

    --bstep|-t          Set minimum step size for bayes haplotype calling.
                        [Int|Default=5]

    --rpen|-r           Set minimum r2 value needed to retain SNPs. This
                        threshold helps to remove false positive SNP call.
                        Defaults to 0.25. [Int|Default=0.25]

  Run Time Options:
    --threads|-p        Set number of threads. [Int|Default=1]

    --quick_mode|-f     Instead of calculating LD for every marker in the
                        window, calculate LD for window_size/10. The next
                        position is the first marker above the R2 threshold.
                        This makes the computation 10-20X faster for the
                        default window size. Time saved increases
                        exponentially as window sizes get larger.
                        [Boolean|Default=False]

    --fast_mode|-q      Rather than use pairwise LD calculations within a
                        window, use LD chaining instead. Starting at marker
                        X, the algorithm searches for the nearest marker, y,
                        with R2 above the threshold. The next iteration
                        starts at marker y and continues until all markers
                        are exhausted. [Boolean|Default=False]

    --filter|-f         Optionally use a file to mask markers. This file
                        should contain the chromosome and position of the
                        marker (example: chr01_1000) on a single line.
                        Useful to cleaning haplotype switch problems or
                        thinning data to specific markers. [String|Optional]

  Misc:
    --help|-h           Print a brief help message and exit.
```

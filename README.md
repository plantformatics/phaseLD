# phaseLD
&copy; Alexandre Marand & Hanain Zhao, 2017


A simple LD based method to quickly phase pseudo-oneway test cross populations. Users should use a two step approach, first identifying the most probable haplotypes, via Bayesian methods (within phase_F1_LD.pl), followed by HMM refinement of called haplotypes with HMM_snp.pl. 

Please cite "" if you use this software.  

## Usage
```
Usage:
    phase_F1_LD.pl --in [OPTIONS]

Options:
    --in    Genotype file [String|REQUIRED]

    --out   Output file name [String|Default=phased.out]

    --win   Set minimum window size for estimating recombination frequencies
            [Int|Default=20]

    --step  Set minimum step size for the sliding window, defaults to a
            quarter of the window size [Int|Default=5]

    --help  Print a brief help message and exit
```

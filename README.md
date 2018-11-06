# ROHan: inference of heterozygosity rates and runs of homozygosity for modern and ancient samples
==========================================================

QUESTIONS :
   gabriel [dot] reno [ at sign ] gmail [dot] com


About
----------------------

ROHan is a Bayesian framework to estimate local rates of heterozygosity, infer runs of homozygosity (ROH) and compute global rates of heterozygosity. ROHan can work on modern and ancient samples with signs of ancient DNA damage.

Description
----------------------

Diploid organisms have minor differences between corresponding pairs of autosomes. Differences of a single nucleotide create heterozygous sites. Very recent inbreeding can lead to large stretches of genomic deserts of heterozygosous sites. Such deserts are called runs of homozygosity (ROH). 

Due to the lack of heterozygosous sites, ROHs can cause an underestimate of global estimates of heterozygosity. Such global estimates of rates of heterozygous sites is useful to infer population genetics parameters such as effective population size (Ne). 

Analysis of ancient samples obtained using ancient DNA (aDNA), in addition to being plagued by potential ROHs due to inbreeding, often have presence of post-mortem molecular damage and can often have low coverage due to a paucity of endogenous material.

ROHan aims at:
1) providing an estimate of local rates of heterozygosity that is robust to low coverage and post-mortem damage. 
2) identifying ROHs using the results from 1). This is done via an hidden markov model whose parameters are optimized using a Markov Chain Monte Carlo
3) using both 1) and 2), provide a global rate of heterozygosity

see LICENSE file for license information.

Downloading:
----------------------

Go to https://github.com/grenaud/rohan and either:

1) Download the ZIP 

or

2) Do a "git clone --depth 1 https://github.com/grenaud/rohan.git"

Installation
----------------------

1) make sure you have "aclocal", "cmake", "libtool", "libpng" and "git" installed, check for it by typing " git --version" and "cmake --version".  

For Ubuntu:

     sudo apt-get install autotools-dev
     sudo apt-get install git
     sudo apt-get install cmake
     sudo apt-get install libtool
     sudo apt-get install libpng-dev

For MacOS, if you have Homebrew (https://brew.sh/) installed: 

    brew install autoconf
    brew install automake
     brew install git
     brew install cmake
     brew install libtool
     brew install libpng

2) make sure you have gcc that supports -std=c++11, gcc version 4.7 or above. Type "gcc -v" and check the version. For both Ubuntu and MacOS, 

3) As the makefile uses "git clone" to download subpackages, please make sure that the computer on which you are installing ROHan has access to the internet. Once done, simply type :
     cd ROHan
     make

For MacOS, if you get the problem: fatal error: 'lzma.h' file not found, this is a problem building htslib with homebrew, please refer to the following htslib page: https://github.com/samtools/htslib/issues/493


4) (optional) Either put the executable in the overall path or add the path to your $PATH environment or add an alias to be able to run "rohan" from any directory.


Quick start
----------------------

For test data, make sure you are connected to the internet and type:
     cd testData/
     make

This will run ROHan on chromosome 21 for: 
* a West African individual from Phase 3 of the 1000G project
* the ancient Loschbour individual (Lazardis et al.). Please be aware that the command for the profile uses no quality score filter: "-minq 0" because the quality scores were artificially decreased to account for damage during calling. We generally do not recommend this as this will inflate substitution rate. 

For these cases, the theta estimate have very large confidence due to the use of a single chromosome. 


Preparing the BAM file
-----------------

1) Make sure you know the transition/transversion ratio for the species/population you are working with. This ratio will be specified via:

    rohan --tstv [TSTV ratio] 

2) Do not apply any filters as mapping quality and base quality are all informative in the model. Duplicate removal is very recommended. Simply use the fail QC flag to remove reads/fragments that have failed basic quality control (e.g. for duplicates). Simply sort and index and provide ROHan the same reference used for mapping. 


3) Is your sample ancient DNA? if not skip to 5). If the sample has aDNA damage, it should be quantified in incorporated in the calculation. The following program can be used:

    bam2prof/bam2prof

To get the best results, take a substantial subsample of your original BAM file and use it in bam2prof. 
a) Find the ideal length for the -length parameter. Try increasing until substition rates level off.
b) Keep increasing  -minq from 0 until the damage levels off.
c) If you have substitutions outside of expected ones (e.g. C->T, G->A) consider using the -mask to filter out polymorphic positions
d) Have a look at the command line in the test data directory.

4) Do you have extensive ancient DNA damage (e.g. 20% at the ends of greater) which could potentially affect the mapping? If so, we recommend filtering for reads in highly mappable regions, see : http://lh3lh3.users.sourceforge.net/snpable.shtml  to create mappability tracks.

5) Create a file with the name of the autosomes, 1 chromosome per line ex:
    chr1
    chr2
    ...

6) Run ROHan, for modern samples run:

     rohan --rohmu 2e-5   -o  [output prefix]  [reference genome]  [BAM file]

For ancient samples run:

     rohan --rohmu 2e-5 --deam5p  deamfile.5p.prof  --deam3p  deamfile.3p.prof    -o  [output prefix]  [reference genome]  [BAM file]

where deamfile.5p.prof and deamfile.3p.prof are the substitution rates due to aDNA damage and not due to sequencing errors. You can add more threads in the calculation using "-t". Be aware that as this is a full likelihood model, the data needs to live in memory and therefore, using multiple threads can result in higher memory usage.

If a mappability track was used, you can add the option: --map to filter wrt to mappable regions.

7) Inspect your results. Refer to [output prefix].het_1_X.pdf and [output prefix].hmm_1_X.pdf, if there are regions labeled as non ROH despite the fact that they clearly show a depression in heterozygosity, re-run using --hmm and a higher value of --rohmu (e.g. for the example above, one could run using 5e-5).

8) For ancient samples, if the result if really unexpected, you can re-run using --tvonly which will limit the estimate to transversions, this needs to be multiplied by (Ts/Tv+1) (ex: if Ts/Tv  = 2.1, multiply the results it by 3.1). This estimate can be an underestimate.


Description of output files
-----------------


|                  File                |                           Contents                               | 
| ------------------------------------ |:----------------------------------------------------------------:| 
[output prefix].hEst.gz                | local estimates of heterozygosity                                | 
[output prefix].het_1_X.pdf            | Plot of the local estimates of heterozygosity 1 of X             | 
[output prefix].hmm_1_X.pdf            | Plot of the posterior decoding for HMM 1 of X                    | 
[output prefix].max.hmmmcmc.gz         | Log of the MCMC using maximum estimates of heterozygosity        | 
[output prefix].max.hmmp.gz            | HMM posterior decoding using maximum estimates of heterozygosity | 
[output prefix].mid.hmmmcmc.gz         | Log of the MCMC using mid. estimates of heterozygosity           | 
[output prefix].mid.hmmp.gz            | HMM posterior decoding using mid. estimates of heterozygosity    | 
[output prefix].min.hmmmcmc.gz         | Log of the MCMC using minimum estimates of heterozygosity        | 
[output prefix].min.hmmp.gz            | HMM posterior decoding using minimum estimates of heterozygosity | 
[output prefix].rginfo.gz              | List of read groups and average coverage                         | 
[output prefix].summary.txt            | Fraction of the genome in ROH and heterozygosity outside of ROH  | 
[output prefix].vcf.gz                 | Posterior probabilities for the genotypes                        | 




FAQ
----------------------



### Why do I have such huge confidence intervals for potential ROH regions?

This is normal, as we use the second derivative for the confidence intervals, the shape of the likelihood function is somewhat flat around 0 which leads to overestimated confidence intervals. 


### I have regions that have clear signs of being runs of homozygosity however they do not get flagged as such, why?

It is possible that you have genuine ROH but if you have a slight overestimate of the rate of local heterozygosity do to him properly calibrated weapon qualities or  base qualities. Therefore such regions will not be flagged as being ROHs. It is advisable to rerun with --rohmu XXX where XXX is the "heterozygosity rate" in ROH regions you are willing to tolerate

### What are the 3 different lines on the HMM plot?

-Green is output using the point estimates. Margenta are the estimates using upper bounds for local het. rates whereas red were computed using lower bounds of het. rates.


### What is the difference between --bed and --map?

- By default, ROHan generates windows of "--size"bp along autosomes, you can override this and use the regions in a bed file. With either option, you can specify individual sites to consider using --map, this is recommended for ancient samples with short fragment size and with some damage.



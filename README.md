inbix
=====

#### **I**nteraction-**N**etwork **BI**onformatics Toolbo**X** (inbix) ####

** This is a pre-release version. **

### Description ###

inbix -- machine learning and epistasis network analysis for high-dimensional data, including GWAS, microarray, RNA-Seq, eQTL, and others -- is a free, open-source, command-line bioinformatics tool, written in C++ and designed to perform a range of large-scale analyses with computational efficiency. The inbix program integrates many different analyses developed by the McKinney [In Silico Lab](insilico.utulsa.edu) for Bioinformatics and Computational Biology at the University of Tulsa.

inbix includes Relief-based and evaporative cooling-based algorithms for feature selection for detecting main effects and interaction effects for case-control and quantitative trait data. Inbix also allows epistasis and expression-epistasis network inference for GWAS and gene expression data; epistasis network centrality analysis; differential co-expression network analysis; interaction QTL (iQTL) network inference and differential-coexpression Variant (dcVar) analysis of eQTL data. 

inbix is built upon the PLINK project developed by Shaun Purcell at the 
Center for Human Genetic Research (CHGR), Massachusetts General Hospital (MGH), 
and the Broad Institute of Harvard & MIT, with the support of others.

For more details, visit the inbix
[website](http://insilico.utulsa.edu/index.php/inbix/)

### Dependencies ###
* [cmake](http:/cmake.org) is required to build the executable from source.

* The libz/zlib compression library is required, but this is installed by default
on most Unix systems.  In MinGW libz is installed via mingw-get.

* LAPACK is a soft dependency: not explicitly required, but highly encouraged
to take advantage of linear algebra routines that have decades of optimization.

* [Boost](http://www.boost.org/) Mostly header extensions to the standard 
C++ libraries.

* [Armadillo](http://arma.sourceforge.net/) Linear algebra header library.

* [GNU Scientific Library](https://www.gnu.org/software/gsl/) Deprecated-
to be replaced by standard c++11 and boost and removed.

* [Doxygen](http://www.stack.nl/~dimitri/doxygen/) to generate documentation 
from source code annotations.

* Finally, OpenMP is required to take advantage of the parallelized epistasis
analysis code.  This is another library typically installed alongside the 
compiler toolchain.

### Compilation Environment and Instructions ###
GNU g++ 6.3.0 successfully compiles the code using the 
cmake 3.7.2 build system.

We have successfully built and run inbix on:

* Linux (64-bit) Debian 9.2, Ubuntu 14.04

To build inbix, change to the build directory. run 'cmake ..' to indicate the
cmake configuration file is in the parent directory, the top level project 
directory. cmake builds a standard 'Makefile'. Finally, run 'make'. 
'make install" will install the build files into the base '/usr/local/'. 

The following commands if successful produce an executable file named 'inbix':

    $ cd build
    $ cmake ..
    $ make
    $ ./inbix --help

To build with documentation:

    $ cd build
    $ cmake .. -DBUILD_DOC=ON
    $ make
    $ sudo make install
    $ ./inbix --help

### Contributors ###
See AUTHORS file.

### References of Methods in inbx ###

Trang T. Le, Ryan J. Urbanowicz, Jason H. Moore, B. A McKinney. “STatistical Inference Relief (STIR) feature selection,” BioRxiv Preprint.

Trang T. Le, W. K. Simmons, M. Misaki, B.C. White, J. Savitz, J. Bodurka, and B. A. McKinney. “Differential privacy-based evaporative cooling feature selection and classification with Relief-F and Random Forests,” Bioinformatics, Volume 33, Issue 18, 15 September 2017, Pages 2906–2913.

B. Rahmani, M. Zimmermann, D. Grill, R. Kennedy, A Oberg, B. C. White, G. A. Poland, B. A. McKinney, “Recursive Indirect-Paths Modularity (RIP-M) for Detecting Community Structure in RNA-Seq Co-Expression Networks,” Frontiers in Genetics, 7:80. doi: 10.3389/fgene.2016.00080.

C. A. Lareau, B.C. White, A.L. Oberg, R.B. Kennedy, G.A. Poland, B.A. McKinney, “An interaction quantitative trait loci (iQTL) tool implicates epistatic functional variants in an apoptosis pathway in smallpox vaccine eQTL data,” Genes and Immunity (Nature Publishing). 17:244–250; doi:10.1038/gene.2016.15. 2016.

C. A. Lareau, B. C. White, Courtney G. Montgomery and B. A. McKinney, “dcVar: A Method for Identifying Common Variants that Modulate Differential Correlation Structures in Gene Expression Data,” Frontiers in Genetics. 6:312. doi: 10.3389/fgene.2015.00312. 2015.

C. Lareau, B.C. White, A.L. Oberg, B.A. McKinney, “Differential co-expression network centrality and machine learning feature selection for identifying susceptibility hubs in networks with scale-free structure,” BMC Biodata Mining 8:5. 2015. 

B.A. McKinney, B.C. White, D.E. Grill, P.W. Li, R.B. Kennedy, G.A. Poland, A.L. Oberg. “ReliefSeq: A gene-wise adaptive-k nearest-neighbor feature selection tool for finding gene-gene interactions and main effects in mRNA-Seq gene expression data,” PLoS ONE 8(12):e81527. 2013. doi:10.1371/journal.pone.0081527.

B.A. McKinney, J.E. Crowe, Jr., J. Guo, and D. Tian, “Capturing the spectrum of interaction effects in genetic association studies by simulated evaporative cooling network analysis,” PLoS Genetics. 5(3): e1000432. doi:10.1371/journal.pgen.1000432; 2009.

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) 
"PLINK: a toolset for whole-genome association and population-based 
linkage analysis."" American Journal of Human Genetics, 81.

[In Silico Lab Publications](http://insilico.utulsa.edu/index.php/publications/) 
for other algorithms that make inbix.

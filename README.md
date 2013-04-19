PLINK
========

#### Whole genome association analysis toolset ####

### Description ###
PLINK is a free, open-source whole genome association analysis toolset,
designed to perform a range of basic, large-scale analyses in a
computationally efficient manner.

The focus of PLINK is purely on analysis of genotype/phenotype data, so there
is no support for steps prior to this (e.g. study design and planning,
generating genotype or CNV calls from raw data). Through integration with
gPLINK and Haploview, there is some support for the subsequent visualization,
annotation and storage of results.

PLINK (one syllable) is being developed by Shaun Purcell at the Center for
Human Genetic Research (CHGR), Massachusetts General Hospital (MGH), and the
Broad Institute of Harvard & MIT, with the support of others.

For more details, visit the PLINK
[website](http://pngu.mgh.harvard.edu/~purcell/plink/)

The latest stable source release of PLINK (1.07) has been modified by our 
research group, [In Silico](http://insilico.utulsa.edu).  All of the 
modifications are licensed GPLv2, same as PLINK.  Our changes,
also highlighted in the ChangeLog file, are:

* Addition of a complete Autotools build system (`./configure && make &&
make install`)
* Creation of PLINK library (libplink) that exposes functionality that
may be imported by third-party software
* OpenMP support in epistasis calculation
* Numeric attribute support (currently only available in library)
* zlib is now a required library

### Dependencies ###
* The libz/zlib compression library is required, but this is installed by default
on most Unix systems.  In MinGW libz is installed via mingw-get.

* LAPACK is a soft dependency:  not explicitly required, but highly encouraged
to take advantage of linear algebra routines that have decades of optimization.

* Finally, OpenMP is required to take advantage of the parallelized epistasis
analysis code.  This is another library typically installed alongside the 
compiler toolchain.

### Compilation Environment and Instructions ###
To compile this code, a GNU toolchain and suitable environment are required.
GNU g++ has been used to successfully compile the code.

We have successfully built and run PLINK on:

* Linux (64-bit Ubuntu) (gcc-4.6)
* Mac (10.6 - 10.7) (gcc-4.2.1)
* Windows 7 (32-bit) using the [MinGW](http://www.mingw.org) compiler system
  (gcc-4.6)

To build PLINK, first run the bootstrap script

    ./bootstrap.sh

This calls  autoreconf and generates the configure script.  From this point, a 
standard

    ./configure && make && sudo make install

will generate the Makefile, compile and link the code, and copy the objects to
the installation directory (default of `/usr/local`).  As is convention, headers
are installed in `$PREFIX/include`, binary in `$PREFIX/bin`, and the library in
`$PREFIX/lib`.

The PLINK library (libplink) can be used to use PLINK functionality in
third-party code.  See our Encore project for a good example of this.

The resulting binary plink/plink.exe will run with the same command-line
options as the official binary version of PLINK that can be downloaded from 
the website.

### Contributors ###
See AUTHORS file.

### References ###
Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) 
PLINK: a toolset for whole-genome association and population-based 
linkage analysis. American Journal of Human Genetics, 81.


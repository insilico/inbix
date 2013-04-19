inbix
=====

#### Insilico Bioinformatics ####

### Description ###
inbix is a free, open-source bioinformatics tool, designed to perform a 
range of large-scale analyses in a computationally efficient manner. The inbix
program integrates many different analyses developed by the McKinney Insilico
Lab for Bioinformatics and Computational Biology at the University of Tulsa.

inbix is built upon the PLINK project developed by Shaun Purcell at the 
Center for Human Genetic Research (CHGR), Massachusetts General Hospital (MGH), 
and the Broad Institute of Harvard & MIT, with the support of others.

For more details, visit the inbix
[website](http://insilico.utulsa.edu/inbix/)

### Dependencies ###
* The libz/zlib compression library is required, but this is installed by default
on most Unix systems.  In MinGW libz is installed via mingw-get.

* LAPACK is a soft dependency:  not explicitly required, but highly encouraged
to take advantage of linear algebra routines that have decades of optimization.

* Finally, OpenMP is required to take advantage of the parallelized epistasis
analysis code.  This is another library typically installed alongside the 
compiler toolchain.

### Compilation Environment and Instructions ###
GNU g++ has been used to successfully compile the code.

We have successfully built and run inbix on:

* Linux (64-bit)

To build inbix, run make, which invokes the Makefile

    $ make

### Contributors ###
See AUTHORS file.

### References ###
Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) 
inbix: a toolset for whole-genome association and population-based 
linkage analysis. American Journal of Human Genetics, 81.

inbix
=====

#### In Silico Bioinformatics ####

** This is a pre-release version. **

### Description ###
inbix is a free, open-source bioinformatics tool, designed to perform a 
range of large-scale analyses in a computationally efficient manner. The inbix
program integrates many different analyses developed by the McKinney [In Silico
Lab](insilico.utulsa.edu) for Bioinformatics and Computational Biology at the 
University of Tulsa.

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

### References ###
Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) 
PLINK: a toolset for whole-genome association and population-based 
linkage analysis. American Journal of Human Genetics, 81.

[In Silico Lab Publications](http://insilico.utulsa.edu/index.php/publications/) 
document the algorithms added to PLINK to make inbix.

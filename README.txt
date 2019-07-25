PSRSALSA - A Suite of ALgorithms for Statistical Analysis of pulsar data

See http://www.jb.man.ac.uk/research/pulsar/Resources/psrsalsa.html for some instructions and a tutorial.

DOWNLOAD/UPDATE
----------------------------

Download the latest version from https://github.com/weltevrede/psrsalsa/ which can be done on the command line with: git clone https://github.com/weltevrede/psrsalsa.git
This will make a directory psrsalsa/

To update existing code you can run: git pull

When the code is updated, do a "make clean" and run "make" to compile the code again.


INSTALLATION
----------------------------

In principle running "make" in the directory where this README.txt
file lives should be enough. However, since PSRSALSA depends on other
libraries, the Makefile might require editing. These dependencies
include:
  - cpgplot
  - cfitsio
  - fftw3
  - gsl
See below for a description of those dependencies.

If you have a Mac OSX you might want to install homebrew (see
https://brew.sh/) and then run the following commands:
brew install gcc
brew install cfitsio
brew install fftw
brew install gsl
brew install pgplot

The entries in the Makefile that might need editing are:

* The LIBS variable.

  The following libraries are directly used by the code: libcfitsio, libcpgplot, libgsl and libfftw3f.

  However, they might depend on other libraries themselves and which
  ones is system dependent.

* The LIBDIRS variable.

  If libraries are located in a place the compiler does not look by
  default, they can be defined with this variable. By default, it is
  adding the location which is set by the PGPLOT_DIR environment
  variable.

* The INCDIRS variable.

  If header files are located in a place the compiler does not look by
  default. Again, by default it is adding the location which is set by
  the PGPLOT_DIR environment variable.

* GSLFLAGS variable

  The library version of the GSL library you use can be specified with
  the -DGSL_VERSION_NUMBER=XXX. Changing this number enables/disables
  different pieces of the code, such that the code is compilable
  agains a range of versions. The functionality of the code might be
  restricted depending on the version of GSL you're compiling the code
  against.


DEPENDENCIES
----------------------------

* cfitsio:

Website = http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html

Probably your linux distribution allows you to download this package.

Version 3.21 is known to work.

* fftw3:

Website = http://www.fftw.org/

Probably your linux distribution allows you to download this
package. Make sure to compile/obtain the float version of the library
as well: libfftw3f.

Version 3.2.1 is known to work.

* gsl:

Website = http://www.gnu.org/software/gsl/

Probably your linux distribution allows you to download this
package. You need libgsl and libgslcblas.

Version 1.15 is known to work

* pgplot:

Website = http://www.astro.caltech.edu/~tjp/pgplot/

Probably your linux distribution allows you to download this
package. This library probably depend on other libraries, such as
libgfortran, libX11 and libpng.

Version 5.2.2 is known to work

see README_pgplot.txt how to compile pgplot.


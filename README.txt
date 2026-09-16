PSRSALSA - A Suite of ALgorithms for Statistical Analysis of pulsar data

See http://www.jb.man.ac.uk/research/pulsar/Resources/psrsalsa.html for some instructions and a tutorial.

Further information about how to use the psrsalsa python wrapper can
be found in pywrapper/manual.pdf.

DOWNLOAD/UPDATE
----------------------------

Download the latest version from https://github.com/weltevrede/psrsalsa/ which can be done on the command line with: git clone https://github.com/weltevrede/psrsalsa.git
This will make a directory psrsalsa/

To update existing code you can run: git pull

INSTALLATION
----------------------------

In principle running "make" in the directory where this README.txt
file lives should be enough. However, since PSRSALSA depends on other
libraries (see the section dependancies below), compilation might well
require some further actions.

When the code is updated, do a "make clean" and run "make" to compile
the code again. should compile everything. But see the below how to
solve problems during compilation because of missing dependencies, or
dependencies in non-standard locations.

There is a python wrapper for psrsalsa, which allows you to run
psrsalsa operations from python. To avoid attempting to build the
python wrapper, use "make nopy" rather than "make".

See the next section how solve problems during compilation when
dependancies are not found.

ENVIRONMENT VARIABLES AFFECTING COMPILATION
-------------------------------------------

Not much should need to be done if your system has the required
libraries installed via a package manager (see the following
dependencies section). However, it is common for at least the location
of pgplot to be not automatically be recognized, as this needs manual
installation (see README_pgplot.txt).

If you have dependencies in locations that are not automatically
found, you might want to set one or mor eof the following environment
variables. These environment variables only affect compilation, so it
is not necessary to set these permanently.

PGPLOT_DIR: This is the default directory that is used to find the
pgplot libraries and include files. If this is not the case, one can
use the PSRSALSA_EXTRA_INCDIRS and PSRSALSA_EXTRA_LIBDIRS
environment variables to point the Makefile where these can be
found.

PSRSALSA_EXTRA_LIBDIRS: These are -L options for the compiler
indicating where libraries can be found that are otherwise not found
automatically. Example in bash: 
export PSRSALSA_EXTRA_LIBDIRS=-L/packages/pulsarsoft/lib/ -L/morepackages/lib/
So for example, when pgplot can't be located during compilation,
make sure to define the directory where libpgplot.a and/or
libpgplot.so, as well as libcpgplot.a and/or libcpgplot.so can be
found.

PSRSALSA_EXTRA_INCDIRS: These are -I options for the compiler
indicating where include files can be found that are otherwise not
found automatically. Example in bash:
export PSRSALSA_EXTRA_INCDIRS=-I/packages/pulsarsoft/include/ -I/morepackages/include/
So for example, when pgplot can't be located during compilation,
make sure to define the directory where cpgplot.h can be found.

PSRSALSA_GSL_VERSION: If unset, the gsl version is automatically
determined using pkg-config. So in general this environment variable
is not needed. To avoid this, one can set PSRSALSA_GSL_VERSION. The
psrsalsa code need to know the GSL version number, as things the
syntax for some things depend on the version number. Here version
115 would correspond to 1.15, and 203 to 2.3, etc. So version =
(100*major version + minor version). Example in bash: export
PSRSALSA_GSL_VERSION=203


DEPENDENCIES
----------------------------

The following libraries are used by psrsalsa:
  cfitsio
  cpgplot
  gsl
  fftw3
  fftw3f

These libraries in turn dependent on further libraries. Which ones
is system dependent, but probably at least includes the following
dependencies because of cpgplot:

 - libpgplot   (fortran library which is called by the cpgplot wrapper)
 - libgfortran (if gfortran was used to compile pgplot)
 - libX11      (if used when compiling pgplot)
 - libpng      (if used when compiling pgplot)

If you have pkg-config you can check if the required dependencies
are installed by running the following commands (it should give some
compiler command-line options, rather than an error that the package
is not found). Note that pgplot always needs manual installation.

  pkg-config --libs cfitsio
  pkg-config --libs gsl
  pkg-config --libs fftw3
  pkg-config --libs fftw3f

If you have Ubuntu, and if you have missing packages, you may wish
to run one or more of the following commands:

  sudo apt install git
  sudo apt install gcc
  sudo apt install gfortran
  sudo apt install make
  sudo apt install libcfitsio-dev
  sudo apt install libgsl-dev
  sudo apt install libfftw3-dev
  
To compile the python wrapper, there are further dependencies. On
Ubuntu, the following commands are at least required to cover all
dependencies.
                                                                 
  sudo apt install python3                   
  sudo apt install python3-dev               
  sudo apt install python3-setuptools        
  sudo apt install python3-numpy             

More information about the dependencies can be found below.

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


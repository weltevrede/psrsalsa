# On MAC OS X:
#  - Expected the following packages to be installed:
#    brew install gcc
#    brew install cfitsio
#    brew install fftw
#    brew install gsl
#    brew install pgplot


#Define directories where include files can be found.
INCDIRS = -Isrc/lib/ -Isrc/slalib/

#Define directories where libraries can be found.
LIBDIRS = -Lsrc/lib/ -Lsrc/slalib 

# Detect Linux vs Mac OS X and set some variables different depending
# on OS. These are system dependent settings, so you may have to
# change this to get it to work. In particular, it defines the
# libraries on which psrsalsa is dependent, which are cfitsio,
# cpgplot, gsl and fftw3f. Furthermore, libcpgplot requires:
# - libpgplot
# - libgfortran (if used gfortran to compile pgplot)
# - libX11 (if used when compiling pgplot)
# - libpng  (if used when compiling pgplot)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	OSFLAG = LINUX
#Define the libraries used by the code in a linux environment
	LIBS = -lm -lcfitsio -lcpgplot -lpgplot -lpng -lX11 -lgsl -lgslcblas -lfftw3f -lgfortran
#Define directories where libraries can be found if installed in unusual places.
	LIBDIRS += 
#Define some linux-specific directories where some required include
#files can be found.  You can probably use something like: locate
#fitsio.h to find the location of for instance the fitsio include
#files if not found automatically by the compiler. The following was
#required on the linux system of the author.
	INCDIRS += -I/local/scratch/wltvrede/puma1soft/trunk/src/Soft/cfitsio/include/
#You need to define the GSL version you're using during
#compilation. Version 1.15 would be 115, 2.3 would be 203,
#etc. (100*major version + minor version).
	GSLFLAGS = -DGSL_VERSION_NUMBER=115
else
ifeq ($(UNAME_S),Darwin)
	OSFLAG = OSX
#Define the libraries used by the code in an OS X environment
	LIBS = -lm -lcfitsio -lcpgplot -lgsl -lfftw3f
#Define directories where libraries can be found if installed in
#unusual places. When installing missing packages with brew, it should
#not be necessary to add things here.
	LIBDIRS += 
#Define some OS X-specific directories where some required include
#files can be found. When installing missing packages with brew, it
#should not be necessary to add things here.
	INCDIRS += 
#You need to define the GSL version you're using during
#compilation. Version 1.15 would be 115, 2.3 would be 203,
#etc. (100*major version + minor version).
	GSLFLAGS = -DGSL_VERSION_NUMBER=203
endif
endif


#If the PGPLOT_DIR environment variable is defined, add it to the path
#where libraries and include files can be found.
ifdef PGPLOT_DIR
	LIBDIRS += -L$(PGPLOT_DIR) -L/local/scratch/wltvrede/puma1soft/trunk/src/Soft/cfitsio/lib/
	INCDIRS += -I$(PGPLOT_DIR)
endif

# Define the C compiler to be used
CC = gcc

# Some flags to pass on to C compiler
# Could do -O3 option to make code faster
CFLAGS = -Wall -g


# Define the fortran compiler to be used (for slalib which is included in the source directory)
F77 = gfortran

# Some flags to pass on to the fortran compiler (for slalib)
FFLAGS = -fno-underscoring -O

#####################################################################
# After this point, it should not be necessary to edit the Makefile #
#####################################################################



#The slalib wrapper library to be generated
SLALIBTARGET = src/slalib/libsla_wrap.a 

#The source files used to make the library, i.e. all .f and .c files in src/slalib/
SLALIBSRC = $(wildcard src/slalib/*.f) $(wildcard src/slalib/*.c)

#These are the object files to be generated from SLALIBSRC
SLALIBOBJ = $(SLALIBSRC:.f=.o) $(SLALIBSRC:.c=.o)



#The name of the library to be generated
LIBTARGET = src/lib/libpsrsalsa.a

#The source files used to make the library, i.e. all .c files in src/lib/
PSRSALSALIBSRC = $(wildcard src/lib/*.c)

#These are the object files to be generated from PSRSALSALIBSRC
PSRSALSALIBOBJ = $(PSRSALSALIBSRC:.c=.o)



#Make a list of executables to be generated.
#Take the source code files in the prog directory, strip the .c extensions and replace src/prog with bin
#This means this variable should be something like:
#EXECUTABLES = bin/pspec bin/pspecFig ....
EXECUTABLESSRC = $(wildcard src/prog/*.c)
EXECUTABLES_TMP = $(subst .c,,$(EXECUTABLESSRC))
EXECUTABLES = $(subst src/prog,bin,$(EXECUTABLES_TMP))


#Below are the rules how to do the compilation

#These are the primary targets: the name of the libraries and executables to be generated
all: startmessage slalibmessage $(SLALIBTARGET) psrsalsalibmessage $(LIBTARGET) psrsalsaprogmessage $(EXECUTABLES)
	@echo ""
	@echo "Finished Makefile. If compiled without errors, the executables should be in the directory bin/"
	@echo ""

startmessage:
ifndef OSFLAG
	@echo "Cannot detect operating system. It might not be Linux or Mac OS X"
	@echo "Compilation process aborted"
	@echo ""
	@false
else
	@echo ""
	@echo "Detected operating system $(OSFLAG)"
	@echo ""
endif

slalibmessage:
	@echo ""
	@echo "======================================================================="
	@echo "Starting compiling slalib. This requires the compilers $(F77) and $(CC)"
ifeq ($(OSFLAG),OSX)
	@echo "On a Mac this can be installed together with gcc using"
	@echo "brew install gcc"
	@echo "See internet how to install the homebrew package manager"
endif
	@echo "======================================================================="
	@echo ""

#This is the rule of how to make the objects to go into the slalib library
src/slalib/%.o:src/slalib/%.f
	$(F77) $(FFLAGS) -c -o $@ $<

#This is the rule of how to make the slalib wrapper object
src/slalib/%.o:src/slalib/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

psrsalsalibmessage:
	@echo ""
	@echo "=================================================================================================="
	@echo "Starting compiling the psrsalsa library. This requires the cfitsio, fftw, gsl and pgplot libraries"
	@echo "You might get errors related to include files not found, which either means that the related"
	@echo "isn't installed, or that the INCDIRS variable in Makefile needs to be adjusted to let the "
	@echo "compiler know where they can be found"
ifeq ($(OSFLAG),OSX)
	@echo "On a Mac this can be installed using"
	@echo "brew install cfitsio"
	@echo "brew install fftw"
	@echo "brew install gsl"
	@echo "brew install pgplot"
	@echo "See internet how to install the homebrew package manager"
endif
	@echo ""
	@echo "gsl related errors can occur if GSL_VERSION_NUMBER is not set correctly in the Makefile"
	@echo "=================================================================================================="
	@echo ""

#This is the rule of how to make the objects to go into the library
src/lib/%.o:src/lib/%.c $(SLALIBTARGET)
	$(CC) $(INCDIRS) $(CFLAGS) $(GSLFLAGS) -c -o $@ $<

#This is the rule of how to make the slalib library from the object files
$(SLALIBTARGET): $(SLALIBOBJ)
	ar rcs $(SLALIBTARGET) $(SLALIBOBJ)

#This is the rule of how to make the library from the object files
$(LIBTARGET): $(PSRSALSALIBOBJ) $(SLALIBTARGET)
	ar rcs $(LIBTARGET) $(PSRSALSALIBOBJ)

psrsalsaprogmessage:
	@echo ""
	@echo "=================================================================================================="
	@echo "Starting compiling the executables, which require the cfitsio, fftw, gsl and pgplot libraries"
	@echo "You might get errors related to libraries not found, which either means that the related"
	@echo "isn't installed, or that the LIBDIRS variable in Makefile needs to be adjusted to let the "
	@echo "compiler know where they can be found"
ifeq ($(OSFLAG),OSX)
	@echo "On a Mac this can be installed using"
	@echo "brew install cfitsio"
	@echo "brew install fftw"
	@echo "brew install gsl"
	@echo "brew install pgplot"
	@echo "See internet how to install the homebrew package manager"
endif
	@echo ""
	@echo "gsl related errors can occur if GSL_VERSION_NUMBER is not set correctly in the Makefile"
	@echo "=================================================================================================="
	@echo ""

bin/%: src/prog/%.c $(LIBTARGET) $(SLALIBTARGET)
	$(CC) $(INCDIRS) $(CFLAGS) $(GSLFLAGS) $(LIBDIRS) $< -lpsrsalsa -lsla_wrap $(LIBS) -o $@

#This is the rule of how to clean up things, so everything can be compiled from scratch
clean:
	rm -rf src/slalib/*.o src/lib/*.o $(SLALIBTARGET) $(LIBTARGET) $(EXECUTABLES)

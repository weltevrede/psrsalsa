# Define the C compiler to be used
CC = gcc

# Some flags to pass on to C compiler
# Could do -O3 option to make code faster
CFLAGS = -Wall -g

# Define the fortran compiler to be used (for slalib which is included in the source directory)
F77 = gfortran

# Some flags to pass on to the fortran compiler (for slalib)
FFLAGS = -fno-underscoring -O

#Define libraries on which psrsalsa is dependent, which are cfitsio, cpgplot, gsl and fftw3f.
#
#Furthermore, libcpgplot requires:
# - libpgplot
# - libgfortran (if used gfortran to compile pgplot)
# - libX11 (if used when compiling pgplot)
# - libpng  (if used when compiling pgplot)
LIBS = -lm -lcfitsio -lcpgplot -lpgplot -lpng -lX11 -lgsl -lgslcblas -lfftw3f -lgfortran

#Define directories where libraries can be found.
LIBDIRS = -L$(PGPLOT_DIR) -L/local/scratch/wltvrede/puma1soft/trunk/src/Soft/cfitsio/lib/

#Define directories where include files can be found.
#You can probably use something like: locate fitsio.h
#to find the location of for instance the fitsio include files
INCDIRS = -Isrc/lib/ -I$(PGPLOT_DIR) -I/local/scratch/wltvrede/puma1soft/trunk/src/Soft/cfitsio/include/

#You may want to define the GSL version you're using during compilation. Version 1.15 would be 115 etc (100*major version + minor version).
GSLFLAGS = -DGSL_VERSION_NUMBER=115

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
all: $(SLALIBTARGET) $(LIBTARGET) $(EXECUTABLES)
	$(SHELL) message.csh

#This is the rule of how to make the objects to go into the slalib library
src/slalib/%.o:src/slalib/%.f
	$(F77) $(FFLAGS) -c -o $@ $<

#This is the rule of how to make the slalib wrapper object
src/slalib/%.o:src/slalib/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

#This is the rule of how to make the objects to go into the library
src/lib/%.o:src/lib/%.c $(SLALIBTARGET)
	$(CC) $(INCDIRS) -I src/slalib/ $(CFLAGS) $(GSLFLAGS) -c -o $@ $<

#This is the rule of how to make the slalib library from the object files
$(SLALIBTARGET): $(SLALIBOBJ)
	ar rcs $(SLALIBTARGET) $(SLALIBOBJ)

#This is the rule of how to make the library from the object files
$(LIBTARGET): $(PSRSALSALIBOBJ) $(SLALIBTARGET)
	ar rcs $(LIBTARGET) $(PSRSALSALIBOBJ)

bin/%: src/prog/%.c $(LIBTARGET) $(SLALIBTARGET)
	$(CC) $(INCDIRS) -I src/lib/ $(CFLAGS) $(GSLFLAGS) $(LIBDIRS) -L src/lib/ -L src/slalib $< -lpsrsalsa -lsla_wrap $(LIBS) -o $@

#This is the rule of how to clean up things, so everything can be compiled from scratch
clean:
	rm -rf src/slalib/*.o src/lib/*.o $(SLALIBTARGET) $(LIBTARGET) $(EXECUTABLES)

# NOTE: Makefile starts with some explanation of ENVIRONMENT VARIABLES
# and DEPENDENCIES, followed by a relatively short section that might
# need edditing to account for the SPECIFICT OF THE CONFIGURATION of
# your OS, followed by a final (largest) bit that should in general be
# left alone.

#####################################################################
# Environment variables that affect the Makefile                    #
#####################################################################

# Not much needs to be done if your system is set-up such that
# pkg-config recognizes the following packages:
#
# gsl
# cfitsio
# fftw3
# fftw3f
#
# If you have dependencies in locations that are not automatically
# found, you might want to set the following environment variables:
#
# PGPLOT_DIR: This is the directory where (hopefully) both the pgplot
# libraries and include files can be found. If not, one can use the
# PSRSALSA_EXTRA_INCDIRS and PSRSALSA_EXTRA_LIBDIRS environment
# variables to point to where these can be found.
#
# PSRSALSA_EXTRA_LIBDIRS: These are -L options for the compiler
# indicating where libraries can be found that are otherwise not
# found automatically. Example in bash:
# export PSRSALSA_EXTRA_LIBDIRS=-L/packages/pulsarsoft/lib/ -L/morepackages/lib/
#
# PSRSALSA_EXTRA_INCDIRS: These are -I options for the compiler
# indicating where include files can be found that are otherwise not
# found automatically. Example in bash:
# export PSRSALSA_EXTRA_INCDIRS=-I/packages/pulsarsoft/include/ -I/morepackages/include/
#
# PSRSALSA_GSL_VERSION: If not set, the gsl version is automatically
# determined using pkg-config. To avoid this, one can set
# PSRSALSA_GSL_VERSION. The psrsalsa code need to know the GSL version number, as things the
# syntax for some things depend on the version number. Here version 115 would correspond to
# 1.15, and 203 to 2.3, etc. So version = (100*major version + minor
# version). Example in bash:
# export PSRSALSA_GSL_VERSION=203

#####################################################################
# Dependencies                                                      #
#####################################################################

# The following libraries are used by psrsalsa:
#  - cfitsio,
#  - cpgplot
#  - gsl
#  - fftw3f

# These libraries in turn dependent on further libraries. Which ones
# is system dependent, but probably at least includes the following
# dependencies of cpgplot:
#
#  - libpgplot   (fortran library which is called by the cpgplot wrapper)
#  - libgfortran (if gfortran was used to compile pgplot)
#  - libX11      (if used when compiling pgplot)
#  - libpng      (if used when compiling pgplot)

#####################################################################
# You might need to edit the following section                      #
#####################################################################

# This defines a list of packages that are recognized by pkg-config,
# which then is used to obtain relevant compiler flags. If one of
# these causes errors, it either means you don't have the package
# installed, or it is a custom installation that will required the
# environment variables (see above) to be set to get the compilation
# to work.
PKG_CONFIG_LIST = gsl cfitsio fftw3 fftw3f

# Define the C compiler to be used
CC = gcc

# Define flags to pass on to C compiler.
# -Wall    - generate more warnings
# -g       - enable debug information
# -O3      - could enable it to make execution faster, but it will make debugging difficult
CFLAGS = -Wall -g -fPIC
#CFLAGS += -O3

# Define the fortran compiler to be used (for slalib which is included in the source directory)
F77 = gfortran

# Some flags to pass on to the fortran compiler (for slalib)
FFLAGS = -fno-underscoring -O -fPIC

# Define other libraries that are needed, which is at least pgplot.
USERLIBS = -lcpgplot -lpgplot

# Detect Linux vs Mac OS X, and set can set some variables differently
# depending on the OS.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	OSFLAG = LINUX
# Your pgplot probably at least depends on libgfortran
	USERLIBS += -lgfortran
# It might also depend on libraries such as:
	USERLIBS += -lpng -lX11
# If you have a custom installation of cfitsio:
#	USERLIBS += -lcfitsio
# cfitsio might also need the following:
#       USERLIBS += -lcurl -lz
# If you have a custom installation of gsl:
#	USERLIBS += -lgsl -lgslcblas
# If you have a custom installation of fftw:
#	USERLIBS += -lfftw3 -lfftw3f
else
ifeq ($(UNAME_S),Darwin)
	OSFLAG = OSX
#Define the libraries used by the code in an OS X environment
#	USERLIBS += 
endif
endif



#####################################################################
# After this point, it should not be necessary to edit the Makefile #
#####################################################################

# Libraries on which psrsalsa depends. There are defined via the
# earlier defined PKG_CONFIG_LIST and USERLIBS, as well as the m
# (maths) library used directly by psrsalsa.
LIBS = `pkg-config --libs $(PKG_CONFIG_LIST)` $(USERLIBS) -lm

# Define the psrsalsa directories where the psrsalsa library source and include files can be found
LIBDIRS = -Lsrc/lib/ -Lsrc/slalib/
INCDIRS = -Isrc/lib/ -Isrc/slalib/


#Define directories where libraries can be found, and initialise it
#with user defined -L compiler commands if they are defined in the environment.
ifdef PSRSALSA_EXTRA_LIBDIRS
	LIBDIRS += $(PSRSALSA_EXTRA_LIBDIRS)
endif

#Define directories where include files can be found, and initialise
#it with user defined -I compiler commands if they are defined in the
#environment.
ifdef PSRSALSA_EXTRA_INCDIRS
	INCDIRS += $(PSRSALSA_EXTRA_INCDIRS)
endif

#If the PGPLOT_DIR environment variable is defined, add it to the path
#as it (hopefully) point towards the pgplot libraries and include
#files can be found. If not, one can use the PSRSALSA_EXTRA_INCDIRS
#and PSRSALSA_EXTRA_LIBDIRS environment variables to point to where
#these can be found.
ifdef PGPLOT_DIR
	LIBDIRS += -L$(PGPLOT_DIR)
	INCDIRS += -I$(PGPLOT_DIR)
endif

# These are the header file dependencies of psrsalsa
LIB_INCLUDE_FILES = src/lib/psrsalsa_defines.h src/lib/psrsalsa.h src/lib/psrsalsa_typedefs.h

# Define a double space for indentation of messages
space := $(empty)  $(empty)
empty :=

$(info The following pkg-config packages are defined in Makefile: $(PKG_CONFIG_LIST))
PKG_CONFIG_LIBS := $(shell pkg-config --libs $(PKG_CONFIG_LIST))
$(info $(space)pkg-config -libs: $(PKG_CONFIG_LIBS))
PKG_CONFIG_CFLAGS := $(shell pkg-config --cflags $(PKG_CONFIG_LIST))
$(info $(space)pkg-config -cflags: $(PKG_CONFIG_CFLAGS))

# Add the pkg-config cflags 
CFLAGS += `pkg-config --cflags $(PKG_CONFIG_LIST)`   

#The psrsalsa code need to know the GSL version number, as things the
#syntax for some things depend on the version number. This is passed
#on to psrsalsa code (using the GSLFLAGS variable) with the option
#-DGSL_VERSION_NUMBER=version. Here version 115 would correspond to
#1.15, and 203 to 2.3, etc. So version = (100*major version + minor
#version). If the PSRSALSA_GSL_VERSION environment variable is set,
#use that. Otherwise, use pkg-config to try to automatically set this
#flag.

# Detect whether PSRSALSA_GSL_VERSION is set in the environment or via make command line
ifeq ($(origin PSRSALSA_GSL_VERSION), undefined)

$(info )
$(info PSRSALSA_GSL_VERSION environment variable not set)
$(info $(space)Try to determine appropriate version using pkg-config)
# Obtain GSL version string "major.minor.patch"
gsl_version_str := $(shell pkg-config --modversion gsl)
$(info $(space)GSL version according to pkg-config: $(gsl_version_str))
# Extract components
gsl_major := $(word 1, $(subst ., ,$(gsl_version_str)))
gsl_minor := $(word 2, $(subst ., ,$(gsl_version_str)))
gsl_patch := $(or $(word 3, $(subst ., ,$(gsl_version_str))),0)
# Compute numeric version = 100*major + minor
PSRSALSA_GSL_VERSION := $(shell echo $$(( $(gsl_major) * 100 + $(gsl_minor) )))
$(info $(space)100*major + minor = $(PSRSALSA_GSL_VERSION))

else

$(info Using PSRSALSA_GSL_VERSION from environment: $(PSRSALSA_GSL_VERSION))

endif

# Create GSLFLAGS
GSLFLAGS := -DGSL_VERSION_NUMBER=$(PSRSALSA_GSL_VERSION)
#$(info $(space)GSL version number: $(PSRSALSA_GSL_VERSION))
$(info $(space)Flags: $(GSLFLAGS))


#The slalib wrapper library to be generated
SLALIBTARGET = src/slalib/libsla_wrap.a 

#The source files used to make the library, i.e. all .f and .c files in src/slalib/
SLALIBSRC = $(wildcard src/slalib/*.f) $(wildcard src/slalib/*.c)

#These are the object files to be generated from SLALIBSRC
SLALIBOBJ = $(SLALIBSRC:.f=.o) $(SLALIBSRC:.c=.o)









#The name of the library to be generated
LIBTARGET = src/lib/libpsrsalsa_release_make.a

#The source files used to make the library, i.e. all .c files in src/lib/
PSRSALSALIBSRC = $(wildcard src/lib/*.c)

#These are the object files to be generated from PSRSALSALIBSRC
PSRSALSALIBOBJ = $(PSRSALSALIBSRC:.c=.o)


PSRSALSALIBS = -lpsrsalsa_release_make -lsla_wrap


#Make a list of executables to be generated.
#Take the source code files in the prog directory, strip the .c extensions and replace src/prog with bin
#This means this variable should be something like:
#EXECUTABLES = bin/pspec bin/pspecFig ....
EXECUTABLESSRC = $(wildcard src/prog/*.c)
EXECUTABLES_TMP = $(subst .c,,$(EXECUTABLESSRC))
EXECUTABLES = $(subst src/prog,bin,$(EXECUTABLES_TMP))


define psrsalsa_extra_libdirs_warning
@if [ -z "$(PSRSALSA_EXTRA_LIBDIRS)" ]; then \
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"; \
	echo "WARNING: You might want to set the PSRSALSA_EXTRA_LIBDIRS environment variable to"; \
	echo "         define places where libraries can be found that are otherwise not found"; \
	echo "         automatically. Example in bash:"; \
	echo "         export PSRSALSA_EXTRA_LIBDIRS=-L$(HOME)/pulsarsoft/lib/"; \
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"; \
fi
@if [ -z "$(PSRSALSA_EXTRA_INCDIRS)" ]; then \
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"; \
	echo "WARNING: You might want to set the PSRSALSA_EXTRA_INCDIRS environment variable to"; \
	echo "         define places where include files can be found that are otherwise not found"; \
	echo "         automatically. Example in bash:"; \
	echo "         export PSRSALSA_EXTRA_INCDIRS=-I$(HOME)/pulsarsoft/include/"; \
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"; \
fi
endef


#Below are the rules how to do the compilation

EXTERNAL_LIB_TARGETS = $(SLALIBTARGET)


#These are the primary targets: the name of the libraries and executables to be generated
all: startmessage slalibmessage $(EXTERNAL_LIB_TARGETS) psrsalsalibmessage $(LIBTARGET) psrsalsaprogmessage $(EXECUTABLES)
	@echo "" 
	@echo "Finished Makefile. If compiled without errors, the executables should be in the directory bin/"
	@echo "In your .bashrc file you may wish to add the lines:"
	@echo ""
	@echo "export PSRSALSAPATH=`pwd`"
	@echo "export PATH=$$""PSRSALSAPATH/bin/:""$$""PATH"



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
#	$(call psrsalsa_extra_libdirs_warning)


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
#	$(call psrsalsa_extra_libdirs_warning)

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
#	$(call psrsalsa_extra_libdirs_warning)

#This is the rule of how to make the objects to go into the library
src/lib/%.o:src/lib/%.c $(EXTERNAL_LIB_TARGETS) $(LIB_INCLUDE_FILES)
	$(CC) $(INCDIRS) $(CFLAGS) $(GSLFLAGS) -c -o $@ $<

#This is the rule of how to make the slalib library from the object files
$(SLALIBTARGET): $(SLALIBOBJ)
	ar rcs $(SLALIBTARGET) $(SLALIBOBJ)


#This is the rule of how to make the library from the object files
$(LIBTARGET): $(PSRSALSALIBOBJ) $(EXTERNAL_LIB_TARGETS)
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
#	$(call psrsalsa_extra_libdirs_warning)

bin/%: src/prog/%.c $(LIBTARGET) $(EXTERNAL_LIB_TARGETS)
	$(CC) $(INCDIRS) $(CFLAGS) $(GSLFLAGS) $(LIBDIRS) $< $(PSRSALSALIBS) $(LIBS) -o $@


#This is the rule of how to clean up things, so everything can be compiled from scratch
clean:
	rm -f src/slalib/*.o src/lib/*.o $(EXTERNAL_LIB_TARGETS) $(LIBTARGET) $(EXECUTABLES) \




#include "revision_autogenerated.psrsalsainfo"


# On MAC OS X:
#  - Expected the following packages to be installed:
#    brew install gcc
#    brew install cfitsio
#    brew install fftw
#    brew install gsl
#    brew install pgplot

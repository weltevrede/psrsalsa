Compilation of PGPLOT
-------------------------------------

Download source code from: ftp://ftp.astro.caltech.edu/pub/pgplot/pgplot5.2.tar.gz

Extract (pgplot subdirectory will be made): tar -xvf pgplot5.2.tar.gz

cd pgplot/sys_linux/
cp g77_gcc_aout.conf custom.conf

Edit custom.conf:
  Change FCOMPL="g77" -> FCOMPL="gfortran"

Go back to pgplot directory: cd ../

Edit drivers.list if you want to include for instance the following devices:
  - /PNG /TPNG /PS /VPS /CPS /VCPS /XSERVE

Generate makefile: ./makemake . linux custom

Compile: make  (see below for solutions for possible errors)

When successful compiled pgplot, compile the cpgplot wrapper library:
   make cpg

You make have to do the following as well:
   make grfont.dat

The following files should be there:
 libpgplot.a
 libcpgplot.a
 cpgplot.h

Check if PGPLOT_DIR and PGPLOT_FONT are already set to something in your login scripts (if it is pointing to a different location than where you installed your own version it probably should probably still be fine). If this is the only installation of pgplot add the following lines to your .bashrc file:
    export PGPLOT_DIR=YOUR_INSTALL_DIRECTORY
    export PGPLOT_FONT=YOUR_INSTALL_DIRECTORY/grfont.dat

Therefore if you run:
  ls $PGPLOT_DIR
you should see something like:
cpgdemo
grfont.dat
pgxwin_server
rgb.txt

If you run:
  ls $PGPLOT_FONT
you should see the file grfont.dat

If you do not get labels/numbers in pgplot graphs after compilation of
your software, copy the grfont.dat from a location where PGPLOT does
work to your $PGPLOT_FONT location.

When compiling psrsalsa, you probably need to set the environment
variables PSRSALSA_EXTRA_LIBDIRS and PSRSALSA_EXTRA_INCDIRS to make it
find your pgplot installation. See README.txt for more info.


Compilation errors
--------------------------------------

# Problem with png support:
#
make: *** No rule to make target `png.h', needed by `pndriv.o'.  Stop.

Solution:
You could run: perl -pi -e 's/^pndriv\.o :/# pndriv\.o :/' makefile
  (this comments out "pndriv.o : ./png.h ./pngconf.h ./zlib.h ./zconf.h" in the makefile which might solve the problem)
  Try to run make again:
     - make clean
     - make


# Problem with function definitions. If you get errors like:
cpgwnad.c:6:3: error: too many arguments to function pgwnad_; expected 0, have 4
    6 |   pgwnad_(&x1, &x2, &y1, &y2);
      |   ^~~~~~~ ~~~
cpgwnad.c:2:13: note: declared here
    2 | extern void pgwnad_();
      |             ^~~~~~~
  
Solution:
Edit makefile and add  -std=gnu11 to the CFLAGC definition:
CFLAGC=-Wall -fPIC -DPG_PPU -g -O0 -I. -fPIC -std=gnu11


Conflicting pgplot installations
--------------------------------------

The following is in genral not required. But if you want to make sure
that your pgplot compilation is picked up by the compiler when compiling
psrsalsa, you can always rename the pgplot libraries. For example

mv libpgplot.a libpgplot_psrsalsa.a
mv libcpgplot.a libcpgplot_psrsalsa.a

Now in the Makefile of you can edit the USERLIBS variable to make use of these psrsalsa specific libraries:

USERLIBS = -lcpgplot_psrsalsa -lpgplot_psrsalsa





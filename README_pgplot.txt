To compile pgplot:

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

The following files should be there:
 libpgplot.a
 libcpgplot.a
 cpgplot.h

If you want to make sure that this library is picked up by the compiler when compiling psrsalsa, you can always rename these libraries. For example

mv libpgplot.a libpgplot_psrsalsa.a
mv libcpgplot.a libcpgplot_psrsalsa.a

Now in the Makefile of psrsalsa make sure to link with those libraries rather than the default file names:

Change: -lcpgplot -lpgplot   => -lcpgplot_psrsalsa -lpgplot_psrsalsa
Change: -L $(PGPLOT_DIR)     => -L YOUR_INSTALL_DIRECTORY
Change: -I $(PGPLOT_DIR)     => -I YOUR_INSTALL_DIRECTORY

Check if PGPLOT_DIR and PGPLOT_FONT are already set to something in login scripts (if it is pointing to a different location than where you installed your own version it probably should be fine). If this is the only installation of pgplot add the following lines to your .cshrc file (assuming you're running tcsh):
    setenv PGPLOT_DIR YOUR_INSTALL_DIRECTORY
    setenv PGPLOT_FONT YOUR_INSTALL_DIRECTORY/grfont.dat

Possible errors during compilation:

make: *** No rule to make target `png.h', needed by `pndriv.o'.  Stop.
  You could run: perl -pi -e 's/^pndriv\.o :/# pndriv\.o :/' makefile
  (this comments out "pndriv.o : ./png.h ./pngconf.h ./zlib.h ./zconf.h" in the makefile which might solve the problem)
  Try to run make again:
     - make clean
     - make

Possible errors during runtime:

If you do not get labels/numbers in pgplot after compilation of your software, copy the grfont.dat from a location where PGPLOT does work to your $PGPLOT_FONT location.


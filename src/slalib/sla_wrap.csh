echo "rm *.o libsla_wrap.a" > makefile.csh
ls -1 *.f | awk '{print "gfortran -fno-underscoring -O -c "$1}' >> makefile.csh
echo "gcc -Wall -c sla_wrap.c" >> makefile.csh
echo "ar ruv libsla_wrap.a *.o" >> makefile.csh

csh -v makefile.csh

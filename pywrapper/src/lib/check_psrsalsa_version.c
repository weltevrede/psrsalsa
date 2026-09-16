/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <string.h>
#include <Python.h>
#ifdef HAVE_NUMPY
#include <numpy/ndarrayobject.h>
#endif
int check_psrsalsa_version(void)
{
  char *libversion;
  libversion = get_psrsalsa_library_packaging_version();
  if(strcmp(libversion, PSRSALSA_PACKAGING_VERSION) != 0) {
    printerror(0, "ERROR check_psrsalsa_version: The psrsalsa version used to compile the psrsalsa python wrapper and library loaded at runtime are different. Expect serious problems, such as segfaults and unexpected behaviour. None of the obtained results can be trusted.");
    printerror(0, "PSRSALSA version: %s (library)", libversion);
    printerror(0, "PSRSALSA version: %s (wrapper)\n", PSRSALSA_PACKAGING_VERSION);
    return 0;
  }
#ifdef HAVE_NUMPY
  if(NPY_VERSION != PyArray_GetNDArrayCVersion()) {
    printerror(0, "Numpy version used by PSRsalsa inferface during compilation is not equal to that during runtime: %X != %X\n", NPY_VERSION, PyArray_GetNDArrayCVersion());
    return 0;
  }
#endif
  return 1;
}

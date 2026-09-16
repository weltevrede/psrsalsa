/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include "psrsalsa.h"
#include <Python.h>
#include "lib/check_psrsalsa_version.c"
static PyObject *method_ppol_from_IQUV(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule, *capsule2;
  static char *kwlist[] = {"archive", "onpulse", "sigma", "norm", "extended", "long_offset", "pa_offset", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  int extended = 0;
  float sigma = 3.0;
  int norm = 1;
  float long_offset = 0.0;
  float pa_offset = 0.0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|fppffpp", kwlist, &capsule, &capsule2, &sigma, &norm, &extended, &long_offset, &pa_offset, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *original;
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  pulselongitude_regions_definition *onpulse;
  onpulse = NULL;
  onpulse = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule2, "pulselongitude_regions_definition");
  if(onpulse == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int spstat = 0;
  int sigmaI = 0;
  int correctLbias = 0;
  int correctPbias = 0;
  float correctQV = 1.0;
  float correctV = 1.0;
  int nolongitudes = 0;
  datafile_definition *rms_file = NULL;
  float rebin_factor = 1.0;
  int onpulseonly = 0;
  if(make_paswing_fromIQUV(original, extended, spstat, sigma, sigmaI, *onpulse, norm, correctLbias, correctPbias, correctQV, correctV, nolongitudes, long_offset, pa_offset, rms_file, rebin_factor, onpulseonly, verbose) == 0) {
    PyErr_SetString(PyExc_TypeError, "Calculating polarization products failed");
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef polMethods[] = {
  {"ppol_from_IQUV", (PyCFunction) method_ppol_from_IQUV, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
    {NULL, NULL, 0, NULL}
};
static struct PyModuleDef pol_module = {
    PyModuleDef_HEAD_INIT,
    "pol",
    "Python interface for PSRSALSA",
    -1,
    polMethods
};
PyMODINIT_FUNC PyInit_pol(void)
{
  Py_Initialize();
#ifdef HAVE_NUMPY
#if NUMPY_MAJOR_VERSION >= 2
  if (PyArray_ImportNumPyAPI() < 0) {
    return NULL;
  }
#else
  import_array();
#endif
#endif
  check_psrsalsa_version();
  return PyModule_Create(&pol_module);
}

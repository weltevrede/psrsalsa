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
#include "lib/psrsalsa_array_to_python_2Dlist.c"
static PyObject *method_foldP3(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "foldp3", "nrp3bins", "onpulse", "offset", "refine", "cyclesperblock", "noSmooth", "smoothWidth", "slope", "no_numpy_out", "verbose", "debug", NULL};
  PyObject *capsule = NULL;
  PyObject *capsule2 = NULL;
  int nrp3bins = 0;
  int refine = 1;
  int cyclesperblock = 1;
  int noSmooth = 0;
  float foldp3 = 0.0;
  float smoothWidth = -1.0;
  float slope = 0.0;
  float offset = 0.0;
  int no_numpy = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Ofi|Ofiipffppp", kwlist, &capsule, &foldp3, &nrp3bins, &capsule2, &offset, &refine, &cyclesperblock, &noSmooth, &smoothWidth, &slope, &no_numpy, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *archive;
  archive = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(archive == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  pulselongitude_regions_definition *onpulse;
  onpulse = NULL;
  if(capsule2 != NULL) {
    onpulse = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule2, "pulselongitude_regions_definition");
    if(onpulse == NULL) {
      PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
      return NULL;
    }
  }
  if(archive->NrPols != 1 || archive->NrFreqChan != 1) {
    PyErr_SetString(PyExc_TypeError, "Expected data to have one polarization channel and one frequency channel");
    return NULL;
  }
  float *map;
  long datasize;
  datasize = archive->NrBins*nrp3bins;
  map = malloc(datasize*sizeof(float));
  if(map == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  slope = slope*360.0/(float)archive->NrBins;
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  char *offsetFileName = NULL;
  int readwriteOffsets = 0;
  int oldcorrelation = 0;
  if(foldP3(archive->data, archive->NrSubints, archive->NrBins, map, nrp3bins, foldp3, refine, cyclesperblock, noSmooth, smoothWidth, slope, offset, onpulse, offsetFileName, readwriteOffsets
     , oldcorrelation, verbose) != 1) {
    PyErr_SetString(PyExc_TypeError, "P3 folding failed");
    return NULL;
  }
  PyObject *ret;
  ret = psrsalsa_array_to_python_2Dlist(map, NULL, nrp3bins, archive->NrBins, 0, no_numpy);
  if(ret == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  free(map);
  return ret;
}
static PyMethodDef foldMethods[] = {
  {"foldP3", (PyCFunction) method_foldP3, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
{NULL, NULL, 0, NULL}
};
static struct PyModuleDef fold_module = {
    PyModuleDef_HEAD_INIT,
    "fold",
    "Python interface for PSRSALSA",
    -1,
    foldMethods
};
PyMODINIT_FUNC PyInit_fold(void)
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
  return PyModule_Create(&fold_module);
}

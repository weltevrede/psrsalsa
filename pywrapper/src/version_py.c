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
static PyObject *method_check_version_mismatch(PyObject *self, PyObject *args, PyObject *kwargs)
{
  int ret;
  ret = check_psrsalsa_version();
  if(ret == 0) {
    Py_INCREF(Py_False);
    return Py_False;
  }
  Py_INCREF(Py_True);
  return Py_True;
}
static PyObject *method_show_libversions(PyObject *self, PyObject *args, PyObject *kwargs)
{
  showlibraryversioninformation(stdout);
  #ifdef HAVE_NUMPY
  printf("Numpy version used by PSRsalsa inferface during compilation:  %d.%d %X %X %X\n", NUMPY_MAJOR_VERSION, NUMPY_MINOR_VERSION, NPY_VERSION, NPY_ABI_VERSION, NPY_FEATURE_VERSION);
  printf("Numpy version used by PSRsalsa inferface during runtime:      %X\n", PyArray_GetNDArrayCVersion());
  #endif
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_get_version(PyObject *self, PyObject *args, PyObject *kwargs)
{
  return PyUnicode_FromString(PSRSALSA_PACKAGING_VERSION);
}
static PyObject *method_numpy_support(PyObject *self, PyObject *args, PyObject *kwargs)
{
  #ifdef HAVE_NUMPY
  Py_INCREF(Py_True);
  return Py_True;
  #else
  Py_INCREF(Py_False);
  return Py_False;
  #endif
}
static PyMethodDef psrsalsaMethods[] = {
  {"check_version_mismatch", (PyCFunction) method_check_version_mismatch, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"numpy_support", (PyCFunction) method_numpy_support, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"show_libversions", (PyCFunction) method_show_libversions, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_version", (PyCFunction) method_get_version, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {NULL, NULL, 0, NULL}
};
static struct PyModuleDef psrsalsa_module = {
    PyModuleDef_HEAD_INIT,
    "version",
    "Python interface for PSRSALSA",
    -1,
    psrsalsaMethods
};
PyMODINIT_FUNC PyInit_version(void)
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
  return PyModule_Create(&psrsalsa_module);
}

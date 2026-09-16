/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <Python.h>
#ifdef HAVE_NUMPY
#include <numpy/ndarrayobject.h>
#endif
void destruct_numpy_array_memory(PyObject *capsule)
{
  double *array;
  array = (double *)PyCapsule_GetPointer(capsule, "numpy_double_array");
  if(array == NULL) {
    printerror(0, "Error destruct_numpy_array_memory: Function argument doesn't appear to be a pulsar data object");
  }else {
    free(array);
  }
}
static PyObject *psrsalsa_array_to_python_list(float *array_float, double *array_double, size_t length, int no_numpy)
{
  PyObject *list;
  int donumpy;
#ifdef HAVE_NUMPY
  donumpy = 1;
#else
  donumpy = 0;
#endif
  if(no_numpy != 0) {
    donumpy = 0;
  }
#ifdef HAVE_NUMPY
  if(donumpy != 0) {
    PyArray_Descr *numpy_descr;
    numpy_descr = PyArray_DescrFromType(NPY_DOUBLE);
    int nr_dimensions = 1;
    npy_intp dims[1];
    dims[0] = length;
    double *numpy_mem;
    numpy_mem = malloc(length*sizeof(double));
    if(numpy_mem == NULL) {
      printerror(0, "psrsalsa_array_to_python_list: Memory allocation error.");
      return NULL;
    }
    if(array_double == NULL && array_float != NULL) {
      size_t i;
      for(i = 0; i < length; i++) {
 numpy_mem[i] = array_float[i];
      }
    }
    if(array_double != NULL && array_float == NULL) {
      memcpy(numpy_mem, array_double, length*sizeof(double));
    }
    int flags = NPY_ARRAY_WRITEABLE;
    list = PyArray_NewFromDescr(&PyArray_Type, numpy_descr, nr_dimensions, dims, NULL, numpy_mem, flags, NULL);
    PyObject *capsule;
    capsule = PyCapsule_New(numpy_mem, "numpy_double_array", &destruct_numpy_array_memory);
    if(PyArray_SetBaseObject((PyArrayObject *)list, capsule) != 0) {
      printerror(0, "psrsalsa_array_to_python_list: Call to PyArray_SetBaseObject failed.");
      return NULL;
    }
  }
#endif
  if(donumpy == 0) {
    Py_ssize_t py_len, i;
    py_len = length;
    list = PyList_New(py_len);
    if(list == NULL) {
      printerror(0, "psrsalsa_array_to_python_list: Creation of python list failed.");
      return NULL;
    }
    int ret;
    for(i = 0; i < py_len; i++) {
      if(array_float != NULL) {
 ret = PyList_SetItem(list, i, PyFloat_FromDouble(array_float[i]));
      }else if(array_double != NULL) {
 ret = PyList_SetItem(list, i, PyFloat_FromDouble(array_double[i]));
      }
      if(ret != 0) {
 printerror(0, "psrsalsa_array_to_python_list:  Filling of python list failed.");
 return NULL;
      }
    }
  }
  return list;
}

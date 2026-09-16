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
void destruct_numpy_2Darray_memory(PyObject *capsule)
{
  double *array;
  array = (double *)PyCapsule_GetPointer(capsule, "numpy_double_2Darray");
  if(array == NULL) {
    printerror(0, "Error destruct_numpy_2Darray_memory: Function argument doesn't appear to be a pulsar data object");
  }else {
    free(array);
  }
}
static PyObject *psrsalsa_array_to_python_2Dlist(float *array_float, double *array_double, size_t length_y, size_t length_x, int transpose, int no_numpy)
{
  PyObject *list2D;
  int donumpy;
#ifdef HAVE_NUMPY
  donumpy = 1;
#else
  donumpy = 0;
#endif
  if(no_numpy != 0) {
    donumpy = 0;
  }
  Py_ssize_t py_lenx, py_leny, x, y;
  if(transpose == 0) {
    py_lenx = length_x;
    py_leny = length_y;
  }else {
    py_lenx = length_y;
    py_leny = length_x;
  }
#ifdef HAVE_NUMPY
  Py_ssize_t i;
  if(donumpy != 0) {
    PyArray_Descr *numpy_descr;
    numpy_descr = PyArray_DescrFromType(NPY_DOUBLE);
    int nr_dimensions = 2;
    npy_intp dims[2];
    dims[0] = py_leny;
    dims[1] = py_lenx;
    double *numpy_mem;
    numpy_mem = malloc(py_lenx*py_leny*sizeof(double));
    if(numpy_mem == NULL) {
      printerror(0, "psrsalsa_array_to_python_2Dlist: Memory allocation error.");
      return NULL;
    }
    if(array_double == NULL && array_float != NULL) {
      if(transpose == 0) {
 for(i = 0; i < (Py_ssize_t)(length_x*length_y); i++) {
   numpy_mem[i] = array_float[i];
 }
      }else {
 for(y = 0; y < (Py_ssize_t)length_y; y++) {
   for(x = 0; x < (Py_ssize_t)length_x; x++) {
     numpy_mem[x*length_y+y] = array_float[y*length_x+x];
   }
 }
      }
    }
    if(array_double != NULL && array_float == NULL) {
      if(transpose == 0) {
 memcpy(numpy_mem, array_double, length_x*length_y*sizeof(double));
      }else {
 for(y = 0; y < (Py_ssize_t)length_y; y++) {
   for(x = 0; x < (Py_ssize_t)length_x; x++) {
     numpy_mem[x*length_y+y] = array_double[y*length_x+x];
   }
 }
      }
    }
    int flags = NPY_ARRAY_WRITEABLE;
    list2D = PyArray_NewFromDescr(&PyArray_Type, numpy_descr, nr_dimensions, dims, NULL, numpy_mem, flags, NULL);
    PyObject *capsule;
    capsule = PyCapsule_New(numpy_mem, "numpy_double_2Darray", &destruct_numpy_2Darray_memory);
    if(PyArray_SetBaseObject((PyArrayObject *)list2D, capsule) != 0) {
      printerror(0, "psrsalsa_array_to_python_2Dlist: Call to PyArray_SetBaseObject failed.");
      return NULL;
    }
  }
#endif
  if(donumpy == 0) {
    list2D = PyList_New(py_leny);
    if(list2D == NULL) {
      printerror(0, "psrsalsa_array_to_python_2Dlist: Creation of python list failed.");
      return NULL;
    }
    int ret;
    for(y = 0; y < py_leny; y++) {
      PyObject *list;
      list = PyList_New(py_lenx);
      if(list == NULL) {
 return NULL;
      }
      ret = PyList_SetItem(list2D, y, list);
      if(ret != 0) {
 printerror(0, "psrsalsa_array_to_python_2Dlist: Setting of python list item failed.");
 return NULL;
      }
      if(transpose == 0) {
 for(x = 0; x < py_lenx; x++) {
   if(array_double == NULL && array_float != NULL) {
     ret = PyList_SetItem(list, x, PyFloat_FromDouble(array_float[y*length_x + x]));
   }
   if(array_double != NULL && array_float == NULL) {
     ret = PyList_SetItem(list, x, PyFloat_FromDouble(array_double[y*length_x + x]));
   }
   if(ret != 0) {
     printerror(0, "psrsalsa_array_to_python_2Dlist: Setting of python list item failed.");
     return NULL;
   }
 }
      }else {
 for(x = 0; x < py_lenx; x++) {
   if(array_double == NULL && array_float != NULL) {
     ret = PyList_SetItem(list, x, PyFloat_FromDouble(array_float[x*length_x + y]));
   }
   if(array_double != NULL && array_float == NULL) {
     ret = PyList_SetItem(list, x, PyFloat_FromDouble(array_double[x*length_x + y]));
   }
   if(ret != 0) {
     printerror(0, "psrsalsa_array_to_python_2Dlist: Setting of python list item failed.");
     return NULL;
   }
 }
      }
    }
  }
  return list2D;
}

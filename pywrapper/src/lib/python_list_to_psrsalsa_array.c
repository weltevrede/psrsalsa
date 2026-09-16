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
static void *python_list_to_psrsalsa_array(PyObject *list, long expected_array_length, float *inplace_float, double *inplace_double, int return_data_type, long *length, int *error_status)
{
  Py_ssize_t py_len, i;
  *error_status = 0;
  *length = 0;
  int list_is_numpy = 0;
#ifdef HAVE_NUMPY
  list_is_numpy = PyArray_Check(list);
  PyObject *list_npy_c;
  list_npy_c = NULL;
  if(list_is_numpy) {
    list_npy_c = PyArray_FROM_OTF(list, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if(list_npy_c == NULL) {
      PyErr_SetString(PyExc_ValueError, "python_list_to_psrsalsa_array: Cannot obtain a C-style array from the provided numpy array.");
      return NULL;
    }
  }
#endif
  if(list_is_numpy == 0) {
    py_len = PyList_Size(list);
#ifdef HAVE_NUMPY
  }else {
    int nrdimensions;
    nrdimensions = PyArray_NDIM((PyArrayObject *)list_npy_c);
    if(nrdimensions != 1) {
      PyErr_SetString(PyExc_ValueError, "python_list_to_psrsalsa_array: Input numpy array with data is expected to be a 1D list.");
      Py_DECREF(list_npy_c);
      return NULL;
    }
    npy_intp *list_array_lenghts;
    list_array_lenghts = PyArray_DIMS((PyArrayObject *)list_npy_c);
    py_len = list_array_lenghts[0];
#endif
  }
  if(py_len == -1) {
#ifdef HAVE_NUMPY
    if(list_is_numpy) {
      PyErr_SetString(PyExc_ValueError, "python_list_to_psrsalsa_array: Unexpected array length -1 of numpy array.");
      Py_DECREF(list_npy_c);
    }else {
      PyErr_SetString(PyExc_ValueError, "python_list_to_psrsalsa_array: Input variable was expected to be a list. But the variable is neither a numpy array nor a python list?");
    }
#else
    PyErr_SetString(PyExc_ValueError, "python_list_to_psrsalsa_array: expected a python list of finite length as the input. If it is a numpy array (numpy support was disabled in the PSRSalsa python interface during compilation), convert it to a python list first using numpy_array_variable_name.tolist().");
#endif
    return NULL;
  }
  if(py_len == 0) {
    PyErr_SetString(PyExc_ValueError, "python_list_to_psrsalsa_array: expected a list of non-zero length as the input.");
#ifdef HAVE_NUMPY
    if(list_is_numpy) {
      Py_DECREF(list_npy_c);
    }
#endif
    return NULL;
  }
  if(expected_array_length > 0) {
    if(py_len != expected_array_length) {
#ifdef HAVE_NUMPY
      if(list_is_numpy) {
 Py_DECREF(list_npy_c);
      }
#endif
      *length = py_len;
      *error_status = 1;
      return NULL;
    }
  }
  float *float_array;
  double *double_array;
  float_array = NULL;
  double_array = NULL;
  if(inplace_float != NULL) {
    float_array = inplace_float;
  }else if(inplace_double != NULL) {
    double_array = inplace_double;
  }else {
    if(return_data_type != sizeof(float) && return_data_type != sizeof(double)) {
      PyErr_SetString(PyExc_ValueError, "python_list_to_psrsalsa_array: expected output data type should either be floats or doubles.");
#ifdef HAVE_NUMPY
      if(list_is_numpy) {
 Py_DECREF(list_npy_c);
      }
#endif
      return NULL;
    }
    int alloc_error;
    alloc_error = 0;
    if(return_data_type == sizeof(float)) {
      float_array = malloc(py_len*return_data_type);
      if(float_array == NULL) {
 alloc_error = 1;
      }
    }else if(return_data_type == sizeof(double)) {
      double_array = malloc(py_len*return_data_type);
      if(double_array == NULL) {
 alloc_error = 1;
      }
    }
    if(alloc_error) {
      PyErr_SetString(PyExc_MemoryError, "python_list_to_psrsalsa_array: Memory allocation error");
#ifdef HAVE_NUMPY
      if(list_is_numpy) {
 Py_DECREF(list_npy_c);
      }
#endif
      return NULL;
    }
  }
  if(list_is_numpy == 0) {
    PyObject *item;
    for(i = 0; i < py_len; i++) {
      item = PyList_GetItem(list, i);
      if(item == NULL) {
 PyErr_SetString(PyExc_ValueError, "python_list_to_psrsalsa_array: Error obtaining item from python list");
 return NULL;
      }
      if(float_array != NULL) {
 float_array[i] = PyFloat_AsDouble(item);
      }else {
 double_array[i] = PyFloat_AsDouble(item);
      }
      if(PyErr_Occurred() != NULL) {
 PyErr_SetString(PyExc_ValueError, "python_list_to_psrsalsa_array: Error converting value from python list to double. A reason could be that the list is a nested list rather than a 1D list.");
 return NULL;
      }
    }
#ifdef HAVE_NUMPY
  }else {
    double *numpy_ptr;
    numpy_ptr = (double *)PyArray_DATA((PyArrayObject *)list_npy_c);
    if(float_array != NULL) {
      for(i = 0; i < py_len; i++) {
 float_array[i] = numpy_ptr[i];
      }
    }else if(double_array != NULL) {
      memcpy(double_array, numpy_ptr, py_len*sizeof(double));
    }
#endif
  }
#ifdef HAVE_NUMPY
  if(list_is_numpy) {
    Py_DECREF(list_npy_c);
  }
#endif
  *length = py_len;
  if(float_array != NULL) {
    return float_array;
  }else {
    return double_array;
  }
}

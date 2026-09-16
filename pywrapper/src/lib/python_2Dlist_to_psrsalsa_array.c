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
static void *python_2Dlist_to_psrsalsa_array(PyObject *matrix, int transpose, long expected_array_length_y, long expected_array_length_x, float *inplace_float, double *inplace_double, int return_data_type, long *length_x, long *length_y, int *error_status)
{
  Py_ssize_t py_len_x, py_len_y, x, y;
  *error_status = 0;
  *length_x = 0;
  *length_y = 0;
  int matrix_is_numpy = 0;
#ifdef HAVE_NUMPY
  matrix_is_numpy = PyArray_Check(matrix);
  PyObject *matrix_npy_c;
  matrix_npy_c = NULL;
  if(matrix_is_numpy) {
    matrix_npy_c = PyArray_FROM_OTF(matrix, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if(matrix_npy_c == NULL) {
      PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: Cannot obtain a C-style array from the provided numpy array.");
      return NULL;
    }
  }
#endif
  PyObject **rows_python_list;
  rows_python_list = NULL;
  if(matrix_is_numpy == 0) {
    py_len_y = PyList_Size(matrix);
    if(py_len_y < 0) {
#ifdef HAVE_NUMPY
      PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: Expected a 2D list (matrix) as input, but input cannot be interpreted as neither a python list nor a numpy array.");
#else
      PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: Expected a 2D list (matrix) as input, but input cannot be interpreted as a python list. If it is a numpy array (numpy support was disabled in the PSRSalsa python interface during compilation), convert it to a python list first using numpy_array_variable_name.tolist().");
#endif
      return NULL;
    }
    if(py_len_y == 0) {
      PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: Expected a 2D list (matrix) as input, but input list appears to have zero length.");
      return NULL;
    }
    rows_python_list = malloc(py_len_y*sizeof(PyObject *));
    if(rows_python_list == NULL) {
      PyErr_SetString(PyExc_MemoryError, "python_2Dlist_to_psrsalsa_array: Memory allocation error");
      return NULL;
    }
    for(y = 0; y < py_len_y; y++) {
      rows_python_list[y] = PyList_GetItem(matrix, y);
      if(rows_python_list[y] == NULL) {
 free(rows_python_list);
 PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: Error obtaining item from list");
 return NULL;
      }
      x = PyList_Size(rows_python_list[y]);
      if(y == 0) {
 py_len_x = x;
 if(py_len_x < 0) {
   free(rows_python_list);
   PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: expected a 2D Python list, but input list has at least one entry that cannot be interpreted as a python list, so it doesn't appear to be a 2D list.");
   return NULL;
 }
 if(py_len_x == 0) {
   free(rows_python_list);
   PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: expected a 2D Python list with non-zero dimensions as the input");
   return NULL;
 }
      }else {
 if(x != py_len_x) {
   free(rows_python_list);
   PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: expected a 2D Python list which is rectangular");
   return NULL;
 }
      }
    }
#ifdef HAVE_NUMPY
  }else {
    int nrdimensions;
    nrdimensions = PyArray_NDIM((PyArrayObject *)matrix_npy_c);
    if(nrdimensions != 2) {
      PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: Input numpy array with data is expected to be a 2D array.");
      Py_DECREF(matrix_npy_c);
      return NULL;
    }
    npy_intp *matrix_array_lenghts;
    matrix_array_lenghts = PyArray_DIMS((PyArrayObject *)matrix_npy_c);
    py_len_y = matrix_array_lenghts[0];
    py_len_x = matrix_array_lenghts[1];
    if(py_len_x <= 0 || py_len_y <= 0) {
      PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: expected a 2D numpy array with non-zero dimensions as the input");
      Py_DECREF(matrix_npy_c);
      return NULL;
    }
#endif
  }
  if(expected_array_length_x > 0) {
    if(py_len_x != expected_array_length_x) {
      *error_status = 1;
    }
  }
  if(expected_array_length_y > 0) {
    if(py_len_y != expected_array_length_y) {
      *error_status = 1;
    }
  }
  if(*error_status == 1) {
#ifdef HAVE_NUMPY
    if(matrix_is_numpy) {
      Py_DECREF(matrix_npy_c);
    }
#endif
    *length_x = py_len_x;
    *length_y = py_len_x;
    return NULL;
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
      PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: expected output data type should either be floats or doubles.");
#ifdef HAVE_NUMPY
      if(matrix_is_numpy) {
 Py_DECREF(matrix_npy_c);
      }
#endif
      return NULL;
    }
    int alloc_error;
    alloc_error = 0;
    if(return_data_type == sizeof(float)) {
      float_array = malloc(py_len_x*py_len_y*return_data_type);
      if(float_array == NULL) {
 alloc_error = 1;
      }
    }else if(return_data_type == sizeof(double)) {
      double_array = malloc(py_len_x*py_len_y*return_data_type);
      if(double_array == NULL) {
 alloc_error = 1;
      }
    }
    if(alloc_error) {
      PyErr_SetString(PyExc_MemoryError, "python_2Dlist_to_psrsalsa_array: Memory allocation error");
#ifdef HAVE_NUMPY
      if(matrix_is_numpy) {
 Py_DECREF(matrix_npy_c);
      }
#endif
      return NULL;
    }
  }
  if(matrix_is_numpy == 0) {
    PyObject *item;
    double value;
    for(y = 0; y < py_len_y; y++) {
      for(x = 0; x < py_len_x; x++) {
 item = PyList_GetItem(rows_python_list[y], x);
 if(item == NULL) {
   free(rows_python_list);
   if(inplace_float != NULL) {
     free(float_array);
   }
   if(inplace_double != NULL) {
     free(double_array);
   }
   PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: Error obtaining item from list");
   return NULL;
 }
 value = PyFloat_AsDouble(item);
 if(PyErr_Occurred() != NULL) {
   free(rows_python_list);
   if(inplace_float != NULL) {
     free(float_array);
   }
   if(inplace_double != NULL) {
     free(double_array);
   }
   PyErr_SetString(PyExc_ValueError, "python_2Dlist_to_psrsalsa_array: Error converting value");
   return NULL;
 }
 if(inplace_float != NULL || return_data_type == sizeof(float)) {
   if(transpose == 0) {
     float_array[y*py_len_x+x] = value;
   }else {
     float_array[x*py_len_y+y] = value;
   }
 }else {
   if(transpose == 0) {
     double_array[y*py_len_x+x] = value;
   }else {
     double_array[x*py_len_y+y] = value;
   }
 }
      }
    }
#ifdef HAVE_NUMPY
  }else {
    double *numpy_ptr;
    numpy_ptr = (double *)PyArray_DATA((PyArrayObject *)matrix_npy_c);
    for(y = 0; y < py_len_y; y++) {
      if(inplace_float != NULL || return_data_type == sizeof(float) || transpose != 0) {
 for(x = 0; x < py_len_x; x++) {
   double value;
   value = numpy_ptr[y*py_len_x+x];
   if(inplace_float != NULL || return_data_type == sizeof(float)) {
     if(transpose == 0) {
       float_array[y*py_len_x+x] = value;
     }else {
       float_array[x*py_len_y+y] = value;
     }
   }else {
     double_array[x*py_len_y+y] = value;
   }
 }
      }else {
 memcpy(&(double_array[y*py_len_x]), &(numpy_ptr[y*py_len_x]), py_len_x*sizeof(double));
      }
    }
#endif
  }
#ifdef HAVE_NUMPY
  if(matrix_is_numpy) {
    Py_DECREF(matrix_npy_c);
  }
#endif
  *length_x = py_len_x;
  *length_y = py_len_y;
  if(matrix_is_numpy == 0) {
    free(rows_python_list);
  }
  if(float_array != NULL) {
    return float_array;
  }else {
    return double_array;
  }
}

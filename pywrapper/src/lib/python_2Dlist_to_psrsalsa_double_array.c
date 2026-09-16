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
static double *python_2Dlist_to_psrsalsa_double_array(PyObject *matrix, unsigned long *nrx, unsigned long *nry, int transpose)
{
  Py_ssize_t py_lenx, py_leny, x, y;
  *nrx = 0;
  *nry = 0;
  py_leny = PyList_Size(matrix);
  if(py_leny <= 0) {
    printerror(0, "python_2Dlist_to_psrsalsa_double_array: expected a 2D Python list with non-zero dimensions as the input. If it is a numpy array, make sure to convert it to a Python list first with variable.tolist().");
    return NULL;
  }
  PyObject **rows;
  rows = malloc(py_leny*sizeof(PyObject *));
  if(rows == NULL) {
    printerror(0, "python_2Dlist_to_psrsalsa_double_array: Memory allocation error");
    return NULL;
  }
  for(y = 0; y < py_leny; y++) {
    rows[y] = PyList_GetItem(matrix, y);
    if(rows[y] == NULL) {
      free(rows);
      printerror(0, "python_list_to_psrsalsa_double_array: Error obtaining item");
      return NULL;
    }
    x = PyList_Size(rows[y]);
    if(y == 0) {
      py_lenx = x;
      if(py_lenx <= 0) {
 free(rows);
 printerror(0, "python_2Dlist_to_psrsalsa_double_array: expected a 2D Python list with non-zero dimensions as the input");
 return NULL;
      }
    }else {
      if(x != py_lenx) {
 free(rows);
 printerror(0, "python_2Dlist_to_psrsalsa_double_array: expected a 2D Python list which is rectangular");
 return NULL;
      }
    }
  }
  double *array;
  array = malloc(py_lenx*py_leny*sizeof(double));
  if(array == NULL) {
    free(rows);
    printerror(0, "python_2Dlist_to_psrsalsa_double_array: Memory allocation error");
    return NULL;
  }
  PyObject *item;
  double value;
  for(y = 0; y < py_leny; y++) {
    for(x = 0; x < py_lenx; x++) {
      item = PyList_GetItem(rows[y], x);
      if(item == NULL) {
 free(rows);
 free(array);
 printerror(0, "python_2Dlist_to_psrsalsa_double_array: Error obtaining item");
 return NULL;
      }
      value = PyFloat_AsDouble(item);
      if(PyErr_Occurred() != NULL) {
 free(rows);
 free(array);
 printerror(0, "python_2Dlist_to_psrsalsa_double_array: Error converting value");
 return NULL;
      }
      if(transpose == 0) {
 array[y*py_lenx+x] = value;
      }else {
 array[x*py_leny+y] = value;
      }
    }
  }
  *nrx = py_lenx;
  *nry = py_leny;
  free(rows);
  return array;
}

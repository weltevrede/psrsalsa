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
static int python_nonrectangular_2Dlist_to_psrsalsa_int_array(PyObject *matrix, int ***intlist, int **nrx, int *nry)
{
  Py_ssize_t py_leny, x, y;
  *nrx = NULL;
  *nry = 0;
  *intlist = NULL;
  py_leny = PyList_Size(matrix);
  if(py_leny <= 0) {
    printerror(0, "python_nonrectangular_2Dlist_to_psrsalsa_int_array: expected a 2D Python list with non-zero dimensions as the input. If it is a numpy array, make sure to convert it to a Python list first with variable.tolist().");
    return 1;
  }
  *nry = py_leny;
  *nrx = malloc(py_leny*sizeof(unsigned long *));
  PyObject **rows;
  rows = malloc(py_leny*sizeof(PyObject *));
  if(*nrx == NULL || rows == NULL) {
    printerror(0, "python_nonrectangular_2Dlist_to_psrsalsa_int_array: Memory allocation error 1");
    return 1;
  }
  for(y = 0; y < py_leny; y++) {
    rows[y] = PyList_GetItem(matrix, y);
    if(rows[y] == NULL) {
      free(rows);
      free(*nrx);
      *nrx = NULL;
      printerror(0, "python_list_to_psrsalsa_float_array: Error obtaining item");
      return 1;
    }
    (*nrx)[y] = PyList_Size(rows[y]);
    if((*nrx)[y] <= 0) {
      free(rows);
      free(*nrx);
      *nrx = NULL;
      printerror(0, "python_nonrectangular_2Dlist_to_psrsalsa_int_array: expected a 2D Python list with non-zero dimensions as the input");
      return 1;
    }
  }
  *intlist = malloc(py_leny*sizeof(int *));
  if(*intlist == NULL) {
    free(rows);
    free(*nrx);
    *nrx = NULL;
    printerror(0, "python_nonrectangular_2Dlist_to_psrsalsa_int_array: Memory allocation error 2");
    return 1;
  }
  for(y = 0; y < py_leny; y++) {
    (*intlist)[y] = malloc((*nrx)[y]*sizeof(int));
    if((*intlist)[y] == NULL) {
      free(rows);
      free(*nrx);
      *nrx = NULL;
      printerror(0, "python_nonrectangular_2Dlist_to_psrsalsa_int_array: Memory allocation error 3");
      return 1;
    }
  }
  PyObject *item;
  int value;
  for(y = 0; y < py_leny; y++) {
    for(x = 0; x < (*nrx)[y]; x++) {
      item = PyList_GetItem(rows[y], x);
      if(item == NULL) {
 free(rows);
 free(*nrx);
 *nrx = NULL;
 int y2;
 for(y2 = 0; y2 < py_leny; y2++) {
   free((*intlist)[y2]);
 }
 free(*intlist);
 *intlist = NULL;
 printerror(0, "python_nonrectangular_2Dlist_to_psrsalsa_int_array: Error obtaining item");
 return 1;
      }
      value = PyLong_AsLong(item);
      if(PyErr_Occurred() != NULL) {
 free(rows);
 free(*nrx);
 *nrx = NULL;
 int y2;
 for(y2 = 0; y2 < py_leny; y2++) {
   free((*intlist)[y2]);
 }
 free(*intlist);
 *intlist = NULL;
 printerror(0, "python_nonrectangular_2Dlist_to_psrsalsa_int_array: Error converting value");
 return 1;
      }
      (*intlist)[y][x] = value;
    }
  }
  return 0;
}

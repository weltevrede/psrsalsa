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
static double *python_list_to_psrsalsa_double_array(PyObject *list)
{
  Py_ssize_t py_len, i;
  py_len = PyList_Size(list);
  if(py_len <= 0) {
    printerror(0, "python_list_to_psrsalsa_double_array: expected a list of finite length as the input. If it is a numpy array, make sure to convert it to a Python list first with variable.tolist().");
    return NULL;
  }
  double *array;
  array = malloc(py_len*sizeof(double));
  if(array == NULL) {
    printerror(0, "python_list_to_psrsalsa_double_array: Memory allocation error");
    return NULL;
  }
  PyObject *item;
  for(i = 0; i < py_len; i++) {
    item = PyList_GetItem(list, i);
    if(item == NULL) {
      printerror(0, "python_list_to_psrsalsa_double_array: Error obtaining item");
      return NULL;
    }
    array[i] = PyFloat_AsDouble(item);
    if(PyErr_Occurred() != NULL) {
      printerror(0, "python_list_to_psrsalsa_double_array: Error converting value");
      return NULL;
    }
  }
  return array;
}

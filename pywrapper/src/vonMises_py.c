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
#include "lib/destruct_vonMises_collection_definition.c"
#include "lib/psrsalsa_array_to_python_list.c"
#include "lib/python_list_to_psrsalsa_array.c"
#include "lib/python_nonrectangular_2Dlist_to_psrsalsa_int_array.c"
static PyObject *method_create_empty(PyObject *self, PyObject *args)
{
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)malloc(sizeof(vonMises_collection_definition));
  if(model == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  vonMises_collection_initialise(model);
  return PyCapsule_New(model, "vonMises_collection_definition", &destruct_vonMises_collection_definition);
}
static PyObject *method_readVonMisesModel(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"filename", "verbose", "debug", NULL};
  char *filename;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s|pp", kwlist, &filename, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)malloc(sizeof(vonMises_collection_definition));
  if(model == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int ret;
  ret = readVonMisesModel(filename, model, verbose);
  if(ret != 0) {
    PyErr_SetString(PyExc_FileNotFoundError, "Reading von Mises collection file failed");
    return NULL;
  }
  return PyCapsule_New(model, "vonMises_collection_definition", &destruct_vonMises_collection_definition);
}
static PyObject *method_add_component(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  static char *kwlist[] = {"model", "centre", "concentration", "height", "verbose", "debug", NULL};
  double centre = 0;
  double concentration = 0;
  double height = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oddd|pp", kwlist, &capsule, &centre, &concentration, &height, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int ret;
  ret = vonMises_collection_add_component(model, centre, concentration, height, verbose);
  if(ret != 0) {
    PyErr_SetString(PyExc_FileNotFoundError, "Adding component to von Mises collection failed");
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_change_component(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  static char *kwlist[] = {"model", "compnr", "centre", "concentration", "height", "verbose", "debug", NULL};
  int compnr = 0;
  double centre = 0;
  double concentration = 0;
  double height = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oiddd|pp", kwlist, &capsule, &compnr, &centre, &concentration, &height, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int ret;
  ret = vonMises_collection_change_component(model, compnr, centre, concentration, height, verbose);
  if(ret != 0) {
    PyErr_SetString(PyExc_FileNotFoundError, "Changing component to von Mises collection failed");
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_print_model(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  static char *kwlist[] = {"model", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|pp", kwlist, &capsule, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  vonMises_collection_print_components(*model, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_create_profile(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  static char *kwlist[] = {"model", "nrbins", "shift", "norm", "no_numpy_out", "verbose", "debug", NULL};
  int nrbins = 0;
  double shift = 0;
  int norm = 0;
  int no_numpy = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi|dpppp", kwlist, &capsule, &nrbins, &shift, &norm, &no_numpy, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  float *profile;
  profile = (float *)malloc(nrbins*sizeof(float));
  if(profile == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  calcVonMisesProfile(model, nrbins, profile, shift, norm);
  PyObject *ret;
  ret = psrsalsa_array_to_python_list(profile, NULL, nrbins, no_numpy);
  if(ret == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  return ret;
}
static PyObject *method_correlate_profile(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"profile", "model", "verbose", "debug", NULL};
  PyObject *capsule, *profile_list;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|pp", kwlist, &profile_list, &capsule, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  float *profile;
  long nrBins;
  int error_status;
  profile = python_list_to_psrsalsa_array(profile_list, -1, NULL, NULL, sizeof(float), &nrBins, &error_status);
  if(profile == NULL) {
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  float shift;
  shift = correlateVonMisesFunction(model, nrBins, profile, verbose);
  free(profile);
  return PyFloat_FromDouble(shift);
}
static PyObject *method_align_archive(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "model", "locked", "verbose", "debug", NULL};
  PyObject *capsule, *archive;
  int verbose_flag = 0;
  int debug_flag = 0;
  int locked = 1;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|ppp", kwlist, &archive, &capsule, &locked, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(archive, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int doalign;
  if(locked == 0) {
    doalign = 2;
  }else {
    doalign = 1;
  }
  float shift;
  int tscr_complete = 0;
  if(apply_doalign(datafile, doalign, model, NULL, &shift, tscr_complete, verbose) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Error aligning the data");
    return NULL;
  }
  if(doalign == 1) {
    if(preprocess_fftshift(*datafile, -1, -shift, 0, 0.0, verbose) == 0) {
      PyErr_SetString(PyExc_RuntimeError, "Error applying alignment");
      return NULL;
    }
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_refine_model(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "model", "no_neg", "fixamp", "fixwidth", "fixphase", "fixrelamp", "fixrelphase", "fixrelamp_2Dlist", "pre_align", "fit_baseline", "suppress_error", "verbose", "debug", NULL};
  PyObject *capsule, *archive, *fixrelamp_2dlist;
  fixrelamp_2dlist = NULL;
  int verbose_flag = 0;
  int debug_flag = 0;
  int suppress_error = 0;
  int fitbaseline = 0;
  int avoid_neg_components = 0;
  int fixamp = 0;
  int fixwidth = 0;
  int fixphase = 0;
  int fixrelamp = 0;
  int fixrelphase = 0;
  int pre_align = 1;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|ppppppOppppp", kwlist, &archive, &capsule, &avoid_neg_components, &fixamp, &fixwidth, &fixphase, &fixrelamp, &fixrelphase, &fixrelamp_2dlist, &pre_align, &fitbaseline, &suppress_error, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  int **refine_fixrelamp_list = NULL;
  int *refine_fixrelamp_list_nrcomps = NULL;
  int refine_fixrelamp_listlength = 0;
  if(fixrelamp_2dlist != NULL) {
    if(python_nonrectangular_2Dlist_to_psrsalsa_int_array(fixrelamp_2dlist, &refine_fixrelamp_list, &refine_fixrelamp_list_nrcomps, &refine_fixrelamp_listlength) != 0) {
      PyErr_SetString(PyExc_TypeError, "Error interpreting fixrelamp_2Dlist");
      return NULL;
    }
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(archive, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  float baseline;
  int no_pre_align = 1;
  if(pre_align) {
    no_pre_align = 0;
  }
  int ret;
  ret = fitvonmises_refine_model(*datafile, model, fitbaseline, avoid_neg_components, &baseline, fixamp, fixwidth, fixphase, fixrelamp, refine_fixrelamp_listlength, refine_fixrelamp_list_nrcomps, refine_fixrelamp_list, fixrelphase, no_pre_align, verbose);
  if(suppress_error == 0) {
    if(ret == 1) {
      PyErr_SetString(PyExc_RuntimeError, "Maximum number of fit itterations exceeded");
      return NULL;
    }else if(ret != 0) {
      PyErr_SetString(PyExc_RuntimeError, "Fitting failed");
      return NULL;
    }
  }else {
    if(ret == 1) {
      printwarning(debug_flag, "WARNING fitvonmises_refine_model: Maximum number of fit itterations exceeded");
    }else if(ret != 0) {
      printwarning(debug_flag, "WARNING fitvonmises_refine_model: Fitting failed");
    }
  }
  return PyLong_FromLong(ret);
}
static PyObject *method_clone(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"model", "verbose", "debug", NULL};
  PyObject *capsule;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|pp", kwlist, &capsule, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  vonMises_collection_definition *clone;
  clone = (vonMises_collection_definition *)malloc(sizeof(vonMises_collection_definition));
  if(clone == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  if(vonMises_collection_clone(*model, clone, verbose) != 0) {
    PyErr_SetString(PyExc_MemoryError, "Cloning of model failed");
    return NULL;
  }
  return PyCapsule_New(clone, "vonMises_collection_definition", &destruct_vonMises_collection_definition);
}
static PyObject *method_get_phase(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"model", "compnr", "verbose", "debug", NULL};
  PyObject *capsule;
  int compnr;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi|pp", kwlist, &capsule, &compnr, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(compnr < 0 || compnr >= model->nrcomponents) {
    PyErr_SetString(PyExc_ValueError, "Requested component nr (counting from zero) is invalid");
    return NULL;
  }
  double centre;
  centre = model->centre[compnr];
  return PyFloat_FromDouble(centre);
}
static PyObject *method_get_concentration(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"model", "compnr", "verbose", "debug", NULL};
  PyObject *capsule;
  int compnr;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi|pp", kwlist, &capsule, &compnr, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(compnr < 0 || compnr >= model->nrcomponents) {
    PyErr_SetString(PyExc_ValueError, "Requested component nr (counting from zero) is invalid");
    return NULL;
  }
  double concentration;
  concentration = model->concentration[compnr];
  return PyFloat_FromDouble(concentration);
}
static PyObject *method_get_amp(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"model", "compnr", "verbose", "debug", NULL};
  PyObject *capsule;
  int compnr;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi|pp", kwlist, &capsule, &compnr, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(compnr < 0 || compnr >= model->nrcomponents) {
    PyErr_SetString(PyExc_ValueError, "Requested component nr (counting from zero) is invalid");
    return NULL;
  }
  double height;
  height = model->height[compnr];
  return PyFloat_FromDouble(height);
}
static PyObject *method_get_integral(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"model", "compnr", "verbose", "debug", NULL};
  PyObject *capsule;
  int compnr = -1;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ipp", kwlist, &capsule, &compnr, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(compnr >= model->nrcomponents) {
    PyErr_SetString(PyExc_ValueError, "Requested component number (should be counting from zero) is invalid");
    return NULL;
  }
  double area;
  if(compnr < 0) {
    area = integrateVonMisesFunction(model);
  }else {
    area = integrateVonMisesFunction2(model->concentration[compnr], model->height[compnr]);
  }
  return PyFloat_FromDouble(area);
}
static PyObject *method_get_nrcomp(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"model", "verbose", "debug", NULL};
  PyObject *capsule;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|pp", kwlist, &capsule, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  vonMises_collection_definition *model;
  model = (vonMises_collection_definition *)PyCapsule_GetPointer(capsule, "vonMises_collection_definition");
  if(model == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a von Mises component collection object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  return PyLong_FromLong(model->nrcomponents);
}
static PyMethodDef vonMisesMethods[] = {
  {"load_model", (PyCFunction) method_readVonMisesModel, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"add_component", (PyCFunction) method_add_component, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"change_component", (PyCFunction) method_change_component, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"clone", (PyCFunction) method_clone, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"create_empty", method_create_empty, METH_VARARGS, "Python interface for psrsalsa library function"},
  {"get_amp", (PyCFunction) method_get_amp, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_concentration", (PyCFunction) method_get_concentration, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_integral", (PyCFunction) method_get_integral, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_nrcomp", (PyCFunction) method_get_nrcomp, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_phase", (PyCFunction) method_get_phase, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"align_archive", (PyCFunction) method_align_archive, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"print_model", (PyCFunction) method_print_model, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"create_profile", (PyCFunction) method_create_profile, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"correlate_profile", (PyCFunction) method_correlate_profile, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"refine_model", (PyCFunction) method_refine_model, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
{NULL, NULL, 0, NULL}
};
static struct PyModuleDef vonMises_module = {
    PyModuleDef_HEAD_INIT,
    "vonMises",
    "Python interface for PSRSALSA",
    -1,
    vonMisesMethods
};
PyMODINIT_FUNC PyInit_vonMises(void)
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
  return PyModule_Create(&vonMises_module);
}

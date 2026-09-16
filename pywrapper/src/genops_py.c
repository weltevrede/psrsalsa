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
#include "lib/destruct_datafile_definition.c"
#include "lib/destruct_pulselongitude_regions_definition.c"
#include "lib/python_list_to_psrsalsa_array.c"
#include "lib/psrsalsa_int_array_to_python_list.c"
static PyObject *method_make_profile(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original, *new;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "stokesI", "median", "verbose", "debug", NULL};
  int stokesI = 0;
  int median = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|pppp", kwlist, &capsule, &stokesI, &median, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  new = malloc(sizeof(datafile_definition));
  if(new == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  if(preprocess_make_profile(*original, new, stokesI, median, verbose) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Making profile failed");
    return NULL;
  }
  return PyCapsule_New(new, "datafile_definition", &destruct_datafile_definition);
}
static PyObject *method_rebin(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original, *new;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "nrbins",
      "verbose", "debug", NULL};
  long nrbins = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  int interp = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Ol|ppp", kwlist, &capsule, &nrbins,
      &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  new = malloc(sizeof(datafile_definition));
  if(new == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  if(preprocess_rebin(*original, new, nrbins, interp, verbose) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Rebinning failed");
    return NULL;
  }
  return PyCapsule_New(new, "datafile_definition", &destruct_datafile_definition);
}
static PyObject *method_shuffle(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original, *new;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "fixseed", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  int fixseed = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ppp", kwlist, &capsule, &fixseed, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  new = malloc(sizeof(datafile_definition));
  if(new == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  if(preprocess_shuffle(*original, new, fixseed, verbose) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Shuffling failed");
    return NULL;
  }
  return PyCapsule_New(new, "datafile_definition", &destruct_datafile_definition);
}
static PyObject *method_tscr(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original, *new;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "nrsubints", "complete_only", "verbose", "debug", NULL};
  long nrsubints = -1;
  int complete = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|lppp", kwlist, &capsule, &nrsubints, &complete, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  if(nrsubints == -1) {
    nrsubints = original->NrSubints;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  new = malloc(sizeof(datafile_definition));
  if(new == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  if(preprocess_addsuccessivepulses(*original, new, nrsubints, complete, verbose) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Adding subints failed");
    return NULL;
  }
  return PyCapsule_New(new, "datafile_definition", &destruct_datafile_definition);
}
static PyObject *method_fscr(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original, *new;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "nrfreqs", "complete_only", "verbose", "debug", NULL};
  long nrfreqs = -1;
  int complete = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|lppp", kwlist, &capsule, &nrfreqs, &complete, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  if(nrfreqs == -1) {
    nrfreqs = original->NrFreqChan;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  new = malloc(sizeof(datafile_definition));
  if(new == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  if(preprocess_addsuccessiveFreqChans(*original, new, nrfreqs, complete, NULL, verbose) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Adding frequency channels failed");
    return NULL;
  }
  return PyCapsule_New(new, "datafile_definition", &destruct_datafile_definition);
}
static PyObject *method_subtract_profile(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original, *profile, *new;
  PyObject *capsule, *capsule2;
  static char *kwlist[] = {"archive", "profile", "mode", "inplace", "verbose", "debug", NULL};
  int mode = 1;
  int inplace = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|ippp", kwlist, &capsule, &capsule2, &mode, &inplace, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  profile = (datafile_definition *)PyCapsule_GetPointer(capsule2, "datafile_definition");
  if(profile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  new = NULL;
  if(inplace == 0) {
    new = malloc(sizeof(datafile_definition));
    if(new == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
      return NULL;
    }
  }
  if(preprocess_subtractprof(original, *profile, mode, new, inplace, verbose) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Subtracting profile failed");
    return NULL;
  }
  PyObject *ret;
  if(inplace != 0) {
    Py_INCREF(Py_None);
    ret = Py_None;
  }else {
    ret = PyCapsule_New(new, "datafile_definition", &destruct_datafile_definition);
  }
  return ret;
}
static PyObject *method_polselect(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original, *new;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "polnr", "verbose", "debug", NULL};
  long polnr = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Ol|pp", kwlist, &capsule, &polnr, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  new = malloc(sizeof(datafile_definition));
  if(new == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  if(preprocess_polselect(*original, new, polnr, verbose) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Polarisation selection failed");
    return NULL;
  }
  return PyCapsule_New(new, "datafile_definition", &destruct_datafile_definition);
}
static PyObject *method_norm(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "normvalue", "globalnorm", "verbose", "debug", NULL};
  float normvalue = 1;
  int verbose_flag = 0;
  int debug_flag = 0;
  int global = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Of|ppp", kwlist, &capsule, &normvalue, &global, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  pulselongitude_regions_definition *onpulse = NULL;
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int normavrg = 0;
  if(preprocess_norm(*original, normvalue, onpulse, global, normavrg, verbose) != 1) {
    PyErr_SetString(PyExc_RuntimeError, "Normalising data failed");
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_stokes(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|pp", kwlist, &capsule, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(preprocess_stokes(original, verbose) != 1) {
    PyErr_SetString(PyExc_RuntimeError, "Conversion to Stokes parameters failed");
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_dedisperse(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "undo", "old_reference_freq", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  int undo = 0;
  double old_reference_freq = sqrt(-1);
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|pdpp", kwlist, &capsule, &undo, &old_reference_freq, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int update = 0;
  if(old_reference_freq != sqrt(-1)) {
    update = 1;
  }
  int ret;
  ret = preprocess_dedisperse(original, undo, update, old_reference_freq, verbose);
  if(ret == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Dedispersion failed");
    return NULL;
  }
  return PyLong_FromLong(ret);
}
static PyObject *method_invariant_interval(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|pp", kwlist, &capsule, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(preprocess_invariant_interval(original, verbose) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Conversion to invariant interval failed");
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_deFarad(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "undo", "old_reference_freq", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  int undo = 0;
  double old_reference_freq = sqrt(-1);
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|pdpp", kwlist, &capsule, &undo, &old_reference_freq, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int update = 0;
  if(old_reference_freq != sqrt(-1)) {
    update = 1;
  }
  int ret;
  ret = preprocess_deFaraday(original, undo, update, old_reference_freq, NULL, verbose);
  if(ret == 0) {
    PyErr_SetString(PyExc_RuntimeError, "De-Faraday rotation failed");
    return NULL;
  }
  return PyLong_FromLong(ret);
}
static PyObject *method_debase(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original;
  PyObject *capsule, *capsule2;
  static char *kwlist[] = {"archive", "onpulse", "slope", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  int slope = 0;
  capsule2 = NULL;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|Oppp", kwlist, &capsule, &capsule2, &slope, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  pulselongitude_regions_definition *regions;
  regions = NULL;
  if(capsule2 != NULL) {
    regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule2, "pulselongitude_regions_definition");
    if(regions == NULL) {
      PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
      return NULL;
    }
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  float **baseline = NULL;
  if(preprocess_debase(original, regions, baseline, slope, verbose) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Removing baseline failed");
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_rotate(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "phaseshift", "subint", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  double phaseshift = 0.0;
  long subint = -1;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od|lpp", kwlist, &capsule, &phaseshift, &subint, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(preprocess_fftshift(*original, subint, phaseshift, 0, 0.0, verbose) != 1) {
    PyErr_SetString(PyExc_RuntimeError, "Rotating data failed");
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_scale(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *original;
  PyObject *capsule;
  static char *kwlist[] = {"archive", "scale", "offset", "polnr", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  double scale = 0.0;
  double offset = 0.0;
  int polnr = -1;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od|dipp", kwlist, &capsule, &scale, &offset, &polnr, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  original = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(original == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(preprocess_scale(*original, scale, offset, polnr, verbose) != 1) {
    PyErr_SetString(PyExc_RuntimeError, "Scaling data failed");
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_regions_select(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *profileI;
  static char *kwlist[] = {"profileI", "onlyOne", "powerTwo", "even", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  int onlyOne = 0;
  int powerTwo = 0;
  int evenNumber = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ppppp", kwlist, &profileI, &onlyOne, &powerTwo, &evenNumber, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  int error_status;
  long nrBins;
  float *profileI_C;
  profileI_C = python_list_to_psrsalsa_array(profileI, -1, NULL, NULL, sizeof(float), &nrBins, &error_status);
  if(profileI_C == NULL) {
    return NULL;
  }
  pgplot_options_definition pgplot;
  pgplot_clear_options(&pgplot);
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  pulselongitude_regions_definition *regions;
  regions = malloc(sizeof(pulselongitude_regions_definition));
  if(regions == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  if(initPulselongitudeRegion(regions, verbose) == 0) {
    PyErr_SetString(PyExc_MemoryError, "Pulse longitude region initialisation failed");
    return NULL;
  }
  selectRegions(profileI_C, nrBins, &pgplot, onlyOne, powerTwo, evenNumber, regions, verbose);
  free(profileI_C);
  return PyCapsule_New(regions, "pulselongitude_regions_definition", &destruct_pulselongitude_regions_definition);
}
static PyObject *method_regions_print(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  pulselongitude_regions_definition *regions;
  static char *kwlist[] = {"regions", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule, "pulselongitude_regions_definition");
  if(regions == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  printRegions(regions, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_regions_int_to_frac(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  pulselongitude_regions_definition *regions;
  static char *kwlist[] = {"regions", "scale", "offset", NULL};
  float scale, offset;
  scale = 1;
  offset = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Off|", kwlist, &capsule, &scale, &offset)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule, "pulselongitude_regions_definition");
  if(regions == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  region_int_to_frac(regions, scale, offset);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_regions_frac_to_int(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  pulselongitude_regions_definition *regions;
  static char *kwlist[] = {"regions", "scale", "offset", NULL};
  float scale, offset;
  scale = 1;
  offset = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Off|", kwlist, &capsule, &scale, &offset)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule, "pulselongitude_regions_definition");
  if(regions == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  region_frac_to_int(regions, scale, offset);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_regions_get_nr_regions(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"regions", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  pulselongitude_regions_definition *regions;
  regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule, "pulselongitude_regions_definition");
  if(regions == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  PyObject *ret;
  ret = PyLong_FromLong(regions->nrRegions);
  return ret;
}
static PyObject *method_regions_check_bin(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"regions", "bin", "regionnr", "verbose", "debug", NULL};
  PyObject *capsule;
  long bin;
  long whichregion = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Ol|lpp", kwlist, &capsule, &bin, &whichregion, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  pulselongitude_regions_definition *regions;
  regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule, "pulselongitude_regions_definition");
  if(regions == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int result;
  result = checkRegions(bin, regions, whichregion, verbose);
  PyObject *ret;
  ret = PyLong_FromLong(result);
  return ret;
}
static PyObject *method_regions_get_binrange(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"regions", "regionnr", NULL};
  PyObject *capsule;
  int regionnr = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi", kwlist, &capsule, &regionnr)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  pulselongitude_regions_definition *regions;
  regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule, "pulselongitude_regions_definition");
  if(regions == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  PyObject *ret;
  int foundinfo = 0;
  if(regionnr < regions->nrRegions && regionnr >= 0) {
    if(regions->bins_defined[regionnr]) {
      int intlist[2];
      intlist[0] = regions->left_bin[regionnr];
      intlist[1] = regions->right_bin[regionnr];
      ret = psrsalsa_int_array_to_python_list(intlist, 2);
      if(ret == NULL) {
 PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
 return NULL;
      }
      foundinfo = 1;
    }
  }
  if(foundinfo == 0) {
    ret = PyList_New(0);
  }
  return ret;
}
static PyObject *method_regions_count_nrbins(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  pulselongitude_regions_definition *regions;
  static char *kwlist[] = {"regions", "debug", NULL};
  int debug_flag;
  debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|p", kwlist, &capsule, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule, "pulselongitude_regions_definition");
  if(regions == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  int nrbins;
  nrbins = count_nrbins_in_Regions(regions, debug_flag);
  return PyLong_FromLong(nrbins);
}
static PyObject *method_regions_add_binrange(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  pulselongitude_regions_definition *regions;
  static char *kwlist[] = {"regions", "left", "right", "debug", NULL};
  int debug_flag;
  int left, right;
  left = 0;
  right = 0;
  debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oii|p", kwlist, &capsule, &left, &right, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule, "pulselongitude_regions_definition");
  if(regions == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  regions->left_bin[regions->nrRegions] = left;
  regions->right_bin[regions->nrRegions] = right;
  regions->bins_defined[regions->nrRegions] = 1;
  (regions->nrRegions)++;
  if(regions->nrRegions == MAX_pulselongitude_regions) {
    PyErr_SetString(PyExc_RuntimeError, "To many regions selected.");
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_regions_create_empty(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "|pp", kwlist, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  pulselongitude_regions_definition *regions;
  regions = malloc(sizeof(pulselongitude_regions_definition));
  if(regions == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  if(initPulselongitudeRegion(regions, verbose) == 0) {
    PyErr_SetString(PyExc_MemoryError, "Pulse longitude region initialisation failed");
    return NULL;
  }
  return PyCapsule_New(regions, "pulselongitude_regions_definition", &destruct_pulselongitude_regions_definition);
}
static PyObject *method_add_archives(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *archive_list;
  static char *kwlist[] = {"archive_list", "nsub", "stokesI", "poladd", "freqadd",
    "verbose", "debug",
    NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  int poladd = 0;
  int freqadd = 0;
  int onlyI = 0;
  long sumNsub = 1;
  char *args_format;
  args_format = "O|lppppp";
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, args_format, kwlist, &archive_list, &sumNsub, &onlyI, &poladd, &freqadd
      , &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  Py_ssize_t nrinputfiles;
  nrinputfiles = PyList_Size(archive_list);
  if(nrinputfiles <= 1) {
    PyErr_SetString(PyExc_ValueError, "Expected at least two input archives.");
    return NULL;
  }
  datafile_definition **fin;
  fin = malloc(nrinputfiles*sizeof(datafile_definition *));
  if(fin == NULL) {
    PyErr_SetString(PyExc_ValueError, "Memory allocation error");
    return NULL;
  }
  PyObject *item;
  Py_ssize_t i;
  for(i = 0; i < nrinputfiles; i++) {
    datafile_definition *single_archive;
    item = PyList_GetItem(archive_list, i);
    if(item == NULL) {
      PyErr_SetString(PyExc_ValueError, "Error obtaining item from list of archives");
      return NULL;
    }
    single_archive = (datafile_definition *)PyCapsule_GetPointer(item, "datafile_definition");
    if(single_archive == NULL) {
      PyErr_SetString(PyExc_TypeError, "Archive list item doesn't appear to be a pulsar data object");
      return NULL;
    }
    fin[i] = single_archive;
  }
  datafile_definition *fout;
  fout = malloc(sizeof(datafile_definition));
  if(fout == NULL) {
    PyErr_SetString(PyExc_ValueError, "Memory allocation error");
    return NULL;
  }
  int oformat = MEMORY_format;
  char *output_fname = NULL;
  int noinput = 1;
  char *PlotDevice = NULL;
  float *shiftedProfile = NULL;
  float *Iprofile_firstfile = NULL;
  int circularShift = 1;
  int shift = 0;
  int memsave = 0;
  int argc = 0;
  char **argv = NULL;
  psrsalsaApplication *application = NULL;
  int history_cmd_only = 0;
  int noclosePSRData = 1;
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int extraquiet;
  if(verbose_flag == 0) {
    verbose.nocounters = 1;
    extraquiet = 1;
  }else {
    verbose.nocounters = 0;
    extraquiet = 0;
  }
  if(padd_check_parameters(poladd, freqadd, sumNsub, circularShift, noinput, shift, memsave,
      verbose) == 0) {
      PyErr_SetString(PyExc_ValueError, "Provided parameters are invalid");
      return NULL;
    }
  if(padd_do_stuff(fin, nrinputfiles, fout, oformat, output_fname,
     poladd, freqadd, noinput, PlotDevice, shiftedProfile, Iprofile_firstfile, onlyI, circularShift, shift, sumNsub,
     memsave, argc, argv, application, history_cmd_only, noclosePSRData, extraquiet, verbose) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "Adding archive list failed");
    return NULL;
  }
  free(fin);
  return PyCapsule_New(fout, "datafile_definition", &destruct_datafile_definition);
}
static PyMethodDef genopsMethods[] = {
  {"add_archives", (PyCFunction) method_add_archives, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"fscr", (PyCFunction) method_fscr, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"make_profile", (PyCFunction) method_make_profile, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"polselect", (PyCFunction) method_polselect, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"rebin", (PyCFunction) method_rebin, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"shuffle", (PyCFunction) method_shuffle, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"tscr", (PyCFunction) method_tscr, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"debase", (PyCFunction) method_debase, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"dedisperse", (PyCFunction) method_dedisperse, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"deFarad", (PyCFunction) method_deFarad, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"invariant_interval", (PyCFunction) method_invariant_interval, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"norm", (PyCFunction) method_norm, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"rotate", (PyCFunction) method_rotate, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"scale", (PyCFunction) method_scale, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"stokes", (PyCFunction) method_stokes, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"subtract_profile", (PyCFunction) method_subtract_profile, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"regions_add_binrange", (PyCFunction) method_regions_add_binrange, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"regions_check_bin", (PyCFunction) method_regions_check_bin, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"regions_count_nrbins", (PyCFunction) method_regions_count_nrbins, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"regions_create_empty", (PyCFunction) method_regions_create_empty, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"regions_get_binrange", (PyCFunction) method_regions_get_binrange, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"regions_get_nr_regions", (PyCFunction) method_regions_get_nr_regions, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"regions_frac_to_int", (PyCFunction) method_regions_frac_to_int, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"regions_int_to_frac", (PyCFunction) method_regions_int_to_frac, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"regions_print", (PyCFunction) method_regions_print, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"regions_select", (PyCFunction) method_regions_select, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
{NULL, NULL, 0, NULL}
};
static struct PyModuleDef genops_module = {
    PyModuleDef_HEAD_INIT,
    "genops",
    "Python interface for PSRSALSA",
    -1,
    genopsMethods
};
PyMODINIT_FUNC PyInit_genops(void)
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
  return PyModule_Create(&genops_module);
}

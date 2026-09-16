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
#include "lib/psrsalsa_array_to_python_list.c"
#include "lib/psrsalsa_array_to_python_2Dlist.c"
#include "lib/python_list_to_psrsalsa_array.c"
#include "lib/python_2Dlist_to_psrsalsa_array.c"
static PyObject *method_io_psrdata_load(PyObject *self, PyObject *args, PyObject *kwargs)
{
  char *filename;
  int ret;
  static char *kwlist[] = {"filename", "format", "nowarnings", "headeronly", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  int nowarnings = 0;
  int headeronly = 0;
  char *format_str = NULL;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s|sippp", kwlist, &filename, &format_str, &nowarnings, &headeronly, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = malloc(sizeof(datafile_definition));
  if(datafile == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  int format = 0;
  if(format_str != NULL) {
    format = parsePSRDataFormats(format_str);
    if(format == 0) {
      PyErr_SetString(PyExc_ValueError, "The specified format cannot be interpretted as a pulsar data format");
      return NULL;
    }
  }
  int enable_write = 0;
  int read_in_memory;
  if(headeronly == 0) {
    read_in_memory = 1;
  }else {
    read_in_memory = 0;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  ret = openPSRData(datafile, filename, format, enable_write, read_in_memory, nowarnings, -1, verbose);
  if(ret == 0) {
    PyErr_SetString(PyExc_FileNotFoundError, "Cannot open file");
    return NULL;
  }
  if(headeronly != 0) {
    int noscales_option = 1;
    if(readHeaderPSRData(datafile, noscales_option, 0, -1, verbose) == 0) {
      PyErr_SetString(PyExc_RuntimeError, "Error reading header");
      return NULL;
    }
    int perserve_header = 1;
    int perserve_data = 0;
    if(closePSRData(datafile, perserve_header, perserve_data, verbose) != 0) {
      PyErr_SetString(PyExc_RuntimeError, "Cannot close file");
      return NULL;
    }
  }
  return PyCapsule_New(datafile, "datafile_definition", &destruct_datafile_definition);
}
static PyObject *method_io_psrdata_save(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  char *filename;
  int ret;
  static char *kwlist[] = {"archive", "filename", "format", "argv", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  char *format_str = NULL;
  PyObject *argv_list = NULL;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|sOpp", kwlist, &capsule, &filename, &format_str, &argv_list, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  int format = 0;
  if(format_str != NULL) {
    format = parsePSRDataFormats(format_str);
  }else {
    format = parsePSRDataFormats("psrfits");
  }
  if(format == 0) {
    PyErr_SetString(PyExc_ValueError, "The specified format cannot be interpretted as a pulsar data format");
    return NULL;
  }
  int enable_write = 1;
  int nowarnings = 0;
  int read_in_memory = 0;
  datafile_definition fout;
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  cleanPSRData(&fout, verbose);
  if(copy_params_PSRData(*datafile, &fout, verbose) == 0) {
    PyErr_SetString(PyExc_FileNotFoundError, "Making a copy of the header parameters failed");
    return NULL;
  }
  ret = openPSRData(&fout, filename, format, enable_write, read_in_memory, nowarnings, -1, verbose);
  if(ret == 0) {
    PyErr_SetString(PyExc_FileNotFoundError, "Cannot open file for writing");
    return NULL;
  }
  Py_ssize_t argc = 0;
  char **argv = NULL;
  if(argv_list != NULL) {
    argc = PyList_GET_SIZE(argv_list);
  }
  argv = malloc((argc+1)*sizeof(char *));
  if(argv == NULL) {
    PyErr_SetString(PyExc_FileNotFoundError, "Memory allocation error");
    return NULL;
  }
  char *first_str = "psrsala_py";
  argv[0] = malloc(strlen(first_str)+1);
  if(argv[0] == NULL) {
    PyErr_SetString(PyExc_FileNotFoundError, "Memory allocation error");
    return NULL;
  }
  strcpy(argv[0], first_str);
  if(argv_list != NULL) {
    Py_ssize_t i;
    for(i = 0; i < argc; i++) {
      PyObject *string;
      string = PyList_GetItem(argv_list, i);
      PyObject *bytesobject;
      bytesobject = PyUnicode_AsASCIIString(string);
      char *cstring;
      cstring = PyBytes_AsString(bytesobject);
      if(cstring == NULL) {
 PyErr_SetString(PyExc_FileNotFoundError, "String conversion failed related to argv argument");
 return NULL;
      }
      argv[i+1] = malloc(strlen(cstring)+1);
      if(argv[i+1] == NULL) {
 PyErr_SetString(PyExc_FileNotFoundError, "Memory allocation error");
 return NULL;
      }
      strcpy(argv[i+1], cstring);
    }
  }
  int cmd_only = 0;
  char *notes = NULL;
  if(writeHeaderPSRData(&fout, argc+1, argv, cmd_only, notes, verbose) != 1) {
    PyErr_SetString(PyExc_FileNotFoundError, "Cannot write header");
    return NULL;
  }
  if(writePSRData(&fout, datafile->data, verbose) != 1) {
    PyErr_SetString(PyExc_FileNotFoundError, "Cannot write data");
    return NULL;
  }
  int perserve_header = 0;
  int perserve_data = 1;
  closePSRData(&fout, perserve_header, perserve_data, verbose);
  int i;
  for(i = 0; i < argc; i++) {
    free(argv[i]);
  }
  free(argv);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_create_empty(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  static char *kwlist[] = {"nrsubints", "nrbins", "nrfreqchans", "nrpols", "header", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  long nrsubints, nrbins;
  long nrfreqs = 1;
  long nrpols = 1;
  capsule = NULL;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "ll|llOpp", kwlist, &nrsubints, &nrbins, &nrfreqs, &nrpols, &capsule, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  datafile_definition *datafile;
  datafile = malloc(sizeof(datafile_definition));
  if(datafile == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  cleanPSRData(datafile, verbose);
  if(capsule != NULL) {
    datafile_definition *datafile_hdr;
    datafile_hdr = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
    if(datafile_hdr == NULL) {
      PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
      return NULL;
    }
    if(copy_params_PSRData(*datafile_hdr, datafile, verbose) == 0) {
      PyErr_SetString(PyExc_RuntimeError, "Copying of header parameters failed");
      return NULL;
    }
    if(datafile->tsampMode == TSAMPMODE_LONGITUDELIST && datafile->tsamp_list != NULL) {
      free(datafile->tsamp_list);
      datafile->tsamp_list = NULL;
      datafile->tsampMode = TSAMPMODE_UNKNOWN;
    }
    if(datafile->tsub_list != NULL) {
      free(datafile->tsub_list);
    }
    datafile->tsub_list = (double *)malloc(sizeof(double));
    if(datafile->tsub_list == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
      return NULL;
    }
    datafile->tsub_list[0] = 0;
    datafile->tsubMode = TSUBMODE_UNKNOWN;
    if(datafile->freqlabel_list != NULL) {
      free(datafile->freqlabel_list);
      datafile->freqlabel_list = NULL;
    }
    datafile->freqMode = FREQMODE_UNKNOWN;
  }
  datafile->data = calloc(nrsubints*nrbins*nrfreqs*nrpols, sizeof(float));
  if(datafile == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  datafile->format = MEMORY_format;
  datafile->NrSubints = nrsubints;
  datafile->NrBins = nrbins;
  datafile->NrFreqChan = nrfreqs;
  datafile->NrPols = nrpols;
  return PyCapsule_New(datafile, "datafile_definition", &destruct_datafile_definition);
}
static PyObject *method_io_printPSRDataFormats(PyObject *self, PyObject *args)
{
  printPSRDataFormats(stdout, 0);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_printHeaderPSRData(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  int update = 0;
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  printHeaderPSRData(*datafile, update, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_PSRData_get_data(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "no_numpy_out", NULL};
  PyObject *capsule;
  int no_numpy = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|p", kwlist, &capsule, &no_numpy)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret;
  size_t datasize;
  datasize = datafile->NrSubints*datafile->NrBins*datafile->NrFreqChan*datafile->NrPols;
  ret = psrsalsa_array_to_python_list(datafile->data, NULL, datasize, no_numpy);
  if(ret == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  return ret;
}
static PyObject *method_io_PSRData_get_pulse_data(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  static char *kwlist[] = {"archive", "subint", "polnr", "freqnr", "no_numpy_out", "verbose", "debug", NULL};
  long subint = 0;
  int polnr = 0;
  int freqnr = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  int no_numpy = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Olii|ppp", kwlist, &capsule, &subint, &polnr, &freqnr, &no_numpy, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  if(subint < 0 || subint >= datafile->NrSubints || freqnr < 0 || freqnr >= datafile->NrFreqChan || polnr < 0 || polnr >= datafile->NrPols) {
    if(subint < 0 || subint >= datafile->NrSubints) {
      printerror(0, "Requested subint %ld unexpected. Total available is %ld", subint, datafile->NrSubints);
    }
    if(freqnr < 0 || freqnr >= datafile->NrFreqChan) {
      printerror(0, "Requested frequency channel %ld unexpected. Total available is %ld", freqnr, datafile->NrFreqChan);
    }
    if(polnr < 0 || polnr >= datafile->NrPols) {
      printerror(0, "Requested polarization channel %ld unexpected. Total available is %ld", polnr, datafile->NrPols);
    }
    PyErr_SetString(PyExc_ValueError, "Selected pulse is outside data range");
    return NULL;
  }
  int binnr = 0;
  float *pulse_ptr;
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(get_pointer_PulsePSRData(datafile, subint, polnr, freqnr, binnr, &pulse_ptr, verbose) != 1) {
    PyErr_SetString(PyExc_TypeError, "Obtaining data failed");
    return NULL;
  }
  PyObject *ret;
  size_t datasize;
  datasize = datafile->NrBins;
  ret = psrsalsa_array_to_python_list(pulse_ptr, NULL, datasize, no_numpy);
  if(ret == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  return ret;
}
static PyObject *method_io_PSRData_set_pulse_data(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  PyObject *pulsedata;
  static char *kwlist[] = {"archive", "pulsedata", "subint", "polnr", "freqnr", "verbose", "debug", NULL};
  long subint = 0;
  int polnr = 0;
  int freqnr = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OOlii|pp", kwlist, &capsule, &pulsedata, &subint, &polnr, &freqnr, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  if(subint < 0 || subint >= datafile->NrSubints || freqnr < 0 || freqnr >= datafile->NrFreqChan || polnr < 0 || polnr >= datafile->NrPols) {
    if(subint < 0 || subint >= datafile->NrSubints) {
      printerror(0, "Requested subint %ld unexpected. Total available is %ld", subint, datafile->NrSubints);
    }
    if(freqnr < 0 || freqnr >= datafile->NrFreqChan) {
      printerror(0, "Requested frequency channel %ld unexpected. Total available is %ld", freqnr, datafile->NrFreqChan);
    }
    if(polnr < 0 || polnr >= datafile->NrPols) {
      printerror(0, "Requested polarization channel %ld unexpected. Total available is %ld", polnr, datafile->NrPols);
    }
    PyErr_SetString(PyExc_ValueError, "Selected pulse is outside data range");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int binnr = 0;
  float *pulse_ptr;
  if(get_pointer_PulsePSRData(datafile, subint, polnr, freqnr, binnr, &pulse_ptr, verbose) != 1) {
    PyErr_SetString(PyExc_TypeError, "Setting data failed");
    return NULL;
  }
  int error_status;
  long nrbins;
  if(python_list_to_psrsalsa_array(pulsedata, datafile->NrBins, pulse_ptr, NULL, 0, &nrbins, &error_status) == NULL) {
    if(error_status == 1) {
      printerror(0, "ERROR method_io_PSRData_set_pulse_data: %ld != %ld.", nrbins, datafile->NrBins);
      PyErr_SetString(PyExc_ValueError, "Provided list with data has a different length than the number of bins in the archive.");
    }
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_PSRData_get_data_as_matrix(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "order", "freqchan", "polchan", "subint", "phasebin", "no_numpy_out", "verbose", "debug", NULL};
  PyObject *capsule;
  char *order;
  long freqchan = -1;
  long polchan = -1;
  long subint = -1;
  long phasebin = -1;
  int verbose_flag = 0;
  int debug_flag = 0;
  int no_numpy = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|llllppp", kwlist, &capsule, &order, &freqchan, &polchan, &subint, &phasebin, &no_numpy, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  datafile_definition *datafile_orig;
  datafile_orig = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile_orig == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  datafile_definition *datafile_use;
  datafile_definition datafile_clone;
  if(polchan >= 0 || freqchan >= 0 || subint >= 0 || phasebin >= 0) {
    if(make_clone(*datafile_orig, &datafile_clone, verbose) == 0) {
      PyErr_SetString(PyExc_RuntimeError, "Making a copy of the data failed.");
      return NULL;
    }
    datafile_use = &datafile_clone;
  }else {
    datafile_use = datafile_orig;
  }
  if(strlen(order) != 2) {
    PyErr_SetString(PyExc_TypeError, "The order variable should be two characters long.");
    return NULL;
  }
  int dim_length[2];
  int transpose = 0;
  if(strcasecmp(order, "SB") == 0 || strcasecmp(order, "BS") == 0) {
    dim_length[0] = datafile_use->NrSubints;
    dim_length[1] = datafile_use->NrBins;
    if(polchan >= 0) {
      datafile_definition clone;
      if(preprocess_polselect(*datafile_use, &clone, polchan, verbose) == 0) {
 PyErr_SetString(PyExc_RuntimeError, "Selecting polarization channel failed.");
 return NULL;
      }
      swap_orig_clone(datafile_use, &clone, verbose);
    }
    if(freqchan >= 0) {
      datafile_definition clone;
      if(preprocess_channelselect(*datafile_use, &clone, freqchan, verbose) == 0) {
 PyErr_SetString(PyExc_RuntimeError, "Selecting frequency channel failed.");
 return NULL;
      }
      swap_orig_clone(datafile_use, &clone, verbose);
    }
    if(subint >= 0) {
      PyErr_SetString(PyExc_RuntimeError, "The subint is an output dimension, so no single subint should be selected.");
      return NULL;
    }
    if(phasebin >= 0) {
      PyErr_SetString(PyExc_RuntimeError, "The phase bin is an output dimension, so no single phase bin should be selected.");
      return NULL;
    }
    if(datafile_use->NrFreqChan != 1) {
      PyErr_SetString(PyExc_RuntimeError, "The input archive should only have 1 frequency channel, or use the freqchan option of this function.");
      return NULL;
    }
    if(datafile_use->NrPols != 1) {
      PyErr_SetString(PyExc_RuntimeError, "The input archive should only have 1 polarization, or use the polchan option of this function.");
      return NULL;
    }
    if(strcasecmp(order, "SB") == 0) {
      transpose = 0;
    }else {
      transpose = 1;
    }
  }else if(strcasecmp(order, "FB") == 0|| strcasecmp(order, "BF") == 0) {
    dim_length[0] = datafile_use->NrFreqChan;
    dim_length[1] = datafile_use->NrBins;
    if(polchan >= 0) {
      datafile_definition clone;
      if(preprocess_polselect(*datafile_use, &clone, polchan, verbose) == 0) {
 PyErr_SetString(PyExc_RuntimeError, "Selecting polarization channel failed.");
 return NULL;
      }
      swap_orig_clone(datafile_use, &clone, verbose);
    }
    if(subint >= 0) {
      datafile_definition clone;
      if(preprocess_pulsesselect(*datafile_use, &clone, subint, 1, verbose) == 0) {
 PyErr_SetString(PyExc_RuntimeError, "Selecting subint failed.");
 return NULL;
      }
      swap_orig_clone(datafile_use, &clone, verbose);
    }
    if(freqchan >= 0) {
      PyErr_SetString(PyExc_RuntimeError, "The frequency channel is an output dimension, so no single frequency channel should be selected.");
      return NULL;
    }
    if(phasebin >= 0) {
      PyErr_SetString(PyExc_RuntimeError, "The phase bin is an output dimension, so no single phase bin should be selected.");
      return NULL;
    }
    if(datafile_use->NrSubints != 1) {
      PyErr_SetString(PyExc_RuntimeError, "The input archive should only have 1 subint, or use the subint option of this function.");
      return NULL;
    }
    if(datafile_use->NrPols != 1) {
      PyErr_SetString(PyExc_RuntimeError, "The input archive should only have 1 polarization, or use the polchan option of this function.");
      return NULL;
    }
    if(strcasecmp(order, "FB") == 0) {
      transpose = 0;
    }else {
      transpose = 1;
    }
  }else if(strcasecmp(order, "FS") == 0|| strcasecmp(order, "SF") == 0) {
    dim_length[0] = datafile_use->NrSubints;
    dim_length[1] = datafile_use->NrFreqChan;
    if(polchan >= 0) {
      datafile_definition clone;
      if(preprocess_polselect(*datafile_use, &clone, polchan, verbose) == 0) {
 PyErr_SetString(PyExc_RuntimeError, "Selecting polarization channel failed.");
 return NULL;
      }
      swap_orig_clone(datafile_use, &clone, verbose);
    }
    if(phasebin >= 0) {
      datafile_definition clone;
      if(preprocess_phasebinselect(*datafile_use, &clone, phasebin, verbose) == 0) {
 PyErr_SetString(PyExc_RuntimeError, "Selecting phase bin failed.");
 return NULL;
      }
      swap_orig_clone(datafile_use, &clone, verbose);
    }
    if(subint >= 0) {
      PyErr_SetString(PyExc_RuntimeError, "The subint is an output dimension, so no single subint should be selected.");
      return NULL;
    }
    if(freqchan >= 0) {
      PyErr_SetString(PyExc_RuntimeError, "The frequency channel is an output dimension, so no single frequency channel should be selected.");
      return NULL;
    }
    if(datafile_use->NrBins != 1) {
      PyErr_SetString(PyExc_RuntimeError, "The input archive should only have 1 phase bin, or use the phasebin option of this function.");
      return NULL;
    }
    if(datafile_use->NrPols != 1) {
      PyErr_SetString(PyExc_RuntimeError, "The input archive should only have 1 polarization, or use the polchan option of this function.");
      return NULL;
    }
    if(strcasecmp(order, "SF") == 0) {
      transpose = 0;
    }else {
      transpose = 1;
    }
  }else {
    PyErr_SetString(PyExc_TypeError, "The order variable is not recognized or supported.");
    return NULL;
  }
  PyObject *ret;
  ret = psrsalsa_array_to_python_2Dlist(datafile_use->data, NULL, dim_length[0], dim_length[1], transpose, no_numpy);
  if(ret == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  if(polchan >= 0 || freqchan >= 0 || subint >= 0 || phasebin >= 0) {
    if(closePSRData(&datafile_clone, 0, 0, verbose) != 0) {
      PyErr_SetString(PyExc_RuntimeError, "Freeing cloned data failed.");
      return NULL;
    }
  }
  return ret;
}
static PyObject *method_io_make_clone(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "verbose", "debug", NULL};
  PyObject *capsule;
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|pp", kwlist, &capsule, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  datafile_definition *clone;
  clone = malloc(sizeof(datafile_definition));
  if(clone == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  if(make_clone(*datafile, clone, verbose) != 1) {
    PyErr_SetString(PyExc_TypeError, "Making a clone of the input data failed.");
    return NULL;
  }
  return PyCapsule_New(clone, "datafile_definition", &destruct_datafile_definition);
}
static PyObject *method_io_get_gentype(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyLong_FromLong(datafile->gentype);
}
static PyObject *method_io_set_gentype(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "gentype", NULL};
  PyObject *capsule;
  int gentype;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi", kwlist, &capsule, &gentype)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  datafile->gentype = gentype;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_format(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyLong_FromLong(datafile->format);
}
static PyObject *method_io_set_format(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "format", NULL};
  PyObject *capsule;
  int format;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi", kwlist, &capsule, &format)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  datafile->format = format;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_nrsubints(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyLong_FromLong(datafile->NrSubints);
}
static PyObject *method_io_get_nrbins(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyLong_FromLong(datafile->NrBins);
}
static PyObject *method_io_get_nrpols(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyLong_FromLong(datafile->NrPols);
}
static PyObject *method_io_get_nrfreqchans(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyLong_FromLong(datafile->NrFreqChan);
}
static PyObject *method_io_get_mjd_start(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyFloat_FromDouble(datafile->mjd_start);
}
static PyObject *method_io_set_mjd_start(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "mjd", NULL};
  PyObject *capsule;
  double mjd;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &capsule, &mjd)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  datafile->mjd_start = mjd;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_ra(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyFloat_FromDouble(datafile->ra);
}
static PyObject *method_io_set_ra(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "ra", NULL};
  PyObject *capsule;
  double ra;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &capsule, &ra)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  datafile->ra = ra;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_dec(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyFloat_FromDouble(datafile->dec);
}
static PyObject *method_io_set_dec(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "dec", NULL};
  PyObject *capsule;
  double dec;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &capsule, &dec)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  datafile->dec = dec;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_xrange(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret;
  if(datafile->xrangeset != 0) {
    int no_numpy = 1;
    ret = psrsalsa_array_to_python_list(datafile->xrange, NULL, 2, no_numpy);
    if(ret == NULL) {
      PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
      return NULL;
    }
  }else {
    ret = PyList_New(0);
  }
  return ret;
}
static PyObject *method_io_get_yrange(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret;
  if(datafile->yrangeset != 0) {
    int no_numpy = 1;
    ret = psrsalsa_array_to_python_list(datafile->yrange, NULL, 2, no_numpy);
    if(ret == NULL) {
      PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
      return NULL;
    }
  }else {
    ret = PyList_New(0);
  }
  return ret;
}
static PyObject *method_io_get_feedtype(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret = NULL;
  switch(datafile->feedtype) {
  case FEEDTYPE_LINEAR:
    ret = PyUnicode_FromString("Linear");
    break;
  case FEEDTYPE_CIRCULAR:
    ret = PyUnicode_FromString("Circular");
    break;
  default:
    ret = PyUnicode_FromString("Unknown");
    break;
  }
  return ret;
}
static PyObject *method_io_set_feedtype(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "feedtype", NULL};
  PyObject *capsule;
  char *string;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Os", kwlist, &capsule, &string)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  if(strcasecmp(string, "Linear") == 0) {
    datafile->feedtype = FEEDTYPE_LINEAR;
  }else if(strcasecmp(string, "Circular") == 0) {
    datafile->feedtype = FEEDTYPE_CIRCULAR;
  }else {
    if(strcasecmp(string, "Unknown") != 0) {
      printwarning(0, "WARNING method_io_set_feedtype: provided string %s is not recognized. Defaulting to \"Unknown\".", string);
    }
    datafile->feedtype = FEEDTYPE_UNKNOWN;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_backend(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret;
  if(datafile->instrument != NULL) {
    ret = PyUnicode_FromString(datafile->instrument);
  }else {
    ret = PyUnicode_FromString("");
  }
  return ret;
}
static PyObject *method_io_set_backend(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "backend", NULL};
  PyObject *capsule;
  char *string;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Os", kwlist, &capsule, &string)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  set_instrument_PSRData(datafile, string, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_psrname(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret;
  if(datafile->psrname != NULL) {
    ret = PyUnicode_FromString(datafile->psrname);
  }else {
    ret = PyUnicode_FromString("");
  }
  return ret;
}
static PyObject *method_io_set_psrname(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "psrname", NULL};
  PyObject *capsule;
  char *string;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Os", kwlist, &capsule, &string)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  set_psrname_PSRData(datafile, string, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_observatory(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret;
  if(datafile->observatory != NULL) {
    ret = PyUnicode_FromString(datafile->observatory);
  }else {
    ret = PyUnicode_FromString("");
  }
  return ret;
}
static PyObject *method_io_set_observatory(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "observatory", NULL};
  PyObject *capsule;
  char *string;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Os", kwlist, &capsule, &string)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  set_observatory_PSRData(datafile, string, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_scanID(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret;
  if(datafile->scanID != NULL) {
    ret = PyUnicode_FromString(datafile->scanID);
  }else {
    ret = PyUnicode_FromString("");
  }
  return ret;
}
static PyObject *method_io_set_scanID(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "scanID", NULL};
  PyObject *capsule;
  char *string;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Os", kwlist, &capsule, &string)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  set_scanID_PSRData(datafile, string, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_observer(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret;
  if(datafile->observer != NULL) {
    ret = PyUnicode_FromString(datafile->observer);
  }else {
    ret = PyUnicode_FromString("");
  }
  return ret;
}
static PyObject *method_io_set_observer(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "observer", NULL};
  PyObject *capsule;
  char *string;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Os", kwlist, &capsule, &string)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  set_observer_PSRData(datafile, string, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_projectID(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret;
  if(datafile->projectID != NULL) {
    ret = PyUnicode_FromString(datafile->projectID);
  }else {
    ret = PyUnicode_FromString("");
  }
  return ret;
}
static PyObject *method_io_set_projectID(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "projectID", NULL};
  PyObject *capsule;
  char *string;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Os", kwlist, &capsule, &string)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  set_projectID_PSRData(datafile, string, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_institute(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret;
  if(datafile->institute != NULL) {
    ret = PyUnicode_FromString(datafile->institute);
  }else {
    ret = PyUnicode_FromString("");
  }
  return ret;
}
static PyObject *method_io_set_institute(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "institute", NULL};
  PyObject *capsule;
  char *string;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Os", kwlist, &capsule, &string)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  set_institute_PSRData(datafile, string, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_centre_freq(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  double freq;
  freq = get_centre_frequency(*datafile, verbose);
  return PyFloat_FromDouble(freq);
}
static PyObject *method_io_set_centre_freq(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "freq", NULL};
  PyObject *capsule;
  double freq;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &capsule, &freq)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  set_centre_frequency(datafile, freq, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_bw(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  double bw;
  bw = get_bandwidth(*datafile, verbose);
  return PyFloat_FromDouble(bw);
}
static PyObject *method_io_set_bw(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "bw", NULL};
  PyObject *capsule;
  double bw;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &capsule, &bw)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  set_bandwidth(datafile, bw, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_dedisp_state(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyLong_FromLong(datafile->isDeDisp);
}
static PyObject *method_io_set_dedisp_state(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "state", NULL};
  PyObject *capsule;
  int state;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi", kwlist, &capsule, &state)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  if(state < -1 || state > 1) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a valid state: -1, 0 or 1 was expected.");
    return NULL;
  }else {
    datafile->isDeDisp = state;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_defarad_state(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyLong_FromLong(datafile->isDeFarad);
}
static PyObject *method_io_set_defarad_state(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "state", NULL};
  PyObject *capsule;
  int state;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi", kwlist, &capsule, &state)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  if(state < -1 || state > 1) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a valid state: -1, 0 or 1 was expected.");
    return NULL;
  }else {
    datafile->isDeFarad = state;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_depar_state(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyLong_FromLong(datafile->isDePar);
}
static PyObject *method_io_set_depar_state(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "state", NULL};
  PyObject *capsule;
  int state;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi", kwlist, &capsule, &state)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  if(state < -1 || state > 1) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a valid state: -1, 0 or 1 was expected.");
    return NULL;
  }else {
    datafile->isDePar = state;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_debase_state(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyLong_FromLong(datafile->isDebase);
}
static PyObject *method_io_set_debase_state(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "state", NULL};
  PyObject *capsule;
  int state;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi", kwlist, &capsule, &state)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  if(state < -1 || state > 1) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a valid state: -1, 0 or 1 was expected.");
    return NULL;
  }else {
    datafile->isDebase = state;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_dm(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyFloat_FromDouble(datafile->dm);
}
static PyObject *method_io_set_dm(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "dm", NULL};
  PyObject *capsule;
  double dm;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &capsule, &dm)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  datafile->dm = dm;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_rm(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyFloat_FromDouble(datafile->rm);
}
static PyObject *method_io_set_rm(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "rm", NULL};
  PyObject *capsule;
  double rm;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &capsule, &rm)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  datafile->rm = rm;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_reference_freq(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  return PyFloat_FromDouble(datafile->freq_ref);
}
static PyObject *method_io_set_reference_freq(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "reference_freq", NULL};
  PyObject *capsule;
  double value;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &capsule, &value)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  datafile->freq_ref = value;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_poltype(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  PyObject *ret = NULL;
  switch(datafile->poltype) {
  case POLTYPE_STOKES:
    ret = PyUnicode_FromString("Stokes");
    break;
  case POLTYPE_COHERENCY:
    ret = PyUnicode_FromString("Coherency");
    break;
  case POLTYPE_ILVPAdPA:
    ret = PyUnicode_FromString("ILVPAdPA");
    break;
  case POLTYPE_PAdPA:
    ret = PyUnicode_FromString("PAdPA");
    break;
  case POLTYPE_ILVPAdPATEldEl:
    ret = PyUnicode_FromString("ILVPAdPATEldEl");
    break;
  default:
    ret = PyUnicode_FromString("Unknown");
    break;
  }
  return ret;
}
static PyObject *method_io_set_poltype(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "poltype", NULL};
  PyObject *capsule;
  char *string;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Os", kwlist, &capsule, &string)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  if(strcasecmp(string, "Stokes") == 0) {
    datafile->poltype = POLTYPE_STOKES;
  }else if(strcasecmp(string, "Coherency") == 0) {
    datafile->poltype = POLTYPE_COHERENCY;
  }else if(strcasecmp(string, "ILVPAdPA") == 0) {
    datafile->poltype = POLTYPE_ILVPAdPA;
  }else if(strcasecmp(string, "PAdPA") == 0) {
    datafile->poltype = POLTYPE_PAdPA;
  }else if(strcasecmp(string, "ILVPAdPATEldEl") == 0) {
    datafile->poltype = POLTYPE_ILVPAdPATEldEl;
  }else {
    if(strcasecmp(string, "Unknown") != 0) {
      printwarning(0, "WARNING method_io_set_poltype: provided string %s is not recognized. Defaulting to \"Unknown\".", string);
    }
    datafile->poltype = POLTYPE_UNKNOWN;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_duration(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  double tobs;
  tobs = get_tobs(*datafile, verbose);
  return PyFloat_FromDouble(tobs);
}
static PyObject *method_io_get_subint_durations(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "no_numpy_out", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  int no_numpy = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|p", kwlist, &capsule, &no_numpy)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  double *list_C;
  list_C = malloc(datafile->NrSubints*sizeof(double));
  if(list_C == NULL) {
    PyErr_SetString(PyExc_TypeError, "Memory allocation error");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  int i;
  for(i = 0; i < datafile->NrSubints; i++) {
    list_C[i] = get_tsub(*datafile, i, verbose);
  }
  PyObject *ret;
  ret = psrsalsa_array_to_python_list(NULL, list_C, datafile->NrSubints, no_numpy);
  if(ret == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  free(list_C);
  return ret;
}
static PyObject *method_io_get_tsamp(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  double tsamp;
  tsamp = get_tsamp(*datafile, 0, verbose);
  return PyFloat_FromDouble(tsamp);
}
static PyObject *method_io_set_tsamp(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "tsamp", NULL};
  PyObject *capsule;
  double tsamp;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &capsule, &tsamp)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  set_tsamp(datafile, tsamp, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_fold_period(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", NULL};
  PyObject *capsule;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &capsule)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  double period;
  get_period(*datafile, 0, &period, verbose);
  return PyFloat_FromDouble(period);
}
static PyObject *method_io_set_fold_period(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "period", NULL};
  PyObject *capsule;
  double period;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Od", kwlist, &capsule, &period)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  set_period(datafile, period, verbose);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_get_subint_freqs(PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"archive", "subint", "no_numpy_out", NULL};
  datafile_definition *datafile;
  PyObject *capsule;
  long subint;
  int no_numpy = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Ol|p", kwlist, &capsule, &subint, &no_numpy)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  if(subint < 0 || subint >= datafile->NrSubints) {
    PyErr_SetString(PyExc_TypeError, "Requested subint doesn't appear to exist");
    return NULL;
  }
  double *list_C;
  list_C = malloc(datafile->NrFreqChan*sizeof(double));
  if(list_C == NULL) {
    PyErr_SetString(PyExc_TypeError, "Memory allocation error");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  int i;
  for(i = 0; i < datafile->NrFreqChan; i++) {
    list_C[i] = get_weighted_channel_freq(*datafile, subint, i, verbose);
  }
  PyObject *ret;
  ret = psrsalsa_array_to_python_list(NULL, list_C, datafile->NrFreqChan, no_numpy);
  if(ret == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  free(list_C);
  return ret;
}
static PyObject *method_io_set_freq_table(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  PyObject *freqtable;
  static char *kwlist[] = {"archive", "freqtable", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|pp", kwlist, &capsule, &freqtable, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(datafile->freqlabel_list != NULL) {
    free(datafile->freqlabel_list);
  }
  datafile->freqlabel_list = malloc(datafile->NrSubints*datafile->NrFreqChan*sizeof(double));
  if(datafile->freqlabel_list == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  datafile->freqMode = FREQMODE_FREQTABLE;
  int error_status;
  double *array;
  long nrx, nry;
  array = python_2Dlist_to_psrsalsa_array(freqtable, 0, datafile->NrSubints, datafile->NrFreqChan, NULL, NULL, sizeof(double), &nrx, &nry, &error_status);
  if(array == NULL) {
    if(error_status == 1) {
      if(nry != datafile->NrSubints) {
 PyErr_SetString(PyExc_ValueError, "Provided list with data has a different length than the number of subints in the archive.");
 return NULL;
      }
      if(nrx != datafile->NrFreqChan) {
 PyErr_SetString(PyExc_ValueError, "Provided list with data has a different length than the number of bins in the archive.");
 return NULL;
      }
    }
    return NULL;
  }
  long chan, subint;
  for(subint = 0; subint < datafile->NrSubints; subint++) {
    for(chan = 0; chan < datafile->NrFreqChan; chan++) {
      if(set_weighted_channel_freq(datafile, subint, chan, array[subint*datafile->NrFreqChan+chan], verbose) == 0) {
 PyErr_SetString(PyExc_ValueError, "Cannot set frequency in header of archive.");
 return NULL;
      }
    }
  }
  free(array);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *method_io_set_tsub_table(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyObject *capsule;
  PyObject *tsubtable;
  static char *kwlist[] = {"archive", "tsubtable", "verbose", "debug", NULL};
  int verbose_flag = 0;
  int debug_flag = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|pp", kwlist, &capsule, &tsubtable, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile_definition *datafile;
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(datafile->tsub_list != NULL) {
    free(datafile->tsub_list);
  }
  datafile->tsub_list = malloc(datafile->NrSubints*sizeof(double));
  if(datafile->tsub_list == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  datafile->tsubMode = TSUBMODE_TSUBLIST;
  long nrSubints;
  int error_status;
  if(python_list_to_psrsalsa_array(tsubtable, datafile->NrSubints, NULL, datafile->tsub_list, 0, &nrSubints, &error_status) == NULL) {
    if(error_status == 1) {
      printerror(0, "ERROR method_io_set_tsub_table: List length %ld != %ld", nrSubints, datafile->NrSubints);
      PyErr_SetString(PyExc_ValueError, "Provided list with data has a different length than the number of subints in the archive.");
    }
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef ioMethods[] = {
  {"psrdata_load", (PyCFunction) method_io_psrdata_load, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"psrdata_save", (PyCFunction) method_io_psrdata_save, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"printPSRDataFormats", method_io_printPSRDataFormats, METH_VARARGS, "Python interface for psrsalsa library function"},
  {"get_data", (PyCFunction) method_io_PSRData_get_data, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_data_as_matrix", (PyCFunction) method_io_PSRData_get_data_as_matrix, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_pulse_data", (PyCFunction) method_io_PSRData_get_pulse_data, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_pulse_data", (PyCFunction) method_io_PSRData_set_pulse_data, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"clone", (PyCFunction) method_io_make_clone, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"create_empty", (PyCFunction) method_io_create_empty, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"printheader", (PyCFunction) method_io_printHeaderPSRData, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_backend", (PyCFunction) method_io_get_backend, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_backend", (PyCFunction) method_io_set_backend, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_bw", (PyCFunction) method_io_get_bw, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_bw", (PyCFunction) method_io_set_bw, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_centre_freq", (PyCFunction) method_io_get_centre_freq, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_centre_freq", (PyCFunction) method_io_set_centre_freq, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_debase_state", (PyCFunction) method_io_get_debase_state, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_debase_state", (PyCFunction) method_io_set_debase_state, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_dec", (PyCFunction) method_io_get_dec, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_dec", (PyCFunction) method_io_set_dec, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_dedisp_state", (PyCFunction) method_io_get_dedisp_state, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_dedisp_state", (PyCFunction) method_io_set_dedisp_state, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_defarad_state", (PyCFunction) method_io_get_defarad_state, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_defarad_state", (PyCFunction) method_io_set_defarad_state, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_depar_state", (PyCFunction) method_io_get_depar_state, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_depar_state", (PyCFunction) method_io_set_depar_state, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_dm", (PyCFunction) method_io_get_dm, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_dm", (PyCFunction) method_io_set_dm, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_duration", (PyCFunction) method_io_get_duration, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_feedtype", (PyCFunction) method_io_get_feedtype, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_feedtype", (PyCFunction) method_io_set_feedtype, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_fold_period", (PyCFunction) method_io_get_fold_period, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_fold_period", (PyCFunction) method_io_set_fold_period, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_format", (PyCFunction) method_io_get_format, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_format", (PyCFunction) method_io_set_format, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_freq_table", (PyCFunction) method_io_set_freq_table, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_gentype", (PyCFunction) method_io_get_gentype, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_gentype", (PyCFunction) method_io_set_gentype, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_institute", (PyCFunction) method_io_get_institute, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_institute", (PyCFunction) method_io_set_institute, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_mjd_start", (PyCFunction) method_io_get_mjd_start, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_mjd_start", (PyCFunction) method_io_set_mjd_start, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_nrbins", (PyCFunction) method_io_get_nrbins, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_nrfreqchans", (PyCFunction) method_io_get_nrfreqchans, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_nrpols", (PyCFunction) method_io_get_nrpols, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_nrsubints", (PyCFunction) method_io_get_nrsubints, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_observatory", (PyCFunction) method_io_get_observatory, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_observatory", (PyCFunction) method_io_set_observatory, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_observer", (PyCFunction) method_io_get_observer, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_observer", (PyCFunction) method_io_set_observer, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_poltype", (PyCFunction) method_io_get_poltype, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_poltype", (PyCFunction) method_io_set_poltype, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_projectID", (PyCFunction) method_io_get_projectID, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_projectID", (PyCFunction) method_io_set_projectID, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_psrname", (PyCFunction) method_io_get_psrname, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_psrname", (PyCFunction) method_io_set_psrname, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_ra", (PyCFunction) method_io_get_ra, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_ra", (PyCFunction) method_io_set_ra, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_reference_freq", (PyCFunction) method_io_get_reference_freq, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_reference_freq", (PyCFunction) method_io_set_reference_freq, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_rm", (PyCFunction) method_io_get_rm, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_rm", (PyCFunction) method_io_set_rm, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_scanID", (PyCFunction) method_io_get_scanID, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_scanID", (PyCFunction) method_io_set_scanID, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_subint_durations", (PyCFunction) method_io_get_subint_durations, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_subint_freqs", (PyCFunction) method_io_get_subint_freqs, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_tsamp", (PyCFunction) method_io_get_tsamp, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_tsamp", (PyCFunction) method_io_set_tsamp, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"set_tsub_table", (PyCFunction) method_io_set_tsub_table, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_xrange", (PyCFunction) method_io_get_xrange, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"get_yrange", (PyCFunction) method_io_get_yrange, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {NULL, NULL, 0, NULL}
};
static struct PyModuleDef io_module = {
    PyModuleDef_HEAD_INIT,
    "io",
    "Python interface for PSRSALSA",
    -1,
    ioMethods
};
PyMODINIT_FUNC PyInit_io(void)
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
  return PyModule_Create(&io_module);
}

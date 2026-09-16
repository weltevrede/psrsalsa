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
#include "lib/psrsalsa_array_to_python_list.c"
#include "lib/psrsalsa_array_to_python_2Dlist.c"
static PyObject *method_lrfs(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *datafile;
  PyObject *capsule, *capsule2;
  unsigned long fft_size;
  static char *kwlist[] = {"archive", "fft_size", "modindex", "stddev", "phaseTrack", "subpulseAmp", "phaseTrackOffsets", "onpulse", "fmin", "fmax", "no_numpy_out", "verbose", "debug", NULL};
  int calcPhaseTrack = 0;
  int calcsubpulseAmplitude = 0;
  int phaseTrackOffsets = 0;
  float freq_min = -1;
  float freq_max = -1;
  int getmodindex = 0;
  int getstddev = 0;
  int no_numpy = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  capsule2 = NULL;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "Ok|pppppOffppp", kwlist, &capsule, &fft_size, &getmodindex, &getstddev, &calcPhaseTrack, &calcsubpulseAmplitude, &phaseTrackOffsets, &capsule2, &freq_min, &freq_max, &no_numpy, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  pulselongitude_regions_definition *regions = NULL;
  if(capsule2 != NULL) {
    regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule2, "pulselongitude_regions_definition");
    if(regions == NULL) {
      PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
      return NULL;
    }
  }
  if(calcPhaseTrack != 0 || calcsubpulseAmplitude != 0) {
    if(freq_min < 0 || freq_max < 0 || freq_max > 0.5 || freq_max > 0.5) {
      PyErr_SetString(PyExc_ValueError, "To calculate subpulse phase / amplitude, please set the fmin and fmax keywords");
      return NULL;
    }
    if(regions == NULL) {
      PyErr_SetString(PyExc_ValueError, "To calculate subpulse phase / amplitude, please set the onpulse keyword");
      return NULL;
    }
  }
  float *lrfs;
  size_t datasize_lrfs = datafile->NrBins * (fft_size/2+1);
  lrfs = malloc(datasize_lrfs*sizeof(float));
  if(lrfs == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  float *phase_track = NULL;
  size_t datasize_phase_track = datafile->NrBins;
  if(calcPhaseTrack != 0) {
    phase_track = malloc(datasize_phase_track*sizeof(float));
    if(phase_track == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
      return NULL;
    }
  }
  float *subpulseAmplitude = NULL;
  size_t datasize_subpulseAmplitude = datafile->NrBins;
  if(calcsubpulseAmplitude != 0) {
    subpulseAmplitude = malloc(datasize_subpulseAmplitude*sizeof(float));
    if(subpulseAmplitude == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
      return NULL;
    }
  }
  float *phase_track_phases = NULL;
  size_t datasize_phaseTrackOffsets = 0;
  if(phaseTrackOffsets != 0) {
    datasize_phaseTrackOffsets = datafile->NrSubints/fft_size;
    phase_track_phases = malloc(datasize_phaseTrackOffsets*sizeof(float));
    if(phase_track_phases == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
      return NULL;
    }
  }
  float *modindex, *stddev, *profileI;
  modindex = NULL;
  stddev = NULL;
  profileI = NULL;
  if(getmodindex != 0 || getstddev != 0) {
    modindex = malloc(datasize_phase_track*sizeof(float));
    stddev = malloc(datasize_phase_track*sizeof(float));
    profileI = malloc(datasize_phase_track*sizeof(float));
    if(modindex == NULL || stddev == NULL || profileI == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
      return NULL;
    }
  }
  int subtractDC = 1;
  float *avrg_offpulse_lrfs_power = NULL;
  int track_only_first_region = 0;
  int mask_freqs = 0;
  int inverseFFT = 0;
  int argc = 0;
  char **argv = NULL;
  float var_rms;
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  if(calcLRFS(datafile->data, datafile->NrSubints, datafile->NrBins, fft_size, lrfs, subtractDC, avrg_offpulse_lrfs_power, phase_track, phase_track_phases, calcPhaseTrack, freq_min, freq_max, track_only_first_region, subpulseAmplitude, calcsubpulseAmplitude, mask_freqs, inverseFFT, regions, &var_rms, argc, argv, verbose) == 0) {
    PyErr_SetString(PyExc_FileNotFoundError, "Cannot calculate LRFS");
    return NULL;
  }
  if(getmodindex != 0 || getstddev != 0) {
    long junk_int;
    junk_int = datafile->NrSubints/fft_size;
    junk_int *= fft_size;
    verbose_definition noverbose;
    cleanVerboseState(&noverbose);
    read_partprofilePSRData(*datafile, profileI, NULL, 0, 0, junk_int, noverbose);
    calcModindex(lrfs, profileI, datafile->NrBins, fft_size, junk_int, stddev, NULL, modindex, NULL, NULL, var_rms, NULL, NULL, 0, verbose);
  }
  int nritems = 0;
  PyObject *lrfs_py;
  lrfs_py = psrsalsa_array_to_python_2Dlist(lrfs, NULL, (fft_size/2+1), datafile->NrBins, 0, no_numpy);
  if(lrfs_py == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  nritems++;
  free(lrfs);
  PyObject *phase_track_py = NULL;
  if(calcPhaseTrack != 0) {
    phase_track_py = psrsalsa_array_to_python_list(phase_track, NULL, datasize_phase_track, no_numpy);
    if(phase_track_py == NULL) {
      PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
      return NULL;
    }
    nritems++;
    free(phase_track);
  }
  PyObject *subpulseAmplitude_py = NULL;
  if(calcsubpulseAmplitude != 0) {
    subpulseAmplitude_py = psrsalsa_array_to_python_list(subpulseAmplitude, NULL, datasize_subpulseAmplitude, no_numpy);
    if(subpulseAmplitude_py == NULL) {
      PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
      return NULL;
    }
    nritems++;
    free(subpulseAmplitude);
  }
  PyObject *phase_track_phases_py = NULL;
  if(phaseTrackOffsets != 0) {
    phase_track_phases_py = psrsalsa_array_to_python_list(phase_track_phases, NULL, datasize_phaseTrackOffsets, no_numpy);
    if(phase_track_phases_py == NULL) {
      PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
      return NULL;
    }
    nritems++;
    free(phase_track_phases);
  }
  PyObject *modindex_py = NULL;
  if(getmodindex != 0) {
    modindex_py = psrsalsa_array_to_python_list(modindex, NULL, datasize_phase_track, no_numpy);
    if(modindex_py == NULL) {
      PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
      return NULL;
    }
    nritems++;
    free(modindex);
  }
  PyObject *stddev_py = NULL;
  if(getstddev != 0) {
    stddev_py = psrsalsa_array_to_python_list(stddev, NULL, datasize_phase_track, no_numpy);
    if(stddev_py == NULL) {
      PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
      return NULL;
    }
    nritems++;
    free(stddev);
  }
  if(getmodindex != 0 || getstddev != 0) {
    free(profileI);
  }
  PyObject *ret;
  ret = PyTuple_New(nritems);
  if(ret == NULL) {
    PyErr_SetString(PyExc_ValueError, "Cannot create return tuple");
    return NULL;
  }
  int itempos = 0;
  if(PyTuple_SetItem(ret, itempos++, lrfs_py) != 0) {
    return NULL;
  }
  if(calcPhaseTrack != 0) {
    if(PyTuple_SetItem(ret, itempos++, phase_track_py) != 0) {
      return NULL;
    }
  }
  if(calcsubpulseAmplitude != 0) {
    if(PyTuple_SetItem(ret, itempos++, subpulseAmplitude_py) != 0) {
      return NULL;
    }
  }
  if(phaseTrackOffsets != 0) {
    if(PyTuple_SetItem(ret, itempos++, phase_track_phases_py) != 0) {
      return NULL;
    }
  }
  if(getmodindex != 0) {
    if(PyTuple_SetItem(ret, itempos++, modindex_py) != 0) {
      return NULL;
    }
  }
  if(getstddev != 0) {
    if(PyTuple_SetItem(ret, itempos++, stddev_py) != 0) {
      return NULL;
    }
  }
  return ret;
}
static PyObject *method_twodfs(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *datafile;
  PyObject *capsule, *capsule2;
  unsigned long fft_size;
  static char *kwlist[] = {"archive", "fft_size", "onpulse", "regionnr", "flip", "noisesub",
    "subtractDC", "no_numpy_out", "verbose", "debug", NULL};
  capsule2 = NULL;
  int regionnr = 0;
  int noisesub = 1;
  int spectype = 0;
  int flip = 1;
  int no_numpy = 0;
  int verbose_flag = 0;
  int debug_flag = 0;
  int subtractDC = 1;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OkO|ipiipppp", kwlist, &capsule, &fft_size, &capsule2, &regionnr, &flip, &noisesub,
      &subtractDC, &no_numpy, &verbose_flag, &debug_flag)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  pulselongitude_regions_definition *regions = NULL;
  regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule2, "pulselongitude_regions_definition");
  if(regions == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  unsigned long nr_onpulse_bins = 0;
  if(regionnr < regions->nrRegions && regionnr >= 0) {
    if(regions->bins_defined[regionnr]) {
      nr_onpulse_bins = regions->right_bin[regionnr] - regions->left_bin[regionnr] + 1;
    }
  }
  float *twodfs;
  size_t datasize = nr_onpulse_bins * (fft_size/2+1);
  twodfs = malloc(datasize*sizeof(float));
  if(twodfs == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
    return NULL;
  }
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.verbose = verbose_flag;
  verbose.debug = debug_flag;
  int allow_cat_offpulse = 0;
  float *var_rms = NULL;
  if(calc2DFS(datafile->data, datafile->NrSubints, datafile->NrBins, fft_size, twodfs, regions, regionnr, allow_cat_offpulse, noisesub, var_rms, spectype, subtractDC, verbose) != 1) {
    PyErr_SetString(PyExc_FileNotFoundError, "Cannot calculate 2DFS");
    return NULL;
  }
  if(flip != 0) {
    float *twodfs_new;
    datasize = nr_onpulse_bins * (fft_size/2+1);
    twodfs_new = malloc(datasize*sizeof(float));
    if(twodfs_new == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory");
      return NULL;
    }
    long unsigned x, y, index, index_new;
    for(y = 0; y < fft_size/2+1; y++) {
      for(x = 0; x < nr_onpulse_bins; x++) {
 index = nr_onpulse_bins*y+x;
 index_new = nr_onpulse_bins*y+(nr_onpulse_bins-1-x);
 twodfs_new[index_new] = twodfs[index];
      }
    }
    memcpy(twodfs, twodfs_new, datasize*sizeof(float));
    free(twodfs_new);
  }
  PyObject *twodfs_py;
  twodfs_py = psrsalsa_array_to_python_2Dlist(twodfs, NULL, fft_size/2+1, nr_onpulse_bins, 0, no_numpy);
  if(twodfs_py == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  free(twodfs);
  return twodfs_py;
}
static PyObject *method_twodfs_xrange(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *datafile;
  PyObject *capsule, *capsule2;
  static char *kwlist[] = {"archive", "onpulse", "regionnr", "flip", NULL};
  capsule2 = NULL;
  int regionnr = 0;
  int flip = 1;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|ip", kwlist, &capsule, &capsule2, &regionnr, &flip)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  pulselongitude_regions_definition *regions = NULL;
  regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule2, "pulselongitude_regions_definition");
  if(regions == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  if(regionnr >= regions->nrRegions) {
    PyErr_SetString(PyExc_TypeError, "Provided onpulse region selection does not include the requested region number to be used in the calculation.");
    return NULL;
  }
  float xrange[2], sign;
  int index1, index2;
  if(flip == 0) {
    index1 = 0;
    index2 = 1;
    sign = 1;
  }else {
    index1 = 1;
    index2 = 0;
    sign = -1;
  }
  xrange[index1] = sign*(-datafile->NrBins/2.0);
  xrange[index2] = sign*(-datafile->NrBins/2.0 +datafile->NrBins*(regions->right_bin[regionnr]-regions->left_bin[regionnr])/(float)(regions->right_bin[regionnr]-regions->left_bin[regionnr]+1));
  PyObject *xrange_py;
  int no_numpy = 1;
  xrange_py = psrsalsa_array_to_python_list(xrange, NULL, 2, no_numpy);
  if(xrange_py == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  return xrange_py;
}
static PyObject *method_twodfs_mplextent(PyObject *self, PyObject *args, PyObject *kwargs)
{
  datafile_definition *datafile;
  PyObject *capsule, *capsule2;
  unsigned long fft_size;
  static char *kwlist[] = {"archive", "fft_size", "onpulse", "regionnr", "flip", NULL};
  capsule2 = NULL;
  int regionnr = 0;
  int flip = 1;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OkO|ip", kwlist, &capsule, &fft_size, &capsule2, &regionnr, &flip)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  datafile = (datafile_definition *)PyCapsule_GetPointer(capsule, "datafile_definition");
  if(datafile == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulsar data object");
    return NULL;
  }
  pulselongitude_regions_definition *regions = NULL;
  regions = (pulselongitude_regions_definition *)PyCapsule_GetPointer(capsule2, "pulselongitude_regions_definition");
  if(regions == NULL) {
    PyErr_SetString(PyExc_TypeError, "Function argument doesn't appear to be a pulse longitude definition object");
    return NULL;
  }
  if(regionnr >= regions->nrRegions) {
    PyErr_SetString(PyExc_TypeError, "Provided onpulse region selection does not include the requested region number to be used in the calculation.");
    return NULL;
  }
  float xrange[2], extent[4];
  xrange[0] = (-datafile->NrBins/2.0);
  xrange[1] = (-datafile->NrBins/2.0 +datafile->NrBins*(regions->right_bin[regionnr]-regions->left_bin[regionnr])/(float)(regions->right_bin[regionnr]-regions->left_bin[regionnr]+1));
  double xbinsize, ybinsize;
  unsigned long nronpulsebins, fft_size_2;
  nronpulsebins = regions->right_bin[regionnr]-regions->left_bin[regionnr] + 1;
  xbinsize = (xrange[1]-xrange[0])/(float)(nronpulsebins-1);
  fft_size_2 = fft_size/2;
  ybinsize = 0.5/(float)fft_size_2;
  if(flip == 0) {
    extent[0] = xrange[0]-0.5*xbinsize;
    extent[1] = xrange[1]+0.5*xbinsize;
  }else {
    extent[1] = -(xrange[0]-0.5*xbinsize);
    extent[0] = -(xrange[1]+0.5*xbinsize);
  }
  extent[2] = 0-0.5*ybinsize;
  extent[3] = 0.5+0.5*ybinsize;
  PyObject *extent_py;
  int no_numpy = 1;
  extent_py = psrsalsa_array_to_python_list(extent, NULL, 4, no_numpy);
  if(extent_py == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  return extent_py;
}
static PyObject *method_lrfs_mplextent(PyObject *self, PyObject *args, PyObject *kwargs)
{
  unsigned long fft_size;
  static char *kwlist[] = {"fft_size", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "k", kwlist, &fft_size)) {
    PyErr_SetString(PyExc_ValueError, "Cannot interpret arguments of function call");
    return NULL;
  }
  float extent[2];
  double ybinsize;
  unsigned long fft_size_2;
  fft_size_2 = fft_size/2;
  ybinsize = 0.5/(float)fft_size_2;
  extent[0] = 0.0-0.5*ybinsize;
  extent[1] = 0.5+0.5*ybinsize;
  PyObject *extent_py;
  int no_numpy = 1;
  extent_py = psrsalsa_array_to_python_list(extent, NULL, 2, no_numpy);
  if(extent_py == NULL) {
    PyErr_SetString(PyExc_ValueError, "Copy to python list failed");
    return NULL;
  }
  return extent_py;
}
static PyMethodDef specMethods[] = {
  {"lrfs", (PyCFunction) method_lrfs, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"lrfs_mplextent", (PyCFunction) method_lrfs_mplextent, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"twodfs", (PyCFunction) method_twodfs, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"twodfs_xrange", (PyCFunction) method_twodfs_xrange, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {"twodfs_mplextent", (PyCFunction) method_twodfs_mplextent, METH_VARARGS | METH_KEYWORDS, "Python interface for psrsalsa library function"},
  {NULL, NULL, 0, NULL}
};
static struct PyModuleDef spec_module = {
    PyModuleDef_HEAD_INIT,
    "spec",
    "Python interface for PSRSALSA",
    -1,
    specMethods
};
PyMODINIT_FUNC PyInit_spec(void)
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
  return PyModule_Create(&spec_module);
}

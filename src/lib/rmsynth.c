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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "psrsalsa.h"
int rmSynthesis(datafile_definition data, float rm_low, float rm_high, float **rmsynth_array, int nrrmsteps, pulselongitude_regions_definition *onpulse, verbose_definition verbose)
{
  long f, n, b;
  int i, ok;
  float rm, *argument_cos_lookup_table, *argument_sin_lookup_table, pc_real, pc_imag, exp_real, exp_imag, spectrum_real, spectrum_imag;
  pulselongitude_regions_definition onpulse_converted;
  if(onpulse != NULL) {
    if(initPulselongitudeRegion(&onpulse_converted, verbose) == 0) {
      printerror(verbose.debug, "ERROR rmSynthesis: Initialising onpulse region failed");
      return 0;
    }
    copyPulselongitudeRegion(*onpulse, &onpulse_converted);
    region_frac_to_int(&onpulse_converted, data.NrBins, 0);
  }
  if(data.isDeFarad == 1) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
    }
    printerror(verbose.debug, "ERROR rmSynthesis: Data in %s is already de-Faraday rotated, so cannot measure RM.", data.filename);
    return 0;
  }
  if(data.isDeFarad == -1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING rmSynthesis (%s): unknown de-Faraday rotation state. Assume data is not yet de-Faraday rotated.", data.filename);
    data.isDeFarad = 0;
  }
  if(data.NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR rmSynthesis (%s): Expected four polarization channels.", data.filename);
    return 0;
  }
  if(data.NrSubints != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR rmSynthesis (%s): Expected all subintegrations to be summed.", data.filename);
    return 0;
  }
  if(data.poltype != POLTYPE_STOKES) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR rmSynthesis (%s): Expected Stokes parameters, but poltype = %d != %d.", data.filename, data.poltype, POLTYPE_STOKES);
    return 0;
  }
  if(data.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR rmSynthesis (%s): only works if data is loaded into memory.", data.filename);
    return 0;
  }
  if(data.poltype == POLTYPE_ILVPAdPA || data.poltype == POLTYPE_PAdPA || data.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR rmSynthesis (%s): Cannot handle PA data.", data.filename);
    return 0;
  }
  if(data.isDebase == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR rmSynthesis (%s): Please remove baseline first, i.e. use pmod -debase.", data.filename);
    return 0;
  }else if(data.isDebase != 1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING rmSynthesis (%s): Unknown baseline state. It is assumed the baseline has already removed from the data.", data.filename);
  }
  if(data.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING rmSynthesis (%s): The frequency channels are possiby non-uniformly separated. The effect on the result is not necessarily understood.", data.filename);
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Determining RM synthesis spectrum from RM %f to %f in %d steps.\n", rm_low, rm_high, nrrmsteps);
  }
  if(*rmsynth_array == NULL) {
    *rmsynth_array = (float *)malloc(2*nrrmsteps*data.NrBins*sizeof(float));
    if(*rmsynth_array == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR rmSynthesis (%s): Cannot allocate memory.", data.filename);
      return 0;
    }
  }
  argument_cos_lookup_table = (float *)malloc(nrrmsteps*data.NrFreqChan*sizeof(float));
  argument_sin_lookup_table = (float *)malloc(nrrmsteps*data.NrFreqChan*sizeof(float));
  if(argument_cos_lookup_table == NULL || argument_sin_lookup_table == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR rmSynthesis (%s): Cannot allocate memory.", data.filename);
    return 0;
  }
  float wavelength2, argument;
  for(f = 0; f < data.NrFreqChan; f++) {
    wavelength2 = 299.792458/get_weighted_channel_freq(data, 0, f, verbose);
    wavelength2 *= wavelength2;
    for(n = 0; n < nrrmsteps; n++) {
      if(nrrmsteps == 1)
 rm = 0.5*(rm_low+rm_high);
      else
 rm = rm_low+n*(rm_high-rm_low)/(float)(nrrmsteps-1);
      argument = -2.0*rm*wavelength2;
      argument_cos_lookup_table[n*data.NrFreqChan+f] = cos(argument);
      argument_sin_lookup_table[n*data.NrFreqChan+f] = sin(argument);
    }
  }
  for(n = 0; n < nrrmsteps; n++) {
    for(b = 0; b < data.NrBins; b++) {
      (*rmsynth_array)[2*(n*data.NrBins+b)] = 0;
      (*rmsynth_array)[2*(n*data.NrBins+b)+1] = 0;
      ok = 0;
      if(onpulse == NULL) {
 ok = 1;
      }else {
 if(checkRegions(b, &onpulse_converted, 0, verbose) != 0 || onpulse_converted.nrRegions == 0) {
   ok = 1;
 }
      }
      if(ok) {
 spectrum_real = 0;
 spectrum_imag = 0;
 for(f = 0; f < data.NrFreqChan; f++) {
   pc_real = data.data[data.NrBins*(1+data.NrPols*f)+b];
   pc_imag = data.data[data.NrBins*(2+data.NrPols*f)+b];
   exp_real = argument_cos_lookup_table[n*data.NrFreqChan+f];
   exp_imag = argument_sin_lookup_table[n*data.NrFreqChan+f];
   spectrum_real += pc_real*exp_real - pc_imag*exp_imag;
   spectrum_imag += pc_real*exp_imag + pc_imag*exp_real;
 }
 (*rmsynth_array)[2*(n*data.NrBins+b)] = sqrt(spectrum_real*spectrum_real+spectrum_imag*spectrum_imag);
 (*rmsynth_array)[2*(n*data.NrBins+b)+1] = atan2(spectrum_imag, spectrum_real);
      }
    }
  }
  free(argument_cos_lookup_table);
  free(argument_sin_lookup_table);
  if(onpulse != NULL) {
    freePulselongitudeRegion(&onpulse_converted);
  }
  return 1;
}
void collapseRMSynthesisArray(float *rmsynth_array, int nrrmsteps, int nrBins, pulselongitude_regions_definition onpulse, float *singlespectrum, verbose_definition verbose)
{
  int b, n;
  region_frac_to_int(&onpulse, nrBins, 0);
  for(n = 0; n < nrrmsteps; n++) {
    singlespectrum[n] = 0;
    for(b = 0; b < nrBins; b++) {
      if(checkRegions(b, &onpulse, 0, verbose) != 0 || onpulse.nrRegions == 0) {
 singlespectrum[n] += rmsynth_array[2*(n*nrBins+b)];
      }
    }
  }
}
int rmSynthesis_instrument_responds(int nrFreqChan, double chanbw, double cfreq0, double rm_low, double rm_high, float **rmsynth_responds, double **rmsynth_responds_double, int usedouble, int nrrmsteps, double rmshift, int callmulti, verbose_definition verbose)
{
  long f, n;
  int i;
  double freq, rm, argument, exp_real, exp_imag, spectrum_real, spectrum_imag, norm;
  static double *wavelength2_lookup_table = NULL;
  static double *rm_lookup_table = NULL;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Determining RM synthesis instrument responds from RM %lf to %lf in %d steps.\n", rm_low, rm_high, nrrmsteps);
  }
  if(usedouble == 0) {
    if(*rmsynth_responds == NULL) {
      *rmsynth_responds = (float *)malloc(nrrmsteps*sizeof(float));
      if(*rmsynth_responds == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR rmSynthesis_instrument_responds: Cannot allocate memory.");
 return 0;
      }
    }
  }else {
    if(*rmsynth_responds_double == NULL) {
      *rmsynth_responds_double = (double *)malloc(nrrmsteps*sizeof(double));
      if(*rmsynth_responds_double == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR rmSynthesis_instrument_responds: Cannot allocate memory.");
 return 0;
      }
    }
  }
  if(callmulti == 0 || callmulti == 2) {
    if(wavelength2_lookup_table != NULL) {
      free(wavelength2_lookup_table);
      wavelength2_lookup_table = NULL;
    }
    if(rm_lookup_table != NULL) {
      free(rm_lookup_table);
      rm_lookup_table = NULL;
    }
  }
  if(wavelength2_lookup_table == NULL) {
    wavelength2_lookup_table = (double *)malloc(nrFreqChan*sizeof(double));
    if(wavelength2_lookup_table == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR rmSynthesis_instrument_responds: Cannot allocate memory.");
      return 0;
    }
    for(f = 0; f < nrFreqChan; f++) {
      freq = cfreq0+chanbw*f;
      wavelength2_lookup_table[f] = 299.792458/freq;
      wavelength2_lookup_table[f] *= wavelength2_lookup_table[f];
    }
  }
  if(rm_lookup_table == NULL) {
    rm_lookup_table = (double *)malloc(nrrmsteps*sizeof(double));
    if(rm_lookup_table == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR rmSynthesis_instrument_responds: Cannot allocate memory.");
      return 0;
    }
    for(n = 0; n < nrrmsteps; n++) {
      if(nrrmsteps == 1)
 rm_lookup_table[n] = 0.5*(rm_low+rm_high);
      else
 rm_lookup_table[n] = rm_low+n*(rm_high-rm_low)/(double)(nrrmsteps-1);
    }
  }
  norm = 1.0/(double)nrFreqChan;
  for(n = 0; n < nrrmsteps; n++) {
    rm = rm_lookup_table[n] - rmshift;
    spectrum_real = 0;
    spectrum_imag = 0;
    for(f = 0; f < nrFreqChan; f++) {
      argument = -2.0*rm*wavelength2_lookup_table[f];
      exp_real = cos(argument);
      exp_imag = sin(argument);
      spectrum_real += exp_real;
      spectrum_imag += exp_imag;
    }
    if(usedouble)
      (*rmsynth_responds_double)[n] = norm * sqrt(spectrum_real*spectrum_real+spectrum_imag*spectrum_imag);
    else
      (*rmsynth_responds)[n] = norm * sqrt(spectrum_real*spectrum_real+spectrum_imag*spectrum_imag);
  }
  if(callmulti == 0 || callmulti == 2) {
    if(wavelength2_lookup_table != NULL) {
      free(wavelength2_lookup_table);
      wavelength2_lookup_table = NULL;
    }
    if(rm_lookup_table != NULL) {
      free(rm_lookup_table);
      rm_lookup_table = NULL;
    }
  }
  return 1;
}
struct {
  float *singlespectrum;
  double rmmin, rmmax, chanbw, cfreq0;
  int nrrmsteps;
  int nrFreqChan;
  double *rmsynth_responds;
  int callmulti;
}internal_internal_fitInstrumentalResponds;
double funk_fitInstrumentalResponds(double *xfit)
{
  int i;
  double chi2, y, dy;
  verbose_definition noverbose;
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;
  rmSynthesis_instrument_responds(internal_internal_fitInstrumentalResponds.nrFreqChan, internal_internal_fitInstrumentalResponds.chanbw, internal_internal_fitInstrumentalResponds.cfreq0, internal_internal_fitInstrumentalResponds.rmmin, internal_internal_fitInstrumentalResponds.rmmax, NULL, &(internal_internal_fitInstrumentalResponds.rmsynth_responds), 1, internal_internal_fitInstrumentalResponds.nrrmsteps, xfit[0], internal_internal_fitInstrumentalResponds.callmulti, noverbose);
  internal_internal_fitInstrumentalResponds.callmulti = 1;
  chi2 = 0;
  for(i = 0; i < internal_internal_fitInstrumentalResponds.nrrmsteps; i++) {
    y = internal_internal_fitInstrumentalResponds.rmsynth_responds[i];
    y *= xfit[1];
    y += xfit[2];
    dy = y-internal_internal_fitInstrumentalResponds.singlespectrum[i];
    chi2 += dy*dy;
  }
  return chi2;
}
int rmSynthesis_fitInstrumentalResponds(float *singlespectrum, float rmmin, float rmmax, int nrrmsteps, float *rm, float *offset, float *scale, int nrFreqChan, float chanbw, float cfreq0, float ftol, verbose_definition verbose)
{
  double xstart[3], dx[3], xfit[3], chi2;
  int fixed[3], nfunk, i, ret;
  internal_internal_fitInstrumentalResponds.singlespectrum = singlespectrum;
  internal_internal_fitInstrumentalResponds.rmmin = rmmin;
  internal_internal_fitInstrumentalResponds.rmmax = rmmax;
  internal_internal_fitInstrumentalResponds.nrrmsteps = nrrmsteps;
  internal_internal_fitInstrumentalResponds.chanbw = chanbw;
  internal_internal_fitInstrumentalResponds.cfreq0 = cfreq0;
  internal_internal_fitInstrumentalResponds.nrFreqChan = nrFreqChan;
  internal_internal_fitInstrumentalResponds.rmsynth_responds = NULL;
  internal_internal_fitInstrumentalResponds.callmulti = 2;
  xstart[0] = (rmmin+rmmax)*0.5;
  xstart[1] = 0;
  xstart[2] = 0;
  for(i = 0; i < nrrmsteps; i++) {
    if(singlespectrum[i] > xstart[1])
      xstart[1] = singlespectrum[i];
  }
  dx[0] = (rmmax-rmmin)*0.05;
  dx[1] = 0.05*xstart[1];
  dx[2] = dx[1];
  fixed[0] = 0;
  fixed[1] = 0;
  fixed[2] = 0;
  ret = doAmoeba_d(0, xstart, dx, fixed, xfit, &chi2, 3, funk_fitInstrumentalResponds, ftol, &nfunk, verbose.verbose, 0, 0, NULL, NULL);
  if(ret == 0) {
    *rm = xfit[0];
    *scale = xfit[1];
    *offset = xfit[2];
    return 0;
  }else {
    if(ret == 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR rmSynthesis_fitInstrumentalResponds: Maximum number of itterations exceeded (maybe try smaller ftol)");
    }else if(ret == 2) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR rmSynthesis_fitInstrumentalResponds: Memory allocation error");
    }else if(ret == 3) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR rmSynthesis_fitInstrumentalResponds: Nr of parameters to fit for should be at least 1");
    }else if(ret == 4) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR rmSynthesis_fitInstrumentalResponds: Specified fit algorithm is not available");
    }
    return 1;
  }
}

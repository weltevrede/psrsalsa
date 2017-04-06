/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include "psrsalsa.h"


void print_fftw_version_used(FILE *stream)
{
  fprintf(stream, "%s (library)", fftwf_version);
}

int rotateSinglepulse(float *data, int npts, float epsilon, verbose_definition verbose)
{
  int i, npts2;
  float fac, dtheta;
  fftwf_complex *dataFFT;
  fftwf_plan plan1, plan2;

  npts2 = npts/2+1;
  dataFFT = (fftwf_complex *)fftwf_malloc(npts2*sizeof(fftwf_complex));
  if(dataFFT == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR rotateSinglepulse: fftwf_malloc failed.");
    return 0;
  }
  plan1 = fftwf_plan_dft_r2c_1d(npts, data, dataFFT, FFTW_ESTIMATE);
  plan2 = fftwf_plan_dft_c2r_1d(npts, dataFFT, data, FFTW_ESTIMATE);


  fftwf_execute(plan1);



  fac = 1.0/(float)npts;
  dtheta = -2.0*M_PI*epsilon/(float)npts;
  for (i=0; i < npts2; i++) {
    dataFFT[i] *= fac*(cos(i*dtheta) + I*sin(i*dtheta));
  }


  fftwf_execute(plan2);

  fftwf_destroy_plan(plan1);
  fftwf_destroy_plan(plan2);
  fftwf_free(dataFFT);
  return 1;
}







int crosscorrelation_fft(float *data1, float *data2, int ndata, float *cc, verbose_definition verbose)
{
  int i, npts2;
  float fac;
  fftwf_complex *dataFFT1, *dataFFT2;
  fftwf_plan plan1, plan2, plan3;

  npts2 = ndata/2+1;
  dataFFT1 = (fftwf_complex *)fftwf_malloc(npts2*sizeof(fftwf_complex));
  dataFFT2 = (fftwf_complex *)fftwf_malloc(npts2*sizeof(fftwf_complex));
  if(dataFFT1 == NULL || dataFFT2 == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR crosscorrelation_fft: fftwf_malloc failed.");
    return 0;
  }
  plan1 = fftwf_plan_dft_r2c_1d(ndata, data1, dataFFT1, FFTW_ESTIMATE);
  plan2 = fftwf_plan_dft_r2c_1d(ndata, data2, dataFFT2, FFTW_ESTIMATE);
  plan3 = fftwf_plan_dft_c2r_1d(ndata, dataFFT1, cc, FFTW_ESTIMATE);


  fftwf_execute(plan1);
  fftwf_execute(plan2);



  fac = 1.0/(float)ndata;
  for (i=0; i < npts2; i++) {
    dataFFT1[i] *= fac*conj(dataFFT2[i]);
  }


  fftwf_execute(plan3);






  fftwf_destroy_plan(plan1);
  fftwf_destroy_plan(plan2);
  fftwf_destroy_plan(plan3);
  fftwf_free(dataFFT1);
  fftwf_free(dataFFT2);
  return 1;
}





int crosscorrelation_fft_padding_cclength(int ndata, int extrazeropad)
{
  int i, ndata_padded;
  i = (int) (log10(1.0 * (ndata+extrazeropad))/log10(2.0));
  ndata_padded = pow(2.0,(i+1));

  if(ndata_padded/2 == ndata)
    ndata_padded = ndata;
  return ndata_padded;
}
int crosscorrelation_fft_padding(float *data1, float *data2, int ndata, int extrazeropad, float **cc, int *cclength, verbose_definition verbose)
{
  float *padded1, *padded2;
  int i, ndata_padded;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Calculating cross correlation of data with length %d\n", ndata);
    if(extrazeropad != 0) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Padding at least %d zero's after data\n", extrazeropad);
    }
  }
  ndata_padded = crosscorrelation_fft_padding_cclength(ndata, extrazeropad);
  if(ndata_padded == ndata) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Data length is already a power of two\n");
    }
    padded1 = data1;
    padded2 = data2;
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Zero-pad it to %d points\n", ndata_padded);
    }
    padded1 = calloc(ndata_padded, sizeof(float));
    padded2 = calloc(ndata_padded, sizeof(float));
    if(padded1 == NULL || padded2 == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR crosscorrelation_fft_padding: Memory allocation error");
      return 0;
    }
    memcpy(padded1, data1, ndata*sizeof(float));
    memcpy(padded2, data2, ndata*sizeof(float));
  }
  *cc = (float *)malloc(ndata_padded*sizeof(float));
  *cclength = ndata_padded;
  if(*cc == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR crosscorrelation_fft_padding: Memory allocation error");
    return 0;
  }
  if(crosscorrelation_fft(padded1, padded2, ndata_padded, *cc, verbose) == 0)
    return 0;
  if(ndata_padded != ndata) {
    free(padded1);
    free(padded2);
  }
  return 1;
}

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
#include <string.h>
#include "psrsalsa.h"


#define NrBoxCarWidths 29
int BoxCars[NrBoxCarWidths] = {1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,40,50,60,75,100,250,500,1000,1500,2000,2500,3000};
float integratePulseEnergy(float *pulse, int bin1, int bin2, float baseline, int squared)
{
  int i;
  float E;
  E = 0;
  for(i = bin1; i <= bin2; i++) {
    if(squared != 0)
      E += (pulse[i]-baseline)*(pulse[i]-baseline);
    else
      E += (pulse[i]-baseline);
  }
  return E;
}
void offpulseStats(float *pulse, int nrBins, float *baseline, float *rms, pulselongitude_regions_definition *onpulse, int nodebase, verbose_definition verbose)
{
  int i, NrOffpulseBins;
  float E;
  E = 0;
  NrOffpulseBins = 0;
  for(i = 0; i < nrBins; i++) {
    if(checkRegions(i, onpulse, 0, verbose) == 0) {
      NrOffpulseBins++;
      E += pulse[i];
    }
  }
  *baseline = E/(float)(NrOffpulseBins);
  if(nodebase)
    *baseline = 0;
  *rms = 0;
  for(i = 0; i < nrBins; i++) {
    if(checkRegions(i, onpulse, 0, verbose) == 0) {
      *rms += (pulse[i]-(*baseline))*(pulse[i]-(*baseline));
    }
  }
  *rms = sqrt(*rms);
  *rms /= sqrt(NrOffpulseBins);
}
void boxcarFindpeak_core(int width, float *pulse, int nrBins, int *bin, int *pulsewidth, float *snrbest, float *E_best, int posOrNeg, int squared, float baseline, float rms, int *firsttime, int *allowedWidths, verbose_definition verbose)
{
  float E, snr;
  int b, ok;
  long nrtrials;
  if(width < nrBins) {
    if(verbose.verbose)
      printf("boxcarFindpeak: Try out width %d: ", width);
    nrtrials = 0;
    for(b = 0; b < nrBins-width; b++) {
      ok = 1;
      if(allowedWidths[b] < b+width-1)
 ok = 0;
      if(ok == 1) {
 nrtrials++;
 E = integratePulseEnergy(pulse, b, b+(width-1), baseline, squared);
 snr = E/(rms*sqrt(width));
 ok = 0;
 if(posOrNeg == 0) {
   if(snr > *snrbest || *firsttime == 1)
     ok = 1;
 }else {
   if(fabs(snr) > *snrbest || *firsttime == 1) {
     ok = 1;
     snr = fabs(snr);
   }
 }
 if(ok) {
   *firsttime = 0;
   *bin = b;
   *pulsewidth = width;
   *snrbest = snr;
   *E_best = E;
 }
      }
    }
    if(verbose.verbose)
      printf("%ld trials (not necessarily independent).\n", nrtrials);
  }
}
int boxcarFindpeak(float *pulse, int nrBins, pulselongitude_regions_definition *onpulse, int *bin, int *pulsewidth, float *snrbest, float *E_best, int squared, int posOrNeg, int allwidths, int refine, int maxwidth, int only_onpulse, int nodebase, verbose_definition verbose)
{
  float baseline, rms;
  int w, w1, w2, dw, width, NrWidths, firsttime, *allowedWidths, b, b2;
  pulselongitude_regions_definition onpulse_search;
  if(initPulselongitudeRegion(&onpulse_search, verbose) == 0) {
    printerror(verbose.debug, "ERROR boxcarFindpeak: Initialising onpulse region failed.");
    return 0;
  }
  offpulseStats(pulse, nrBins, &baseline, &rms, onpulse, nodebase, verbose);
  *snrbest = 0;
  firsttime = 1;
  if(allwidths) {
    NrWidths = nrBins;
  }else {
    NrWidths = NrBoxCarWidths;
  }
  onpulse_search.nrRegions = 1;
  onpulse_search.left_bin[0] = 0;
  onpulse_search.right_bin[0] = nrBins-1;
  onpulse_search.bins_defined[0] = 1;
  if(onpulse != NULL) {
    copyPulselongitudeRegion(*onpulse, &onpulse_search);
    if(only_onpulse == 0) {
      onpulse_search.nrRegions = 1;
      onpulse_search.left_bin[0] = 0;
      onpulse_search.right_bin[0] = nrBins-1;
      onpulse_search.bins_defined[0] = 1;
    }
    if(only_onpulse == 2)
      onpulse_search.nrRegions = 1;
  }else {
    if(only_onpulse != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR boxcarFindpeak: should define onpulse region when limiting search range.");
      return 0;
    }
  }
  allowedWidths = (int *)malloc(nrBins*sizeof(int));
  if(allowedWidths == 0) {
    fflush(stdout);
    printerror(verbose.debug, "boxcarFindpeak: Cannot allocate memory.");
    return 0;
  }
  for(b = 0; b < nrBins; b++) {
    allowedWidths[b] = -1;
    for(b2 = b; b2 < nrBins; b2++) {
      if(checkRegions(b2, &onpulse_search, 0, verbose) == 0) {
 break;
      }else {
 allowedWidths[b] = b2;
      }
    }
  }
  *pulsewidth = 0;
  for(w = 0; w < NrWidths; w++) {
    if(allwidths)
      width = w;
    else
      width = BoxCars[w];
    if(width > 0) {
      if(maxwidth <= 0 || width <= maxwidth)
 boxcarFindpeak_core(width, pulse, nrBins, bin, pulsewidth, snrbest, E_best, posOrNeg, squared, baseline, rms, &firsttime, allowedWidths, verbose);
    }
  }
  if(refine == 1 && allwidths == 0) {
    for(w = 0; w < NrWidths; w++) {
      if(BoxCars[w] == *pulsewidth)
 break;
    }
    if(w > 0 && w < NrBoxCarWidths-1) {
      w1 = BoxCars[w-1];
      w2 = BoxCars[w+1];
      if(w2 > nrBins)
 w2 = nrBins;
      if(w1 < 1)
 w1 = 1;
      if(w2-w1 > 1) {
 do {
   dw = (w2-w1)/7;
   if(dw < 1)
     dw = 1;
   for(w = w1+dw; w <= w2-dw; ) {
     if(maxwidth <= 0 || w <= maxwidth)
       boxcarFindpeak_core(w, pulse, nrBins, bin, pulsewidth, snrbest, E_best, posOrNeg, squared, baseline, rms, &firsttime, allowedWidths, verbose);
     w += dw;
   }
   w1 = *pulsewidth - dw;
   w2 = *pulsewidth + dw;
   if(w2 > nrBins)
     w2 = nrBins;
   if(w1 < 1)
     w1 = 1;
 }while(dw != 1);
      }
    }
  }
  freePulselongitudeRegion(&onpulse_search);
  free(allowedWidths);
  return 1;
}

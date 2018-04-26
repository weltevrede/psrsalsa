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
#define carousel_interpolation_limit 0.005
#define carousel_from_p3fold_smoothing_multiplier 10
void foldP3_simple(float *data, long nry, long starty, long nrx, float *map, float *nrcounts, int nr_p3_bins, float foldp3, float offset, int noNormalise, int noSmooth, float smoothWidth, float slope, float offset2,
     int debug)
{
  float p3, j_frac, weight, weightnext, dp3;
  int i, j, b, jnext;
  for(i = 0; i < nr_p3_bins; i++) {
    for(b = 0; b < nrx; b++) {
      nrcounts[i*nrx+b] = 0;
      map[i*nrx+b] = 0;
    }
  }
  for(i = starty; i < nry; i++) {
    if(smoothWidth <= 0) {
      for(b = 0; b < nrx; b++) {
 p3 = 360.0*(float)(i+offset)/foldp3-slope*b + offset2;
 p3 = derotate_deg(p3);
 if(p3 == 360.0)
   p3 = 0;
 if(b == 0) {
   if(debug) {
     printf("DEBUG foldP3_simple: pulse=%d (block=%ld ... %ld) folded at P3=%f P, with an offset=%f P and additional offset2=%f deg: Subpulse phase of pulse longitude bin 0 = %f deg\n", i, starty, nry-1, foldp3, offset, offset2, p3);
   }
 }
 j = (p3/360.0)*nr_p3_bins;
 if(noSmooth) {
   j_frac = 0;
 }else {
   j_frac = (p3/360.0)*nr_p3_bins - j;
 }
 weight = 1.0-j_frac;
 if(j >= nr_p3_bins || j < 0)
   fprintf(stderr, "ERROR foldP3_simple: Buffer overflow for i=%d, j=%d (nr_p3_bins=%d) p3=%f\n", i, j, nr_p3_bins, p3);
 map[j*nrx+b] += weight*data[i*nrx+b];
 nrcounts[j*nrx+b] += weight;
 if(noSmooth == 0) {
   jnext = j+1;
   if(jnext == nr_p3_bins)
     jnext = 0;
   weightnext = j_frac;
   map[jnext*nrx+b] += weightnext*data[i*nrx+b];
   nrcounts[jnext*nrx+b] += weightnext;
 }
      }
    }else {
      for(j = 0; j < nr_p3_bins; j++) {
 for(b = 0; b < nrx; b++) {
   p3 = 360.0*(float)(i+offset)/foldp3-slope*b + offset2;
   p3 = derotate_deg(p3);
   if(p3 == 360.0)
     p3 = 0;
   if(b == 0 && j == 0) {
     if(debug) {
       printf("DEBUG foldP3_simple: pulse=%d (block=%ld ... %ld) folded at P3=%f P, with an offset=%f P and additional offset2=%f deg: Subpulse phase of pulse longitude bin 0 = %f deg\n", i, starty, nry-1, foldp3, offset, offset2, p3);
     }
   }
   p3 *= nr_p3_bins/360.0;
   dp3 = fabs(p3-j);
   if(fabs(p3-j+nr_p3_bins) < dp3) {
     dp3 = fabs(p3-j+nr_p3_bins);
   }
   if(fabs(p3-j-nr_p3_bins) < dp3) {
     dp3 = fabs(p3-j-nr_p3_bins);
   }
   weight = exp(-(dp3*dp3/(smoothWidth*smoothWidth)));
   map[j*nrx+b] += weight*data[i*nrx+b];
   nrcounts[j*nrx+b] += weight;
 }
      }
    }
  }
  if(noNormalise == 0) {
    for(i = 0; i < nr_p3_bins; i++) {
      for(b = 0; b < nrx; b++) {
 if(nrcounts[i*nrx+b] > 0) {
   map[i*nrx+b] /= nrcounts[i*nrx+b];
 }
      }
    }
  }
}
int foldP3(float *data, long nry, long nrx, float *map, int nr_p3_bins, float foldp3, int refine, int cyclesperblock, int noSmooth, float smoothWidth, float slope, float subpulse_offset, pulselongitude_regions_definition *onpulse
, verbose_definition verbose)
{
  float *blockmap, correl, maxcorrel, *template, *nrcounts, *nrcounts_block;
  int i, b, offset, *bestoffset, itt, ok;
  long startpulse, pulsesleft, dN, blockcounter;
  if(cyclesperblock < 1) {
    fflush(stdout);
    printerror(verbose.debug, "foldP3: cyclesperblock (%d) makes no sense", cyclesperblock);
    return 0;
  }
  dN = foldp3*cyclesperblock;
    {
      if(verbose.verbose)
 printf("Folding data onto P3=%f using %ld subint blocks of data (refine=%d cyclesperblock=%d smoothWidth=%f slope=%f deg/bin offset=%f)\n", foldp3, dN, refine, cyclesperblock, smoothWidth, slope, subpulse_offset);
    }
  nrcounts = (float *)malloc(nrx*nr_p3_bins*sizeof(float));
  if(nrcounts == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "foldP3: cannot allocate memory");
    return 0;
  }
  bestoffset = (int *)malloc((1+nry/dN)*sizeof(int));
  if(bestoffset == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "foldP3: cannot allocate memory");
    return 0;
  }
  if(refine <= 0) {
    foldP3_simple(data, nry, 0, nrx, map, nrcounts, nr_p3_bins, foldp3, 0, 0, noSmooth, smoothWidth, slope, subpulse_offset,
    0*verbose.debug);
  }else {
    blockmap = (float *)malloc(nr_p3_bins*nrx*sizeof(float));
    if(blockmap == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "foldP3: cannot allocate memory");
      return 0;
    }
    nrcounts_block = (float *)malloc(nrx*nr_p3_bins*sizeof(float));
    if(nrcounts_block == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "foldP3: cannot allocate memory");
      return 0;
    }
    if(refine > 1) {
      template = (float *)malloc(nr_p3_bins*nrx*sizeof(float));
      if(template == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "foldP3: cannot allocate memory");
 return 0;
      }
    }
    for(itt = 0; itt < refine; itt++) {
      for(i = 0; i < nr_p3_bins; i++) {
 for(b = 0; b < nrx; b++) {
   nrcounts[i*nrx+b] = 0;
   map[i*nrx+b] = 0;
 }
      }
      startpulse = 0;
      pulsesleft = nry;
      blockcounter = 0;
      while(pulsesleft > foldp3*cyclesperblock) {
 {
   maxcorrel = 0;
   for(offset = 0; offset < nr_p3_bins; offset++) {
     foldP3_simple(data, startpulse+dN, startpulse, nrx, blockmap, nrcounts_block, nr_p3_bins, foldp3, offset*foldp3/(float)nr_p3_bins, 1, noSmooth, smoothWidth, slope, subpulse_offset,
0*verbose.debug);
     correl = 0;
     for(i = 0; i < nr_p3_bins; i++) {
       for(b = 0; b < nrx; b++) {
  ok = 0;
  if(onpulse == NULL) {
    ok = 1;
  }else if(checkRegions(b, onpulse, 0, verbose)) {
    ok = 1;
  }
  if(ok) {
    if(itt > 0)
      correl += template[i*nrx+b]*template[i*nrx+b]*blockmap[i*nrx+b]*blockmap[i*nrx+b];
    else
      correl += map[i*nrx+b]*map[i*nrx+b]*blockmap[i*nrx+b]*blockmap[i*nrx+b];
  }
       }
     }
     if(correl > maxcorrel || offset == 0) {
       maxcorrel = correl;
       bestoffset[blockcounter] = offset;
     }
   }
 }
 foldP3_simple(data, startpulse+dN, startpulse, nrx, blockmap, nrcounts_block, nr_p3_bins, foldp3, bestoffset[blockcounter]*foldp3/(float)nr_p3_bins, 0, noSmooth, smoothWidth, slope, subpulse_offset,
verbose.debug);
 for(i = 0; i < nr_p3_bins; i++) {
   for(b = 0; b < nrx; b++) {
     nrcounts[i*nrx+b] += nrcounts_block[i*nrx+b];
     map[i*nrx+b] += (float)blockmap[i*nrx+b];
   }
 }
 startpulse += dN;
 pulsesleft -= dN;
 if(verbose.verbose && verbose.nocounters == 0) {
   printf("  Itteration %d/%d, pulse %ld/%ld     \r", itt+1, refine, startpulse, nry);
   fflush(stdout);
 }
 blockcounter++;
      }
      if(refine > 1) {
 for(i = 0; i < nr_p3_bins; i++) {
   for(b = 0; b < nrx; b++) {
     template[i*nrx+b] = (float)map[i*nrx+b];
   }
 }
      }
    }
    if(verbose.verbose && verbose.nocounters == 0) {
      printf("\n");
    }
    free(blockmap);
    free(nrcounts_block);
    if(refine > 1) {
      free(template);
    }
  }
  free(nrcounts);
  free(bestoffset);
  if(verbose.verbose)
    printf("  Done\n");
  return 1;
}

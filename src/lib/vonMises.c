/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "psrsalsa.h"






int readVonMisesModel(char *filename, vonMises_def *components, verbose_definition verbose)
{
  int i;
  FILE *fin;
  if(verbose.verbose) printf("Opening %s for reading\n", filename);
  fin = fopen(filename, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readVonMisesModel: Cannot open %s", filename);
    return 0;
  }
  for(components->nrcomponents = 0; components->nrcomponents < maxNrVonMisesComponents; (components->nrcomponents)++) {
    i = fscanf(fin, "%lf %lf %lf", &(components->centre[components->nrcomponents]), &(components->concentration[components->nrcomponents]), &(components->height[components->nrcomponents]));
    if(i != 3) {
      break;
    }else if(verbose.verbose) {
      printf("  component %d: phase=%e concentration=%e amplitude=%e\n", components->nrcomponents, components->centre[components->nrcomponents], components->concentration[components->nrcomponents], components->height[components->nrcomponents]);
    }
  }
  if(verbose.verbose) printf("Closing %s\n", filename);
  fclose(fin);
  if(components->nrcomponents == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readVonMisesModel: No components found in %s", filename);
    return 0;
  }
  return 1;
}


double calcVonMisesFunction2(double centre, double concentration, double height, double phase, double shift)
{
  double y;
  y = exp((cos(2.0*M_PI*(phase-centre-shift))-1.0)*concentration) * height;
  return y;
}



double calcVonMisesFunction(vonMises_def components, double phase, double shift)
{
  int n;
  double y;
  y = 0;
  if(components.nrcomponents > 0) {
    for(n = 0; n < components.nrcomponents; n++) {
      y += calcVonMisesFunction2(components.centre[n], components.concentration[n], components.height[n], phase, shift);
    }
  }
  return y;
}



void calcVonMisesProfile(vonMises_def components, int nrbins, float *profile, double shift, int normalize)
{
  int i;
  double x, Imax = -1;

  for(i = 0; i < nrbins; i++) {
    profile[i] = 0;
  }
  for(i = 0; i < nrbins; i++) {
    x = i/(double)nrbins;
    profile[i] = calcVonMisesFunction(components, x, shift);
  }
  if(normalize) {
    for(i = 0; i < nrbins; i++) {
      if(profile[i] > Imax)
 Imax = profile[i];
    }
    for(i = 0; i < nrbins; i++) {
      profile[i] /= Imax;
    }
  }
}




float correlateVonMisesFunction(vonMises_def components, int nrbins, float *profile, verbose_definition verbose)
{
  float correl_max, *profile2;
  int i;


  int ishift;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Correlating template with profile\n");
    verbose.indent += 2;
  }
  profile2 = (float *)malloc(nrbins*sizeof(float));
  if(profile2 == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR correlateVonMisesFunction: Memory allocation error.");
    return 0;
  }

  calcVonMisesProfile(components, nrbins, profile2, 0, 0);
  find_peak_correlation(profile, profile2, nrbins, 0, 1, 1, &ishift, &correl_max, verbose);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Found a shift of %d bins (%f phase) and max/min correlation is %f.\n", ishift, ishift/(float)nrbins, correl_max);
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("so you could do a shift by %d bins (%f phase) to rotate profile to align with model.\n", -ishift, -ishift/(float)nrbins);
  }
  free(profile2);
  return ishift/(float)nrbins;
}

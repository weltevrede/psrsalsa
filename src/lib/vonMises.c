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
#include <gsl/gsl_sf.h>
#include "psrsalsa.h"
#define AmoebaAlgorithm 0
int readVonMisesModel(char *filename, vonMises_collection_definition *components, verbose_definition verbose)
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
void vonMises_simplify_parameters(vonMises_collection_definition *components)
{
  int i;
  for(i = 0; i < components->nrcomponents; i++) {
    if(components->concentration[i] < 0) {
      components->centre[i] += 0.5;
      if(components->centre[i] >= 1)
 components->centre[i] -= 1.0;
      components->concentration[i] *= -1.0;
      components->height[i] *= exp(2.0*components->concentration[i]);
    }
  }
}
int writeVonMisesModel(char *filename, vonMises_collection_definition *components, verbose_definition verbose)
{
  int n;
  FILE *fout;
  if(filename != NULL) {
    if(verbose.verbose) printf("Opening %s for writing\n", filename);
    fout = fopen(filename, "w");
    if(fout == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeVonMisesModel: Cannot open %s", filename);
      return 0;
    }
  }else {
    fout = stdout;
  }
  for(n = 0; n < components->nrcomponents; n++) {
    fprintf(fout, "%e %e %e\n", components->centre[n], components->concentration[n], components->height[n]);
    if(verbose.verbose) {
      if(n == 0) {
 printf("y = A*exp((cos(2.0*M_PI*(x-x0))-1.0)*concentration)\n");
 printf("  component: x0 concentration H\n");
      }
      printf("  component %d: %e %e %e\n", n, components->centre[n], components->concentration[n], components->height[n]);
    }
  }
  if(filename != NULL) {
    if(verbose.verbose) printf("Closing %s\n", filename);
    fclose(fout);
  }
  return 1;
}
double calcVonMisesFunction2(double centre, double concentration, double height, double phase, double shift)
{
  double y;
  y = exp((cos(2.0*M_PI*(phase-centre-shift))-1.0)*concentration) * height;
  return y;
}
double calcVonMisesFunction2_deriv(double centre, double concentration, double height, double phase, double shift)
{
  double y;
  y = calcVonMisesFunction2(centre, concentration, height, phase, shift);
  y *= -2.0*M_PI*concentration*sin(2.0*M_PI*(phase-centre-shift));
  return y;
}
double integrateVonMisesFunction2(double concentration, double height)
{
  double area;
  area = 2*M_PI*height*gsl_sf_bessel_I0_scaled(concentration);
  return area;
}
double widthVonMisesFunction2(double concentration, double ampfrac)
{
  double value = log(ampfrac)/concentration + 1.0;
  if(value < -1.0 || value >= 1.0) {
    return sqrt(-1);
  }
  return 2.0*acos(value);
}
double calcVonMisesFunction(vonMises_collection_definition *components, double phase, double shift)
{
  int n;
  double y;
  y = 0;
  if(components->nrcomponents > 0) {
    for(n = 0; n < components->nrcomponents; n++) {
      y += calcVonMisesFunction2(components->centre[n], components->concentration[n], components->height[n], phase, shift);
    }
  }
  return y;
}
double calcVonMisesFunction_deriv(vonMises_collection_definition *components, double phase, double shift)
{
  int n;
  double y;
  y = 0;
  if(components->nrcomponents > 0) {
    for(n = 0; n < components->nrcomponents; n++) {
      y += calcVonMisesFunction2_deriv(components->centre[n], components->concentration[n], components->height[n], phase, shift);
    }
  }
  return y;
}
double integrateVonMisesFunction(vonMises_collection_definition *components)
{
  int n;
  double area;
  area = 0;
  if(components->nrcomponents > 0) {
    for(n = 0; n < components->nrcomponents; n++) {
      area += integrateVonMisesFunction2(components->concentration[n], components->height[n]);
    }
  }
  return area;
}
void calcVonMisesProfile(vonMises_collection_definition *components, int nrbins, float *profile, double shift, int normalize)
{
  int i;
  double x, yvalue, Imax = -1;
  for(i = 0; i < nrbins; i++) {
    profile[i] = 0;
  }
  for(i = 0; i < nrbins; i++) {
    x = i/(double)nrbins;
    yvalue = calcVonMisesFunction(components, x, shift);
    profile[i] = yvalue;
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
void calcVonMisesProfile_resid_rms(vonMises_collection_definition *components, int nrbins, float *profile, double shift, double *rms)
{
  int i;
  double x, y;
  *rms = 0;
  for(i = 0; i < nrbins; i++) {
    x = i/(double)nrbins;
    y = profile[i] - calcVonMisesFunction(components, x, shift);
    *rms += y*y;
  }
  *rms /= (double)nrbins;
  *rms = sqrt(*rms);
}
float correlateVonMisesFunction(vonMises_collection_definition *components, int nrbins, float *profile, verbose_definition verbose)
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
  find_peak_correlation(profile, profile2, nrbins, 0, 1, 1, 1, &ishift, &correl_max, verbose);
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
int getboundry_single(float *profile, int nrbins, int istart, float y, int *bin)
{
  int i, slope, found;
  if(profile[istart] < y)
    slope = 1;
  else
    slope = -1;
  found = 0;
  for(i = istart; i < nrbins; i++) {
    if(profile[i] > y && slope == 1) {
      *bin = i-1;
      if(*bin < istart)
 *bin = istart;
      found = 1;
      break;
    }
    if(profile[i] < y && slope == -1) {
      *bin = i+1;
      if(*bin >= nrbins)
 *bin = nrbins;
      found = 1;
      break;
    }
  }
  return found*slope;
}
void find_boundaries(float *profile, int nrbins, float y, pulselongitude_regions_definition *regions)
{
  int ret, lastbin, bin, lastbin_set, firstpulse;
  lastbin = -1;
  firstpulse = 1;
  clearPulselongitudeRegion(regions);
  lastbin_set = 0;
  do {
    ret = getboundry_single(profile, nrbins, lastbin+1, y, &bin);
    if(firstpulse == 1 && ret == -1) {
      regions->left_bin[regions->nrRegions] = 0;
      regions->right_bin[regions->nrRegions] = bin;
      regions->bins_defined[regions->nrRegions] = 1;
      (regions->nrRegions)++;
      lastbin = bin;
      firstpulse = 0;
    }else if(ret == 1) {
      lastbin = bin;
      lastbin_set = 1;
      firstpulse = 0;
    }else if(ret == -1) {
      regions->left_bin[regions->nrRegions] = lastbin;
      regions->right_bin[regions->nrRegions] = bin;
      regions->bins_defined[regions->nrRegions] = 1;
      (regions->nrRegions)++;
      lastbin = bin;
      lastbin_set = 0;
      firstpulse = 0;
    }
    if(regions->nrRegions == MAX_pulselongitude_regions)
      ret = 0;
  }while(ret != 0);
  if(lastbin_set && (regions->nrRegions < MAX_pulselongitude_regions)) {
    regions->left_bin[regions->nrRegions] = lastbin;
    regions->right_bin[regions->nrRegions] = nrbins-1;
    regions->bins_defined[regions->nrRegions] = 1;
    (regions->nrRegions)++;
  }
}
int internal_fitvonmises_nrcomponents;
int internal_fitvonmises_nrbins;
int internal_fitvonmises_fitbaseline;
int internal_fitvonmises_fitdummy;
int internal_fitvonmises_plotallsteps;
float *internal_fitvonmises_profile;
float *internal_fitvonmises_fitprofile;
int internal_fitvonmises_relative_phases = 0;
int internal_fitvonmises_relative_amps = 0;
int internal_fitvonmises_debug = 0;
int internal_fitvonmises_showed_overflow_warning = 0;
int internal_fitvonmises_avoid_neg_components = 0;
pulselongitude_regions_definition *internal_fitvonmises_onpulse = NULL;
float internal_fitvonmises_funk(float *params)
{
  static vonMises_collection_definition *components = NULL;
  int n, i, got_neg_component, got_component_outside_onpulse_region;
  double chi2, chi, y;
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.debug = internal_fitvonmises_debug;
  if(components == NULL) {
    components = malloc(sizeof(vonMises_collection_definition));
    if(components == NULL) {
      fprintf(stderr, "ERROR internal_fitvonmises_funk: Cannot allocate memory\n");
      exit(0);
    }
  }
  got_neg_component = 0;
  got_component_outside_onpulse_region = 0;
  components->nrcomponents = internal_fitvonmises_nrcomponents;
  for(n = 0; n < internal_fitvonmises_nrcomponents; n++) {
    double centre = params[3*n+internal_fitvonmises_fitbaseline];
    if(internal_fitvonmises_relative_phases && n > 0) {
      centre += components->centre[0];
    }
    components->centre[n] = centre;
    components->concentration[n] = params[3*n+1+internal_fitvonmises_fitbaseline];
    components->height[n] = params[3*n+2+internal_fitvonmises_fitbaseline];
    if(components->height[n] <= 0.0) {
      got_neg_component = 1;
    }
    if(internal_fitvonmises_onpulse != NULL) {
      if(internal_fitvonmises_onpulse->nrRegions > 0) {
 centre -= floor(centre);
 centre *= internal_fitvonmises_nrbins;
 if(centre >= internal_fitvonmises_nrbins) {
   centre -= internal_fitvonmises_nrbins;
 }
 if(checkRegions(centre, internal_fitvonmises_onpulse, 0, verbose) == 0) {
   got_component_outside_onpulse_region = 1;
 }
      }
    }
    if(internal_fitvonmises_relative_amps && n > 0) {
      components->height[n] *= params[2+internal_fitvonmises_fitbaseline];
    }
    if(internal_fitvonmises_plotallsteps) {
      ppgsci(4);
      y = calcVonMisesFunction2(components->centre[n], components->concentration[n], components->height[n], 0, 0);
      if(internal_fitvonmises_fitbaseline)
 y += params[0];
      ppgmove(0, y);
      for(i = 1; i < internal_fitvonmises_nrbins; i++) {
 y = calcVonMisesFunction2(components->centre[n], components->concentration[n], components->height[n], i/(double)internal_fitvonmises_nrbins, 0);
 if(internal_fitvonmises_fitbaseline)
   y += params[0];
 ppgdraw(i*360.0/(double)internal_fitvonmises_nrbins, y);
      }
    }
  }
  calcVonMisesProfile(components, internal_fitvonmises_nrbins, internal_fitvonmises_fitprofile, 0, 0);
  chi2 =0;
  for(n = 0; n < internal_fitvonmises_nrbins; n++) {
    int accept_bin = 1;
    if(internal_fitvonmises_onpulse != NULL) {
      if(internal_fitvonmises_onpulse->nrRegions > 0) {
 if(checkRegions(n, internal_fitvonmises_onpulse, 0, verbose) == 0) {
   accept_bin = 0;
 }
      }
    }
    if(accept_bin) {
      chi = (internal_fitvonmises_fitprofile[n] - internal_fitvonmises_profile[n]);
      if(internal_fitvonmises_fitbaseline)
 chi += params[0];
      chi2 += chi*chi;
    }
  }
  if(internal_fitvonmises_fitdummy) {
    n = 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline;
    chi2 += params[n]*params[n];
    if(isnan(params[n]) || isinf(params[n])) {
      fflush(stdout);
      printerror(verbose.debug, "internal_fitvonmises_funk: chi2=%f produced by dummy variable\n", chi2);
    }
  }
  if(isnan(chi2) || isinf(chi2)) {
    int tooneg;
    tooneg = 0;
    for(n = 0; n < internal_fitvonmises_nrcomponents; n++) {
      if(components->concentration[n] < -15 && isinf(components->concentration[n]) == 0) {
 tooneg = 1;
      }
    }
    if(tooneg) {
      if(internal_fitvonmises_showed_overflow_warning == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING internal_fitvonmises_funk: At least one trial resulted in a too negative concentration causing a floating-point overflow. These solution will be ignored. This warning will only be shown once.");
 internal_fitvonmises_showed_overflow_warning = 1;
      }
    }
    if(tooneg == 0 || verbose.verbose || verbose.debug) {
      fflush(stdout);
      printerror(verbose.debug, "internal_fitvonmises_funk: chi2=%f", chi2);
      if(verbose.debug) {
 writeVonMisesModel(NULL, components, verbose);
 if(internal_fitvonmises_fitbaseline) {
   printf("Baseline=%lf\n", params[0]);
 }
 calcVonMisesProfile(components, internal_fitvonmises_nrbins, internal_fitvonmises_fitprofile, 0, 0);
 printf("Fitted function:");
 for(n = 0; n < internal_fitvonmises_nrbins; n++) {
   printf(" %f", internal_fitvonmises_fitprofile[n]);
 }
 printf("\n");
 printf("Data:");
 for(n = 0; n < internal_fitvonmises_nrbins; n++) {
   printf(" %f", internal_fitvonmises_profile[n]);
 }
 printf("\n");
      }
    }
    chi2 = 1e100;
  }
  if((got_neg_component && internal_fitvonmises_avoid_neg_components) || got_component_outside_onpulse_region) {
    chi2 *= 1e10;
  }
  return sqrt(chi2);
}
int fitvonmises_refine_model(datafile_definition psrdata, vonMises_collection_definition *components, int fitbaseline, int avoid_neg_components, float *baseline, int fixamp, int fixwidth, int fixphase, int fixrelamp, int fixrelphase, verbose_definition verbose)
{
  float *xstart, *dx, *xfit;
  float chi2, av;
  int *fixed, nfunk;
  int i, superverbose;
  internal_fitvonmises_debug = verbose.debug;
  internal_fitvonmises_avoid_neg_components = avoid_neg_components;
  internal_fitvonmises_relative_phases = fixrelphase;
  internal_fitvonmises_relative_amps = fixrelamp;
  superverbose = 0;
  if(verbose.verbose) {
    printf("fitvonmises_refine_model: Start refining model with %d components.\n", components->nrcomponents);
  }
  if(fixamp && fixwidth && fixphase) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitvonmises_refine_model: Cannot simultaneously fix the phases, concentrations and amplitudes of the components.");
    return 0;
  }
  if(psrdata.NrFreqChan != 1 || psrdata.NrSubints != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitvonmises_refine_model: This is not a pulse profile.");
    return 0;
  }
  if(psrdata.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitvonmises_refine_model: only works if data is loaded into memory.");
    return 0;
  }
  internal_fitvonmises_nrbins = psrdata.NrBins;
  internal_fitvonmises_fitprofile = malloc(psrdata.NrBins*sizeof(float));
  internal_fitvonmises_profile = psrdata.data;
  if(internal_fitvonmises_fitprofile == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitvonMises_refine_model: Memory allocation error.");
    return 0;
  }
  xstart = malloc((3*maxNrVonMisesComponents+1)*sizeof(float));
  dx = malloc((3*maxNrVonMisesComponents+1)*sizeof(float));
  xfit = malloc((3*maxNrVonMisesComponents+1)*sizeof(float));
  fixed = malloc((3*maxNrVonMisesComponents+1)*sizeof(int));
  if(xstart == NULL || dx == NULL || xfit == NULL || fixed == NULL) {
    printerror(verbose.debug, "ERROR fitvonmises_refine_model: Cannot allocate memory");
    return 0;
  }
  if(fitbaseline)
    internal_fitvonmises_fitbaseline = 1;
  else
    internal_fitvonmises_fitbaseline = 0;
  av = 0;
  for(i = 0; i < internal_fitvonmises_nrbins; i++) {
    av += internal_fitvonmises_profile[i];
  }
  av /= (float)internal_fitvonmises_nrbins;
  if(fixamp == 0) {
    float min, max, value, scale;
    min = max = psrdata.data[0];
    for(i = 0; i < psrdata.NrBins; i++) {
      if(psrdata.data[i] < min)
 min = psrdata.data[i];
      if(psrdata.data[i] > max)
 max = psrdata.data[i];
    }
    scale = max - min;
    for(i = 0; i < psrdata.NrBins; i++) {
      value = calcVonMisesFunction(components, i/(float)psrdata.NrBins, 0);
      if(i == 0 || value > max)
 max = value;
    }
    scale /= max;
    if(verbose.debug)
      printf("Scaling input model amplitudes with factor %e\n", scale);
    for(i = 0; i < components->nrcomponents; i++) {
      components->height[i] *= scale;
    }
  }
  if(fixphase == 0) {
    float value;
    value = correlateVonMisesFunction(components, psrdata.NrBins, psrdata.data, verbose);
    if(verbose.debug)
      printf("Shifting input model %e in phase\n", value);
    for(i = 0; i < components->nrcomponents; i++) {
      components->centre[i] += value;
    }
  }
  if(internal_fitvonmises_fitbaseline) {
    xstart[0] = av;
    dx[0] = av*0.1;
    fixed[0] = 0;
  }
  internal_fitvonmises_nrcomponents = components->nrcomponents;
  int nrfitparameters;
  nrfitparameters = 0;
  for(i = 0; i < components->nrcomponents; i++) {
    xstart[3*i+internal_fitvonmises_fitbaseline] = components->centre[i];
    if(fixrelphase && i > 0)
      xstart[3*i+internal_fitvonmises_fitbaseline] -= components->centre[0];
    xstart[3*i+1+internal_fitvonmises_fitbaseline] = components->concentration[i];
    xstart[3*i+2+internal_fitvonmises_fitbaseline] = components->height[i];
    if(fixrelamp && i > 0)
      xstart[3*i+2+internal_fitvonmises_fitbaseline] /= components->height[0];
    if(verbose.debug) {
      printf("Refine component %d with phi=%f, con=%f, amp=%f\n", i+1, components->centre[i], components->concentration[i], components->height[i]);
    }
    if(fixphase || (fixrelphase && i > 0)) {
      fixed[3*i+internal_fitvonmises_fitbaseline] = 1;
      dx[3*i+internal_fitvonmises_fitbaseline] = 0;
    }else {
      fixed[3*i+internal_fitvonmises_fitbaseline] = 0;
      dx[3*i+internal_fitvonmises_fitbaseline] = 0.01;
      nrfitparameters++;
    }
    if(fixwidth) {
      dx[3*i+1+internal_fitvonmises_fitbaseline] = 0;
      fixed[3*i+1+internal_fitvonmises_fitbaseline] = 1;
    }else {
      dx[3*i+1+internal_fitvonmises_fitbaseline] = 0.1*xstart[3*i+1+internal_fitvonmises_fitbaseline];
      fixed[3*i+1+internal_fitvonmises_fitbaseline] = 0;
      nrfitparameters++;
    }
    if(fixamp || (fixrelamp && i > 0)) {
      fixed[3*i+2+internal_fitvonmises_fitbaseline] = 1;
      dx[3*i+2+internal_fitvonmises_fitbaseline] = 0;
    }else {
      fixed[3*i+2+internal_fitvonmises_fitbaseline] = 0;
      dx[3*i+2+internal_fitvonmises_fitbaseline] = 0.1*xstart[3*i+2+internal_fitvonmises_fitbaseline];
      nrfitparameters++;
    }
  }
  if(nrfitparameters == 1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING fitvonMises_refine_model: Adding a dummy variable to construct a multi-parameter fit as required by the downhill-simplex method.");
    i = 3*components->nrcomponents+internal_fitvonmises_fitbaseline;
    internal_fitvonmises_fitdummy = 1;
    xstart[i] = 0;
    dx[i] = 0;
    fixed[i] = 0;
  }else {
    internal_fitvonmises_fitdummy = 0;
  }
  if(superverbose > 2)
    internal_fitvonmises_plotallsteps = 1;
  else
    internal_fitvonmises_plotallsteps = 0;
  if(verbose.verbose) {
    printf("    Start the fitting with %d params, of which ", 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy);
    int n = 0;
    for(i = 0; i < 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy; i++) {
      if(fixed[i] == 0) {
 n++;
      }
    }
    printf("%d are free parameters.", n);
    if(internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy) {
      printf(" This includes ");
      if(internal_fitvonmises_fitbaseline) {
 printf("a baseline offset");
 if(internal_fitvonmises_fitdummy) {
   printf(" and");
 }
      }
      if(internal_fitvonmises_fitdummy) {
 printf("a dummy variable");
      }
    }
    printf("\n");
    if(verbose.debug) {
      printf("Initial values are:     ");
      for(i = 0; i < 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy; i++) {
 if(fixed[i] == 0) {
   printf("%e ", xstart[i]);
 }
      }
      printf("\n");
      printf("Initial step sizes are: ");
      for(i = 0; i < 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy; i++) {
 if(fixed[i] == 0) {
   printf("%e ", dx[i]);
 }
      }
      printf("\n");
    }
  }
  int ret;
  ret = doAmoeba(AmoebaAlgorithm, xstart, dx, fixed, xfit, &chi2, 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy, &internal_fitvonmises_funk, 1e-4, &nfunk, 0*verbose.verbose, 0, 3.0, NULL, NULL);
  if(ret != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitvonMises_refine_model: Amoeba failed.");
    if(ret == 1) {
      printerror(verbose.debug, "ERROR fitvonMises_refine_model: Maximum number of itterations exceeded.");
    }else if(ret == 2) {
      printerror(verbose.debug, "ERROR fitvonMises_refine_model: Memory allocation error.");
    }else if(ret == 2) {
      printerror(verbose.debug, "ERROR fitvonMises_refine_model: Not enough free parameters to perform fit.");
    }else if(ret == 2) {
      printerror(verbose.debug, "ERROR fitvonMises_refine_model: Algorithm is not available.");
    }else if(ret == 2) {
      printerror(verbose.debug, "ERROR fitvonMises_refine_model: Other fit fail error (such as singular matrix).");
    }else if(ret == 2) {
      printerror(verbose.debug, "ERROR fitvonMises_refine_model: Undefined error code.");
    }
    free(xstart);
    free(dx);
    free(xfit);
    free(fixed);
    return 0;
  }
  if(verbose.verbose) {
    printf("    chi2=%f after %d steps.\n", chi2, nfunk);
  }
  double offset;
  offset = xfit[internal_fitvonmises_fitbaseline]-components->centre[0];
  for(i = 0; i < internal_fitvonmises_nrcomponents; i++) {
    if(fixrelphase) {
      components->centre[i] += offset;
    }else {
      components->centre[i] = xfit[3*i+internal_fitvonmises_fitbaseline];
    }
    components->concentration[i] = xfit[3*i+1+internal_fitvonmises_fitbaseline];
    components->height[i] = xfit[3*i+2+internal_fitvonmises_fitbaseline];
    if(fixrelamp && i > 0)
      components->height[i] *= xfit[2+internal_fitvonmises_fitbaseline];
  }
  if(internal_fitvonmises_fitbaseline) {
    *baseline = xfit[0];
    if(verbose.verbose) {
      printf("    Baseline = %f\n", *baseline);
    }
  }
  if(internal_fitvonmises_fitdummy) {
    i = 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline;
    if(verbose.verbose) {
      printf("    Dummy = %e (should be zero)\n", xfit[i]);
    }
  }
  free(internal_fitvonmises_fitprofile);
  if(verbose.verbose) {
    printf("fitvonmises_refine_model: done.\n");
  }
  internal_fitvonmises_relative_phases = 0;
  internal_fitvonmises_relative_amps = 0;
  vonMises_simplify_parameters(components);
  free(xstart);
  free(dx);
  free(xfit);
  free(fixed);
  return 1;
}
void calcVonMisesProfile_shape_parameter_get_location_peak(double *xmax, double *ymax, vonMises_collection_definition *components, double shift, long nrbins_start, double phaseprecision, verbose_definition verbose)
{
  int cont;
  double x, y;
  long nrbins, bin1, bin2, i;
  nrbins = nrbins_start;
  bin1 = 0;
  bin2 = nrbins-1;
  do {
    for(i = bin1; i < bin2; i++) {
      x = i/(double)nrbins;
      y = calcVonMisesFunction(components, x, shift);
      if(i == 0 || y > *ymax) {
 *ymax = y;
 *xmax = x;
      }
    }
    if(verbose.debug) {
      printf("  Found maximum at phase %lf (max=%e)\n", *xmax, *ymax);
    }
    if(1.0/phaseprecision > nrbins)
      cont = 1;
    else
      cont = 0;
    nrbins *= 10;
    bin1 = (*xmax)*nrbins-20;
    bin2 = (*xmax)*nrbins+20;
  }while(cont);
}
void calcVonMisesProfile_shape_parameter_get_edges(double *xedge1, double *xedge2, double *measurement, double xmax, vonMises_collection_definition *components, double shift, double phaseprecision, double yrequest, long nrbins_start, verbose_definition verbose)
{
  int cont;
  long i, nrbins, bin1, bin2;
  double x, y;
  int direction;
  int found;
  for(direction = -1; direction <= 1; direction += 2) {
    found = 0;
    nrbins = nrbins_start;
    if(direction == -1) {
      bin1 = (xmax+0.5)*nrbins;
      bin2 = (xmax-0.5)*nrbins;
    }else {
      bin1 = (xmax-0.5)*nrbins;
      bin2 = (xmax+0.5)*nrbins;
    }
    do {
 i = bin1;
 do {
   x = i/(double)nrbins;
   y = calcVonMisesFunction(components, x, shift);
   if(found == 0 && y <= yrequest) {
     found = 1;
   }
   if(found == 1 && y >= yrequest) {
     found = 2;
     break;
   }
   i += direction;
   if(direction == 1) {
     if(i <= bin2)
       cont = 1;
     else
       cont = 0;
   }else {
     if(i >= bin2)
       cont = 1;
     else
       cont = 0;
   }
 }while(cont);
 if(found < 2) {
   break;
 }
 if(verbose.debug) {
   printf("  Found edge at phase %lf (y=%e is close to %e)\n", x, y, yrequest);
 }
 if(1.0/phaseprecision > nrbins) {
   cont = 1;
   found = 0;
 }else {
   cont = 0;
 }
 nrbins *= 10;
 if(direction == 1) {
   bin1 = x*nrbins-20;
   bin2 = x*nrbins+20;
 }else {
   bin1 = x*nrbins+20;
   bin2 = x*nrbins-20;
 }
      }while(cont);
      if(verbose.debug && direction == -1) {
 printf("  Found one edge, try to find next\n");
      }
      if(found < 2)
 break;
      if(direction == -1) {
 *xedge2 = x;
      }else if(direction == 1) {
 *xedge1 = x;
      }
    }
    if(found == 0) {
      printwarning(verbose.debug, "WARNING calcVonMisesProfile_shape_parameter: Requested width at intensity level (%e) is 360 deg since the intensity level appears to be always higher.\n", yrequest);
      *measurement = 1;
      return;
    }else if(found == 1) {
      printwarning(verbose.debug, "WARNING calcVonMisesProfile_shape_parameter: Requested width at intensity level (%e) is 0 deg since the intensity level appears to be always lower.\n", yrequest);
      *measurement = 0;
      return;
    }else {
      *measurement = (*xedge2)-(*xedge1);
      return;
    }
}
void calcVonMisesProfile_shape_parameter_search_peak_in_restricted_phase_range(int *found, double *x, double startphase, double endphase, vonMises_collection_definition *components, double shift, double phaseprecision, long nrbins_start, verbose_definition verbose)
{
  int cont;
  long nrbins, bin1, bin2, i;
    *found = 0;
    nrbins = nrbins_start;
    bin1 = startphase*nrbins;
    bin2 = endphase*nrbins;
    do {
      for(i = bin1; i <= bin2; i++) {
 double dydx;
 int oldsign;
 *x = i/(double)nrbins;
 dydx = calcVonMisesFunction_deriv(components, *x, shift);
 if(i == bin1) {
   if(dydx >= 0)
     oldsign = 1;
   else
     oldsign = -1;
 }
 if(dydx <= 0 && oldsign == 1) {
   *found = 1;
   break;
 }
 if(oldsign == -1 && dydx >= 0)
   oldsign = 1;
      }
      if(*found == 0) {
 break;
      }
      if(verbose.debug) {
 printf("  Found secondary maximum at phase %lf\n", *x);
      }
      if(1.0/phaseprecision > nrbins) {
 cont = 1;
 *found = 0;
      }else {
 cont = 0;
      }
      nrbins *= 10;
      bin1 = (*x)*nrbins-20;
      bin2 = (*x)*nrbins+20;
    }while(cont);
}
void calcVonMisesProfile_shape_parameter(vonMises_collection_definition *components, double shift, double phaseprecision, int shapepar, double *shapepar_aux, double *measurement, verbose_definition verbose)
{
  long nrbins_start;
  double x, xmax, ymax;
  nrbins_start = 2048;
  if(verbose.debug) {
    printf("Find shape parameter %d\n", shapepar);
  }
  if(components->nrcomponents <= 0) {
    printerror(verbose.debug, "ERROR calcVonMisesProfile_shape_parameter: Analytic model is undefined, so no shape parameter could be measured.");
    *measurement = sqrt(-1);
    return;
  }
  if(shapepar == SHAPEPAR_AMPRATIO_ATPHASE) {
    double value1, value2;
    value1 = calcVonMisesFunction(components, shapepar_aux[0], shift);
    value2 = calcVonMisesFunction(components, shapepar_aux[1], shift);
    *measurement = value1/value2;
    return;
  }else if(shapepar == SHAPEPAR_COMPAMPRATIO) {
    int comp1, comp2;
    comp1 = round(shapepar_aux[0]);
    comp2 = round(shapepar_aux[1]);
    if(comp1 >= components->nrcomponents || comp2 >= components->nrcomponents) {
      printerror(verbose.debug, "ERROR calcVonMisesProfile_shape_parameter: Analytic model has %d components, but components %d and %d (counting from zero) were specified to calculate the amplitude radio.", components->nrcomponents, comp1, comp2);
      *measurement = sqrt(-1);
      return;
    }
    *measurement = components->height[comp1]/components->height[comp2];
    return;
  }
  calcVonMisesProfile_shape_parameter_get_location_peak(&xmax, &ymax, components, shift, nrbins_start, phaseprecision, verbose);
  if(shapepar == SHAPEPAR_PEAKPHASE) {
    *measurement = xmax;
    return;
  }else if(shapepar == SHAPEPAR_PEAKAMP) {
    *measurement = ymax;
    return;
  }
  double yrequest, xedge1, xedge2;
  int search_edges = 1;
  if(shapepar == SHAPEPAR_W10) {
    yrequest = 0.1*ymax;
  }else if(shapepar == SHAPEPAR_W25) {
    yrequest = 0.25*ymax;
  }else if(shapepar == SHAPEPAR_W50) {
    yrequest = 0.5*ymax;
  }else if(shapepar == SHAPEPAR_W75) {
    yrequest = 0.75*ymax;
  }else if(shapepar == SHAPEPAR_W90) {
    yrequest = 0.9*ymax;
  }else {
    search_edges = 0;
  }
  if(search_edges) {
    calcVonMisesProfile_shape_parameter_get_edges(&xedge1, &xedge2, measurement, xmax, components, shift, phaseprecision, yrequest, nrbins_start, verbose);
    return;
  }
  int search_other_peak = 0;
  double startphase, endphase, amplitude_firstpeak;
  if(shapepar == SHAPEPAR_PEAKSEPPHASE || shapepar == SHAPEPAR_PEAKAMPRATIO || shapepar == SHAPEPAR_PEAKAMPRATIO_RECI) {
    search_other_peak = 1;
    if(shapepar_aux[1] > shapepar_aux[0]) {
      startphase = (xmax+shapepar_aux[0]);
      endphase = (xmax+shapepar_aux[1]);
    }else {
      startphase = (xmax+shapepar_aux[1]);
      endphase = (xmax+shapepar_aux[0]);
    }
  }else if(shapepar == SHAPEPAR_PEAKPHASE_RESTRICTED || shapepar == SHAPEPAR_PEAKAMPRATIO_RESTRICTED) {
    search_other_peak = 1;
    if(shapepar_aux[1] > shapepar_aux[0]) {
      startphase = (shapepar_aux[0]);
      endphase = (shapepar_aux[1]);
    }else {
      startphase = (shapepar_aux[1]);
      endphase = (shapepar_aux[0]);
    }
  }
  if(search_other_peak) {
    int found;
    calcVonMisesProfile_shape_parameter_search_peak_in_restricted_phase_range(&found, &x, startphase, endphase, components, shift, phaseprecision, nrbins_start, verbose);
    if(found == 0) {
      printwarning(verbose.debug, "WARNING calcVonMisesProfile_shape_parameter: No peak found in phase range %lf .. %lf. The phase of the main peak is at %lf. The requested shape parameter is undefined as a consequence.\n", startphase, endphase, xmax);
      *measurement = sqrt(-1);
      return;
    }else {
      if(shapepar == SHAPEPAR_PEAKSEPPHASE) {
 *measurement = x-xmax;
 return;
      }else if(shapepar == SHAPEPAR_PEAKPHASE_RESTRICTED) {
 *measurement = x;
 return;
      }else if(shapepar == SHAPEPAR_PEAKAMPRATIO || shapepar == SHAPEPAR_PEAKAMPRATIO_RECI || shapepar == SHAPEPAR_PEAKAMPRATIO_RESTRICTED) {
 double y;
 y = calcVonMisesFunction(components, x, shift);
 if(shapepar == SHAPEPAR_PEAKAMPRATIO_RECI) {
   *measurement = ymax/y;
   return;
 }else if(shapepar == SHAPEPAR_PEAKAMPRATIO) {
   *measurement = y/ymax;
   return;
 }else {
   amplitude_firstpeak = y;
 }
      }
    }
  }
  search_other_peak = 0;
  if(shapepar == SHAPEPAR_PEAKAMPRATIO_RESTRICTED) {
    search_other_peak = 1;
    if(shapepar_aux[3] > shapepar_aux[2]) {
      startphase = (shapepar_aux[2]);
      endphase = (shapepar_aux[3]);
    }else {
      startphase = (shapepar_aux[3]);
      endphase = (shapepar_aux[2]);
    }
  }
  if(search_other_peak) {
    int found;
    calcVonMisesProfile_shape_parameter_search_peak_in_restricted_phase_range(&found, &x, startphase, endphase, components, shift, phaseprecision, nrbins_start, verbose);
    if(found == 0) {
      printwarning(verbose.debug, "WARNING calcVonMisesProfile_shape_parameter: No peak found in phase range %lf .. %lf. The phase of the main peak is at %lf. The requested shape parameter is undefined as a consequence.\n", startphase, endphase, xmax);
      *measurement = sqrt(-1);
      return;
    }else {
      if(shapepar == SHAPEPAR_PEAKAMPRATIO_RESTRICTED) {
 double y;
 y = calcVonMisesFunction(components, x, shift);
 *measurement = amplitude_firstpeak/y;
 return;
      }
    }
  }
  if(shapepar == SHAPEPAR_W10_MAXAMP || shapepar == SHAPEPAR_W25_MAXAMP || shapepar == SHAPEPAR_W50_MAXAMP || shapepar == SHAPEPAR_W75_MAXAMP || shapepar == SHAPEPAR_W90_MAXAMP) {
    int nmax = 0;
    double amp_max = components->height[0];
    if(components->nrcomponents > 1) {
      int n;
      for(n = 1; n < components->nrcomponents; n++) {
 if(components->height[n] > amp_max) {
   nmax = n;
   amp_max = components->height[n];
 }
      }
    }
    if(shapepar == SHAPEPAR_W10_MAXAMP) {
      *measurement = widthVonMisesFunction2(components->concentration[nmax], 0.1)/(2.0*M_PI);
    }else if(shapepar == SHAPEPAR_W25_MAXAMP) {
      *measurement = widthVonMisesFunction2(components->concentration[nmax], 0.25)/(2.0*M_PI);
    }else if(shapepar == SHAPEPAR_W50_MAXAMP) {
      *measurement = widthVonMisesFunction2(components->concentration[nmax], 0.5)/(2.0*M_PI);
    }else if(shapepar == SHAPEPAR_W75_MAXAMP) {
      *measurement = widthVonMisesFunction2(components->concentration[nmax], 0.75)/(2.0*M_PI);
    }else if(shapepar == SHAPEPAR_W90_MAXAMP) {
      *measurement = widthVonMisesFunction2(components->concentration[nmax], 0.9)/(2.0*M_PI);
    }
    if(*measurement < 0) {
      *measurement = sqrt(-1);
    }
    return;
  }
  printerror(verbose.debug, "ERROR calcVonMisesProfile_shape_parameter: Requested shape parameter is not implemented");
  *measurement = sqrt(-1);
  return;
}
void print_shape_par(FILE *fout, int showdescr, int shapepar, double measurement, double error)
{
  if(showdescr) {
    switch(shapepar) {
    case SHAPEPAR_PEAKPHASE:
    case SHAPEPAR_PEAKPHASE_RESTRICTED:
      fprintf(fout, "Pulse phase of peak: ");
      break;
    case SHAPEPAR_PEAKAMP: fprintf(fout, "Amplitude of peak: "); break;
    case SHAPEPAR_W10: fprintf(fout, "Width at 10%% of maximum amplitude (in phase): "); break;
    case SHAPEPAR_W25: fprintf(fout, "Width at 25%% of maximum amplitude (in phase): "); break;
    case SHAPEPAR_W50: fprintf(fout, "Width at 50%% of maximum amplitude (in phase): "); break;
    case SHAPEPAR_W75: fprintf(fout, "Width at 75%% of maximum amplitude (in phase): "); break;
    case SHAPEPAR_W90: fprintf(fout, "Width at 90%% of maximum amplitude (in phase): "); break;
    case SHAPEPAR_W10_MAXAMP: fprintf(fout, "Width at 10%% of amplitude vonMises tallest component (in phase): "); break;
    case SHAPEPAR_W25_MAXAMP: fprintf(fout, "Width at 25%% of amplitude vonMises tallest component (in phase): "); break;
    case SHAPEPAR_W50_MAXAMP: fprintf(fout, "Width at 50%% of amplitude vonMises tallest component (in phase): "); break;
    case SHAPEPAR_W75_MAXAMP: fprintf(fout, "Width at 75%% of amplitude vonMises tallest component (in phase): "); break;
    case SHAPEPAR_W90_MAXAMP: fprintf(fout, "Width at 90%% of amplitude vonMises tallest component (in phase): "); break;
    case SHAPEPAR_PEAKSEPPHASE: fprintf(fout, "Peak separation (in phase): "); break;
    case SHAPEPAR_PEAKAMPRATIO:
    case SHAPEPAR_PEAKAMPRATIO_RESTRICTED:
      fprintf(fout, "Peak amplitude ratio: ");
      break;
    case SHAPEPAR_COMPAMPRATIO: fprintf(fout, "Component amplitude ratio: "); break;
    case SHAPEPAR_AMPRATIO_ATPHASE: fprintf(fout, "Amplitude ratio at fixed phases: "); break;
    case SHAPEPAR_PEAKAMPRATIO_RECI: fprintf(fout, "Reciprocal peak amplitude ratio: "); break;
    default: printerror(0, "ERROR print_shape_par: Requested shape parameter is not implemented"); break;
    }
  }
  if(shapepar == SHAPEPAR_PEAKPHASE || shapepar == SHAPEPAR_PEAKPHASE_RESTRICTED || shapepar == SHAPEPAR_W10 || shapepar == SHAPEPAR_W25 || shapepar == SHAPEPAR_W50 || shapepar == SHAPEPAR_W75 || shapepar == SHAPEPAR_W90 || shapepar == SHAPEPAR_W10_MAXAMP || shapepar == SHAPEPAR_W25_MAXAMP || shapepar == SHAPEPAR_W50_MAXAMP || shapepar == SHAPEPAR_W75_MAXAMP || shapepar == SHAPEPAR_W90_MAXAMP || shapepar == SHAPEPAR_PEAKSEPPHASE) {
    if(error >= 0) {
      fprintf(fout, "%e +- %e (phase) = %e +- %e (deg)\n", measurement, error, measurement*360.0, error*360.0);
    }else {
      fprintf(fout, "%e (phase) = %e (deg)\n", measurement, measurement*360.0);
    }
  }else {
    if(error >= 0) {
      fprintf(fout, "%e +- %e\n", measurement, error);
    }else {
      fprintf(fout, "%e\n", measurement);
    }
  }
}

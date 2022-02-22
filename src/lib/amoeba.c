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
#include <stdlib.h>
#include "psrsalsa.h"
static int nrparams_internal_psrsalsa;
static int *fixed_internal_psrsalsa;
static float *xstart_internal_psrsalsa;
static float *x_internal_psrsalsa;
static float (*funk_remember_user_function)(float []);
static int algorithm_internal_psrsalsa;
float funk_internal_psrsalsa(float x[])
{
  int i, j;
  if(algorithm_internal_psrsalsa == 1)
    j = 1;
  else
    j = 0;
  for(i = 0; i < nrparams_internal_psrsalsa; i++) {
    if(fixed_internal_psrsalsa[i] == 0) {
      x_internal_psrsalsa[i+1] = x[j++];
    }else {
      x_internal_psrsalsa[i+1] = xstart_internal_psrsalsa[i];
    }
  }
  return funk_remember_user_function(x_internal_psrsalsa+1);
}
int doAmoeba(int algorithm, float *xstart, float *dx, int *fixed, float *xfit, float *yfit, int nrparams, float (*funk)(float []), float ftol, int *nfunk, int verbose, int finderrors, float sigma, float *dplus, float *dmin)
{
  int i, j, nfitparameters, ret;
  extern float amoeba_nmsimplex(float (*objfunc)(float[]), float start[], float dx[], int n, float EPSILON, int *nritterations, float *reachedEpsilon, int verbose);
#ifdef NRAVAIL
  float **p_nr, *y_nr;
  extern int amoeba_nr(float **p, float y[], int ndim, float ftol, float (*funk)(float []), int *nfunk);
#endif
  if(algorithm < 0 || algorithm > 1) {
    fprintf(stderr, "ERROR doAmoeba: Unknown algorithm requested.\n");
    return 4;
  }
#ifndef NRAVAIL
  if(algorithm == 1) {
    fprintf(stderr, "ERROR doAmoeba: Numerical Recipies requested, but not included during compilation.\n");
    return 4;
  }
#endif
  nfitparameters = 0;
  for(i = 0; i < nrparams; i++) {
    if(fixed[i] == 0)
      nfitparameters++;
  }
  if(nfitparameters < 2) {
    fprintf(stderr, "ERROR doAmoeba: Cannot fit for less than 2 parameters (now have %d).\n", nfitparameters);
    return 3;
  }
  fixed_internal_psrsalsa = fixed;
  xstart_internal_psrsalsa = xstart;
  nrparams_internal_psrsalsa = nrparams;
  funk_remember_user_function = funk;
  algorithm_internal_psrsalsa = algorithm;
  if(algorithm == 0) {
    float reachedEpsilon, *xstart_nmsimplex, *dx_nmsimplex;
    xstart_nmsimplex = malloc(nfitparameters*sizeof(float));
    dx_nmsimplex = malloc(nfitparameters*sizeof(float));
    x_internal_psrsalsa = malloc((nrparams+1)*sizeof(float));
    if(xstart_nmsimplex == NULL || dx_nmsimplex == NULL || x_internal_psrsalsa == NULL) {
      fprintf(stderr, "ERROR doAmoeba: Memory allocation error.\n");
      return 2;
    }
    j = 0;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
 xstart_nmsimplex[j] = xstart[i];
 dx_nmsimplex[j] = dx[i];
 j++;
      }
    }
    *yfit = amoeba_nmsimplex(funk_internal_psrsalsa, xstart_nmsimplex, dx_nmsimplex, nfitparameters, ftol, nfunk, &reachedEpsilon, 0);
    j = 0;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
 xfit[i] = xstart_nmsimplex[j];
 j++;
      }else {
 xfit[i] = xstart[i];
      }
    }
    free(xstart_nmsimplex);
    free(dx_nmsimplex);
    free(x_internal_psrsalsa);
    if(reachedEpsilon > ftol)
      return 1;
  }
#ifdef NRAVAIL
  int n;
  if(algorithm == 1) {
    p_nr = matrix(1, nfitparameters+1, 1, nfitparameters);
    y_nr = vector(1, nfitparameters+1);
    x_internal_psrsalsa = vector(1, nrparams+1);
    if(p_nr == NULL || y_nr == NULL || x_internal_psrsalsa == NULL) {
      fprintf(stderr, "ERROR doAmoeba: Memory allocation error.\n");
      return 2;
    }
    j = 1;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
 for(n = 1; n <= nfitparameters+1; n++)
   p_nr[n][j] = xstart[i];
 j++;
      }
    }
    j = 1;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
 p_nr[j+1][j] += dx[i];
 j++;
      }
    }
    for(i = 0; i < nfitparameters+1; i++) {
      y_nr[i+1] = funk_internal_psrsalsa(p_nr[i+1]);
    }
    if(verbose) {
      for(j = 1; j <= nfitparameters; j++) {
 if(j == 1)
   fprintf(stderr, "matrix p: ");
 else
   fprintf(stderr, "          ");
 for(n = 1; n <= nfitparameters+1; n++) {
   fprintf(stderr, "%f ", p_nr[n][j]);
 }
 fprintf(stderr, "\n");
      }
    }
    ret = amoeba_nr(p_nr, y_nr, nfitparameters, ftol, funk_internal_psrsalsa, nfunk);
    if(ret != 0)
      return ret;
    *yfit = funk_internal_psrsalsa(p_nr[1]);
    j = 1;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
 xfit[i] = p_nr[1][j];
 j++;
      }else {
 xfit[i] = xstart[i];
      }
    }
    free_matrix(p_nr,1, nfitparameters+1, 1, nfitparameters);
    free_vector(y_nr, 1, nfitparameters+1);
    free_vector(x_internal_psrsalsa, 1, nfitparameters+1);
  }
#endif
  for(i = 0; i < nrparams; i++) {
    if(dplus != NULL)
      dplus[i] = 0;
    if(dmin != NULL)
      dmin[i] = 0;
  }
  if(finderrors && dplus != NULL && dmin != NULL) {
    if(nfitparameters < 3) {
      fprintf(stderr, "ERROR doAmoeba: Cannot estimate errors if the number of fit parameters is less than 3.\n");
    }else {
      for(i = 0; i < nrparams; i++) {
 ret = find_errors_amoeba(algorithm, dx, fixed, xfit, *yfit, nrparams, funk, ftol, i, &dplus[i], &dmin[i], sigma);
 if(ret != 0) {
   fprintf(stderr, "ERROR doAmoeba: find_errors_amoeba failed with error code %d\n", ret);
   return ret;
 }
      }
    }
  }
  return 0;
}
int find_errors_amoeba(int algorithm, float *dx, int *fixed, float *xfit, float yfit, int nrparams, float (*funk)(float []), float ftol, int paramnr, float *dplus, float *dmin, float sigma)
{
  float *xstartnew, *xfitnew, *dxnew, yfitnew, x0, x1, yold, step, fsign;
  int *fixednew, j, nfunknew, ret, sign;
  long n, n2;
  fixednew = malloc(nrparams*sizeof(int));
  xstartnew = malloc(nrparams*sizeof(float));
  dxnew = malloc(nrparams*sizeof(float));
  xfitnew = malloc(nrparams*sizeof(float));
  if(fixednew == NULL || xstartnew == NULL || xfitnew == NULL) {
    fprintf(stderr, "ERROR find_point_amoeba: Memory allocation error.\n");
    return 2;
  }
  if(fixed[paramnr] == 0) {
    for(j = 0; j < nrparams; j++) {
      fixednew[j] = fixed[j];
      if(j == paramnr)
 fixednew[j] = 1;
    }
    for(sign = 0; sign < 2; sign++) {
      if(sign == 0)
 fsign = 1;
      else
 fsign = -1;
      for(j = 0; j < nrparams; j++) {
 xstartnew[j] = xfit[j];
 dxnew[j] = dx[j];
      }
      x0 = xstartnew[paramnr];
      step = fsign*dx[paramnr];
      yold = yfit;
      n = 0;
      n2 = 0;
      do {
 n++;
 x1 = x0 + step;
 xstartnew[paramnr] = x1;
 dxnew[paramnr] = step;
 do {
   ret = doAmoeba(algorithm, xstartnew, dxnew, fixednew, xfitnew, &yfitnew, nrparams, funk, ftol, &nfunknew, 0, 0, sigma, NULL, NULL);
   if(ret == 3) {
     free(fixednew);
     free(xstartnew);
     free(dxnew);
     free(xfitnew);
     return 3;
   }
   if(ret == 1) {
     fprintf(stderr, "WARNING find_point_amoeba: Adjusting amoeba tollerance to try to converge (for errorbar estimation).\n");
     ftol *= 10;
   }
   if(ret == 0) {
     break;
   }
 }while(ftol < 0.01);
 if(ret == 1) {
   free(fixednew);
   free(xstartnew);
   free(dxnew);
   free(xfitnew);
   return 1;
 }
 if(yfitnew < yfit*(sigma+1.0) && yold < yfit*(sigma+1.0)) {
   x0 = x1;
 }else if(yfitnew > yfit*(sigma+1.0) && yold < yfit*(sigma+1.0)) {
   step = 0.5*step;
 }else if(yfitnew < yfit*(sigma+1.0) && yold > yfit*(sigma+1.0)) {
   step = 0.5*step;
 }else {
   x0 = x1;
 }
 if(n % 100 == 0) {
   if(fabs((yfitnew-yfit*(sigma+1.0))/yfit*(sigma+1.0)) < 0.05) {
     fprintf(stderr, "ERROR find_errors_amoeba: error parameter %d didn't quite converge, but possibly it is close enough: chi %f != %f\n", paramnr, yfitnew, yfit*(sigma+1.0));
     yfitnew = yfit*(sigma+1.0);
   }else {
     x0 = xfit[paramnr] + (x0-xfit[paramnr])*0.987654321;
     step = fsign*dx[paramnr];
     n2 ++;
   }
 }
      }while(fabs((yfitnew-yfit*(sigma+1.0))/yfit*(sigma+1.0)) > 0.001 && n2 < 10);
      if(n2 == 10) {
 fprintf(stderr, "ERROR find_errors_amoeba: converging error parameter %d: chi %f != %f\n", paramnr, yfitnew, yfit*(sigma+1.0));
 x0 = sqrt(-1);
      }
      if(sign == 0)
 *dplus = x0 - xfit[paramnr];
      else
 *dmin = x0 - xfit[paramnr];
    }
  }else {
    *dplus = 0;
    *dmin = 0;
  }
  free(fixednew);
  free(xstartnew);
  free(dxnew);
  free(xfitnew);
  return 0;
}

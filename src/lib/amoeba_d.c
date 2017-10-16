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
#include <math.h>
#include "psrsalsa.h"
static int nrparams_internal_psrsalsa_d;
static int *fixed_internal_psrsalsa_d;
static double *xstart_internal_psrsalsa_d;
static double *x_internal_psrsalsa_d;
static double (*funk_remember_user_function_d)(double []);
static int algorithm_internal_psrsalsa_d;






double funk_internal_psrsalsa_d(double x[])
{
  int i, j;
  if(algorithm_internal_psrsalsa_d == 1)
    j = 1;
  else
    j = 0;
  for(i = 0; i < nrparams_internal_psrsalsa_d; i++) {
    if(fixed_internal_psrsalsa_d[i] == 0) {
      x_internal_psrsalsa_d[i+1] = x[j++];
    }else {
      x_internal_psrsalsa_d[i+1] = xstart_internal_psrsalsa_d[i];
    }
  }
  return funk_remember_user_function_d(x_internal_psrsalsa_d+1);
}
int doAmoeba_d(int algorithm, double *xstart, double *dx, int *fixed, double *xfit, double *yfit, int nrparams, double (*funk)(double []), double ftol, int *nfunk, int verbose, int finderrors, double sigma, double *dplus, double *dmin)
{
  int i, j, nfitparameters, ret;
  extern double amoeba_nmsimplex_d(double (*objfunc)(double[]), double start[], double dx[], int n, double EPSILON, int *nritterations, double *reachedEpsilon, int verbose);
  if(algorithm < 0 || algorithm > 1) {
    fprintf(stderr, "ERROR doAmoeba_d: Unknown algorithm requested.\n");
    return 4;
  }
  nfitparameters = 0;
  for(i = 0; i < nrparams; i++) {
    if(fixed[i] == 0)
      nfitparameters++;
  }
  if(nfitparameters < 2) {
    fprintf(stderr, "ERROR doAmoeba_d: Cannot fit for less than 2 parameters (now have %d).\n", nfitparameters);
    return 3;
  }
  fixed_internal_psrsalsa_d = fixed;
  xstart_internal_psrsalsa_d = xstart;
  nrparams_internal_psrsalsa_d = nrparams;
  funk_remember_user_function_d = funk;
  algorithm_internal_psrsalsa_d = algorithm;
  if(algorithm == 0) {
    double reachedEpsilon, *xstart_nmsimplex_d, *dx_nmsimplex_d;
    xstart_nmsimplex_d = malloc(nfitparameters*sizeof(double));
    dx_nmsimplex_d = malloc(nfitparameters*sizeof(double));
    x_internal_psrsalsa_d = malloc((nrparams+1)*sizeof(double));
    if(xstart_nmsimplex_d == NULL || dx_nmsimplex_d == NULL || x_internal_psrsalsa_d == NULL) {
      fprintf(stderr, "ERROR doAmoeba_d: Memory allocation error.\n");
      return 2;
    }
    j = 0;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
 xstart_nmsimplex_d[j] = xstart[i];
 dx_nmsimplex_d[j] = dx[i];
 j++;
      }
    }
    *yfit = amoeba_nmsimplex_d(funk_internal_psrsalsa_d, xstart_nmsimplex_d, dx_nmsimplex_d, nfitparameters, ftol, nfunk, &reachedEpsilon, 0);
    j = 0;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
 xfit[i] = xstart_nmsimplex_d[j];
 j++;
      }else {
 xfit[i] = xstart[i];
      }
    }
    free(xstart_nmsimplex_d);
    free(dx_nmsimplex_d);
    free(x_internal_psrsalsa_d);
    if(reachedEpsilon > ftol)
      return 1;
  }
  for(i = 0; i < nrparams; i++) {
    if(dplus != NULL)
      dplus[i] = 0;
    if(dmin != NULL)
      dmin[i] = 0;
  }
  if(finderrors && dplus != NULL && dmin != NULL) {
    if(nfitparameters < 3) {
      fprintf(stderr, "ERROR doAmoeba_d: Cannot estimate errors if the number of fit parameters is less than 3.\n");
    }else {
      for(i = 0; i < nrparams; i++) {
 ret = find_errors_amoeba_d(algorithm, dx, fixed, xfit, *yfit, nrparams, funk, ftol, i, &dplus[i], &dmin[i], sigma);
 if(ret != 0) {
   fprintf(stderr, "ERROR doAmoeba_d: find_errors_amoeba_d failed with error code %d\n", ret);
   return ret;
 }
      }
    }
  }
  return 0;
}
int find_errors_amoeba_d(int algorithm, double *dx, int *fixed, double *xfit, double yfit, int nrparams, double (*funk)(double []), double ftol, int paramnr, double *dplus, double *dmin, double sigma)
{
  double *xstartnew, *xfitnew, *dxnew, yfitnew, x0, x1, yold, step, fsign;
  int *fixednew, j, nfunknew, ret, sign;
  long n, n2;
  fixednew = malloc(nrparams*sizeof(int));
  xstartnew = malloc(nrparams*sizeof(double));
  dxnew = malloc(nrparams*sizeof(double));
  xfitnew = malloc(nrparams*sizeof(double));
  if(fixednew == NULL || xstartnew == NULL || xfitnew == NULL) {
    fprintf(stderr, "ERROR find_point_amoeba_d: Memory allocation error.\n");
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
   ret = doAmoeba_d(algorithm, xstartnew, dxnew, fixednew, xfitnew, &yfitnew, nrparams, funk, ftol, &nfunknew, 0, 0, sigma, NULL, NULL);
   if(ret == 3) {
     free(fixednew);
     free(xstartnew);
     free(dxnew);
     free(xfitnew);
     return 3;
   }
   if(ret == 1) {
     fprintf(stderr, "WARNING find_point_amoeba_d: Adjusting amoeba tollerance to try to converge (for errorbar estimation).\n");
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
     fprintf(stderr, "ERROR find_errors_amoeba_d: error parameter %d didn't quite converge, but possibly it is close enough: chi %lf != %lf\n", paramnr, yfitnew, yfit*(sigma+1.0));
     yfitnew = yfit*(sigma+1.0);
   }else {
     x0 = xfit[paramnr] + (x0-xfit[paramnr])*0.987654321;
     step = fsign*dx[paramnr];
     n2 ++;
   }
 }
      }while(fabs((yfitnew-yfit*(sigma+1.0))/yfit*(sigma+1.0)) > 0.001 && n2 < 10);
      if(n2 == 10) {
 fprintf(stderr, "ERROR find_errors_amoeba_d: converging error parameter %d: chi %lf != %lf\n", paramnr, yfitnew, yfit*(sigma+1.0));
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

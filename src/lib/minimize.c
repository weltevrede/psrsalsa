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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include "gsl/gsl_roots.h"
#include "psrsalsa.h"
int minimize_1D_double_refine_borders(int findroot, double (*funk)(double, void *), void *params, int gridsearch, int investigateLocalMinima, double *x_lower, double *x_upper, int debug_verbose)
{
  double x, y, ymin[3], ymin_new[3], dy[3], ymax, x_lower_new, x_upper_new, sign_old, gridsearch_margin, imin_new[3], i_originalgrid, imin1, imin2;
  int imin[3], unset[3], i, minima;
  unset[0] = unset[1] = unset[2] = 1;
  for(i = 0; i < gridsearch; i++) {
    x = *x_lower + i*(*x_upper - *x_lower)/(double)(gridsearch-1);
    y = funk(x, params);
    if(debug_verbose) {
      fflush(stdout); fprintf(stderr, "x=%f y=%f\n", x, y);
    }
    if(findroot) {
      if(i == 0) {
 sign_old = y;
 imin[0] = 0;
      }else {
 if(y < 0 && sign_old > 0) {
   imin[0] = i;
 }else if(y > 0 && sign_old < 0) {
   imin[0] = i;
 }
 if(imin[0] != 0)
   break;
      }
    }else {
      if(unset[0]) {
 imin[0] = i;
 ymin[0] = y;
 unset[0] = 0;
      }else if(unset[1]) {
 imin[1] = i;
 ymin[1] = y;
 unset[1] = 0;
      }else if(unset[2]) {
 imin[2] = i;
 ymin[2] = y;
 unset[2] = 0;
      }else {
 for(minima = 0; minima < 3; minima++) {
   dy[minima] = ymin[minima] - y;
 }
 if(dy[0] > dy[1] && dy[0] > dy[2]) {
   if(dy[0] > 0) {
     imin[0] = i;
     ymin[0] = y;
   }
 }else if(dy[1] > dy[0] && dy[1] > dy[2]) {
   if(dy[1] > 0) {
     imin[1] = i;
     ymin[1] = y;
   }
 }else {
   if(dy[2] > 0) {
     imin[2] = i;
     ymin[2] = y;
   }
 }
      }
      if(i == 0 || y > ymax) {
 ymax = y;
      }
    }
  }
  if(findroot == 0) {
    imin1 = imin[0];
    imin2 = imin[0];
    if(imin[1] < imin1)
      imin1 = imin[1];
    if(imin[2] < imin1)
      imin1 = imin[2];
    if(imin[1] > imin2)
      imin2 = imin[1];
    if(imin[2] > imin2)
      imin2 = imin[2];
    if(imin2-imin1 != 2 && debug_verbose) {
      fflush(stdout); fprintf(stderr, "Found more than one possible minima (range = %.2f .. %.2f)\n", imin1, imin2);
    }
    if(imin2-imin1 > 0.3*gridsearch) {
      if(investigateLocalMinima == 0) {
 imin1 = imin[0];
 imin2 = imin[0];
 if(ymin[1] < ymin[0]) {
   imin1 = imin[1];
   imin2 = imin[1];
 }else if(ymin[2] < ymin[0] && ymin[2] < ymin[1]) {
   imin1 = imin[2];
   imin2 = imin[2];
 }
      }else {
 if(debug_verbose) {
   fflush(stdout); fprintf(stderr, "Refining gridsearch around local minima\n");
 }
 unset[0] = unset[1] = unset[2] = 1;
 for(minima = 0; minima < 3; minima++) {
   for(x = *x_lower + (imin[minima]-1.5)*(*x_upper - *x_lower)/(double)(gridsearch-1); x < *x_lower + (imin[minima]+1.5)*(*x_upper - *x_lower)/(double)(gridsearch-1); x += 0.09*(*x_upper - *x_lower)/(double)(gridsearch-1)) {
     y = funk(x, params);
     if(debug_verbose) {
       fflush(stdout); fprintf(stderr, "x=%f y=%f\n", x, y);
     }
     i_originalgrid = (x - *x_lower)*(double)(gridsearch-1)/(*x_upper - *x_lower);
     if(unset[0] || y < ymin_new[0]) {
       imin_new[0] = i_originalgrid;
       ymin_new[0] = y;
       unset[0] = 0;
     }else if(unset[1] || y < ymin_new[1]) {
       imin_new[1] = i_originalgrid;
       ymin_new[1] = y;
       unset[1] = 0;
     }else if(unset[2] || y < ymin_new[2]) {
       imin_new[2] = i_originalgrid;
       ymin_new[2] = y;
       unset[2] = 0;
     }
   }
 }
 imin1 = imin_new[0];
 imin2 = imin_new[0];
 if(imin_new[1] < imin1)
   imin1 = imin_new[1];
 if(imin_new[2] < imin1)
   imin1 = imin_new[2];
 if(imin_new[1] > imin2)
   imin2 = imin_new[1];
 if(imin_new[2] > imin2)
   imin2 = imin_new[2];
 if(imin2-imin1 != 2 && debug_verbose) {
   fflush(stdout); fprintf(stderr, "Found more than one possible minima (range = %.2f .. %.2f)\n", imin1, imin2);
 }
 if(imin2-imin1 > 0.3*gridsearch) {
   fflush(stdout); fprintf(stderr, "minimize_1D_double_refine_borders: Refining gridsearch around local minima failed\n");
   exit(0);
 }
      }
    }
  }
  gridsearch_margin = 0.1*(double)gridsearch;
  if(gridsearch_margin < 1)
    gridsearch_margin = 1;
  if(findroot) {
    if(imin[0] == 0) {
      return 2;
    }
    x_lower_new = *x_lower + (imin[0]-gridsearch_margin)*(*x_upper - *x_lower)/(double)(gridsearch-1);
    x_upper_new = *x_lower + (imin[0]+gridsearch_margin-1)*(*x_upper - *x_lower)/(double)(gridsearch-1);
    *x_lower = x_lower_new;
    *x_upper = x_upper_new;
  }else {
    x_lower_new = *x_lower + (imin1-gridsearch_margin)*(*x_upper - *x_lower)/(double)(gridsearch-1);
    x_upper_new = *x_lower + (imin2+gridsearch_margin)*(*x_upper - *x_lower)/(double)(gridsearch-1);
    *x_lower = x_lower_new;
    *x_upper = x_upper_new;
  }
  return 0;
}
int minimize_1D_double(int findroot, double (*funk)(double, void *), void *params, double x_lower, double x_upper, int gridsearch, int investigateLocalMinima, int nested, double *x_minimum, int max_iter, double epsabs, double epsrel, int verbose, int debug_verbose)
{
  const gsl_min_fminimizer_type *T_minimizer;
  gsl_min_fminimizer *s_minimizer;
  const gsl_root_fsolver_type *T_root;
  gsl_root_fsolver *s_root;
  int status;
  int iter, i, ret, nest;
  gsl_function F;
  F.function = funk;
  F.params = params;
  iter = 0;
#if GSL_VERSION_NUMBER < 102
  printerror(0, "ERROR minimize_1D_double: Not supported for GSL < 1.2");
  exit(0);
#endif
  if(verbose) {
    for(i = 0; i < verbose - 1; i++)
      printf(" ");
    if(findroot == 0)
      printf("Minimising function between %f and %f\n", x_lower, x_upper);
    else
      printf("Finding root of function between %f and %f\n", x_lower, x_upper);
  }
  if(gridsearch > 1) {
    if(verbose) {
      for(i = 0; i < verbose - 1; i++)
 printf(" ");
      printf("  Refining boundaries with a %d point grid search\n", gridsearch);
    }
    for(nest = 0; nest < 1+nested; nest++) {
      ret = minimize_1D_double_refine_borders(findroot, funk, params, gridsearch, investigateLocalMinima, &x_lower, &x_upper, debug_verbose);
      if(ret != 0) {
 if(ret == 1) {
   if(verbose) {
     for(i = 0; i < verbose - 1; i++)
       printf(" ");
     printf("  Maximum itterations (%d) exceeded while refining borders\n", max_iter);
   }
 }
 return ret;
      }
      if(verbose) {
 for(i = 0; i < verbose - 1; i++)
   printf(" ");
 printf("  Refined boundaries are %f and %f\n", x_lower, x_upper);
      }
    }
  }
  *x_minimum = 0.5*(x_upper+x_lower);
  if(findroot == 0) {
    T_minimizer = gsl_min_fminimizer_brent;
    s_minimizer = gsl_min_fminimizer_alloc(T_minimizer);
    fflush(stdout);
    gsl_min_fminimizer_set(s_minimizer, &F, *x_minimum, x_lower, x_upper);
  }else {
    T_root = gsl_root_fsolver_brent;
    s_root = gsl_root_fsolver_alloc(T_root);
    if(gridsearch == 0) {
      double val1, val2;
      val1 = funk(x_lower, params);
      val2 = funk(x_upper, params);
      if((val1 > 0 && val2 > 0) || (val1 < 0 && val2 < 0)) {
 return 3;
      }
    }
    gsl_root_fsolver_set(s_root, &F, x_lower, x_upper);
  }
  do {
    iter++;
    if(findroot == 0) {
      status = gsl_min_fminimizer_iterate(s_minimizer);
#if GSL_VERSION_NUMBER >= 102
      *x_minimum = gsl_min_fminimizer_x_minimum (s_minimizer);
#endif
      x_lower = gsl_min_fminimizer_x_lower(s_minimizer);
      x_upper = gsl_min_fminimizer_x_upper(s_minimizer);
    }else {
      status = gsl_root_fsolver_iterate(s_root);
      *x_minimum = gsl_root_fsolver_root(s_root);
      x_lower = gsl_root_fsolver_x_lower(s_root);
      x_upper = gsl_root_fsolver_x_upper(s_root);
    }
    status = gsl_min_test_interval (x_lower, x_upper, epsabs, epsrel);
  }while(status == GSL_CONTINUE && iter < max_iter);
  if(findroot)
    gsl_root_fsolver_free(s_root);
  else
    gsl_min_fminimizer_free(s_minimizer);
  if(status == GSL_SUCCESS) {
    if(verbose) {
      for(i = 0; i < verbose - 1; i++)
 printf(" ");
      if(findroot)
 printf("  Root found at %f with error %f\n", *x_minimum, fabs(x_upper - x_lower));
      else
 printf("  Minimum found at %f with error %f\n", *x_minimum, fabs(x_upper - x_lower));
    }
    return 0;
  }else if(iter == max_iter) {
    if(verbose) {
      for(i = 0; i < verbose - 1; i++)
 printf(" ");
      printf("  Maximum itterations (%d) exceeded. Value confined to [%f, %f]\n", max_iter, x_lower, x_upper);
    }
    return 1;
  }else {
    if(verbose) {
      for(i = 0; i < verbose - 1; i++)
 printf(" ");
      printf("  Unknown error during minimization\n");
    }
    return 5;
  }
  return status;
}
double (*internal_funk_provided_by_user)(double *, void *);
int internal_paramnr;
double *internal_xminimum;
double internal_desired_chi2;
double internal_find_1D_error_funk(double x, void *params)
{
  double chi2;
  internal_xminimum[internal_paramnr] = x;
  chi2 = internal_funk_provided_by_user(internal_xminimum, params);
  chi2 -= internal_desired_chi2;
  return chi2;
}
int find_1D_error(double (*funk)(double *, void *), double *xminimum, int paramnr, int nrparameters, double dx, double dxmax, void *params, double sigma, double chi2min, int max_itr, double epsabs, double epsrel, double *errorbar, int verbose)
{
  int ittr, i, ret, debug_verbose;
  double x_lower, x_upper, diff, xval, *xminimum_fiddle;
  debug_verbose = 0;
  dx = fabs(dx);
  if(verbose) {
    for(i = 0; i < verbose - 1; i++)
      printf(" ");
    printf("Finding error at value %f with stepsize %f\n", xminimum[paramnr], dx);
  }
  xminimum_fiddle = malloc(nrparameters*sizeof(double));
  if(xminimum_fiddle == NULL) {
    fflush(stdout);
    fprintf(stderr, "ERROR find_1D_error: Memory allocation error\n");
    return 4;
  }
  memcpy(xminimum_fiddle, xminimum, nrparameters*sizeof(double));
  internal_paramnr = paramnr;
  internal_xminimum = xminimum_fiddle;
  internal_funk_provided_by_user = funk;
  internal_desired_chi2 = chi2min*(1+fabs(sigma));
  if(verbose) {
    for(i = 0; i < verbose - 1; i++)
      printf(" ");
    printf("  Mimimum chi2 = %f, so %f sigma point corresponds to %f\n", chi2min, fabs(sigma), internal_desired_chi2);
  }
  x_lower = xminimum[paramnr];
  x_upper = xminimum[paramnr];
  ittr = 0;
  diff = internal_find_1D_error_funk(x_upper, params);
  if(diff > 0) {
    fflush(stdout);
    fprintf(stderr, "ERROR find_1D_error: Function called with initial parameters outside specified sigma limit: chi2 = %f higher than sigma border (%f).\n", diff, internal_desired_chi2);
    exit(0);
  }
  do {
    if(sigma >= 0) {
      x_upper += dx;
      if(dxmax >= 0) {
 if(fabs(x_upper - xminimum[paramnr]) > dxmax) {
   *errorbar = dxmax;
   free(xminimum_fiddle);
   return 0;
 }
      }
      diff = internal_find_1D_error_funk(x_upper, params);
    }else {
      x_lower -= dx;
      if(dxmax >= 0) {
 if(fabs(x_lower - xminimum[paramnr]) > dxmax) {
   *errorbar = dxmax;
   free(xminimum_fiddle);
   return 0;
 }
      }
      diff = internal_find_1D_error_funk(x_lower, params);
    }
    ittr++;
    if(ittr == max_itr) {
      free(xminimum_fiddle);
      return 1;
    }
  }while(diff < 0);
  if(verbose) {
    for(i = 0; i < verbose - 1; i++)
      printf(" ");
    printf("  Found brackets: [%f %f]\n", x_lower, x_upper);
    verbose += 2;
  }
  ret = minimize_1D_double(1, internal_find_1D_error_funk, params, x_lower, x_upper, 0, 0, 0, &xval, max_itr, epsabs, epsrel, verbose, debug_verbose);
  *errorbar = fabs(xval - xminimum[paramnr]);
  if(verbose) {
    verbose -= 2;
    for(i = 0; i < verbose - 1; i++)
      printf(" ");
    if(ret == 0)
      printf("  Found errorbar: %f\n", *errorbar);
    else if(ret == 1)
      printf("  Maximum nr of itterations reached, did not converge fully\n");
    else if(ret == 2)
      printf("  Did not find root in specified range\n");
    else if(ret == 3)
      printf("  Lower and upper limit do not bracket a root\n");
    else {
      ret = 5;
      printf("  Other unspecified error\n");
    }
  }
  free(xminimum_fiddle);
  return ret;
}

/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_cdf.h>
#include "psrsalsa.h"



long randomUnsignedInt()
{
  time_t seconds;
  struct timeval precisetime;




  time(&seconds);
  gettimeofday(&precisetime,0x0);

  return (long)seconds*(long)precisetime.tv_usec;
}



void randomize_idnum(long *idnum)
{

  *idnum = -randomUnsignedInt();
}
int find_peak_correlation(float *data1, float *data2, int ndata, int zeropad, int circularpad, int duplicate, int *lag, float *correl_max, verbose_definition verbose)
{
  int i, lag_max;
  int npoints;
  float *paddata1, *paddata2, *ans, ans_max, ans_min;
  if(duplicate)
    duplicate = ndata;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("%d points in the data\n", ndata);
    if(zeropad != 0) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Padding at least %d points before and after data\n", zeropad);
    }
    if(duplicate) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Duplicating data to avoid wrap problems enabled.\n");
    }
  }
  i = (int) (log10(1.0 * (ndata+2*zeropad+duplicate))/0.301031);
  npoints = pow(2.0,(i+1));
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Going to zero-pad it to %d points\n", npoints);
  }
  paddata1 = (float *)malloc(npoints*sizeof(float));
  paddata2 = (float *)malloc(npoints*sizeof(float));
  ans = (float *)malloc(2*npoints*sizeof(float));
  if(paddata1 == NULL || paddata2 == NULL || ans == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR find_peak_correlation: Memory allocation error.");
    return 0;
  }
  zeropad = (npoints - duplicate - ndata)/2;
  for(i = 0; i < ndata; i++) {
    paddata1[i+zeropad] = data1[i];
    paddata2[i+zeropad] = data2[i];
    if(duplicate) {
      paddata1[i+ndata+zeropad] = data1[i];
      paddata2[i+ndata+zeropad] = data2[i];
    }
  }
  for(i = 0; i < zeropad; i++) {
    paddata1[i] = data1[i-zeropad+ndata];
    paddata2[i] = data2[i-zeropad+ndata];
  }
  for(i = ndata+duplicate+zeropad; i < npoints; i++) {
    paddata1[i] = data1[i-duplicate-zeropad-ndata];
    paddata2[i] = data2[i-duplicate-zeropad-ndata];
  }
  if(crosscorrelation_fft(paddata1, paddata2, npoints, ans, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR find_peak_correlation: Cross correlation failed.");
    return 0;
  }
  lag_max = 0;
  ans_max = ans[0];
  ans_min = ans[0];
  for(i = 1; i < npoints; i++) {
    if(ans[i] > ans_max) {
      ans_max = ans[i];
      lag_max = i;
    }
    if(ans[i] < ans_min)
      ans_min = ans[i];
  }
  if(lag_max >= npoints/2)
    lag_max -= npoints;
  *correl_max = ans_max/ans_min;
  *lag = lag_max;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Found a lag of %d (correlation %f higher)\n", *lag, *correl_max);
  }
  free(paddata1);
  free(paddata2);
  free(ans);
  return 1;
}
long calculate_bin_number(double x, double dx, double min_x, int centered_at_zero, double extra_phase)
{
  long bin, step, binzero;
  if(min_x < 0)
    step = -min_x/dx+10;
  else
    step = 0;
  bin = ( x+(step+0.5*centered_at_zero-extra_phase)*dx)/dx;
  binzero = (min_x+(step+0.5*centered_at_zero-extra_phase)*dx)/dx;
  bin -= binzero;
  return bin;
}
double calculate_bin_location(long binnr, double dx, double min_x, int centered_at_zero, double extra_phase)
{
  long step, binzero;
  double x;
  if(min_x < 0)
    step = -min_x/dx+10;
  else
    step = 0;
  binzero = (min_x+(step+0.5*centered_at_zero-extra_phase)*dx)/dx;
  binnr = binnr + binzero;
  x = binnr*dx - (step+0.5*centered_at_zero-extra_phase)*dx;
  x += 0.5*dx;
  return x;
}
double calculate_required_bin_width(double x, long binnr, double min_x, int centered_at_zero, double extra_phase, verbose_definition verbose)
{
  long lastbin, ok, timesinloop;
  double dx, offset;
  dx = (x - min_x)/(double)(binnr+0.5);
  timesinloop = 0;
  offset = 0;
  do {
    ok = 1;
    dx = (x - min_x)/(double)(binnr+0.5) + offset;
    lastbin = calculate_bin_number(x, dx, min_x, centered_at_zero, extra_phase);
    if(lastbin < binnr) {
      offset += (x-min_x-0.5*dx)/(double)(binnr+0.5) - dx;
      ok = 0;
    }else if(lastbin > binnr) {
      offset += (x-min_x+0.5*dx)/(double)(binnr+0.5) -dx;
      ok = 0;
    }
    timesinloop++;
    if(timesinloop > 10) {
      printerror(verbose.debug, "ERROR calculate_required_bin_width: Cannot find suitable binsize.\n");
      return dx;
    }
  }while(ok == 0);
  return dx;
}
int set_binning_histogram(double min_x_data, double max_x_data, int rangex_set, double rangex_min, double rangex_max, int nrbins_specified, long nrbins, int centered_at_zero, double extra_phase, double *min_x, double *max_x, double *dx, verbose_definition verbose)
{
  int reset;
  do {
    reset = 0;
    *min_x = min_x_data;
    *max_x = max_x_data;
    if(rangex_set) {
      *min_x = rangex_min;
      *max_x = rangex_max;
    }
    if(nrbins_specified) {
      *dx = calculate_required_bin_width(*max_x, nrbins-1, *min_x, centered_at_zero, extra_phase, verbose);
      if(verbose.verbose)
 fprintf(stdout, "Going to use binsize %e.\n", *dx);
    }
    {
      long firstbin;
      firstbin = calculate_bin_number(*min_x, *dx, *min_x, centered_at_zero, extra_phase);
      if(firstbin != 0) {
 printerror(verbose.debug, "ERROR set_binning_histogram: Expected first bin to be number zero (it is %ld).\n", firstbin);
 return 2;
      }
    }
    long i;
    double diff;
    i = calculate_bin_number(*max_x, *dx, *min_x, centered_at_zero, extra_phase);
    diff = *max_x - calculate_bin_location(i, *dx, *min_x, centered_at_zero, extra_phase);
    diff = diff/(*dx);
    if(verbose.debug) {
      printf("Current set max value (%e) is falling at phase=%e w.r.t. centre of last generated bin.\n", *max_x, diff);
    }
    if(diff > -0.501 && diff < -0.499) {
      if(rangex_set) {
 rangex_max -= 0.5*(*dx);
 printwarning(verbose.debug, "WARNING set_binning_histogram: Adjusting maximum value of the specified range to %e to avoid rounding errors. Going to reset choosen binning.", rangex_max);
 reset = 1;
      }
    }
    if(diff > 0.499 && diff < 0.501) {
      if(rangex_set) {
 rangex_max += 0.5*(*dx);
 printwarning(verbose.debug, "WARNING set_binning_histogram: Adjusting maximum value of the specified range to %e to avoid rounding errors. Going to reset choosen binning.", rangex_max);
 reset = 1;
      }else {
 *max_x += 0.5*(*dx);
 if(nrbins_specified) {
   printwarning(verbose.debug, "WARNING set_binning_histogram: Nr of bins might be different from what was requested to avoid rounding errors.");
 }
      }
    }
    i = calculate_bin_number(*min_x, *dx, *min_x, centered_at_zero, extra_phase);
    diff = *min_x - calculate_bin_location(i, *dx, *min_x, centered_at_zero, extra_phase);
    diff = diff/(*dx);
    if(verbose.debug) {
      printf("Current set min value (%e) is falling at phase=%e w.r.t. centre of first generated bin.\n", *min_x, diff);
    }
    if(diff > -0.501 && diff < -0.499) {
      if(rangex_set) {
 rangex_min -= 0.5*(*dx);
 printwarning(verbose.debug, "WARNING set_binning_histogram: Adjusting minimum value of the specified range to %e to avoid rounding errors. Going to reset choosen binning.", rangex_min);
 reset = 1;
      }else {
 *min_x -= 0.5*(*dx);
 if(nrbins_specified) {
   printwarning(verbose.debug, "WARNING set_binning_histogram: Nr of bins might be different from what was requested to avoid rounding errors.");
 }
      }
    }
    if(diff > 0.499 && diff < 0.501) {
      if(rangex_set) {
 rangex_min += 0.5*(*dx);
 printwarning(verbose.debug, "WARNING set_binning_histogram: Adjusting minimum value of the specified range to %e to avoid rounding errors. Going to reset choosen binning.", rangex_min);
 reset = 1;
      }
    }
    if(reset) {
      printwarning(verbose.debug, "WARNING set_binning_histogram: Re-adjusting choosen binning.\n", rangex_max);
    }
  }while(reset == 1);
  return 0;
}
double kstest_cdf_flat(double x, double min_x, double max_x)
{
  if(x <= min_x)
    return 0;
  if(x >= max_x)
    return 1;
  return (x-min_x)/(max_x-min_x);
}
double kstest_cdf_sin(double x)
{
  if(x <= 0)
    return 0;
  if(x >= 90)
    return 1;
  return 1-cos(x*M_PI/180.0);
}
void kstest(double *data1, long n1, double *data2, long n2, int cdf_type, double (*cdf)(double), double *max_diff, double *prob, verbose_definition verbose)
{
  long i1, i2;
  double effective_n, ks_statistic, sign, cur_term, last_term, coeff;
  int converged;
  gsl_sort(data1, 1, n1);
  if(n2 > 0 && data2 != NULL)
    gsl_sort(data2, 1, n2);
  *max_diff = 0;
  if(n2 > 0 && data2 != NULL) {
    i1 = 0;
    i2 = 0;
    while(i1 < n1 && i2 < n2) {
      double diff;
      if(data1[i1] == data2[i2]) {
 i1++;
 i2++;
      }else if(data1[i1] < data2[i2]) {
 i1++;
      }else {
 i2++;
      }
      diff = fabs(i1/(double)n1 - i2/(double)n2);
      if(diff > *max_diff)
 *max_diff = diff;
    }
    effective_n=n1*n2/(double)(n1+n2);
  }else {
    double cdf_right, cdf_left, cdf_model;
    cdf_left = 0;
    for(i1 = 0; i1 < n1; i1++) {
      cdf_right = (i1+1)/(double)n1;
      if(cdf_type == 0) {
 cdf_model = (*cdf)(data1[i1]);
      }else if(cdf_type == 1) {
 cdf_model = kstest_cdf_flat(data1[i1], data1[0], data1[n1-1]);
 if(i1 == 0) {
   printwarning(verbose.debug, "WARNING kstest: Probability will be overestimated, since the minimum/maximum of the uniform distribution is based on the input values.");
 }
      }else if(cdf_type == 2) {
 cdf_model = kstest_cdf_sin(data1[i1]);
      }else {
 fflush(stdout);
 printerror(verbose.debug, "ERROR kstest: Undefined type of cdf is specified.");
 exit(0);
      }
      if(fabs(cdf_right-cdf_model) > *max_diff)
 *max_diff = fabs(cdf_right-cdf_model);
      if(fabs(cdf_left-cdf_model) > *max_diff)
 *max_diff = fabs(cdf_left-cdf_model);
      cdf_left = cdf_right;
    }
    effective_n=n1;
  }
  if(effective_n < 4) {
    printwarning(verbose.debug, "WARNING kstest: Number of data-points is too low to make use of approximations used in this implementation of the KS-test.");
  }
  effective_n=sqrt(effective_n);
  ks_statistic = (*max_diff)*(effective_n+0.12+0.11/effective_n);
  coeff = -2.0*ks_statistic*ks_statistic;
  *prob = 0;
  sign = 1;
  last_term = 0;
  converged = 0;
  for(i1 = 1; i1 <= 100; i1++) {
    cur_term = sign*2.0*exp(coeff*i1*i1);
    *prob += cur_term;
    if(fabs(cur_term) <= 1e-5*fabs(last_term) || fabs(cur_term) <= 1e-10*(*prob)) {
      converged = 1;
      break;
    }
    last_term = cur_term;
    sign = -sign;
  }
  if(!converged)
    *prob = 1;
  if(verbose.verbose) {
    printf("KS-test statistic max_diff: %lf = %e\n", *max_diff, *max_diff);
    printf("KS-test probability:        %lf = %e\nA small probability means the two sets of points are drawn from a different distribution.\n", *prob, *prob);
    printf("The null hypothesis ");
    if(n2 > 0 && data2 != NULL) {
      printf("(data sets are drawn from the same distribution) ");
    }else {
      printf("(data set is drawn from the specified distribution) ");
    }
    printf("can be rejected at the %.2lf sigma level.\n", gsl_cdf_gaussian_Pinv(0.5*(1.0-*prob)+0.5, 1.0));
  }
}

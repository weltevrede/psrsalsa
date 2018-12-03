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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_version.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include "psrsalsa.h"
typedef struct {
  size_t n;
  double *x;
  double *y;
  double *sigma;
  fitfunc_collection_type *fitfunction;
  size_t iter;
}levmar_internal_fitter_data_def;
static levmar_internal_fitter_data_def levmar_internal_fitter_data;
typedef struct {
  gsl_vector_view func_params;
  const gsl_multifit_fdfsolver_type *solver_type;
  gsl_multifit_function_fdf user_itteration_funcs;
  gsl_multifit_fdfsolver *solver;
  gsl_matrix *covar;
}levmar_internal_gsl_info;
typedef struct {
  double *covar;
  double *alpha;
  double *beta;
  double *alpha_tmp;
  double *beta_tmp;
  double *func_params_laststep;
  double *derivatives;
  double alambda;
  double chisq;
  long nrfitparams;
  levmar_internal_fitter_data_def *data;
}levmar_internal_psrsalsa_info;
void print_gsl_version_used(FILE *stream)
{
  int major, minor;
  major = GSL_VERSION_NUMBER/100.0;
  minor = GSL_VERSION_NUMBER-100*major;
  fprintf(stream, "%s (header) %d.%d (specified during compilation)", GSL_VERSION, major, minor);
}
double evaluate_fitfunc_collection(fitfunc_collection_type *function, double x, verbose_definition verbose)
{
  int i;
  double y;
  y = 0;
  for(i = 0; i < function->nrfuncs; i++) {
    if(function->func[i].type == FUNC_POLYNOMAL) {
      y += function->func[i].value[0]*pow(x, function->func[i].param[0]);
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR evaluate_fitfunction_collection: Unknown funtional type in specified function.");
      exit(0);
    }
  }
  return y;
}
double evaluate_fitfunc_collection_deriv_param(fitfunc_collection_type *function, int functionnr, int paramnr, double x, verbose_definition verbose)
{
  if(function->func[functionnr].type == FUNC_POLYNOMAL) {
    if(paramnr == 0) {
      return pow(x, function->func[functionnr].param[0]);
    }
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR evaluate_fitfunc_collection_deriv_param: Unknown funtional type in specified function.");
    exit(0);
  }
  fflush(stdout);
  printerror(verbose.debug, "ERROR evaluate_fitfunc_collection_deriv_param: Unknown paramer number specified for given funtional type.");
  exit(0);
}
int print_fitfunctions(fitfunc_collection_type *function, int novalue, int showerror, int index, verbose_definition verbose)
{
  int i, j, n;
  if(showerror)
    showerror = 1;
  for(n = 0; n <= showerror; n++) {
    if(n == 0)
      printf("y = ");
    j = 1;
    for(i = 0; i < function->nrfuncs; i++) {
      if(index == i || index < 0) {
 if(function->func[i].type == FUNC_POLYNOMAL) {
   if(novalue || showerror) {
     if(showerror && n == 1) {
       printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j, function->func[i].value[0], function->func[i].error[0], function->func[i].value[0], function->func[i].error[0]);
     }else {
       printf("a%d*x**%lf", j, function->func[i].param[0]);
     }
   }else {
     printf("%lf*x**%lf", function->func[i].value[0], function->func[i].param[0]);
   }
   j += 1;
 }else {
   fflush(stdout);
   printerror(verbose.debug, "ERROR print_fitfunctions: Unrecognized function type.");
   return 0;
 }
 if(i == function->nrfuncs - 1)
   printf("\n");
 else if(n == 0)
   printf(" + ");
      }
    }
  }
  return 1;
}
void countnrparameters_fitfunction(fitfunc_collection_type *function, int *nrfitparameters)
{
  int n, fnr, nrvalues;
  *nrfitparameters = 0;
  for(fnr = 0; fnr < function->nrfuncs; fnr++) {
    if(function->func[fnr].type == FUNC_POLYNOMAL) {
      nrvalues = 1;
    }else {
      fflush(stdout);
      printerror(0, "ERROR countnrparameters_fitfunction: Unknown funtional type.");
      exit(0);
    }
    for(n = 0; n < nrvalues; n++) {
      if(function->func[fnr].fit_flag[n]) {
 *nrfitparameters += 1;
      }
    }
  }
}
int set_fitted_parameters_fitfunc_collection(fitfunc_collection_type *function, double *fitparameters, verbose_definition verbose)
{
  int fnr, nrvalues, valnr, fitparamnr;
  fitparamnr = 0;
  for(fnr = 0; fnr < function->nrfuncs; fnr++) {
    if(function->func[fnr].type == FUNC_POLYNOMAL) {
      nrvalues = 1;
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR set_fitted_parameters_fitfunc_collection: Unknown funtional type.");
      return 0;
    }
    for(valnr = 0; valnr < nrvalues; valnr++) {
      if(function->func[fnr].fit_flag[valnr]) {
 function->func[fnr].value[valnr] = fitparameters[fitparamnr++];
      }
    }
  }
  return 1;
}
int levmar_itteration_calc_function_internal_gsl(const gsl_vector *fitparams, void *data, gsl_vector *f)
{
  int n;
  size_t npts = ((levmar_internal_fitter_data_def *)data)->n;
  double *xdata = ((levmar_internal_fitter_data_def *)data)->x;
  double *ydata = ((levmar_internal_fitter_data_def *)data)->y;
  double *sigma = ((levmar_internal_fitter_data_def *)data)->sigma;
  fitfunc_collection_type *fitfunction_internal = ((levmar_internal_fitter_data_def *)data)->fitfunction;
  double value;
  verbose_definition noverbose;
  cleanVerboseState(&noverbose);
  if(fitparams->stride != 1) {
    printerror(0, "ERROR levmar_itteration_calc_function_internal_gsl: It is assumed the gsl vectors have a stride of 1.");
    exit(0);
  }
  fitfunc_collection_type newfunction;
  memcpy(&newfunction, fitfunction_internal, sizeof(fitfunc_collection_type));
  if(set_fitted_parameters_fitfunc_collection(&newfunction, fitparams->data, noverbose) == 0) {
    printerror(0, "ERROR levmar_itteration_calc_function_internal_gsl: Copying new fit parameters failed.");
    exit(0);
  }
  for(n = 0; n < npts; n++) {
    value = evaluate_fitfunc_collection(&newfunction, xdata[n], noverbose);
    gsl_vector_set(f, n, (value - ydata[n])/sigma[n]);
  }
  return GSL_SUCCESS;
}
int levmar_itteration_calc_deriv_internal_gsl(const gsl_vector *fitparams, void *data, gsl_matrix *J)
{
  int fnr, j, n, nrvalues, paramnr;
  size_t npts = ((levmar_internal_fitter_data_def *)data)->n;
  double *xdata = ((levmar_internal_fitter_data_def *)data)->x;
  double *sigma = ((levmar_internal_fitter_data_def *)data)->sigma;
  fitfunc_collection_type *fitfunction_internal = ((levmar_internal_fitter_data_def *)data)->fitfunction;
  double deriv;
  verbose_definition noverbose;
  cleanVerboseState(&noverbose);
  if(fitparams->stride != 1) {
    printerror(0, "ERROR levmar_itteration_calc_function_internal_gsl: It is assumed the gsl vectors have a stride of 1.");
    exit(0);
  }
  fitfunc_collection_type newfunction;
  memcpy(&newfunction, fitfunction_internal, sizeof(fitfunc_collection_type));
  if(set_fitted_parameters_fitfunc_collection(&newfunction, fitparams->data, noverbose) == 0) {
    printerror(0, "ERROR levmar_itteration_calc_function_internal_gsl: Copying new fit parameters failed.");
    exit(0);
  }
  j = 0;
  for(fnr = 0; fnr < newfunction.nrfuncs; fnr++) {
    if(newfunction.func[fnr].type == FUNC_POLYNOMAL) {
      nrvalues = 1;
    }else {
      fflush(stdout);
      printerror(noverbose.debug, "ERROR levmar_itteration_calc_deriv_internal_gsl: Unknown funtional type.");
      return 0;
    }
    for(paramnr = 0; paramnr < nrvalues; paramnr++) {
      if(newfunction.func[fnr].fit_flag[paramnr]) {
 for(n = 0; n < npts; n++) {
   deriv = evaluate_fitfunc_collection_deriv_param(&newfunction, fnr, paramnr, xdata[n], noverbose);
   gsl_matrix_set(J, n, j, deriv/sigma[n]);
 }
 j++;
      }
    }
  }
  return GSL_SUCCESS;
}
int levmar_itteration_calc_func_and_deriv_internal_gsl(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
{
  levmar_itteration_calc_function_internal_gsl(x, data, f);
  levmar_itteration_calc_deriv_internal_gsl(x, data, J);
  return GSL_SUCCESS;
}
void levmar_itteration_calc_alpha_beta_chi2_internal_psrsalsa(double *params, void *info)
{
  long curDataPoint, nrDataPoints, curparam, curparam2, nrfitparams, fnr, nrvalues, paramnr;
  double *alpha, *beta, *derivatives, *x, *y, *sigma, chi2, ypred, tmp1, tmp2, tmp3, tmp4;
  verbose_definition noverbose;
  cleanVerboseState(&noverbose);
  nrDataPoints = ((levmar_internal_psrsalsa_info *)info)->data->n;
  x = ((levmar_internal_psrsalsa_info *)info)->data->x;
  y = ((levmar_internal_psrsalsa_info *)info)->data->y;
  sigma = ((levmar_internal_psrsalsa_info *)info)->data->sigma;
  alpha = ((levmar_internal_psrsalsa_info *)info)->alpha;
  beta = ((levmar_internal_psrsalsa_info *)info)->beta;
  derivatives = ((levmar_internal_psrsalsa_info *)info)->derivatives;
  nrfitparams = ((levmar_internal_psrsalsa_info *)info)->nrfitparams;
  fitfunc_collection_type newfunction;
  memcpy(&newfunction, levmar_internal_fitter_data.fitfunction, sizeof(fitfunc_collection_type));
  if(set_fitted_parameters_fitfunc_collection(&newfunction, params, noverbose) == 0) {
    printerror(0, "ERROR levmar_itteration_calc_func_and_deriv_internal_nr: Copying new fit parameters failed.");
    exit(0);
  }
  for(curparam = 0; curparam < nrfitparams; curparam++) {
    beta[curparam] = 0;
    for(curparam2 = curparam; curparam2 < nrfitparams; curparam2++) {
      alpha[curparam*nrfitparams+curparam2] = 0;
    }
  }
  chi2 = 0;
  for(curDataPoint = 0; curDataPoint < nrDataPoints; curDataPoint++) {
    ypred = evaluate_fitfunc_collection(&newfunction, x[curDataPoint], noverbose);
    curparam = 0;
    for(fnr = 0; fnr < newfunction.nrfuncs; fnr++) {
      if(newfunction.func[fnr].type == FUNC_POLYNOMAL) {
 nrvalues = 1;
      }else {
 fflush(stdout);
 printerror(noverbose.debug, "ERROR levmar_itteration_calc_alpha_beta_chi2_internal_psrsalsa: Unknown funtional type.");
 exit(0);
      }
      for(paramnr = 0; paramnr < nrvalues; paramnr++) {
 if(newfunction.func[fnr].fit_flag[paramnr]) {
   derivatives[curparam++] = evaluate_fitfunc_collection_deriv_param(&newfunction, fnr, paramnr, x[curDataPoint], noverbose);
 }
      }
    }
    tmp1 = 1.0/(sigma[curDataPoint]*sigma[curDataPoint]);
    tmp2 = (y[curDataPoint] - ypred);
    tmp3 = tmp1*tmp2;
    chi2 += tmp2*tmp3;
    for(curparam = 0; curparam < nrfitparams; curparam++) {
      tmp4 = derivatives[curparam]*tmp3;
      beta[curparam] += tmp4;
      for(curparam2 = curparam; curparam2 < nrfitparams; curparam2++) {
 alpha[curparam*nrfitparams+curparam2] += derivatives[curparam]*derivatives[curparam2]*tmp1;
      }
    }
  }
  for(curparam = 0; curparam < nrfitparams; curparam++) {
    for(curparam2 = curparam+1; curparam2 < nrfitparams; curparam2++) {
      alpha[curparam2*nrfitparams+curparam] = alpha[curparam*nrfitparams+curparam2];
    }
  }
  ((levmar_internal_psrsalsa_info *)info)->chisq = chi2;
}
int initialize_levmar_check_start_and_current_are_the_same(fitfunc_collection_type *function, verbose_definition verbose)
{
  int fnr, nrvalues, paramnr;
  for(fnr = 0; fnr < function->nrfuncs; fnr++) {
    if(function->func[fnr].type == FUNC_POLYNOMAL) {
      nrvalues = 1;
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR levmar_itteration_calc_func_and_deriv_internal_nr: Unknown funtional type.");
      return 0;
    }
    for(paramnr = 0; paramnr < nrvalues; paramnr++) {
      if(function->func[fnr].fit_flag[paramnr]) {
 if(function->func[fnr].value[paramnr] != function->func[fnr].start[paramnr]) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR initialize_levmar_check_start_and_current_are_the_same: Start and current value appear to be different (%lf != %lf)", function->func[fnr].value[paramnr], function->func[fnr].start[paramnr]);
   return 0;
 }
      }
    }
  }
  return 1;
}
int levmar_initialize(fitfunc_collection_type *function, double **func_params, double **func_params_err, int *nrfitparameters, double *datax, double *datay, double *datasigma, long nrdatapoints, levmar_internal_gsl_info *gsl_params, levmar_internal_psrsalsa_info *psrsalsa_params,
#ifdef NRAVAIL
levmar_internal_nr_info *nr_params,
#endif
int algorithm, int mode, int showcovariance, verbose_definition verbose)
{
  int fnr, nrvalues, paramnr, n, n2;
  if(mode == 0) {
    countnrparameters_fitfunction(function, nrfitparameters);
    if(verbose.verbose) printf("  Initializing fitter with %d fit parameters\n", *nrfitparameters);
  }else if(mode == 1) {
    levmar_internal_fitter_data.n = nrdatapoints;
    levmar_internal_fitter_data.x = datax;
    levmar_internal_fitter_data.y = datay;
    levmar_internal_fitter_data.sigma = datasigma;
    levmar_internal_fitter_data.fitfunction = function;
    *func_params = malloc(*nrfitparameters*sizeof(double));
    *func_params_err = malloc(*nrfitparameters*sizeof(double));
    if(*func_params == NULL || *func_params_err == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR levmar_initialize: Cannot allocate memory");
      return 0;
    }
    paramnr = 0;
    for(fnr = 0; fnr < function->nrfuncs; fnr++) {
      if(function->func[fnr].type == FUNC_POLYNOMAL) {
 nrvalues = 1;
      }else {
 fflush(stdout);
 printerror(verbose.debug, "ERROR levmar_initialize: Unknown funtional type.");
 return 0;
      }
      for(n = 0; n < nrvalues; n++) {
 if(function->func[fnr].fit_flag[n]) {
   (*func_params)[paramnr++] = function->func[fnr].value[n];
 }
 function->func[fnr].error[n] = sqrt(-1);
      }
    }
    if(algorithm == 1) {
      gsl_params->func_params = gsl_vector_view_array(*func_params, *nrfitparameters);
      gsl_params->covar = gsl_matrix_alloc(*nrfitparameters, *nrfitparameters);
      gsl_params->user_itteration_funcs.f = &levmar_itteration_calc_function_internal_gsl;
      gsl_params->user_itteration_funcs.df = &levmar_itteration_calc_deriv_internal_gsl;
      gsl_params->user_itteration_funcs.fdf = &levmar_itteration_calc_func_and_deriv_internal_gsl;
      gsl_params->user_itteration_funcs.n = nrdatapoints;
      gsl_params->user_itteration_funcs.p = *nrfitparameters;
      gsl_params->user_itteration_funcs.params = &levmar_internal_fitter_data;
      gsl_params->solver_type = gsl_multifit_fdfsolver_lmsder;
      gsl_params->solver = gsl_multifit_fdfsolver_alloc (gsl_params->solver_type, nrdatapoints, *nrfitparameters);
      gsl_multifit_fdfsolver_set(gsl_params->solver, &(gsl_params->user_itteration_funcs), &(gsl_params->func_params.vector));
    }else if(algorithm == 2) {
      psrsalsa_params->func_params_laststep = malloc(sizeof(double)*(*nrfitparameters));
      psrsalsa_params->beta = malloc(sizeof(double)*(*nrfitparameters));
      psrsalsa_params->derivatives = malloc(sizeof(double)*(*nrfitparameters));
      psrsalsa_params->covar = malloc(sizeof(double)*(*nrfitparameters)*(*nrfitparameters));
      psrsalsa_params->alpha = malloc(sizeof(double)*(*nrfitparameters)*(*nrfitparameters));
      psrsalsa_params->alpha_tmp = malloc(sizeof(double)*(*nrfitparameters)*(*nrfitparameters));
      psrsalsa_params->beta_tmp = malloc(sizeof(double)*(*nrfitparameters));
      if(psrsalsa_params->func_params_laststep == NULL || psrsalsa_params->beta == NULL || psrsalsa_params->derivatives == NULL || psrsalsa_params->alpha == NULL || psrsalsa_params->covar == NULL || psrsalsa_params->beta_tmp == NULL || psrsalsa_params->alpha_tmp == NULL) {
 printerror(verbose.debug, "ERROR levmar_initialize: Memory allocation error");
 return 0;
      }
      psrsalsa_params->nrfitparams = *nrfitparameters;
      psrsalsa_params->data = &levmar_internal_fitter_data;
      psrsalsa_params->alambda = 0.001;
      levmar_itteration_calc_alpha_beta_chi2_internal_psrsalsa(*func_params, psrsalsa_params);
    }else if(algorithm == 3) {
#ifdef NRAVAIL
      nr_params->ia = malloc(sizeof(int)*(*nrfitparameters));
      nr_params->func_params_laststep = malloc(sizeof(double)*(*nrfitparameters));
      if(nr_params->ia == NULL || nr_params->func_params_laststep == NULL) {
 printerror(verbose.debug, "ERROR levmar_initialize: Memory allocation error");
 return 0;
      }
      for(n = 0; n < *nrfitparameters; n++)
 nr_params->ia[n] = 1;
      double **matrix_nr_d(long nrl, long nrh, long ncl, long nch);
      nr_params->covar = matrix_nr_d(1, *nrfitparameters, 1, *nrfitparameters);
      nr_params->alpha = matrix_nr_d(1, *nrfitparameters, 1, *nrfitparameters);
      void mrqmin_nr_d(double x[], double y[], double sig[], int ndata, double a[], int ia[],
     int ma, double **covar, double **alpha, double *chisq,
     void (*funcs)(double, double [], double *, double [], int), double *alamda);
      nr_params->alambda = -1;
      mrqmin_nr_d(levmar_internal_fitter_data.x-1, levmar_internal_fitter_data.y-1, levmar_internal_fitter_data.sigma-1, nrdatapoints, (*func_params)-1, nr_params->ia-1, *nrfitparameters, nr_params->covar, nr_params->alpha, &nr_params->chisq, levmar_itteration_calc_func_and_deriv_internal_nr, &(nr_params->alambda));
#else
      printerror(verbose.debug, "ERROR levmar_initialize: Code is not compiled with NR support");
      return 0;
#endif
    }
  }else if(mode == 2) {
    paramnr = 0;
    for(fnr = 0; fnr < function->nrfuncs; fnr++) {
      if(function->func[fnr].type == FUNC_POLYNOMAL) {
 nrvalues = 1;
      }else {
 fflush(stdout);
 printerror(0, "ERROR levmar_initialize: Unknown funtional type.");
 return 0;
      }
      for(n = 0; n < nrvalues; n++) {
 if(function->func[fnr].fit_flag[n]) {
   function->func[fnr].value[n] = (*func_params)[paramnr];
   function->func[fnr].error[n] = (*func_params_err)[paramnr];
   paramnr++;
 }
      }
    }
    if(algorithm == 1) {
      if(showcovariance) {
 printf("Covariance matrix: \n");
 for(n = 0; n < *nrfitparameters; n++) {
   for(n2 = 0; n2 < *nrfitparameters; n2++) {
     printf("%lf ", gsl_matrix_get(gsl_params->covar,n,n2));
   }
   printf("\n");
 }
      }
      gsl_multifit_fdfsolver_free(gsl_params->solver);
      gsl_matrix_free(gsl_params->covar);
    }else if(algorithm == 2) {
      if(showcovariance) {
 printf("Covariance matrix: \n");
 for(n = 0; n < *nrfitparameters; n++) {
   for(n2 = 0; n2 < *nrfitparameters; n2++) {
     printf("%lf ", psrsalsa_params->covar[n*(*nrfitparameters)+n2]);
   }
   printf("\n");
 }
      }
      free(psrsalsa_params->beta);
      free(psrsalsa_params->derivatives);
      free(psrsalsa_params->alpha);
      free(psrsalsa_params->covar);
      free(psrsalsa_params->func_params_laststep);
      free(psrsalsa_params->alpha_tmp);
      free(psrsalsa_params->beta_tmp);
    }else if(algorithm == 3) {
#ifdef NRAVAIL
      if(showcovariance) {
 printf("Covariance matrix: \n");
 for(n = 0; n < *nrfitparameters; n++) {
   for(n2 = 0; n2 < *nrfitparameters; n2++) {
     printf("%lf ", nr_params->covar[n+1][n2+1]);
   }
   printf("\n");
 }
      }
      void free_matrix_nr_d(double **m, long nrl, long nrh, long ncl, long nch);
      free_matrix_nr_d(nr_params->covar, 1, *nrfitparameters, 1, *nrfitparameters);
      free_matrix_nr_d(nr_params->alpha, 1, *nrfitparameters, 1, *nrfitparameters);
      free(nr_params->ia);
      free(nr_params->func_params_laststep);
#else
      printerror(verbose.debug, "ERROR levmar_initialize: Code is not compiled with NR support");
      return 0;
#endif
    }
    free(*func_params);
    free(*func_params_err);
  }
  return 1;
}
int fit_levmar_internal(int algorithm, int nrdatapoints, int nrfitparameters, double *func_params, double *func_params_err, levmar_internal_gsl_info *gsl_params, levmar_internal_psrsalsa_info *psrsalsa_params,
#ifdef NRAVAIL
   levmar_internal_nr_info *nr_params,
#endif
   double epsabs, double epsrel, int maxiter, verbose_definition verbose)
{
  int status, i, converged;
  long iter;
  double chisq_last_nr;
#ifdef NRAVAIL
  void mrqmin_nr_d(double x[], double y[], double sig[], int ndata, double a[], int ia[],
     int ma, double **covar, double **alpha, double *chisq,
     void (*funcs)(double, double [], double *, double [], int), double *alamda);
#else
  if(algorithm == 3) {
    printerror(verbose.debug, "ERROR fit_levmar_internal: Code is not compiled with NR support");
    return 0;
  }
#endif
  if(algorithm == 2) {
    memcpy(psrsalsa_params->func_params_laststep, func_params, sizeof(double)*nrfitparameters);
    chisq_last_nr = psrsalsa_params->chisq;
#ifdef NRAVAIL
  }else if(algorithm == 3) {
    memcpy(nr_params->func_params_laststep, func_params, sizeof(double)*nrfitparameters);
    chisq_last_nr = nr_params->chisq;
#endif
  }
  iter = 0;
  do {
    iter++;
    levmar_internal_fitter_data.iter += 1;
    if(algorithm == 1) {
      status = gsl_multifit_fdfsolver_iterate(gsl_params->solver);
      if(status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG) {
 int status2;
 status2 = gsl_multifit_test_delta(gsl_params->solver->dx, gsl_params->solver->x, epsabs, epsrel);
 if(status2 == GSL_SUCCESS) {
   status = 0;
 }else {
   if(verbose.debug == 0) {
     printwarning(verbose.debug, "Machine precision reached, but required precision tolerance conditions are not met.");
   }else {
     printwarning(verbose.debug, "Machine precision reached, but required precision tolerance conditions are not met:");
     for(i = 0; i < nrfitparameters; i++) {
       printwarning(verbose.debug, "parameter %d: value=%e stepsize=%e", i+1, gsl_vector_get(gsl_params->solver->x, i), gsl_vector_get(gsl_params->solver->dx, i));
     }
   }
 }
      }
      if(verbose.debug) {
 printf("fit_levmar_internal: itteration %ld - status = %s\n", iter, gsl_strerror(status));
      }
      if(status) {
 break;
      }else {
 status = gsl_multifit_test_delta(gsl_params->solver->dx, gsl_params->solver->x, epsabs, epsrel);
      }
    }else if(algorithm == 2) {
      for(i = 0; i < nrfitparameters*nrfitparameters; i++)
 psrsalsa_params->alpha_tmp[i] = psrsalsa_params->alpha[i];
      for(i = 0; i < nrfitparameters; i++)
 psrsalsa_params->beta_tmp[i] = psrsalsa_params->beta[i];
      for (i = 0; i < nrfitparameters; i++) {
 psrsalsa_params->alpha[i*nrfitparameters+i] *= 1.0+psrsalsa_params->alambda;
      }
      int linalg_solve_matrix_eq_gauss_jordan(double *matrixa, double *matrixb, int n, int m, verbose_definition verbose);
      if(linalg_solve_matrix_eq_gauss_jordan(psrsalsa_params->alpha, psrsalsa_params->beta, nrfitparameters, 1, verbose) != 0) {
 printerror(verbose.debug, "ERROR fit_levmar_internal: Solving matrix equation failed.");
 return 9;
      }
      for (i = 0; i < nrfitparameters; i++)
 func_params[i] += psrsalsa_params->beta[i];
      levmar_itteration_calc_alpha_beta_chi2_internal_psrsalsa(func_params, psrsalsa_params);
      status = GSL_SUCCESS;
      if(chisq_last_nr < psrsalsa_params->chisq) {
 status = GSL_ETOLF;
 memcpy(func_params, psrsalsa_params->func_params_laststep, sizeof(double)*nrfitparameters);
 for(i = 0; i < nrfitparameters*nrfitparameters; i++)
   psrsalsa_params->alpha[i] = psrsalsa_params->alpha_tmp[i];
 for(i = 0; i < nrfitparameters; i++)
   psrsalsa_params->beta[i] = psrsalsa_params->beta_tmp[i];
 psrsalsa_params->alambda *= 10.1;
      }
      converged = 1;
      for(i = 0; i < nrfitparameters; i++) {
 if(fabs(func_params[i] - psrsalsa_params->func_params_laststep[i]) >= 1.2e-16+epsabs+epsrel*fabs(psrsalsa_params->func_params_laststep[i])) {
   converged = 0;
   break;
 }
      }
      if(converged && status == GSL_ETOLF)
 status = GSL_SUCCESS;
      if(status == GSL_ETOLF) {
       printwarning(verbose.debug, "Machine precision reached, but required precision tolerance conditions are not met.");
      }else if(converged == 0) {
 status = GSL_CONTINUE;
      }
      if(chisq_last_nr > psrsalsa_params->chisq) {
 memcpy(psrsalsa_params->func_params_laststep, func_params, sizeof(double)*nrfitparameters);
 chisq_last_nr = psrsalsa_params->chisq;
 psrsalsa_params->alambda *= 0.1;
      }
      if(verbose.debug) {
 printf("fit_levmar_internal: itteration %ld - status = %s\n", iter, gsl_strerror(status));
      }
#ifdef NRAVAIL
    }else if(algorithm == 3) {
      mrqmin_nr_d(levmar_internal_fitter_data.x-1, levmar_internal_fitter_data.y-1, levmar_internal_fitter_data.sigma-1, nrdatapoints, func_params-1, nr_params->ia-1, nrfitparameters, nr_params->covar, nr_params->alpha, &(nr_params->chisq), levmar_itteration_calc_func_and_deriv_internal_nr, &(nr_params->alambda));
      status = GSL_SUCCESS;
      for(i = 0; i < nrfitparameters; i++) {
 if(fabs(func_params[i] - nr_params->func_params_laststep[i]) >= epsabs+epsrel*fabs(nr_params->func_params_laststep[i])) {
   status = GSL_CONTINUE;
   break;
 }
      }
      if(verbose.debug) {
 printf("fit_levmar_internal: itteration %ld - status = %s\n", iter, gsl_strerror(status));
      }
      if(status == GSL_CONTINUE) {
 if(chisq_last_nr <= nr_params->chisq) {
   status = GSL_ETOLF;
   printwarning(verbose.debug, "Machine precision reached, but required precision tolerance conditions are not met.");
 }else {
   memcpy(nr_params->func_params_laststep, func_params, sizeof(double)*nrfitparameters);
   chisq_last_nr = nr_params->chisq;
 }
      }
#endif
    }
  }while(status == GSL_CONTINUE && iter < maxiter);
  if(algorithm == 1) {
#if GSL_VERSION_NUMBER >= 200
    if(verbose.debug) {
      printf("fit_levmar_internal: Entering gsl >= 2.0 specific part of the code\n");
    }
    gsl_matrix *J = gsl_matrix_alloc(nrdatapoints, nrfitparameters);
    gsl_multifit_fdfsolver_jac(gsl_params->solver, J);
    gsl_multifit_covar(J, 0.0, gsl_params->covar);
#else
    if(verbose.debug) {
      printf("fit_levmar_internal: Entering gsl < 2.0 specific part of the code\n");
    }
    gsl_multifit_covar(gsl_params->solver->J, 0.0, gsl_params->covar);
#endif
    levmar_internal_fitter_data.fitfunction->chi2 = gsl_blas_dnrm2(gsl_params->solver->f);
    levmar_internal_fitter_data.fitfunction->chi2_red = (levmar_internal_fitter_data.fitfunction->chi2)/sqrt(nrdatapoints - nrfitparameters);
    levmar_internal_fitter_data.fitfunction->chi2_red *= levmar_internal_fitter_data.fitfunction->chi2_red;
    levmar_internal_fitter_data.fitfunction->chi2 *= levmar_internal_fitter_data.fitfunction->chi2;
  }else if(algorithm == 2) {
    int linalg_solve_matrix_eq_gauss_jordan(double *matrixa, double *matrixb, int n, int m, verbose_definition verbose);
    if(linalg_solve_matrix_eq_gauss_jordan(psrsalsa_params->alpha, psrsalsa_params->beta, nrfitparameters, 1, verbose) != 0) {
      printerror(verbose.debug, "ERROR fit_levmar_internal: Solving matrix equation failed.");
      return 9;
    }
    for(i = 0; i < nrfitparameters*nrfitparameters; i++)
      psrsalsa_params->covar[i] = psrsalsa_params->alpha[i];
    levmar_internal_fitter_data.fitfunction->chi2 = psrsalsa_params->chisq;
    levmar_internal_fitter_data.fitfunction->chi2_red = psrsalsa_params->chisq/(double)(nrdatapoints - nrfitparameters);
#ifdef NRAVAIL
  }else if(algorithm == 3) {
    nr_params->alambda = 0;
    mrqmin_nr_d(levmar_internal_fitter_data.x-1, levmar_internal_fitter_data.y-1, levmar_internal_fitter_data.sigma-1, nrdatapoints, func_params-1, nr_params->ia-1, nrfitparameters, nr_params->covar, nr_params->alpha, &nr_params->chisq, levmar_itteration_calc_func_and_deriv_internal_nr, &nr_params->alambda);
    levmar_internal_fitter_data.fitfunction->chi2 = nr_params->chisq;
    levmar_internal_fitter_data.fitfunction->chi2_red = nr_params->chisq/(double)(nrdatapoints - nrfitparameters);
#endif
  }
  if(verbose.debug) {
    printf("chisq/dof = %g after %ld iterations\n", levmar_internal_fitter_data.fitfunction->chi2_red, iter);
  }
  if(algorithm == 1) {
    for(i = 0; i < nrfitparameters; i++) {
      func_params[i] = gsl_vector_get(gsl_params->solver->x, i);
      func_params_err[i] = sqrt(gsl_matrix_get(gsl_params->covar,i,i));
    }
  }else if(algorithm == 2) {
    for(i = 0; i < nrfitparameters; i++) {
      func_params_err[i] = sqrt(psrsalsa_params->covar[i*nrfitparameters+i]);
    }
#ifdef NRAVAIL
  }else if(algorithm == 3) {
    for(i = 0; i < nrfitparameters; i++) {
      func_params_err[i] = sqrt(nr_params->covar[i+1][i+1]);
    }
#endif
  }
  if(levmar_internal_fitter_data.iter == maxiter)
    return 2;
  if(algorithm == 1) {
    if(status == 0)
      return 0;
    if(status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG) {
      return 3;
    }
    if(status == GSL_ENOPROG) {
      fflush(stdout);
      printerror(verbose.debug, "fit_levmar_internal: Cannot determine suitable trial step. Continue itterating might help");
      return 4;
    }
  }else if(algorithm == 2 || algorithm == 3) {
    return 0;
  }
  fflush(stdout);
  printerror(verbose.debug, "fit_levmar_internal: UNKNOWN RETURN CODE");
  return 10;
}
int fit_levmar(int algorithm, fitfunc_collection_type *function, double *data_x, double *data_y, double *data_sigma, long ndata, int oneatatime, int force_chi2_1, double epsabs, double epsrel, int maxiter, int *status, int showresults, int showcovariance, verbose_definition verbose)
{
  int loopnr, improved, i, j, k, fnr, valuenr, nrfitparameters_one, nrfitparameters, free_data_sigma, nrvalues;
  double old_chi2, sig_scale, *func_params, *func_params_err;
  verbose_definition noverbose;
  levmar_internal_gsl_info gsl_params;
  levmar_internal_psrsalsa_info psrsalsa_params;
#ifndef NRAVAIL
  if(algorithm == 3) {
    printerror(verbose.debug, "ERROR fit_levmar: Code is not compiled with NR support");
    return 0;
  }
#else
  levmar_internal_nr_info nr_params;
#endif
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;
  if(verbose.verbose) {
    printf("Levenberg-Marquardt algorithm:\n");
    if(verbose.debug) {
      printf("  Fit function: ");
      print_fitfunctions(function, 1, 0, -1, verbose);
    }
  }
  if(initialize_levmar_check_start_and_current_are_the_same(function, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fit_levmar: Start and current value appear to be different");
    return 1;
  }
  free_data_sigma = 0;
  if(data_sigma == NULL) {
    data_sigma = malloc(ndata*sizeof(double));
    for(i = 0; i < ndata; i++)
      data_sigma[i] = 1;
    free_data_sigma = 1;
    force_chi2_1 = 1;
  }
  if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
         &nr_params,
#endif
         algorithm, 0, showcovariance, verbose) == 0) {
    return 1;
  }
  levmar_internal_fitter_data.iter = 0;
  if(oneatatime == 0) {
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
    &nr_params,
#endif
    algorithm, 1, showcovariance, verbose) == 0) {
      return 1;
    }
    *status = fit_levmar_internal(algorithm, ndata, nrfitparameters, func_params, func_params_err, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
      &nr_params,
#endif
      epsabs, epsrel, maxiter, verbose);
  }else {
    fitfunc_collection_type *original_function;
    original_function = malloc(sizeof(fitfunc_collection_type));
    if(original_function == NULL) {
      printerror(verbose.debug, "ERROR fit_levmar: Memory allocation error.");
      return 1;
    }
    memcpy(original_function, function, sizeof(fitfunc_collection_type));
    old_chi2 = -1;
    loopnr = 0;
    do {
      improved = 0;
      i = 0;
      for(fnr = 0; fnr < function->nrfuncs; fnr++) {
 if(function->func[fnr].type == FUNC_POLYNOMAL) {
   nrvalues = 1;
 }else {
   fflush(stdout);
   printerror(verbose.debug, "ERROR fit_levmar: Unknown funtional type.");
   return 1;
 }
 for(valuenr = 0; valuenr < nrvalues; valuenr++) {
   if(original_function->func[fnr].fit_flag[valuenr]) {
     for(j = 0; j < function->nrfuncs; j++) {
       for(k = 0; k < MaxNrFitParameters; k++)
  function->func[j].fit_flag[k] = 0;
     }
     function->func[fnr].fit_flag[valuenr] = 1;
     if(verbose.debug) {
       printf("  Loop %d: chi2 = %f", loopnr, old_chi2);
       printf("\n");
     }
     if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters_one, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
     &nr_params,
#endif
     algorithm, 0, showcovariance, noverbose) == 0) {
       free(original_function);
       return 1;
     }
     if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters_one, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
     &nr_params,
#endif
     algorithm, 1, showcovariance, verbose) == 0) {
       free(original_function);
       return 1;
     }
     *status = fit_levmar_internal(algorithm, ndata, nrfitparameters_one, func_params, func_params_err, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
       &nr_params,
#endif
       epsabs, epsrel, maxiter, verbose);
     if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters_one, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
     &nr_params,
#endif
     algorithm, 2, showcovariance, verbose) == 0) {
       free(original_function);
       return 1;
     }
     if(function->chi2 < old_chi2 || old_chi2 < 0) {
       old_chi2 = function->chi2;
       improved = 1;
     }
     loopnr++;
     i++;
   }
 }
      }
    }while(improved);
    for(j = 0; j < function->nrfuncs; j++) {
      for(k = 0; k < MaxNrFitParameters; k++)
 function->func[j].fit_flag[k] = original_function->func[j].fit_flag[k];
    }
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
    &nr_params,
#endif
    algorithm, 0, showcovariance, verbose) == 0) {
      free(original_function);
      return 1;
    }
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
    &nr_params,
#endif
    algorithm, 1, showcovariance, verbose) == 0) {
      free(original_function);
      return 1;
    }
    *status = fit_levmar_internal(algorithm, ndata, nrfitparameters, func_params, func_params_err, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
      &nr_params,
#endif
      epsabs, epsrel, maxiter, verbose);
    free(original_function);
  }
  if(verbose.verbose) {
    printf("  After %ld iterations found chisq/dof = %g\n", levmar_internal_fitter_data.iter, levmar_internal_fitter_data.fitfunction->chi2_red);
  }
  if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
         &nr_params,
#endif
         algorithm, 2, showcovariance, verbose) == 0) {
    return 1;
  }
  if(force_chi2_1) {
    sig_scale = sqrt(function->chi2_red);
    for(i = 0; i < ndata; i++) {
      data_sigma[i] *= sig_scale;
    }
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
    &nr_params,
#endif
    algorithm, 0, showcovariance, verbose) == 0) {
      return 1;
    }
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
    &nr_params,
#endif
    algorithm, 1, showcovariance, verbose) == 0) {
      return 1;
    }
    *status = fit_levmar_internal(algorithm, ndata, nrfitparameters, func_params, func_params_err, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
      &nr_params,
#endif
      epsabs, epsrel, maxiter, verbose);
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params,
#ifdef NRAVAIL
    &nr_params,
#endif
    algorithm, 2, showcovariance, verbose) == 0) {
      return 1;
    }
  }
  if(verbose.verbose || showresults) {
    printf("  Fit function: ");
    print_fitfunctions(function, 0, 1, -1, verbose);
    printf("  ");
    print_fitfunctions(function, 0, 0, -1, verbose);
    printf("\n  chi2 = %f\n  reduced chi2 = %f (%ld points - %d params)\n", function->chi2, function->chi2_red, ndata, nrfitparameters);
  }
  if(free_data_sigma) {
    free(data_sigma);
    data_sigma = NULL;
  }
  return *status;
}

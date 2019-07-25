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
#include <string.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_randist.h"
#include "psrsalsa.h"
double get_FWHM_RMsynthesis(datafile_definition datain, verbose_definition verbose);
int main(int argc, char **argv)
{
  int i, j, index, nrrmsteps, status, parabola, pgplot_mapdevice, pgplot_rmdevice, pgplot_profiledevice;
  int collapse, write_ascii, bootstrap, bootstrap_itt, nrOnpulseBins, not_normal_result_warning;
  long binnr, freqchannelnr, polnr, idnum;
  float rmmin, rmmax, sample;
  float *rmsynth_array, *cmap, *rms_channels;
  float *singlespectrum;
  double *singlespectrum_double, *rmgrid, *sigmagrid, *rmestimate, *rmcurestimate;
  double *rm_av, *rm_square, rmsigma, expectedRMerror, ftol;
  char outputname[MaxFilenameLength];
  FILE *ofile;
  psrsalsaApplication application;
  datafile_definition datain, clone;
  pgplot_options_definition pgplot_options;
  fitfunc_collection_type function;
  datafile_definition profiledata;
  gsl_rng *rand_num_gen;
  const gsl_rng_type *rand_num_gen_type;
  initApplication(&application, "rmsynth", "[options] inputfile");
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_rebin = 1;
  application.switch_rot = 1;
  application.switch_onpulse = 1;
  application.switch_onpulsef = 1;
  application.switch_onpulse2 = 1;
  application.switch_onpulsef2 = 1;
  application.switch_nocounters = 1;
  application.switch_stokes = 1;
  application.oformat = FITS_format;
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_dedisperse = 1;
  application.switch_fixseed = 1;
  nrrmsteps = 0;
  parabola = 0;
  pgplot_mapdevice = 0;
  pgplot_rmdevice = 0;
  pgplot_profiledevice = 0;
  collapse = 1;
  write_ascii = 0;
  bootstrap = 0;
  not_normal_result_warning = 0;
  ftol = -1;
  if(argc < 2) {
    printApplicationHelp(&application);
    fprintf(stdout, "Other options:\n");
    fprintf(stdout, "  -device       pgplot device for the Faraday depth plot.\n");
    fprintf(stdout, "  -device2      pgplot device for the longitude resolved map.\n");
    fprintf(stdout, "  -device3      pgplot device for the profile (used with -bootstrap).\n");
    fprintf(stdout, "  -rm           \"min max steps\" Set the minimum and maximum rm to search over and the nr of rm steps to use.\n");
    fprintf(stdout, "  -parabola     Fit data with parabola rather than estimated reponds.\n");
    fprintf(stdout, "  -long         Get a RM for each bin independently.\n");
    fprintf(stdout, "  -ascii        An ascii file the found RM values is generated.\n");
    fprintf(stdout, "  -bootstrap    Bootstrap the errorbars on the RM this number of times.\n");
    fprintf(stdout, "  -ftol         Set the fractional tolerance of the fitting algorithm used.\n");
    printf("\nPlease use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about RM syntesis method as implemented here can be found in: Ilie et al. 2018, accepted for publication in MNRAS, astro-ph/1811.12831\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
 i = index;
      }else if(strcmp(argv[i], "-rm") == 0) {
 j = sscanf(argv[i+1], "%f %f %d", &rmmin, &rmmax, &nrrmsteps);
 if(j != 3) {
   printerror(application.verbose_state.debug, "Error parsing %s option", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-bootstrap") == 0) {
 j = sscanf(argv[i+1], "%d", &bootstrap);
 if(j != 1) {
   printerror(application.verbose_state.debug, "Error parsing %s option", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-ftol") == 0) {
 j = sscanf(argv[i+1], "%lf", &ftol);
 if(j != 1) {
   printerror(application.verbose_state.debug, "Error parsing %s option", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-parabola") == 0) {
 parabola = 1;
      }else if(strcmp(argv[i], "-long") == 0) {
 collapse = 0;
      }else if(strcmp(argv[i], "-ascii") == 0) {
 write_ascii = 1;
      }else if(strcmp(argv[i], "-device") == 0) {
 pgplot_rmdevice = i+1;
 i++;
      }else if(strcmp(argv[i], "-device2") == 0) {
 pgplot_mapdevice = i+1;
 i++;
      }else if(strcmp(argv[i], "-device3") == 0) {
 pgplot_profiledevice = i+1;
 i++;
      }else {
 if(argv[i][0] == '-') {
   printerror(application.verbose_state.debug, "rmsynth: Unknown option: %s\n\nRun rmsynth without command line arguments to show help", argv[i]);
   terminateApplication(&application);
   return 0;
 }else {
   if(applicationAddFilename(i, application.verbose_state) == 0)
     return 0;
 }
      }
    }
  }
  if(application.dodebase) {
    application.dodebase = 2;
  }
  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "rmsynth: No files specified");
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) > 1) {
    printerror(application.verbose_state.debug, "rmsynth: Does not support more than one input file");
    return 0;
  }
  if(nrrmsteps == 0) {
    printerror(application.verbose_state.debug, "Use the -rm option to specify range of rm's to explore");
    return 0;
  }
  if(nrrmsteps <= 1) {
    printerror(application.verbose_state.debug, "Use the -rm option to specify more than one rm to explore");
    return 0;
  }
  gsl_rng_env_setup();
  rand_num_gen_type = gsl_rng_default;
  rand_num_gen = gsl_rng_alloc(rand_num_gen_type);
  if(application.fixseed)
    idnum = 1;
  else
    randomize_idnum(&idnum);
  gsl_rng_set(rand_num_gen, idnum);
  if(!openPSRData(&datain, argv[argc-1], 0, 0, 1, 0, application.verbose_state))
    return 0;
  if(PSRDataHeader_parse_commandline(&datain, argc, argv, application.verbose_state) == 0)
    return 0;
  if(datain.isDeFarad == -1) {
    printwarning(application.verbose_state.debug, "WARNING rmsynth: De-Faraday rotation state is unknown. Assume the data is not yet de-Faraday rotated.\n");
    datain.isDeFarad = 0;
  }
  if(datain.isDebase == 0) {
    printerror(application.verbose_state.debug, "ERROR rmsynth: Baseline is not subtracted. Use pmod -debase first.\n");
    return 0;
  }else if(datain.isDebase != 1) {
    printwarning(application.verbose_state.debug, "WARNING rmsynth: It is not known if baseline is already subtracted, which might affect the fitting.\n");
  }
  if(check_baseline_subtracted(datain, application.verbose_state) == 0) {
    printwarning(application.verbose_state.debug, "WARNING rmsynth: Baseline does not appear to be subtracted, which might affect the fitting.\n");
  }
  region_frac_to_int(&(application.onpulse), datain.NrBins, 0);
    region_frac_to_int(&(application.onpulse2), datain.NrBins, 0);
  for(i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-header") == 0) {
      printwarning(application.verbose_state.debug, "WARNING: If using the -header option, be aware it applied BEFORE the preprocessing.");
    }
  }
  if(preprocessApplication(&application, &datain) == 0) {
    return 0;
  }
  if(preprocess_dedisperse(&datain, 0, 0, 0, application.verbose_state) == 0)
    return 0;
    if(!preprocess_addsuccessiveFreqChans(datain, &profiledata, datain.NrFreqChan, NULL, application.verbose_state))
      return 0;
  pgplot_clear_options(&pgplot_options);
  if(pgplot_profiledevice)
    strcpy(pgplot_options.viewport.plotDevice, argv[pgplot_profiledevice]);
  else
    strcpy(pgplot_options.viewport.plotDevice, "?");
  strcpy(pgplot_options.box.ylabel, "Stokes I");
  strcpy(pgplot_options.box.title, "Select on-pulse region (for noise calculation only)");
  if(application.onpulse2.nrRegions == 0) {
    selectRegions(profiledata.data, profiledata.NrBins, &pgplot_options, 0, 0, 0, &(application.onpulse2), application.verbose_state);
  }else {
    pgplotGraph1(&pgplot_options, profiledata.data, NULL, NULL, profiledata.NrBins, 0, profiledata.NrBins-1, 0, 0, profiledata.NrBins-1, 0, 0, 0, 1, 0, 1, 0, 1, 1, &(application.onpulse2), -1, application.verbose_state);
  }
    region_int_to_frac(&(application.onpulse2), 1.0/(float)datain.NrBins, 0);
  regionShowNextTimeUse(application.onpulse2, "-onpulse2", "-onpulsef2", stdout);
  if(bootstrap) {
    rms_channels = (float *)calloc(datain.NrFreqChan, sizeof(float));
    rm_av = (double *)calloc(datain.NrBins, sizeof(double));
    rm_square = (double *)calloc(datain.NrBins, sizeof(double));
    if(rm_av == NULL || rm_square == NULL || rms_channels == NULL) {
      printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot allocate memory");
      return 0;
    }
    for(i = 0; i < datain.NrFreqChan; i++) {
 if(read_rmsPSRData(datain, &rms_channels[i], NULL, NULL, &(application.onpulse2), 0, 0, i, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot determine RMS");
   return 0;
 }
    }
  }
  cleanPSRData(&clone, application.verbose_state);
  if(copy_params_PSRData(datain, &clone, application.verbose_state) == 0)
    return 0;
  clone.data = (float *)malloc(datain.NrBins*datain.NrFreqChan*datain.NrPols*sizeof(float));
  if(clone.data == NULL) {
    printerror(application.verbose_state.debug, "ERROR rmsynth: Memory allocation error");
    return 0;
  }
  rmsynth_array = NULL;
  cmap = (float *)malloc(nrrmsteps*datain.NrBins*sizeof(float));
  rmestimate = (double *)calloc(datain.NrBins, sizeof(double));
  rmcurestimate = (double *)calloc(datain.NrBins, sizeof(double));
  singlespectrum = malloc(nrrmsteps*sizeof(float));
  singlespectrum_double = malloc(nrrmsteps*sizeof(double));
  sigmagrid = malloc(nrrmsteps*sizeof(double));
  rmgrid = malloc(nrrmsteps*sizeof(double));
  if(cmap == NULL || rmestimate == NULL || rmcurestimate == NULL || rmgrid == NULL || singlespectrum == NULL || singlespectrum_double == NULL || sigmagrid == NULL) {
    printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot allocate memory");
    return 0;
  }
  for(bootstrap_itt = 0; bootstrap_itt <= bootstrap; bootstrap_itt++) {
    if(bootstrap_itt == 1 && application.verbose_state.verbose)
      printf("Starting bootstrap\n");
    if(application.verbose_state.debug)
      printf("  itteration %d: Making clone of data with added noise\n", bootstrap_itt+1);
    for(freqchannelnr = 0; freqchannelnr < datain.NrFreqChan; freqchannelnr++) {
      for(binnr = 0; binnr < datain.NrBins; binnr++) {
 for(polnr = 0; polnr < datain.NrPols; polnr++) {
   if(readPulsePSRData(&datain, 0, polnr, freqchannelnr, binnr, 1, &sample, application.verbose_state) != 1) {
     printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot read data.");
     return 0;
   }
   if(bootstrap_itt > 0)
     sample += gsl_ran_gaussian(rand_num_gen, rms_channels[freqchannelnr]);
   if(writePulsePSRData(&clone, 0, polnr, freqchannelnr, binnr, 1, &sample, application.verbose_state) != 1) {
     printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot write data.");
     return 0;
   }
 }
      }
    }
    verbose_definition verbose2;
    copyVerboseState(application.verbose_state, &verbose2);
    verbose2.verbose = application.verbose_state.verbose;
    if(bootstrap_itt != 0)
      verbose2.verbose = 0;
    if(application.verbose_state.debug)
      verbose2.verbose = 1;
    if(application.verbose_state.debug)
      printf("  itteration %d: Applying RM synthesis\n", bootstrap_itt+1);
    if(rmSynthesis(clone, rmmin, rmmax, &rmsynth_array, nrrmsteps, &(application.onpulse), verbose2) == 0)
      return 0;
    if((application.verbose_state.verbose && bootstrap_itt == 0) || application.verbose_state.debug)
      printf("RM synthesis calculation done\n");
    if(bootstrap_itt == 0) {
      for(i = 0; i < datain.NrBins; i++) {
 for(j = 0; j < nrrmsteps; j++) {
   cmap[j*datain.NrBins + i] = rmsynth_array[2*(j*datain.NrBins+i)];
 }
      }
      pgplot_clear_options(&pgplot_options);
      if(pgplot_mapdevice)
 strcpy(pgplot_options.viewport.plotDevice, argv[pgplot_mapdevice]);
      else
 strcpy(pgplot_options.viewport.plotDevice, "?");
      pgplot_options.viewport.ysize = 0.7;
      pgplot_options.viewport.dontclose = 1;
      strcpy(pgplot_options.box.xlabel, "Pulse longitude (bin)");
      strcpy(pgplot_options.box.ylabel, "RM (rad/m\\u2\\d)");
      pgplotMap(&pgplot_options, cmap, datain.NrBins, nrrmsteps, 0, datain.NrBins-1, 0-0.5, datain.NrBins-1 + 0.5, rmmin, rmmax, rmmin, rmmax, PPGPLOT_HEAT, 0, 0, 0, NULL, 0, 0, 1.0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, application.verbose_state);
      pgplot_clear_options(&pgplot_options);
      if(pgplot_mapdevice)
 strcpy(pgplot_options.viewport.plotDevice, argv[pgplot_mapdevice]);
      else
 strcpy(pgplot_options.viewport.plotDevice, "?");
      pgplot_options.viewport.dontopen = 1;
      pgplot_options.viewport.noclear = 1;
      pgplot_options.viewport.ysize = 0.25;
      pgplot_options.viewport.dxplot = 0.05;
      pgplot_options.viewport.xsize = 0.75;
      pgplot_options.viewport.dyplot = 0.62;
      strcpy(pgplot_options.box.ylabel, "Pulse profile");
      strcpy(pgplot_options.box.title, "RM synthesis map");
      strcpy(pgplot_options.box.box_xopt, "bcsti");
      if(application.onpulse.nrRegions == 0) {
 selectRegions(profiledata.data, profiledata.NrBins, &pgplot_options, 0, 0, 0, &(application.onpulse), application.verbose_state);
      }else {
 pgplotGraph1(&pgplot_options, profiledata.data, NULL, NULL, profiledata.NrBins, 0, profiledata.NrBins-1, 0, 0, profiledata.NrBins-1, 0, 0, 0, 1, 0, 1, 0, 1, 1, &(application.onpulse), -1, application.verbose_state);
      }
      ppgend();
      region_int_to_frac(&(application.onpulse), 1.0/(float)datain.NrBins, 0);
      regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);
    }
    if(bootstrap_itt == 0) {
      pgplot_clear_options(&pgplot_options);
      if(pgplot_rmdevice)
 strcpy(pgplot_options.viewport.plotDevice, argv[pgplot_rmdevice]);
      else
 strcpy(pgplot_options.viewport.plotDevice, "?");
    }
    for(binnr = 0; binnr < datain.NrBins; binnr++) {
      if(collapse) {
 if(application.verbose_state.debug)
   printf("  itteration %d: Collapsing RM synthesis spectrum\n", bootstrap_itt+1);
 collapseRMSynthesisArray(rmsynth_array, nrrmsteps, datain.NrBins, application.onpulse, singlespectrum, application.verbose_state);
 for(i = 0; i < nrrmsteps; i++) {
   singlespectrum_double[i] = singlespectrum[i];
 }
      }else {
 if(application.verbose_state.debug)
   printf("  itteration %d: Extracting RM synthesis spectrum for bin %ld\n", bootstrap_itt+1, binnr);
 for(i = 0; i < nrrmsteps; i++) {
   singlespectrum[i] = rmsynth_array[2*(i*datain.NrBins+binnr)];
   singlespectrum_double[i] = rmsynth_array[2*(i*datain.NrBins+binnr)];
 }
      }
      if(checkRegions(binnr, &application.onpulse, 0, application.verbose_state) != 0 || application.onpulse.nrRegions == 0 || collapse) {
 if(application.verbose_state.debug)
   printf("  itteration %d: Plotting spectrum for bin %ld\n", bootstrap_itt+1, binnr);
 for(i = 0; i < nrrmsteps; i++) {
   rmgrid[i] = rmmin + (rmmax-rmmin)*i/(float)(nrrmsteps-1);
   sigmagrid[i] = 1.0;
 }
 pgplot_options.viewport.dontclose = 1;
 strcpy(pgplot_options.box.ylabel, "RM synthesis power");
 strcpy(pgplot_options.box.xlabel, "RM (rad/m\\u2\\d)");
 pgplotGraph1(&pgplot_options, singlespectrum, NULL, NULL, nrrmsteps, rmmin, rmmax, 0, rmmin, rmmax, 0, 0, 0, 0, 0, 1, -10, 1, 1, NULL, -1, application.verbose_state);
 pgplot_options.viewport.dontopen = 1;
 if(parabola) {
   if(application.verbose_state.debug)
     printf("  itteration %d: Start fitting of parabola for bin %ld\n", bootstrap_itt+1, binnr);
   function.nrfuncs = 3;
   function.func[0].type = FUNC_POLYNOMAL;
   function.func[0].param[0] = 0;
   function.func[0].start[0] = singlespectrum[0];
   function.func[0].value[0] = singlespectrum[0];
   function.func[0].fit_flag[0] = 1;
   function.func[1].type = FUNC_POLYNOMAL;
   function.func[1].param[0] = 1;
   function.func[1].start[0] = 0.0;
   function.func[1].value[0] = 0.0;
   function.func[1].fit_flag[0] = 1;
   function.func[2].type = FUNC_POLYNOMAL;
   function.func[2].param[0] = 2;
   function.func[2].start[0] = 0.0;
   function.func[2].value[0] = 0.0;
   function.func[2].fit_flag[0] = 1;
   if(application.verbose_state.debug) {
   }
   if(ftol < 0) {
     ftol = 1e-5;
   }
   fit_levmar(1, &function, rmgrid, singlespectrum_double, sigmagrid, nrrmsteps, 0, 1, 0, ftol, 1000, &status, 0, 0, application.verbose_state);
   if(status != 0) {
     printerror(application.verbose_state.debug, "Fitting failed: error %d", status);
     switch(status) {
     case 1: printerror(application.verbose_state.debug, "Memory allocation error"); return 0;
     case 2: printerror(application.verbose_state.debug, "Max number of itterations reached"); return 0;
     case 3: printwarning(application.verbose_state.debug, "WARNING: Machine precision limit reached"); break;
     case 4: printerror(application.verbose_state.debug, "Cannot determine suitable trial step."); return 0;
     case 10: printerror(application.verbose_state.debug, "Shouldn't happen, not documented error code of gsl"); return 0;
     }
   }
   if(application.verbose_state.debug) {
     printf("  itteration %d: Start fitting of parabola for bin %ld done\n", bootstrap_itt+1, binnr);
     print_fitfunctions(&function, 0, 0, -1, application.verbose_state);
   }
   ppgbbuf();
   ppgsci(2);
   for(i = 0; i < nrrmsteps; i++) {
     if(i == 0)
       ppgmove(rmgrid[i], evaluate_fitfunc_collection(&function, rmgrid[i], application.verbose_state));
     else
       ppgdraw(rmgrid[i], evaluate_fitfunc_collection(&function, rmgrid[i], application.verbose_state));
   }
   ppgebuf();
   rmcurestimate[binnr] = -0.5*function.func[1].value[0]/function.func[2].value[0];
 }else {
   float best_rm, best_offset, best_scale, *rmsynth_responds;
   if(datain.freqMode != FREQMODE_UNIFORM) {
     fflush(stdout);
     printwarning(application.verbose_state.debug, "WARNING rmsynth: Frequency range is non-uniform. The instrumental responds is determined using the non-weighted frequencies.");
   }
   if(application.verbose_state.debug)
     printf("  itteration %d: Fitting expected instrumental shape for bin %ld\n", bootstrap_itt+1, binnr);
   double chanbw;
   if(get_channelbandwidth(datain, &chanbw, application.verbose_state) == 0) {
     printerror(application.verbose_state.debug, "ERROR rmsynth (%s): Cannot obtain channel bandwidth.", datain.filename);
     return 0;
   }
   if(ftol < 0) {
     ftol = 0.001;
   }
   if(rmSynthesis_fitInstrumentalResponds(singlespectrum, rmmin, rmmax, nrrmsteps, &best_rm, &best_offset, &best_scale, datain.NrFreqChan, chanbw, get_nonweighted_channel_freq(datain, 0, application.verbose_state), ftol, application.verbose_state) != 0) {
     printerror(application.verbose_state.debug, "ERROR rmsynth: Fitting of instrumental responds shape failed");
     return 0;
   }
   if(application.verbose_state.debug)
     printf("  itteration %d: Generating instrumental shape for bin %ld\n", bootstrap_itt+1, binnr);
   rmsynth_responds = NULL;
   if(rmSynthesis_instrument_responds(datain.NrFreqChan, chanbw, get_nonweighted_channel_freq(datain, 0, application.verbose_state), rmmin, rmmax, &rmsynth_responds, NULL, 0, nrrmsteps, best_rm, 0, application.verbose_state) == 0) {
     return 0;
   }
   rmcurestimate[binnr] = best_rm;
   if(application.verbose_state.debug)
     printf("  itteration %d: Plotting instrumental shape for bin %ld\n", bootstrap_itt+1, binnr);
   ppgbbuf();
   ppgsci(2);
   for(i = 0; i < nrrmsteps; i++) {
     rmsynth_responds[i] *= best_scale;
     rmsynth_responds[i] += best_offset;
     if(i == 0)
       ppgmove(rmgrid[i], rmsynth_responds[i]);
     else
       ppgdraw(rmgrid[i], rmsynth_responds[i]);
   }
   ppgebuf();
   free(rmsynth_responds);
 }
 if(!isnormal(rmcurestimate[binnr])) {
   printwarning(application.verbose_state.debug, "Not a normal number for RM: %f", rmcurestimate[binnr]);
   not_normal_result_warning = 1;
 }
      }else {
 if(application.verbose_state.debug)
   printf("  itteration %d: Ignoring bin %ld\n", bootstrap_itt+1, binnr);
      }
      if(bootstrap_itt > 0) {
 rm_av[binnr] += rmcurestimate[binnr];
 rm_square[binnr] += rmcurestimate[binnr]*rmcurestimate[binnr];
      }else {
 rmestimate[binnr] = rmcurestimate[binnr];
      }
      if(collapse)
 break;
      if(application.verbose_state.debug)
 printf("\n");
      if(application.verbose_state.nocounters == 0 && bootstrap && application.verbose_state.verbose && bootstrap_itt > 0) {
 printf("\r%.2f%%       ", 100.0*((bootstrap_itt-1)*datain.NrBins+binnr)/(float)(datain.NrBins*bootstrap));
 fflush(stdout);
      }
      if(application.verbose_state.debug)
 printf("\n");
    }
  }
  if(application.verbose_state.debug)
    printf("\n");
  if(bootstrap && application.verbose_state.verbose) {
    printf("\rDone                    \n");
  }
  if(application.verbose_state.verbose)
    printf("Creating de-Faraday rotated profile\n");
  if(application.verbose_state.verbose)
    printf("de-Faraday rotate data\n");
  if(collapse) {
    datain.rm = rmestimate[0];
  }else {
    nrOnpulseBins = 0;
    datain.rm = 0;
    for(binnr = 0; binnr < datain.NrBins; binnr++) {
      if(checkRegions(binnr, &application.onpulse, 0, application.verbose_state) != 0 || application.onpulse.nrRegions == 0) {
 datain.rm += rmestimate[binnr];
 nrOnpulseBins++;
      }
    }
    datain.rm /= (double)nrOnpulseBins;
    printwarning(application.verbose_state.debug, "WARNING: To get the analytic errorbar, the average of the longitude-resolved RM was used = %f", datain.rm);
  }
  int ret = preprocess_deFaraday(&datain, 0, 0, 0, NULL, application.verbose_state);
  if(ret == 1) {
    if(fabs(datain.rm) < 1e-3) {
      printwarning(application.verbose_state.debug, "WARNING rmsynth: de-Faraday rotation wasn't applied because of the very low RM. Something appears to be strange.");
      ret = 2;
    }else {
      printerror(application.verbose_state.debug, "de-Faraday rotation failed, cannot determine analytic errorbar");
      return 0;
    }
  }
  if(ret != 2) {
    printerror(application.verbose_state.debug, "de-Faraday rotation failed, cannot determine analytic errorbar");
    return 0;
  }
  if(application.verbose_state.verbose)
    printf("Adding frequency channels\n");
  closePSRData(&clone, 0, application.verbose_state);
  if(!preprocess_addsuccessiveFreqChans(datain, &clone, datain.NrFreqChan, NULL, application.verbose_state))
    return 0;
  if(application.verbose_state.verbose)
    printf("Calculating L from Q and U\n");
    if(make_paswing_fromIQUV(&clone, 0, application.onpulse2, 0, 1, 1, 1, 0, 0.0, 0.0, NULL, 1.0, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "Creation of L profile failed, cannot determine analytic errorbar");
      return 0;
    }
    double sumL, fwhm, s2nL;
    sumL = 0;
    nrOnpulseBins = 0;
    for(binnr = 0; binnr < datain.NrBins; binnr++) {
      if(checkRegions(binnr, &application.onpulse2, 0, application.verbose_state) != 0 || application.onpulse2.nrRegions == 0) {
 sumL += clone.data[datain.NrBins+binnr];
 nrOnpulseBins++;
      }
    }
    printwarning(application.verbose_state.debug, "WARNING: The analytic errorbar is based on the total L (within -onpulse2) as determined using the overall RM (using -onpulse).");
    if(collapse) {
      s2nL = sumL/(sqrt(nrOnpulseBins)*clone.offpulse_rms[1]);
    }else {
      s2nL = (sumL/(double)nrOnpulseBins)/(sqrt(nrOnpulseBins)*clone.offpulse_rms[1]);
      printwarning(application.verbose_state.debug, "WARNING: The analytic errorbar as computed is per definition the same for each pulse longitude bin.");
    }
    if(application.verbose_state.verbose) {
      printf("Total L = %f in %d bins\n", sumL, nrOnpulseBins);
      printf("RMS L   = %f\n", clone.offpulse_rms[1]);
      printf("S/N L   = %f\n", s2nL);
    }
    fwhm = get_FWHM_RMsynthesis(datain, application.verbose_state);
    expectedRMerror = 0.5*fwhm/s2nL;
    if(application.verbose_state.verbose) {
      printf("Expected FWHM RM synthesis peak = %f\n", fwhm);
      printf("Expected errorbar on RM = %f\n", expectedRMerror);
    }
  closePSRData(&clone, 0, application.verbose_state);
  if(collapse || write_ascii == 0) {
    for(binnr = 0; binnr < datain.NrBins; binnr++) {
      if(checkRegions(binnr, &application.onpulse, 0, application.verbose_state) != 0 || application.onpulse.nrRegions == 0 || collapse) {
 printf("RM ");
 if(collapse == 0)
   printf("for bin %ld ", binnr);
 printf("is %lf", rmestimate[binnr]);
 if(bootstrap) {
   rmsigma = sqrt(rm_square[binnr]/(double)bootstrap - rm_av[binnr]*rm_av[binnr]/(double)(bootstrap*bootstrap));
   printf(" +- %f (from bootstrapping)", rmsigma);
   if(!isnormal(rmsigma)) {
     printwarning(application.verbose_state.debug, "Not a normal number for rmsigma: %f (%f %f)", rmsigma, rm_square[binnr], rm_av[binnr]);
     not_normal_result_warning = 2;
   }
 }
 printf(" +- %f (analytic prediction for uniform frequency channel weights)", expectedRMerror);
 printf("\n");
 if(collapse)
   break;
      }
    }
  }
  if(write_ascii) {
    if(change_filename_extension(argv[argc-1], outputname, "RMtable", 1000, application.verbose_state) == 0)
      return 0;
    ofile = fopen(outputname, "w+");
    if(ofile == NULL) {
      printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot open file %s for writing", outputname);
      return 0;
    }
    fprintf(ofile, "#");
    if(collapse == 0)
      fprintf(ofile, "binnr Profile ");
    fprintf(ofile, "RM RMerror(analytic)");
    if(bootstrap)
      fprintf(ofile, " RMerror(bootstrap)");
    fprintf(ofile, "\n");
    for(binnr = 0; binnr < datain.NrBins; binnr++) {
      if(checkRegions(binnr, &application.onpulse, 0, application.verbose_state) != 0 || application.onpulse.nrRegions == 0 || collapse) {
 if(collapse == 0)
   fprintf(ofile, "%ld %f ", binnr, profiledata.data[binnr]);
 fprintf(ofile, "%lf ", rmestimate[binnr]);
 fprintf(ofile, "%lf ", expectedRMerror);
 if(bootstrap) {
   rmsigma = sqrt(rm_square[binnr]/(float)bootstrap - rm_av[binnr]*rm_av[binnr]/(float)(bootstrap*bootstrap));
   fprintf(ofile, " %f", rmsigma);
 }
 fprintf(ofile, "\n");
 if(collapse)
   break;
      }
    }
    fclose(ofile);
  }
  ppgend();
  free(singlespectrum);
  free(singlespectrum_double);
  free(rmgrid);
  free(sigmagrid);
  free(rmsynth_array);
  free(cmap);
  free(rmestimate);
  free(rmcurestimate);
  closePSRData(&datain, 0, application.verbose_state);
  closePSRData(&profiledata, 0, application.verbose_state);
  gsl_rng_free (rand_num_gen);
  if(bootstrap) {
    free(rms_channels);
    free(rm_av);
    free(rm_square);
  }
  if(not_normal_result_warning == 1) {
    printwarning(application.verbose_state.debug, "\nNote that warnings were generated for one or more suspicious values found for the RM, suggesting something is wrong.");
  }else if(not_normal_result_warning == 2) {
    printwarning(application.verbose_state.debug, "\nNote that warnings were generated for one or more suspicious values found for the error on the RM, suggesting something is wrong.");
  }
  terminateApplication(&application);
  return 0;
}
double get_FWHM_RMsynthesis(datafile_definition datain, verbose_definition verbose)
{
  int freqnr;
  double lambda0, lambda, c, f, norm, sigmal2, fwhm, channelbw;
  c = 299792458.0;
  lambda0 = 0;
  norm = 0;
  if(get_channelbandwidth(datain, &channelbw, verbose) == 0) {
    printerror(verbose.debug, "ERROR get_FWHM_RMsynthesis (%s): Cannot obtain channel bandwidth.", datain.filename);
    exit(0);
  }
  for(freqnr = 0; freqnr < datain.NrFreqChan; freqnr++) {
    f = 1e6*get_nonweighted_channel_freq(datain, freqnr, verbose);
    lambda = c/f;
    lambda0 += lambda*lambda;
    norm += 1;
  }
  lambda0 /= norm;
  lambda0 = sqrt(lambda0);
  sigmal2 = 0;
  for(freqnr = 0; freqnr < datain.NrFreqChan; freqnr++) {
    f = 1e6*get_nonweighted_channel_freq(datain, freqnr, verbose);
    lambda = c/f;
    sigmal2 += pow(lambda, 4) - pow(lambda0, 4);
  }
  sigmal2 /= (double)(datain.NrFreqChan-1);
  sigmal2 = sqrt(sigmal2);
  fwhm = 1.0/(sigmal2);
  return fwhm;
}

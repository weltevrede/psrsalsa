/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "psrsalsa.h"
#include "gsl/gsl_randist.h"
int readZapFile(char *zaplistname, int *zapMask, int zapSkipLines, int nrZapCols, int zapColumn, int zapColumn2, int inverseZap, verbose_definition verbose);
void make_blocks(long baseline_length, long blockSize, long nrPulses, long *nrOutputBlocks, int *zapMask, verbose_definition verbose);
int main(int argc, char **argv)
{
  int debase_flag, debase_offset_flag, index, deviceOpened, read_whole_file;
  int zapoption, inverseZap, fzapoption, finverseZap, zapColumn, zapColumn2, nrZapCols, zapSkipLines;
  int blockMode, remove_pulses_flag, prange_set;
  int nrPol, nrBins, NrFreqChan, addnoise_flag, removeOnPulse_flag;
  int selectMoreOnpulseRegions;
  int outputlist, filename;
  int *zapMask, *fzapMask;
  long i, j, k, l, m, tmp_used_pulses;
  long nrOutputBlocks, nrZapped, baseline_length, blockSize, firstPulseToKeep, lastPulseToKeep, nrPulses;
  long dataout2_pulse, pulse_nr_in_output, idnum;
  float *profileI, debase_offset_value, *baseline, *runningBaseline, *runningRMS, *rms, noiseRMS;
  char *filename_ptr, zaplistname[MaxFilenameLength], fzaplistname[MaxFilenameLength];
  char output_suffix[MaxFilenameLength], output_suffix2[MaxFilenameLength], output_name[MaxFilenameLength], output_name2[MaxFilenameLength];
  char txt[1000];
  datafile_definition datain, dataout, dataout2;
  psrsalsaApplication application;
  gsl_rng *rand_num_gen;
  const gsl_rng_type *rand_num_gen_type;
  pgplot_options_definition pgplot_options;
  initApplication(&application, "pmod", "[options] inputfile(s)");
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.oformat = FITS_format;
  application.switch_formatlist = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_onpulse = 1;
  application.switch_onpulsef = 1;
  application.switch_filelist = 1;
  application.switch_device = 1;
  application.switch_conshift= 1;
  application.switch_circshift= 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_rebin = 1;
  application.switch_tscr = 1;
  application.switch_TSCR = 1;
  application.switch_tscr_complete = 1;
  application.switch_fscr = 1;
  application.switch_FSCR = 1;
  application.switch_nocounters = 1;
  application.switch_dedisperse = 1;
  application.switch_deFaraday = 1;
  application.switch_polselect = 1;
  application.switch_stokes = 1;
  application.switch_coherence = 1;
  application.switch_noweights = 1;
  application.switch_useweights = 1;
  application.switch_uniformweights = 1;
  application.switch_scale = 1;
  application.switch_insertparang = 1;
  application.switch_deparang = 1;
  application.switch_changeRefFreq = 1;
  application.switch_history_cmd_only = 1;
  application.switch_fixseed = 1;
  application.switch_templatedata = 1;
  application.switch_template = 1;
  application.switch_align = 1;
  application.switch_shuffle = 1;
  application.switch_rotateStokes = 1;
  application.switch_libversions = 1;
  debase_flag = 0;
  debase_offset_flag = 0;
  read_whole_file = 1;
  zapoption = 0;
  inverseZap = -1;
  finverseZap = -1;
  fzapoption = 0;
  nrZapCols = 1;
  zapColumn = 1;
  zapColumn2 = 0;
  zapSkipLines = 0;
  blockMode = 0;
  remove_pulses_flag = 0;
  prange_set = 0;
  debase_offset_value = 0;
  baseline_length = 0;
  nrOutputBlocks = 0;
  outputlist = 0;
  addnoise_flag = 0;
  noiseRMS = 0;
  selectMoreOnpulseRegions = 0;
  removeOnPulse_flag = 0;
  strcpy(output_suffix, "debase.gg");
  strcpy(output_suffix2, "zapped.gg");
  filename = 0;
  pgplot_clear_options(&pgplot_options);
  if(argc < 2) {
    printf("Program to modify pulsar data in various ways. Usage:\n\n");
    printApplicationHelp(&application);
    printf("Action options:\n\n");
    printf("-debase              Subtract baseline from observation\n");
    printf("-debase_length       Specify number of pulses before and after the pulse to\n");
    printf("                     determine running baseline (default is -debase_length 0).\n");
    printf("                     Note that the number of pulses used is 2N+1\n");
    printf("-debase_value        Subtract the specified value baseline (fixed value) rather\n");
    printf("                     than the determined baseline, i.e. -debase_value 1.234\n");
    printf("-list                List the zap list to a file\n");
    printf("-addnoise rms        Add white noise with RMS rms to the data.\n");
    printf("-onpulse_subst_noise Substitute the onpulse region with white noise based on the\n");
    printf("                     (running) off-pulse rms\n");
    printf("\nData selection options:\n\n");
    printf("-zapfile file     Specify filename with pulse numbers to zap (first pulse is 0).\n");
    printf("                  Expected format: see -format.\n");
    printf("-zapfile_i        or alternatively, specify file with pulse numbers NOT to zap\n");
    printf("-zap              or alternatively, specify first and last pulse to zap, \n");
    printf("                  e.g. -zap \"0 10\" (can specify -zap multiple times)\n");
    printf("-prange           or alternatively, specify first and last pulse to keep\n");
    printf("-format           Specify column number, total number columns and lines to skip\n");
    printf("                  in zap file (default is \"%d %d %d\")\n", zapColumn, nrZapCols, zapSkipLines);
    printf("-fzapfile         Specify file with freq channels to zap (first channel is 0)\n");
    printf("-fzapfile_i       or alternatively, specify file with channels NOT to zap\n");
    printf("-zapfile2         If this flag is specified, the -zapfile or -zapfile_i or\n");
    printf("                  equivalent frequency channel zap options should specify a file\n");
    printf("                  with two columns, being the start and end values of ranges\n");
    printf("                  rather than individual values.\n");
    printf("-fzap             Specify first and last frequency channel to zap\n");
    printf("                  (can specify -fzap multiple times)\n");
    printf("-blocksize        Nr of succesive nonzapped pulses that should be writen out\n");
    printf("-remove           Remove zapped pulses instead of making them zero\n");
    printf("\nOn-pulse definition:\n\n");
    printf("-onpulsegr        Enables selecting more on-pulse regions graphically than\n");
    printf("                  defined by -onpulse\n");
    printf("\nOutput options:\n\n");
    printf("-ext              Specify suffix, default is '%s'\n", output_suffix);
    printf("-output filename  Write output to filename rather than changing the extension\n");
    printf("                  to '%s'.\n", output_suffix);
    printf("-memsave          Don't read the file in as a whole at the start of the program\n");
    printf("\n");
    printCitationInfo();
   terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
 i = index;
      }else if(strcmp(argv[i], "-output") == 0) {
 filename = i+1;
 i++;
      }else if(strcmp(argv[i], "-debase") == 0) {
 debase_flag = 1;
      }else if(strcmp(argv[i], "-list") == 0) {
 outputlist = 1;
      }else if(strcmp(argv[i], "-debase_value") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &debase_offset_value, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 debase_offset_flag = 1;
 debase_flag = 1;
 i++;
      }else if(strcmp(argv[i], "-onpulse_subst_noise") == 0) {
 removeOnPulse_flag = 1;
      }else if(strcmp(argv[i], "-debase_length") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &baseline_length, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 if(baseline_length != 0)
   remove_pulses_flag = 1;
        i++;
      }else if(strcmp(argv[i], "-zapfile") == 0) {
 strcpy(zaplistname, argv[i+1]);
 if(zapoption) {
   printerror(application.verbose_state.debug, "pmod: Can only either use -zapfile or -zapfile_i, and they cannot be used more than once.");
   return 0;
 }
 if(inverseZap == 1) {
   printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
   return 0;
 }
 zapoption = 1;
 inverseZap = 0;
        i++;
      }else if(strcmp(argv[i], "-zapfile_i") == 0) {
 strcpy(zaplistname, argv[i+1]);
 zapoption = 1;
 inverseZap = 1;
 if(inverseZap == 0) {
   printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-zapfile2") == 0) {
 if(zapColumn2 == 0) {
   zapColumn2 = zapColumn+1;
 }
      }else if(strcmp(argv[i], "-fzapfile") == 0) {
 strcpy(fzaplistname, argv[i+1]);
 fzapoption = 1;
 if(finverseZap == 1) {
   printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
   return 0;
 }
 finverseZap = 0;
        i++;
      }else if(strcmp(argv[i], "-fzapfile_i") == 0) {
 strcpy(fzaplistname, argv[i+1]);
 fzapoption = 1;
 if(finverseZap == 0) {
   printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
   return 0;
 }
 finverseZap = 1;
        i++;
      }else if(strcmp(argv[i], "-blocksize") == 0) {
 blockMode = 1;
 remove_pulses_flag = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &blockSize, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-format") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d %d", &zapColumn, &nrZapCols, &zapSkipLines, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-prange") == 0) {
 prange_set = 1;
 if(inverseZap == 0) {
   printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
   return 0;
 }
 inverseZap = 1;
        i++;
      }else if(strcmp(argv[i], "-zap") == 0) {
 if(inverseZap == 1) {
   printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-fzap") == 0) {
 if(finverseZap == 1) {
   printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-addnoise") == 0) {
 addnoise_flag = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &noiseRMS, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-ext") == 0) {
 strcpy(output_suffix, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-remove") == 0) {
 remove_pulses_flag = 1;
      }else if(strcmp(argv[i], "-onpulsegr") == 0) {
 selectMoreOnpulseRegions = 1;
      }else if(strcmp(argv[i], "-memsave") == 0) {
 read_whole_file = 0;
      }else {
 if(argv[i][0] == '-') {
   printerror(application.verbose_state.debug, "pmod: Unknown option: %s\n\nRun pmod without command line arguments to show help", argv[i]);
   terminateApplication(&application);
   return 0;
 }else {
   if(applicationAddFilename(i, application.verbose_state) == 0)
     return 0;
 }
      }
    }
  }
  if(inverseZap == -1)
    inverseZap = 0;
  if(finverseZap == -1)
    finverseZap = 0;
  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pmod: No files specified");
    return 0;
  }
  strcpy(pgplot_options.viewport.plotDevice, application.pgplotdevice);
  pgplot_options.viewport.dontclose = 1;
  gsl_rng_env_setup();
  rand_num_gen_type = gsl_rng_default;
  rand_num_gen = gsl_rng_alloc(rand_num_gen_type);
  if(application.fixseed)
    idnum = 1;
  else
    randomize_idnum(&idnum);
  gsl_rng_set(rand_num_gen, idnum);
  while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
    deviceOpened = 0;
    if(filename == 0) {
      if(change_filename_extension(filename_ptr, output_name, output_suffix, MaxFilenameLength, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pmod: Cannot change extension in output name.");
 return 0;
      }
    }else {
      strcpy(output_name, argv[filename]);
    }
    if(change_filename_extension(filename_ptr, output_name2, output_suffix2, MaxFilenameLength, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR pmod: Cannot change extension in output name.");
      return 0;
    }
    cleanPSRData(&dataout, application.verbose_state);
    cleanPSRData(&dataout2, application.verbose_state);
    if(application.iformat <= 0)
      application.iformat = guessPSRData_format(filename_ptr, 0, application.verbose_state);
    if(isValidPSRDATA_format(application.iformat) == 0) {
      printerror(application.verbose_state.debug, "ERROR pmod: Please specify a valid input format with the -iformat option.");
      closePSRData(&dataout, 0, application.verbose_state);
      closePSRData(&dataout2, 0, application.verbose_state);
      gsl_rng_free(rand_num_gen);
      terminateApplication(&application);
      return 0;
    }
    i = openPSRData(&datain, filename_ptr, application.iformat, 0, read_whole_file, 0, application.verbose_state);
    if(i == 0) {
      printerror(application.verbose_state.debug, "ERROR pmod: Error opening data");
      return 0;
    }
    if(!read_whole_file) {
      if(readHeaderPSRData(&datain, 0, 0, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "pmod: Error reading header");
 return 0;
      }
    }
    if(PSRDataHeader_parse_commandline(&datain, argc, argv, application.verbose_state) == 0)
      return 0;
    region_frac_to_int(&(application.onpulse), datain.NrBins, 0);
    for(i = 1; i < argc; i++) {
      if(strcmp(argv[i], "-header") == 0) {
 fflush(stdout);
 printwarning(application.verbose_state.debug, "WARNING: If using the -header option, be aware it applied BEFORE the preprocessing.");
 break;
      }
    }
    if(preprocessApplication(&application, &datain) == 0) {
      return 0;
    }
    zapMask = (int *)malloc(datain.NrSubints*sizeof(int));
    fzapMask = (int *)malloc(datain.NrFreqChan*sizeof(int));
    if(zapMask == NULL || fzapMask == NULL) {
      printerror(application.verbose_state.debug, "pmod: Memory allocation error.");
      return 0;
    }
    for(i = 0; i < datain.NrSubints; i++)
      zapMask[i] = inverseZap;
    for(i = 0; i < datain.NrFreqChan; i++)
      fzapMask[i] = finverseZap;
    if(argc > 2) {
      for(i = 1; i < argc-1; i++) {
 if(strcmp(argv[i], "-zap") == 0) {
   if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld %ld", &k, &l, NULL) == 0) {
     printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
     return 0;
   }
   for(j = k; j <= l; j++) {
     if(j >= 0 && j < datain.NrSubints)
       zapMask[j] = 1;
   }
   i++;
 }
      }
      for(i = 1; i < argc-1; i++) {
 if(strcmp(argv[i], "-fzap") == 0) {
   if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld %ld", &k, &l, NULL) == 0) {
     printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
     return 0;
   }
   for(j = k; j <= l; j++) {
     if(j >= 0 && j < datain.NrFreqChan) {
       fzapMask[j] = 1;
     }
   }
   i++;
 }
      }
    }
    if(zapoption != 0) {
      if(readZapFile(zaplistname, zapMask, zapSkipLines, nrZapCols, zapColumn, zapColumn2, inverseZap, application.verbose_state) == 0)
 return 0;
      if(application.verbose_state.verbose) {
 j = 0;
 for(i = 0; i < datain.NrSubints; i++) {
   if(zapMask[i] != 0)
     j++;
 }
 printf("Zapfile zapped %ld pulses.\n", j);
      }
    }
    if(fzapoption != 0) {
      if(readZapFile(fzaplistname, fzapMask, zapSkipLines, nrZapCols, zapColumn, zapColumn2, finverseZap, application.verbose_state) == 0)
 return 0;
      if(application.verbose_state.verbose) {
 j = 0;
 for(i = 0; i < datain.NrFreqChan; i++) {
   if(fzapMask[i] != 0)
     j++;
 }
 printf("Zapfile zapped %ld/%ld frequency channels.\n", j, datain.NrFreqChan);
 if(j == 0)
   application.fzapMask = NULL;
 else
   application.fzapMask = fzapMask;
      }
    }
    if(prange_set) {
      for(i = 1; i < argc-1; i++) {
 if(strcmp(argv[i], "-prange") == 0) {
   if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld %ld", &firstPulseToKeep, &lastPulseToKeep, NULL) == 0) {
     printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
     return 0;
   }
   if(lastPulseToKeep >= datain.NrSubints) {
     printerror(application.verbose_state.debug, "ERROR pmod: Invalid -prange option. Last subint to keep %ld >= %ld", lastPulseToKeep, datain.NrSubints);
     return 0;
   }
   if(firstPulseToKeep < 0) {
     printerror(application.verbose_state.debug, "ERROR pmod: Invalid -prange option. First subint to keep %ld < 0", firstPulseToKeep);
     return 0;
   }
   for(j = firstPulseToKeep; j <= lastPulseToKeep; j++)
     zapMask[j] = 0;
   i++;
 }
      }
    }
    nrPol = datain.NrPols;
    nrBins = datain.NrBins;
    nrPulses= datain.NrSubints;
    NrFreqChan = datain.NrFreqChan;
    if(application.verbose_state.verbose) printf("Input data contains %d bins, %ld pulses, %d polarizations and %d frequencies.\n", nrBins, nrPulses, nrPol, NrFreqChan);
    profileI = (float *)malloc(nrBins*sizeof(float));
    rms = (float *)malloc(nrPulses*nrPol*sizeof(float));
    runningBaseline = (float *)malloc(nrPulses*nrPol*sizeof(float));
    baseline = (float *)malloc(nrPulses*nrPol*sizeof(float));
    runningRMS = (float *)malloc(nrPulses*nrPol*sizeof(float));
    if(profileI == NULL || rms == NULL || runningBaseline == NULL || baseline == NULL || runningRMS == NULL
       ) {
      printerror(application.verbose_state.debug, "ERROR pmod: Memory allocation error.");
      return 0;
    }
    if(debase_flag || removeOnPulse_flag
       ) {
      if(read_profilePSRData(datain, profileI, zapMask, 0, application.verbose_state) != 1) {
 printerror(application.verbose_state.debug, "Error pmod: Cannot form pulse profile");
 return 0;
      }
    }
      if((debase_flag && debase_offset_flag == 0) || removeOnPulse_flag
  ) {
 if(selectMoreOnpulseRegions || application.onpulse.nrRegions == 0) {
   pgplot_options.viewport.dontopen = deviceOpened;
   strcpy(pgplot_options.box.xlabel, "Bin");
   strcpy(pgplot_options.box.ylabel, "Intensity");
   strcpy(pgplot_options.box.title, "Select onpulse region(s) (which will be ignored from baseline calculation) ");
   strcat(pgplot_options.box.title, datain.psrname);
   selectRegions(profileI, nrBins, &pgplot_options, 0, 0, 0, &(application.onpulse), application.verbose_state);
   deviceOpened = 1;
 }else {
   if(application.verbose_state.debug) {
     pgplot_options.viewport.dontopen = deviceOpened;
     strcpy(pgplot_options.box.xlabel, "Bin");
     strcpy(pgplot_options.box.ylabel, "I");
     strcpy(pgplot_options.box.title, "Profile");
     pgplotGraph1(&pgplot_options, profileI, NULL, NULL, nrBins, 0, nrBins, 0, 0, nrBins, 0, 0, 0, 1, 0, 0, 1, 1, &(application.onpulse), application.verbose_state);
     deviceOpened = 1;
     printf("Press return to continue\n");
     scanf("%c", txt);
   }
 }
      }
    region_int_to_frac(&(application.onpulse), 1.0/(float)datain.NrBins, 0);
    regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);
    if(baseline_length) {
      for(i = 0; i < baseline_length; i++)
 zapMask[i] = 0;
      for(i = nrPulses-baseline_length; i < nrPulses; i++)
 zapMask[i] = 0;
    }
    if(blockMode != 0)
      make_blocks(baseline_length, blockSize, nrPulses, &nrOutputBlocks, zapMask, application.verbose_state);
    if(application.verbose_state.debug) {
      if(baseline_length)
 printf("DEBUG: Zapping %ld subints because using a running baseline\n", 2*baseline_length);
    }
    nrZapped = 0;
    int prev_zap_state = -1;
    long prev_pulse_nr = 0;
    for(i = baseline_length; i < nrPulses-baseline_length; i++) {
      if(zapMask[i]) {
 nrZapped++;
      }
      if(application.verbose_state.debug) {
 if(prev_zap_state != zapMask[i]) {
   if(zapMask[i]) {
     if(i != baseline_length)
       printf("%ld = %ld subints\n", i-1, i-prev_pulse_nr);
     printf("DEBUG: Zapping subint %ld - ", i);
     prev_pulse_nr = i;
   }else {
     if(i != baseline_length)
       printf("%ld = %ld subints\n", i-1, i-prev_pulse_nr);
     printf("DEBUG: Keeping subint %ld - ", i);
     prev_pulse_nr = i;
   }
 }
 prev_zap_state = zapMask[i];
      }
    }
    if(application.verbose_state.debug) {
      printf("%ld = %ld subints\n", nrPulses-baseline_length-1, nrPulses-baseline_length-prev_pulse_nr);
    }
    if(application.verbose_state.verbose) printf("Total number of zapped pulses: %ld\n", nrZapped);
    if(outputlist) {
      FILE *fout_zap;
      fout_zap = fopen("zaplist.txt", "w");
      if(fout_zap == NULL) {
 printerror(application.verbose_state.debug, "pmod: Cannot open zaplist.txt");
 return 0;
      }
      for(i = 0; i < nrPulses; i++) {
 if(zapMask[i] != 0) {
   fprintf(fout_zap, "%ld\n", i);
 }
      }
      fclose(fout_zap);
      printf("Zap list is written to zaplist.txt\n");
    }
    copy_params_PSRData(datain, &dataout, application.verbose_state);
    copy_params_PSRData(datain, &dataout2, application.verbose_state);
    if(debase_flag) {
      dataout.isDebase = 1;
      dataout2.isDebase = 1;
    }
    dataout.NrSubints = datain.NrSubints - 2*baseline_length;
    if(remove_pulses_flag)
      dataout.NrSubints -= nrZapped;
    if(application.verbose_state.verbose) printf("Output data contains %ld bins, %ld pulses, %ld polarizations and %ld frequencies.\n", dataout.NrBins, dataout.NrSubints, dataout.NrPols, dataout.NrFreqChan);
    if(openPSRData(&dataout, output_name, application.oformat, 1, 0, 0, application.verbose_state) == 0) {
      printf("Cannot open %s\n\n", output_name);
      return 0;
    }
    if(writeHeaderPSRData(&dataout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
      printf("Cannot write header to %s\n\n", output_name);
      return 0;
    }
    if(zapoption
       ) {
      dataout2.NrBins = nrBins;
      dataout2.NrSubints = nrZapped;
      if(dataout2.NrSubints > 0) {
 if(openPSRData(&dataout2, output_name2, application.oformat, 1, 0, 0, application.verbose_state) == 0) {
   printf("Cannot open %s\n\n", output_name);
   return 0;
 }
 if(writeHeaderPSRData(&dataout2, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
   printf("Cannot write header to %s\n\n", output_name);
   return 0;
 }
      }
    }
      for(k=0; k < nrPol; k++) {
 for(l = 0; l < NrFreqChan; l++) {
   if(l == 0) {
     if(application.verbose_state.verbose) {
       printf("Processing polarization channel %ld (of the %d)\n", k+1, nrPol);
     }
   }
   if(debase_flag || removeOnPulse_flag) {
     read_rmsPSRData(datain, &rms[datain.NrSubints*k], &runningBaseline[datain.NrSubints*k], zapMask, &(application.onpulse), 0, k, l, application.verbose_state);
     if(l == 0 && application.verbose_state.verbose) {
       if(datain.NrFreqChan > 1)
  printf("  Preparing baseline for first frequency channel\n");
       else
  printf("  Preparing baseline\n");
     }
     for(i = 0; i < baseline_length; i++) {
       baseline[i+datain.NrSubints*k] = runningBaseline[i+datain.NrSubints*k];
       runningRMS[i+datain.NrSubints*k] = rms[i+datain.NrSubints*k];
     }
     for(i = nrPulses-baseline_length; i < nrPulses; i++) {
       baseline[i+datain.NrSubints*k] = runningBaseline[i+datain.NrSubints*k];
       runningRMS[i+datain.NrSubints*k] = rms[i+datain.NrSubints*k];
     }
     int itteration, ok;
     float rms_old;
     rms_old = 0;
     for(i = baseline_length; i < nrPulses-baseline_length; i++) {
       for(itteration = 0; itteration < 2; itteration++) {
  baseline[i+datain.NrSubints*k] = 0;
  runningRMS[i+datain.NrSubints*k] = 0;
  tmp_used_pulses = 0;
  for(j = -baseline_length; j <= baseline_length; j++) {
    if(runningBaseline[i+j+datain.NrSubints*k] != 0.0) {
        ok = 1;
      if(ok) {
        baseline[i+datain.NrSubints*k] += runningBaseline[i+j+datain.NrSubints*k];
        runningRMS[i+datain.NrSubints*k] += rms[i+j+datain.NrSubints*k];
        rms_old += runningBaseline[i+j+datain.NrSubints*k]*runningBaseline[i+j+datain.NrSubints*k];
        tmp_used_pulses++;
      }
    }
  }
  if(tmp_used_pulses > 0) {
    baseline[i+datain.NrSubints*k] /= (float)(tmp_used_pulses);
    runningRMS[i+datain.NrSubints*k] /= (float)(tmp_used_pulses);
    rms_old /= (float)(tmp_used_pulses);
    rms_old -= baseline[i+datain.NrSubints*k]*baseline[i+datain.NrSubints*k];
  }
    break;
       }
       m = i/10;
       if((10*m == i) && (m > 0) && application.verbose_state.nocounters == 0 && l == 0 && application.verbose_state.verbose) {
  printf("    %.1f%%     \r",(100.0*(i-baseline_length))/(float)(nrPulses-2.0*baseline_length));
  fflush(stdout);
       }
     }
     if(l == 0 && application.verbose_state.verbose)
       printf("    Done                      \n");
   }
   if(debase_flag && l == 0) {
     if(application.verbose_state.verbose)
       printf("  Subtracting baseline\n");
   }
   dataout2_pulse = 0;
   pulse_nr_in_output = 0;
   for(i = 0; i < nrPulses; i++) {
     readPulsePSRData(&datain, i, k, l, 0, nrBins, profileI, application.verbose_state);
     if(debase_flag && debase_offset_flag == 0) {
       for(j = 0; j < nrBins; j++) {
  if(runningBaseline[i+datain.NrSubints*k] != 0.0) {
    profileI[j] -= baseline[i+datain.NrSubints*k];
  }
       }
     }
     if(debase_flag && debase_offset_flag) {
       for(j = 0; j < nrBins; j++) {
  if(runningBaseline[i+datain.NrSubints*k] != 0.0) {
    profileI[j] -= debase_offset_value;
  }
       }
     }
     if(removeOnPulse_flag) {
       for(j = 0; j < nrBins; j++) {
  if(checkRegions(j, &(application.onpulse), 0, application.verbose_state) != 0) {
    profileI[j] = gsl_ran_gaussian(rand_num_gen, runningRMS[i+datain.NrSubints*k]);
  }
       }
     }
     if(addnoise_flag) {
       for(j = 0; j < nrBins; j++) {
  profileI[j] += gsl_ran_gaussian(rand_num_gen, noiseRMS);
       }
     }
     if(zapoption
        ) {
       if(zapMask[i] != 0) {
  if(writePulsePSRData(&dataout2, dataout2_pulse, k, l, 0, nrBins, profileI, application.verbose_state) != 1) {
    printerror(application.verbose_state.debug, "pmod: Writing data failed.");
    return 0;
  }
  dataout2_pulse += 1;
       }
     }
     if(zapMask[i] != 0) {
       for(j = 0; j < nrBins; j++)
  profileI[j] = 0;
       baseline[i] = 0;
     }
     if(fzapMask[l] != 0) {
       for(j = 0; j < nrBins; j++)
  profileI[j] = 0;
       baseline[i] = 0;
     }
     if(i >= baseline_length && i < nrPulses-baseline_length) {
       if(remove_pulses_flag == 0 || zapMask[i] == 0) {
  if(remove_pulses_flag == 0) {
    pulse_nr_in_output = i;
  }
    if(writePulsePSRData(&dataout, pulse_nr_in_output, k, l, 0, nrBins, profileI, application.verbose_state) != 1) {
      printerror(application.verbose_state.debug, "pmod: Writing data failed.");
      return 0;
    }
  pulse_nr_in_output++;
       }
     }
     m = i/10;
     if((10*m == i) && (m > 0) && application.verbose_state.nocounters == 0 && l != 0 && application.verbose_state.verbose) {
       printf("   %.1f%%     \r",(100.0*(i+nrPulses*l))/(float)(nrPulses*NrFreqChan));
       fflush(stdout);
     }
   }
   if(NrFreqChan > 1 && l == 0 && application.verbose_state.verbose)
     printf("  Processing other frequency channels\n");
 }
 if(application.verbose_state.verbose)
   printf("    Done                             \n");
 if(k == 0 && (debase_flag && debase_offset_flag==0) && nrPulses-2*baseline_length > 1) {
   strcpy(txt, "Baseline ");
   strcat(txt, datain.psrname);
   pgplot_options.viewport.dontopen = deviceOpened;
   strcpy(pgplot_options.box.xlabel, "Pulse number");
   strcpy(pgplot_options.box.ylabel, "Baseline");
   strcpy(pgplot_options.box.title, txt);
   pgplotGraph1(&pgplot_options, baseline+baseline_length, NULL, NULL, nrPulses-2*baseline_length, baseline_length, nrPulses-baseline_length, 0, baseline_length, nrPulses-baseline_length, 0, 0, 0, 0, 0, 0, 1, 1, NULL, application.verbose_state);
   deviceOpened = 1;
 }
      }
    free(profileI);
    free(baseline);
    free(rms);
    free(runningBaseline);
    free(runningRMS);
    free(zapMask);
    free(fzapMask);
    closePSRData(&datain, 0, application.verbose_state);
    closePSRData(&dataout, 0, application.verbose_state);
    closePSRData(&dataout2, 0, application.verbose_state);
    if(deviceOpened) {
      ppgclos();
      ppgend();
    }
    printf("Writing of %s done\n", output_name);
  }
  gsl_rng_free(rand_num_gen);
  terminateApplication(&application);
  return 0;
}
int readZapFile(char *zaplistname, int *zapMask, int zapSkipLines, int nrZapCols, int zapColumn, int zapColumn2, int inverseZap, verbose_definition verbose)
{
  FILE *fin_zap;
  long i, i2, j;
  char txt[100];
  if(zapColumn2 > 0) {
    if(zapColumn2 <= zapColumn) {
      printerror(verbose.debug, "pmod: column %d was expected to be larger than %d.", zapColumn2, zapColumn);
      return 0;
    }
  }
  fin_zap = fopen(zaplistname, "r");
  if(fin_zap == NULL) {
    printerror(verbose.debug, "pmod: Cannot open %s", zaplistname);
    return 0;
  }else {
    fprintf(stdout, "Opened zapfile '%s'\n", zaplistname);
  }
  if(zapSkipLines > 0) {
    printf("Skipping %d lines\n", zapSkipLines);
    for(i = 0; i < zapSkipLines; i++) {
      do {
 j = fgetc(fin_zap);
      }while(j !='\n');
    }
  }
  do {
    if(zapColumn > 1) {
      for(i = 0; i < zapColumn-1; i++) {
 j = fscanf(fin_zap, "%s", txt);
 if(j != 1)
   break;
      }
    }
    j = fscanf(fin_zap, "%ld", &i);
    if(zapColumn2 > 0 && j == 1) {
      if(zapColumn2 - zapColumn > 1) {
 for(i = zapColumn+1; i < zapColumn2; i++) {
   j = fscanf(fin_zap, "%s", txt);
   if(j != 1)
     break;
 }
      }
      j = fscanf(fin_zap, "%ld", &i2);
    }
    if(j == 1) {
      if(verbose.debug) {
 printf("DEBUG: read in value %ld", i);
 if(zapColumn2 > 0)
   printf(" and %ld", i2);
 printf(" (inverseZap = %d)\n", inverseZap);
      }
      if(zapColumn2 <= 0) {
 if(inverseZap == 0) {
   zapMask[i] = 1;
 }else {
   zapMask[i] = 0;
 }
      }else {
 long loopnr;
 for(loopnr = i; loopnr <= i2; loopnr++) {
   if(inverseZap == 0) {
     zapMask[loopnr] = 1;
   }else {
     zapMask[loopnr] = 0;
   }
 }
      }
    }
    long nskip;
    if(zapColumn2 > 0)
      nskip = nrZapCols - zapColumn2;
    else
      nskip = nrZapCols - zapColumn;
    if(nskip > 0) {
      for(i = 0; i < nrZapCols - zapColumn; i++) {
 j = fscanf(fin_zap, "%s", txt);
 if(j != 1)
   break;
      }
    }
  }while(j == 1);
  fclose(fin_zap);
  return 1;
}
void make_blocks(long baseline_length, long blockSize, long nrPulses, long *nrOutputBlocks, int *zapMask, verbose_definition verbose)
{
  long i, j, k, nrZapped;
  nrZapped = 0;
  i = baseline_length;
  do {
    k = 1;
    for(j = i; j < i+blockSize; j++) {
      if(j == nrPulses-baseline_length) {
 k = -j;
 break;
      }
      if(zapMask[j]) {
 k = -j;
 break;
      }
    }
    if(k <= 0) {
      for(j = i; j <= -k; j++) {
 if(j < nrPulses) {
   zapMask[j] = 1;
   nrZapped++;
 }
      }
      i = -k+1;
    }else {
      i += blockSize;
      nrOutputBlocks++;
    }
  }while(i < nrPulses-baseline_length);
  if(verbose.verbose) printf("%ld pulses zapped in total to make integer number of %ld pulse blocks.\n", nrZapped, blockSize);
}

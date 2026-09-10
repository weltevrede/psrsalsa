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
#define _LARGEFILE_SOURCE 1
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_sort.h>
#include "psrsalsa.h"
int padd_check_parameters(long poladd, long freqadd, long sumNsub, int circularShift, int noinput, int shift, int memsave,
     verbose_definition verbose);
int main(int argc, char **argv)
{
  char output_fname[1000], PlotDevice[100], *inputname;
  int circularShift, noinput, onlyI, memsave, currentfilenumber, dummy_int;
  int shift, poladd, freqadd, ignore_dedisp_check;
  long i, nrinputfiles, sumNsub;
  float *Iprofile_firstfile, *shiftedProfile;
  datafile_definition **fin;
  datafile_definition fout;
  psrsalsaApplication application;
  initApplication(&application, "padd", "[options] inputfiles");
  application.switch_blocksize = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_filelist = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.switch_formatlist = 1;
  application.switch_fscr = 1;
  application.switch_FSCR = 1;
  application.switch_nocounters = 1;
  application.switch_changeRefFreq = 1;
  application.switch_noweights = 1;
  application.switch_useweights = 1;
  application.switch_uniformweights = 1;
  application.switch_history_cmd_only = 1;
  application.oformat = FITS_format;
  strcpy(PlotDevice, "?");
  onlyI = 0;
  strcpy(output_fname, "addedfile.gg");
  circularShift = 0;
  noinput = 0;
  shift = 0;
  memsave = 0;
  sumNsub = 1;
  poladd = 0;
  freqadd = 0;
  ignore_dedisp_check = 0;
  if(argc < 2) {
    printf("Program to add data files together. Usage:\n\n");
    printApplicationHelp(&application);
    printf("Optional options:\n");
    printf("-w                      Output name. Default is \"%s\"\n", output_fname);
    printf("-I                      Only process the first polarization channel\n");
    printf("-c                      Turn on circular shifting (so that no subint is lost).\n");
    printf("                        The first and last subintegrations are spilling over\n");
    printf("                        into each other.\n");
    printf("-n                      No graphical input, circular shifting by this number of\n");
    printf("                        bins, so this option implies -c. So -n 0 results in a\n");
    printf("                        simple concatenation of the input files.\n");
    printf("-nsub                   This number of subintegrations are summed before being\n");
    printf("                        written out\n");
    printf("-poladd                 The input files are to be interpretted as separate\n");
    printf("                        polarization channels (option implies -n 0).\n");
    printf("-freqadd                The input files are to be interpretted as separate\n");
    printf("                        frequency bands (option implies -n 0).\n");
    printf("-memsave                Only one full input data-set exists in memory at a time,\n");
    printf("                        but every input file will be opened twice.\n");
    printf("\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      dummy_int = i;
      if(processCommandLine(&application, argc, argv, &dummy_int)) {
 i = dummy_int;
      }else if(strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "-W") == 0) {
 strcpy(output_fname,argv[i+1]);
        i++;
      }else if(strcmp(argv[i], "-memsave") == 0) {
 memsave = 1;
      }else if(strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "-C") == 0) {
 circularShift = 1;
      }else if(strcmp(argv[i], "-I") == 0) {
 onlyI = 1;
      }else if(strcmp(argv[i], "-poladd") == 0) {
 poladd = 1;
 circularShift = 1;
 noinput = 1;
 shift = 0;
      }else if(strcmp(argv[i], "-freqadd") == 0) {
 freqadd = 1;
 circularShift = 1;
 noinput = 1;
 shift = 0;
      }else if(strcmp(argv[i], "-no_dedisp_check") == 0) {
 ignore_dedisp_check = 1;
      }else if(strcmp(argv[i], "-n") == 0) {
 circularShift = 1;
 noinput = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &shift, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR padd: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-nsub") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &sumNsub, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR padd: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 if(sumNsub < 1) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR padd: Cannot parse option %s, expected one integer number > 1", argv[i]);
   return 0;
 }
 i++;
      }else {
 if(argv[i][0] == '-') {
   printerror(application.verbose_state.debug, "ERROR padd: Unknown option %s. Run padd without command-line options to get help.", argv[i]);
   terminateApplication(&application);
   return 0;
 }else {
   if(applicationAddFilename(i, application.verbose_state) == 0)
     return 0;
 }
      }
    }
  }
  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  nrinputfiles = numberInApplicationFilenameList(&application, argv, application.verbose_state);
  if(nrinputfiles < 2) {
    printerror(application.verbose_state.debug, "ERROR padd: Need at least two input files");
    return 0;
  }
  if(padd_check_parameters(poladd, freqadd, sumNsub, circularShift, noinput, shift, memsave,
      application.verbose_state) == 0)
    {
      return 0;
    }
  fin = malloc(nrinputfiles*sizeof(datafile_definition *));
  if(fin == NULL) {
    printerror(application.verbose_state.debug, "ERROR padd: Memory allocation error");
    return 0;
  }
  for(i = 0; i < nrinputfiles; i++) {
    fin[i] = malloc(sizeof(datafile_definition));
    if(fin[i] == NULL) {
      printerror(application.verbose_state.debug, "ERROR padd: Memory allocation error");
      return 0;
    }
  }
  currentfilenumber = 0;
  while((inputname = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
    verbose_definition verbose2;
    copyVerboseState(application.verbose_state, &verbose2);
    verbose2.indent = application.verbose_state.indent + 2;
    if(memsave == 0 || (currentfilenumber == 0 && noinput == 0)) {
      if(currentfilenumber == 0) {
 printf("Read in input files:\n");
      }
      if(openPSRData(fin[currentfilenumber], inputname, application.iformat, 0, 1, 0, application.obsnr, verbose2) == 0) {
 printerror(application.verbose_state.debug, "ERROR padd: Cannot open %s\n", inputname);
 return 0;
      }
      if(currentfilenumber == 0) {
 for(i = 1; i < argc; i++) {
   if(strcmp(argv[i], "-header") == 0) {
     fflush(stdout);
     printwarning(application.verbose_state.debug, "WARNING: If using the -header option, be aware it applied BEFORE the preprocessing.");
   }
 }
      }
      if(preprocessApplication(&application, fin[currentfilenumber]) == 0) {
 printerror(application.verbose_state.debug, "ERROR padd: preprocess option failed on file %s\n", inputname);
 return 0;
      }
    }else {
      if(currentfilenumber == 0)
 printf("Read in headers of input files:\n");
      if(openPSRData(fin[currentfilenumber], inputname, application.iformat, 0, 0, 0, application.obsnr, verbose2) == 0) {
 printerror(application.verbose_state.debug, "ERROR padd: Cannot open %s\n", inputname);
 return 0;
      }
      if(readHeaderPSRData(fin[currentfilenumber], 1, 0, application.obsnr, verbose2) == 0) {
 printerror(application.verbose_state.debug, "ERROR padd: Cannot read header of file %s\n", inputname);
 return 0;
      }
    }
    if(currentfilenumber == 0 && noinput == 0) {
      Iprofile_firstfile = (float *)malloc(fin[0]->NrPols*fin[0]->NrBins*sizeof(float));
      shiftedProfile = (float *)malloc(fin[0]->NrPols*fin[0]->NrBins*sizeof(float));
      if(Iprofile_firstfile == NULL || shiftedProfile == NULL) {
 printerror(application.verbose_state.debug, "ERROR padd: Cannot allocate memory.");
 return 0;
      }
      if(read_profilePSRData(*fin[0], Iprofile_firstfile, NULL, 0, application.verbose_state) != 1) {
 printerror(application.verbose_state.debug, "ERROR padd: Reading pulse profile of first input file failed.");
 return 0;
      }
    }
    if(memsave) {
      if(closePSRData(fin[currentfilenumber], 1, 0, application.verbose_state) != 0) {
 printerror(application.verbose_state.debug, "ERROR padd: Closing file %s failed\n", inputname);
 return 0;
      }
    }
    if(fin[currentfilenumber]->NrBins != fin[0]->NrBins) {
      printerror(application.verbose_state.debug, "ERROR padd: Number of pulse longitude bins not equal in input files.");
      return 0;
    }
    if(fin[0]->NrPols != fin[currentfilenumber]->NrPols && onlyI == 0) {
      printerror(application.verbose_state.debug, "ERROR padd: Number of polarization channels are different. Maybe you want to use the -I option?");
      return 0;
    }
    if(freqadd == 0) {
      if(fin[currentfilenumber]->NrFreqChan != fin[0]->NrFreqChan) {
 printerror(application.verbose_state.debug, "ERROR padd: Number of frequency channels not equal in input files.");
 return 0;
      }
    }
    if(poladd || freqadd) {
      if(fin[0]->NrSubints != fin[currentfilenumber]->NrSubints) {
 printerror(application.verbose_state.debug, "ERROR padd: Number of subints are different. This is not allowed with the -poladd or -freqadd options.");
 return 0;
      }
    }
    if(poladd) {
      if(fin[currentfilenumber]->NrPols != 1) {
 printerror(application.verbose_state.debug, "ERROR padd: Number of polarization channels in input files should be 1 if using the -poladd option.");
 return 0;
      }
    }
    if(ignore_dedisp_check == 0) {
      if(fin[currentfilenumber]->isDeDisp != fin[0]->isDeDisp) {
 printerror(application.verbose_state.debug, "ERROR padd: Dedispersion state is not equal in input files.");
 return 0;
      }
    }
    if(fin[currentfilenumber]->isDeFarad != fin[0]->isDeFarad) {
      printerror(application.verbose_state.debug, "ERROR padd: De-Faraday rotation state is not equal in input files.");
      return 0;
    }
    if(fin[currentfilenumber]->isDePar != fin[0]->isDePar) {
      printerror(application.verbose_state.debug, "ERROR padd: Parallactic angle state is not equal in input files.");
      return 0;
    }
    currentfilenumber++;
  }
  printf("Reading of input files done\n");
  int noclosePSRData = 0;
  padd_do_stuff(fin, nrinputfiles, &fout, application.oformat, output_fname,
  poladd, freqadd, noinput, PlotDevice, shiftedProfile, Iprofile_firstfile, onlyI, circularShift, shift, sumNsub,
  memsave, argc, argv, &application, application.history_cmd_only, noclosePSRData, 0, application.verbose_state);
  for(i = 0; i < nrinputfiles; i++) {
    free(fin[i]);
  }
  free(fin);
  closePSRData(&fout, 0, 0, application.verbose_state);
  if(noinput == 0)
    ppgend();
  if(noinput == 0) {
    free(shiftedProfile);
    free(Iprofile_firstfile);
  }
  terminateApplication(&application);
  return 0;
}

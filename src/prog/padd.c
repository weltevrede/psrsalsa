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
#include "psrsalsa.h"

#define MaxNrOnpulseRegions 10



int NrOnpulseRegions, OnPulseRegion[MaxNrOnpulseRegions][2];

void PlotProfile(int NrBins, float *Ipulse, char *xlabel, char *ylabel, char *title, int Highlight, int color, int clearPage);
void ShiftProfile(int shift, int NrBins, float *Iprofile, float *outputProfile);

int main(int argc, char **argv)
{
  char output_fname[1000], PlotDevice[100], singlechar, *inputname;
  int deviceID, circularShift, noinput, onlyI, memsave, currentfilenumber, dummy_int;
  int shift, bin1, bin2, subintWritten, poladd;
  long i, j, pol, fchan, nsub, binnr, nrinputfiles, subintslost, currentOutputSubint, sumNsub, curNrInsubint;
  float *Iprofile, *Iprofile_firstfile, *shiftedProfile, *subint, x, y, *float_ptr, *float_ptr2;
  datafile_definition **fin;
  datafile_definition fout, clone;


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



  application.oformat = FITS_format;

  strcpy(PlotDevice, "?");
  onlyI = 0;
  strcpy(output_fname, "addedfile.gg");
  NrOnpulseRegions = 0;
  circularShift = 0;
  noinput = 0;
  shift = 0;
  memsave = 0;
  sumNsub = 1;
  poladd = 0;
  x = y = singlechar = 0;

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


  if(sumNsub != 1 && poladd) {
    printerror(application.verbose_state.debug, "ERROR padd: Cannot use the -nsub and -poladd flag simultaneously.");
    return 0;
  }
  if(circularShift != 1 && poladd) {
    printerror(application.verbose_state.debug, "ERROR padd: You must use the -c option together with -poladd.");
    return 0;
  }
  if((noinput != 1 || shift != 0) && poladd) {
    printerror(application.verbose_state.debug, "ERROR padd: You must use the -n 0 option together with -poladd.");
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
      if(currentfilenumber == 0)
 printf("Read in input files:\n");
      if(openPSRData(fin[currentfilenumber], inputname, application.iformat, 0, 1, 0, verbose2) == 0) {
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
      if(openPSRData(fin[currentfilenumber], inputname, application.iformat, 0, 0, 0, verbose2) == 0) {
 printerror(application.verbose_state.debug, "ERROR padd: Cannot open %s\n", inputname);
 return 0;
      }
      if(readHeaderPSRData(fin[currentfilenumber], 1, 0, verbose2) == 0) {
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
      if(closePSRData(fin[currentfilenumber], 1, application.verbose_state) != 0) {
 printerror(application.verbose_state.debug, "ERROR padd: Closing file %s failed\n", inputname);
 return 0;
      }
    }


    if(fin[currentfilenumber]->NrFreqChan != fin[0]->NrFreqChan) {
      printerror(application.verbose_state.debug, "ERROR padd: Nr of frequency channel are not equal in input files.");
      return 0;
    }
    if(fin[currentfilenumber]->NrBins != fin[0]->NrBins) {
      printerror(application.verbose_state.debug, "ERROR padd: Nr of bins in not equal in input files.");
      return 0;
    }
    if(fin[0]->NrPols != fin[currentfilenumber]->NrPols && onlyI == 0) {
      printerror(application.verbose_state.debug, "ERROR padd: Number of polarization channels are different. Maybe you want to use the -I option?");
      return 0;
    }
    if(poladd) {
      if(fin[0]->NrSubints != fin[currentfilenumber]->NrSubints) {
 printerror(application.verbose_state.debug, "ERROR padd: Number of subints are different. This is not allowed with the -poladd option.");
 return 0;
      }
      if(fin[currentfilenumber]->NrPols != 1) {
 printerror(application.verbose_state.debug, "ERROR padd: Number of polarization in input files should be 1 if using the -poladd option.");
 return 0;
      }
    }
    currentfilenumber++;
  }
  printf("Reading of input files done\n");



  cleanPSRData(&fout, application.verbose_state);
  copy_params_PSRData(*(fin[0]), &fout, application.verbose_state);
  if(onlyI)
    fout.NrPols = 1;
  fout.format = application.oformat;
  fout.NrSubints = 0;

  subintslost = 0;
  if(poladd) {
    fout.NrPols = nrinputfiles;
    fout.NrSubints = fin[0]->NrSubints;
  }else {
    for(i = 0; i < nrinputfiles; i++) {
      fout.NrSubints += fin[i]->NrSubints;


      if(circularShift == 0) {
 fout.NrSubints -= 1;
 subintslost += 1;

      }



    }
  }
  printf("\nInput data contains %ld subints, %ld phase bins %ld frequency channels and %ld polarizations.\n", fout.NrSubints+subintslost, fout.NrBins, fout.NrFreqChan, fout.NrPols);
  printf("%ld subints are lost because of the alignment of input data", subintslost);
  if(subintslost > 0)
    printf(" (consider using -c option)");
  fout.tsubMode = TSUBMODE_TSUBLIST;
  if(fout.tsub_list != NULL)
    free(fout.tsub_list);
  fout.tsub_list = (double *)malloc(fout.NrSubints*sizeof(double));
  if(fout.tsub_list == NULL) {
    printerror(application.verbose_state.debug, "ERROR padd: Memory allocation error");
    return 0;
  }

  currentOutputSubint = 0;
  fout.tsub_list[0] = 0;
  subintWritten = 0;
  curNrInsubint = 0;
  if(poladd) {
    for(nsub = 0; nsub < fin[0]->NrSubints; nsub++) {
      fout.tsub_list[nsub] = get_tsub(*(fin[0]), nsub, application.verbose_state);
    }
  }else {
    for(i = 0; i < nrinputfiles; i++) {
      for(nsub = 0; nsub < fin[i]->NrSubints; nsub++) {

 if(currentOutputSubint < fout.NrSubints)
   fout.tsub_list[currentOutputSubint] += get_tsub(*(fin[i]), nsub, application.verbose_state);
 if(sumNsub == 1) {
   subintWritten = 1;
 }else if(curNrInsubint == sumNsub - 1) {
   subintWritten = 1;
 }
 curNrInsubint++;
 if(subintWritten) {
   subintWritten = 0;
   curNrInsubint = 0;
   currentOutputSubint++;
   if(currentOutputSubint < fout.NrSubints)
     fout.tsub_list[currentOutputSubint] = 0;
 }
      }
    }
  }

  dummy_int = fout.NrSubints % sumNsub;
  fout.NrSubints = fout.NrSubints/sumNsub;
  printf("\nOutput data will contain %ld subints", fout.NrSubints);
  if(sumNsub > 1)
    printf(" after summing every %ld input subints", sumNsub);
  if(dummy_int)
    printf(" (%d input subints lost because of incomplete last subint)", dummy_int);
  printf("\n\n");


  if(fout.gentype == GENTYPE_PULSESTACK && sumNsub != 1) {
    if(fout.NrSubints != 1)
      fout.gentype = GENTYPE_SUBINTEGRATIONS;
    else
      fout.gentype = GENTYPE_PROFILE;
  }



  if(openPSRData(&fout, output_fname, fout.format, 1, 0, 0, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR padd: Cannot open %s", output_fname);
    return 0;
  }
  if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) != 1) {
    printerror(application.verbose_state.debug, "ERROR padd: Cannot write header to %s", output_fname);
    return 0;
  }



  Iprofile = (float *)malloc(fout.NrPols*fout.NrBins*sizeof(float));
  if(Iprofile == NULL) {
    printerror(application.verbose_state.debug, "ERROR padd: Cannot allocate memory.");
    return 0;
  }
  if(sumNsub > 1) {
    subint = (float *)malloc(fout.NrPols*fout.NrBins*fout.NrFreqChan*sizeof(float));
    if(subint == NULL) {
      printerror(application.verbose_state.debug, "ERROR padd: Cannot allocate memory.");
      return 0;
    }
  }


  if(noinput == 0) {
    deviceID = ppgopen(PlotDevice);
    ppgask(0);
    ppgslw(1);
  }

  currentfilenumber = 0;
  currentOutputSubint = 0;
  curNrInsubint = 0;
  subintWritten = 0;
  rewindFilenameList(&application);
  while((inputname = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {




    if(memsave) {
      closePSRData(fin[currentfilenumber], 0, application.verbose_state);
      if(openPSRData(fin[currentfilenumber], inputname, application.iformat, 0, 1, 0, application.verbose_state) == 0) {
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
    }



    if(shift >= fout.NrBins)
      shift -= fout.NrBins;
    if(shift < 0)
      shift += fout.NrBins;


    if(noinput == 0 && currentfilenumber != 0) {
      if(read_profilePSRData(*fin[currentfilenumber], Iprofile, NULL, 0, application.verbose_state) != 1) {
 printerror(application.verbose_state.debug, "ERROR padd: Reading pulse profile failed.");
 return 0;
      }
      do {

 ShiftProfile(shift, fout.NrBins, Iprofile, shiftedProfile);
 PlotProfile(fout.NrBins, Iprofile_firstfile, "Bin number", "Intensity", "Click to shift profile, press s to stop", 0, 1, 1);
 PlotProfile(fout.NrBins, shiftedProfile, "", "", "", 0, 2, 0);


 j = 0;
 do {
   if(j == 0)
     ppgband(0, 0, 0.0, 0.0, &x, &y, &singlechar);
   else
     ppgband(4, 0, bin1, 0.0, &x, &y, &singlechar);
   if(singlechar == 65) {
     if(j == 0)
       bin1 = x;
     else
       bin2 = x;
     j++;
   }else if(singlechar == 115 || singlechar ==83 ) {
     j = 10;
   }
 }while(j < 2);
 if(j < 10) {
   shift += bin2 - bin1;
   if(shift >= fout.NrBins)
     shift -= fout.NrBins;
   if(shift < 0)
     shift += fout.NrBins;
 }
      }while(j < 10);
    }


    if(shift != 0) {
      if(continuous_shift(*fin[currentfilenumber], &clone, shift, circularShift, "padd", MEMORY_format, 0, NULL, application.verbose_state, application.verbose_state.debug) != 1) {
 printerror(application.verbose_state.debug, "ERROR padd: circular shift failed.");
      }
      swap_orig_clone(fin[currentfilenumber], &clone, application.verbose_state);
    }

    long nrPulsesInCurFile;
    nrPulsesInCurFile = fin[currentfilenumber]->NrSubints;
    if(shift == 0 && circularShift == 0)
      nrPulsesInCurFile -= 1;




    for(nsub = 0; nsub < nrPulsesInCurFile; nsub++) {
      int nrpolsinloop;
      nrpolsinloop = fout.NrPols;
      if(poladd)
 nrpolsinloop = 1;
      for(pol = 0; pol < nrpolsinloop; pol++) {
 for(fchan = 0; fchan < fout.NrFreqChan; fchan++) {
   if(readPulsePSRData(fin[currentfilenumber], nsub, pol, fchan, 0, fin[currentfilenumber]->NrBins, Iprofile, application.verbose_state) != 1) {
     printerror(application.verbose_state.debug, "ERROR padd: Read error");
     return 0;
   }
   if(sumNsub == 1) {

     if(poladd == 0) {
       if(writePulsePSRData(&fout, currentOutputSubint, pol, fchan, 0, fout.NrBins, Iprofile, application.verbose_state) != 1) {
  printerror(application.verbose_state.debug, "ERROR padd: Write error");
  return 0;
       }
     }else {
       if(writePulsePSRData(&fout, nsub, currentfilenumber, fchan, 0, fout.NrBins, Iprofile, application.verbose_state) != 1) {
  printerror(application.verbose_state.debug, "ERROR padd: Write error");
  return 0;
       }
     }
     subintWritten = 1;
   }else {
     float_ptr = &subint[fout.NrBins*(pol+fout.NrPols*fchan)];
     float_ptr2 = Iprofile;
     if(curNrInsubint == 0) {
       for(binnr = 0; binnr < fout.NrBins; binnr++) {
  *float_ptr = *float_ptr2;
  float_ptr++;
  float_ptr2++;
       }
     }else {
       for(binnr = 0; binnr < fout.NrBins; binnr++) {
  *float_ptr += *float_ptr2;
  float_ptr++;
  float_ptr2++;
       }
     }
     if(curNrInsubint == sumNsub - 1) {
       if(writePulsePSRData(&fout, currentOutputSubint, pol, fchan, 0, fout.NrBins, float_ptr-fout.NrBins, application.verbose_state) != 1) {
  printerror(application.verbose_state.debug, "ERROR padd: Write error");
  return 0;
       }
       subintWritten = 1;
     }
   }
 }
      }
      curNrInsubint++;
      if(subintWritten) {
 subintWritten = 0;
 curNrInsubint = 0;
 currentOutputSubint++;
      }
      if(application.verbose_state.nocounters == 0) {
 printf("Processing input file %d: %.1f%%     \r", currentfilenumber+1, (100.0*(nsub+1))/(float)(fin[currentfilenumber]->NrSubints));
 fflush(stdout);
      }
    }
    if(application.verbose_state.nocounters == 0) {
      printf("Processing file %d is done.                                \n", currentfilenumber+1);
    }

    closePSRData(fin[currentfilenumber], 0, application.verbose_state);
    currentfilenumber++;
  }

  closePSRData(&fout, 0, application.verbose_state);
  if(noinput == 0)
    ppgend();

  free(Iprofile);
  if(noinput == 0) {
    free(shiftedProfile);
    free(Iprofile_firstfile);
  }
  if(sumNsub > 1)
    free(subint);

  for(i = 0; i < nrinputfiles; i++) {
    free(fin[i]);
  }
  free(fin);
  terminateApplication(&application);
  return 0;
}


int CheckOnPulse(int bin, int NrRegions, int Regions[MaxNrOnpulseRegions][2])
{
  int i;
  for(i = 0; i < NrRegions; i++) {
    if(bin >= Regions[i][0] && bin <= Regions[i][1])
      return i+1;
  }
  return 0;
}

void PlotProfile(int NrBins, float *Ipulse, char *xlabel, char *ylabel, char *title, int Highlight, int color, int clearPage)
{
  long j;
  float ymin, ymax;

  ymin = ymax = Ipulse[0];
  for(j = 1; j < NrBins; j++) {
    if(Ipulse[j] > ymax)
      ymax = Ipulse[j];
    if(Ipulse[j] < ymin)
      ymin = Ipulse[j];
  }

  if(clearPage) {
    ppgpage();
    ppgsci(1);
    ppgsvp(0.1, 0.9, 0.1, 0.9);
    ppgswin(0,NrBins-1,-0.1,1.1);
    ppgbox("bcnsti",0.0,0,"bcnti",0.0,0);
    ppglab(xlabel, ylabel, title);
  }
  ppgsci(color);
  ppgmove(0, Ipulse[0]/ymax);
  for(j = 1; j < NrBins; j++) {
    if(Highlight != 0)
      ppgsci(color+CheckOnPulse(j,NrOnpulseRegions,OnPulseRegion));
    else
      ppgsci(color);
    ppgdraw(j, Ipulse[j]/ymax);
  }
  ppgsci(1);
}


void ShiftProfile(int shift, int NrBins, float *Iprofile, float *outputProfile)
{
  int b, b2;
  for(b = 0; b < NrBins; b++) {
    b2 = b+shift;
    if(b2 < 0)
      b2 += NrBins;
    if(b2 >= NrBins)
      b2 -= NrBins;
    outputProfile[b2] = Iprofile[b];
  }
}

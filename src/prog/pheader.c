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
psrsalsaApplication application;
int main(int argc, char **argv)
{
  int index, c_index, nrwords, iformat_initial_value, didhistory, didweights, didfreqlist, nohead, noweights, show_linenumber;
  int maxfilenamelength, maxobservatorylength, maxgentypelength, maxscanidlength, maxobserverlength, maxprojectidlength, maxinstrumentlength, maxfileformatlength, j;
  int showfootnotes, footnote_length, footnote_length2, footnote_search, footnote_parang, precision;
  long i;
  char *filename_ptr, cmd[5000];
  datafile_definition *datain;
  initApplication(&application, "pheader", "[options] inputfile(s)");
  c_index = 0;
  application.switch_formatlist = 1;
  application.switch_iformat = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_filelist = 1;
  nohead = 0;
  noweights = 1;
  showfootnotes = 1;
  footnote_length = 0;
  footnote_length2 = 0;
  footnote_search = 0;
  footnote_parang = 0;
  show_linenumber = 0;
  precision = 0;
  if(argc < 2) {
    printf("Program to show the header information of pulsar data. Usage:\n\n");
    printApplicationHelp(&application);
    printf("-c            Specify things to show (or run in verbose to show all)\n");
    printf("              Example: -c \"nbin nfreq\".\n");
    printf("-H            show list of things that can be specified with -c\n");
    printf("-linenr       Print line numbers\n");
    printf("-nohead       Do not print a header at the top of the table\n");
    printf("-nofootnotes  Do not print at the bottom of the table\n");
    printf("-precision d  Add d decimal places to floating point numbers\n");
    printf("\n");
    printCitationInfo();
    return 0;
  }else if(argc >= 2) {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
 i = index;
      }else if(strcmp(argv[i], "-nohead") == 0) {
 nohead = 1;
      }else if(strcmp(argv[i], "-nofootnotes") == 0) {
 showfootnotes = 0;
      }else if(strcmp(argv[i], "-precision") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &precision, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pheader: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-H") == 0) {
 printf("backend       Name of the backend\n");
 printf("bw            Bandwidth\n");
 printf("chbw          Channel bandwidth\n");
 printf("dedisp        De-dispersion state\n");
 printf("defarad       De-Faraday rotation state\n");
 printf("depar         Parallactic angle correction state\n");
 printf("dec           Declination\n");
 printf("dm            Dispersion measure\n");
 printf("dt            Time interval per phase bin\n");
 printf("format        Format identifier of the observation\n");
 printf("freq          Centre frequency\n");
 printf("freqlist      Show the (weighted) frequency for each channel/subint\n");
 printf("gentype       Identifies the type of data stored in file.\n");
 printf("hist          Show history (if supported in file)\n");
 printf("lat           Latitude of the telescope (GRS80 derived from ITRF)\n");
 printf("length        Duration of the observation (if defined).\n");
 printf("length2       Number bins times number subints times sampling time, so this time\n");
 printf("              can be less than the duration of the observation for folded data.\n");
 printf("long          Longitude of the telescope (GRS80 derived from ITRF)\n");
 printf("mjd           Start MJD of the observation\n");
 printf("name          Name of the pulsar\n");
 printf("nbin          Number of phase bins\n");
 printf("nbits         Number of bits per sample\n");
 printf("nchan         Number of frequency channels\n");
 printf("npol          Number of polarization channels\n");
 printf("nsub          Number of subints (or pulses in single-pulse data)\n");
 printf("observatory   Name of the observatory\n");
 printf("observer      Observer identifier string\n");
 printf("p0            Period\n");
 printf("parang        Show (derived) parallactic angle at midpoint observation\n");
 printf("parang1       Show (derived) parallactic angle for first subint\n");
 printf("poltype       Polarization type (Stokes, coherency, etc.)\n");
 printf("projectid     Project identifier string\n");
 printf("ra            Right ascension\n");
 printf("reffreq       Reference frequency (used for dedispersion etc)\n");
 printf("rm            Rotation measure\n");
 printf("scanid        Scan identifier string\n");
 printf("tsub          Average subintegration time, and a list of individual tsub times\n");
 printf("weights       Show weights (if supported in file)\n");
 return 0;
      }else if(strcmp(argv[i], "-c") == 0) {
 if(c_index != 0) {
   printerror(application.verbose_state.verbose, "Cannot use -c multiple times. Use the notation -c \"param1 param2\" instead.");
   return 0;
 }
 c_index = ++i;
      }else if(strcmp(argv[i], "-linenr") == 0) {
 show_linenumber = 1;
      }else {
 if(argv[i][0] == '-') {
   printerror(application.verbose_state.verbose, "Unknown option: %s\n\nRun pheader without command line arguments to show help", argv[i]);
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
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.verbose, "pheader: No files specified");
    return 0;
  }
  if(c_index != 0) {
    pickWordFromString(argv[c_index], 1, &nrwords, 0, ' ', application.verbose_state);
    for(i = 0; i < nrwords; i++) {
      sscanf(pickWordFromString(argv[c_index], i+1, &nrwords, 0, ' ', application.verbose_state), "%s", cmd);
      if(strcasecmp(cmd, "p0") == 0 || strcasecmp(cmd, "period") == 0) {
      }else if(strcasecmp(cmd, "freq") == 0 || strcasecmp(cmd, "cfreq") == 0) {
      }else if(strcasecmp(cmd, "freqlist") == 0) {
      }else if(strcasecmp(cmd, "reffreq") == 0) {
      }else if(strcasecmp(cmd, "npulses") == 0 || strcasecmp(cmd, "nsub") == 0 || strcasecmp(cmd, "nsubint") == 0 || strcasecmp(cmd, "nsubints") == 0 || strcasecmp(cmd, "nrsubint") == 0 || strcasecmp(cmd, "nrsubints") == 0 || strcasecmp(cmd, "nrpulses") == 0) {
      }else if(strcasecmp(cmd, "nbin") == 0 || strcasecmp(cmd, "nbins") == 0) {
      }else if(strcasecmp(cmd, "npol") == 0) {
      }else if(strcasecmp(cmd, "nchan") == 0 || strcasecmp(cmd, "nfreq") == 0) {
      }else if(strcasecmp(cmd, "nbits") == 0) {
      }else if(strcasecmp(cmd, "poltype") == 0) {
      }else if(strcasecmp(cmd, "dt") == 0 || strcasecmp(cmd, "tsamp") == 0) {
      }else if(strcasecmp(cmd, "length") == 0 || strcasecmp(cmd, "dur") == 0 || strcasecmp(cmd, "tobs") == 0) {
      }else if(strcasecmp(cmd, "length2") == 0) {
      }else if(strcasecmp(cmd, "tsub") == 0 || strcasecmp(cmd, "tsubint") == 0 || strcasecmp(cmd, "t_sub") == 0) {
      }else if(strcasecmp(cmd, "name") == 0) {
      }else if(strcasecmp(cmd, "bw") == 0) {
      }else if(strcasecmp(cmd, "chbw") == 0 || strcasecmp(cmd, "chanbw") == 0) {
      }else if(strcasecmp(cmd, "dm") == 0) {
      }else if(strcasecmp(cmd, "rm") == 0) {
      }else if(strcasecmp(cmd, "ra") == 0) {
      }else if(strcasecmp(cmd, "dec") == 0) {
      }else if(strcasecmp(cmd, "mjd") == 0) {
      }else if(strcasecmp(cmd, "format") == 0) {
      }else if(strcasecmp(cmd, "hist") == 0) {
      }else if(strcasecmp(cmd, "weights") == 0) {
 noweights = 0;
      }else if(strcasecmp(cmd, "observatory") == 0) {
      }else if(strcasecmp(cmd, "gentype") == 0) {
      }else if(strcasecmp(cmd, "long") == 0) {
      }else if(strcasecmp(cmd, "lat") == 0) {
      }else if(strcasecmp(cmd, "backend") == 0) {
      }else if(strcasecmp(cmd, "scanid") == 0) {
      }else if(strcasecmp(cmd, "observer") == 0 || strcasecmp(cmd, "observers") == 0) {
      }else if(strcasecmp(cmd, "project") == 0 || strcasecmp(cmd,"projectid") == 0 || strcasecmp(cmd,"projid") == 0) {
      }else if(strcasecmp(cmd, "parang") == 0) {
      }else if(strcasecmp(cmd, "parang1") == 0) {
      }else if(strcasecmp(cmd, "dedisp") == 0) {
      }else if(strcasecmp(cmd, "defarad") == 0) {
      }else if(strcasecmp(cmd, "depar") == 0) {
      }else {
 printerror(application.verbose_state.verbose, "\npheader: Unknown header parameter (%s), specify -H for a list", cmd);
 return 0;
      }
    }
  }
  datain = malloc(numberInApplicationFilenameList(&application, argv, application.verbose_state)*sizeof(datafile_definition));
  if(datain == NULL) {
    printerror(application.verbose_state.verbose, "pheader: Cannot allocate memory");
    return 0;
  }
  iformat_initial_value = application.iformat;
  i = 0;
  while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
    application.iformat = iformat_initial_value;
    if(application.iformat <= 0) {
      application.iformat = guessPSRData_format(filename_ptr, 0, application.verbose_state);
      if(application.iformat == -2 || application.iformat == -3)
 return 0;
    }
    if(isValidPSRDATA_format(application.iformat) == 0) {
      printerror(application.verbose_state.debug, "ERROR pheader: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", filename_ptr);
      free(datain);
      return 0;
    }
    if(openPSRData(&datain[i], filename_ptr, application.iformat, 0, 0, 0, application.verbose_state) == 0) {
      printerror(application.verbose_state.verbose, "pheader: Error opening data");
      return 0;
    }
    verbose_definition verbose2;
    copyVerboseState(application.verbose_state, &verbose2);
    if(c_index == 0)
      verbose2.verbose = 1;
    if(readHeaderPSRData(&datain[i], noweights, 0, verbose2) == 0) {
      printerror(application.verbose_state.verbose, "pheader: Error reading header");
      return 0;
    }
    if(c_index) {
      verbose_definition noverbose;
      cleanVerboseState(&noverbose);
      didhistory = 0;
      pickWordFromString(argv[c_index], 1, &nrwords, 0, ' ', application.verbose_state);
      for(j = 0; j < nrwords; j++) {
 sscanf(pickWordFromString(argv[c_index], j+1, &nrwords, 0, ' ', application.verbose_state), "%s", cmd);
 if(strcasecmp(cmd, "hist") == 0) {
   printf("History for %s\n", filename_ptr);
   showHistory(datain[i], noverbose);
   didhistory = 1;
 }
      }
      didweights = 0;
      didfreqlist = 0;
      pickWordFromString(argv[c_index], 1, &nrwords, 0, ' ', application.verbose_state);
      for(j = 0; j < nrwords; j++) {
 sscanf(pickWordFromString(argv[c_index], j+1, &nrwords, 0, ' ', application.verbose_state), "%s", cmd);
 if(strcasecmp(cmd, "weights") == 0) {
   printf("Weights for %s\n", filename_ptr);
   if(datain[i].weights == NULL) {
     printf("  no weights defined in file\n");
   }else {
     long nsub, nfreq;
     for(nsub = 0; nsub < datain[i].NrSubints; nsub++) {
       for(nfreq = 0; nfreq < datain[i].NrFreqChan; nfreq++) {
  printf("  subint %04ld channel %04ld = %lf MHz: %f\n", nsub, nfreq, get_weighted_channel_freq(datain[i], nsub, nfreq, application.verbose_state), datain[i].weights[nsub*datain[i].NrFreqChan+nfreq]);
       }
     }
   }
   didweights = 1;
 }else if(strcasecmp(cmd, "freqlist") == 0) {
   printf("Frequencies for %s\n", filename_ptr);
   long nsub, nfreq;
   for(nsub = 0; nsub < datain[i].NrSubints; nsub++) {
     for(nfreq = 0; nfreq < datain[i].NrFreqChan; nfreq++) {
       printf("  subint %04ld channel %04ld = %lf MHz\n", nsub, nfreq, get_weighted_channel_freq(datain[i], nsub, nfreq, application.verbose_state));
     }
   }
   didfreqlist = 1;
 }
      }
      pickWordFromString(argv[c_index], 1, &nrwords, 0, ' ', application.verbose_state);
      for(j = 0; j < nrwords; j++) {
 sscanf(pickWordFromString(argv[c_index], j+1, &nrwords, 0, ' ', application.verbose_state), "%s", cmd);
 if(strcasecmp(cmd, "tsub") == 0 || strcasecmp(cmd, "tsubint") == 0 || strcasecmp(cmd, "t_sub") == 0) {
   printf("tsub for %s: ", filename_ptr);
   long nsub;
   for(nsub = 0; nsub < datain[i].NrSubints; nsub++) {
     if(nsub != 0)
       printf(",");
     printf("%.*lf", 1+precision, get_tsub(datain[i], nsub, application.verbose_state));
   }
   printf(" sec\n");
 }
      }
    }
    closePSRData(&datain[i], 1, application.verbose_state);
    i++;
  }
  if(c_index == 0) {
    for(i = 0; i < numberInApplicationFilenameList(&application, argv, application.verbose_state); i++) {
      closePSRData(&datain[i], 0, application.verbose_state);
    }
    free(datain);
    terminateApplication(&application);
    return 0;
  }
  if(nrwords == didhistory + didweights + didfreqlist) {
    for(i = 0; i < numberInApplicationFilenameList(&application, argv, application.verbose_state); i++) {
      closePSRData(&datain[i], 0, application.verbose_state);
    }
    free(datain);
    terminateApplication(&application);
    return 0;
  }
  if(c_index) {
    pickWordFromString(argv[c_index], 1, &nrwords, 0, ' ', application.verbose_state);
    if(nrwords == 0) {
      printerror(application.verbose_state.verbose, "pheader: No parameters in -c option?");
      return 0;
    }
  }else {
    printerror(application.verbose_state.verbose, "pheader: Nothing to do, specify -c");
    return 0;
  }
  maxfilenamelength = strlen("filename");
  for(index = 0; index < numberInApplicationFilenameList(&application, argv, application.verbose_state); index++) {
    if(strlen(datain[index].filename) > maxfilenamelength)
      maxfilenamelength = strlen(datain[index].filename);
  }
  maxobservatorylength = strlen("observatory");
  for(index = 0; index < numberInApplicationFilenameList(&application, argv, application.verbose_state); index++) {
    if(strlen(datain[index].observatory) > maxobservatorylength)
      maxobservatorylength = strlen(datain[index].observatory);
  }
  maxgentypelength = strlen("gentype");
  for(index = 0; index < numberInApplicationFilenameList(&application, argv, application.verbose_state); index++) {
    if(strlen(returnGenType_str(datain[index].gentype)) > maxgentypelength)
      maxgentypelength = strlen(returnGenType_str(datain[index].gentype));
  }
  maxfileformatlength = strlen("form.");
  for(index = 0; index < numberInApplicationFilenameList(&application, argv, application.verbose_state); index++) {
    if(strlen(returnFileFormat_str(datain[index].format)) > maxfileformatlength)
      maxfileformatlength = strlen(returnFileFormat_str(datain[index].format));
  }
  maxinstrumentlength = strlen("backend");
  for(index = 0; index < numberInApplicationFilenameList(&application, argv, application.verbose_state); index++) {
    if(strlen(datain[index].instrument) > maxinstrumentlength)
      maxinstrumentlength = strlen(datain[index].instrument);
  }
  maxscanidlength = strlen("scanid");
  for(index = 0; index < numberInApplicationFilenameList(&application, argv, application.verbose_state); index++) {
    if(strlen(datain[index].scanID) > maxscanidlength)
      maxscanidlength = strlen(datain[index].scanID);
  }
  maxobserverlength = strlen("observer");
  for(index = 0; index < numberInApplicationFilenameList(&application, argv, application.verbose_state); index++) {
    if(strlen(datain[index].observer) > maxobserverlength)
      maxobserverlength = strlen(datain[index].observer);
  }
  maxprojectidlength = strlen("projectID");
  for(index = 0; index < numberInApplicationFilenameList(&application, argv, application.verbose_state); index++) {
    if(strlen(datain[index].projectID) > maxprojectidlength)
      maxprojectidlength = strlen(datain[index].projectID);
  }
  if(nohead == 0) {
    if(show_linenumber)
      printf("line ");
    printf("%s", "filename");
    for(j = 0; j < maxfilenamelength - strlen("filename"); j++)
      printf(" ");
    for(i = 0; i < nrwords; i++) {
      sscanf(pickWordFromString(argv[c_index], i+1, &nrwords, 0, ' ', application.verbose_state), "%s", cmd);
      int extra_precision;
      extra_precision = 0;
      if(strcasecmp(cmd, "p0") == 0 || strcasecmp(cmd, "period") == 0) {
 printf("  period");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "freq") == 0 || strcasecmp(cmd, "cfreq") == 0) {
 printf(" Frequency");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "reffreq") == 0) {
 printf("  Ref freq");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "npulses") == 0 || strcasecmp(cmd, "nsub") == 0 || strcasecmp(cmd, "nsubint") == 0 || strcasecmp(cmd, "nsubints") == 0 || strcasecmp(cmd, "nrsubint") == 0 || strcasecmp(cmd, "nrsubints") == 0 || strcasecmp(cmd, "nrpulses") == 0) {
 printf("      nsub");
      }else if(strcasecmp(cmd, "nbin") == 0 || strcasecmp(cmd, "nbins") == 0) {
 printf("      nbin");
      }else if(strcasecmp(cmd, "npol") == 0) {
 printf(" npol");
      }else if(strcasecmp(cmd, "nchan") == 0 || strcasecmp(cmd, "nfreq") == 0) {
 printf(" nchan");
      }else if(strcasecmp(cmd, "nbits") == 0) {
 printf(" nbits");
      }else if(strcasecmp(cmd, "poltype") == 0) {
 printf(" poltype");
      }else if(strcasecmp(cmd, "dt") == 0 || strcasecmp(cmd, "tsamp") == 0) {
 printf(" samp time");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "length") == 0 || strcasecmp(cmd, "dur") == 0 || strcasecmp(cmd, "tobs") == 0) {
 printf("   length");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "length2") == 0) {
 printf("  length2");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "tsub") == 0 || strcasecmp(cmd, "tsubint") == 0 || strcasecmp(cmd, "t_sub") == 0) {
 printf("  <tsub>");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "name") == 0) {
 printf("           name");
      }else if(strcasecmp(cmd, "bw") == 0) {
 printf("   bandw");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "chbw") == 0 || strcasecmp(cmd, "chanbw") == 0) {
 printf("   chan bw");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "dm") == 0) {
 printf("       DM");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "rm") == 0) {
 printf("       RM");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "ra") == 0) {
 printf("       RA");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "dec") == 0) {
 printf("     DEC");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "long") == 0) {
 printf("    long");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "lat") == 0) {
 printf("     lat");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "mjd") == 0) {
 printf("       MJD");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "dedisp") == 0) {
 printf(" dedisp");
      }else if(strcasecmp(cmd, "defarad") == 0) {
 printf(" defarad");
      }else if(strcasecmp(cmd, "depar") == 0) {
 printf(" depar");
      }else if(strcasecmp(cmd, "hist") == 0) {
      }else if(strcasecmp(cmd, "weights") == 0) {
      }else if(strcasecmp(cmd, "freqlist") == 0) {
      }else if(strcasecmp(cmd, "observatory") == 0) {
 printf(" observatory");
 for(j = 0; j < maxobservatorylength - strlen("observatory"); j++)
   printf(" ");
      }else if(strcasecmp(cmd, "gentype") == 0) {
 printf(" gentype");
 for(j = 0; j < maxgentypelength - strlen("gentype"); j++)
   printf(" ");
      }else if(strcasecmp(cmd, "format") == 0) {
 printf(" form.");
 for(j = 0; j < maxfileformatlength - strlen("form."); j++)
   printf(" ");
      }else if(strcasecmp(cmd, "backend") == 0) {
 printf(" backend");
 for(j = 0; j < maxinstrumentlength - strlen("backend"); j++)
   printf(" ");
      }else if(strcasecmp(cmd, "scanid") == 0) {
 printf(" scanid");
 for(j = 0; j < maxscanidlength - strlen("scanid"); j++)
   printf(" ");
      }else if(strcasecmp(cmd, "observer") == 0 || strcasecmp(cmd, "observers") == 0) {
 printf(" observer");
 for(j = 0; j < maxobserverlength - strlen("observer"); j++)
   printf(" ");
      }else if(strcasecmp(cmd, "project") == 0 || strcasecmp(cmd,"projectid") == 0 || strcasecmp(cmd,"projid") == 0) {
 printf(" projectID");
 for(j = 0; j < maxprojectidlength - strlen("projectID"); j++)
   printf(" ");
      }else if(strcasecmp(cmd, "parang") == 0) {
 printf(" parallactic angle (mid)");
 extra_precision = 1;
      }else if(strcasecmp(cmd, "parang1") == 0) {
 printf(" parallactic angle (start)");
 extra_precision = 1;
      }else {
 printerror(application.verbose_state.verbose, "\npheader: Unknown header parameter (%s), specify -H for a list", cmd);
 return 0;
      }
      if(extra_precision) {
 if(precision > 0) {
   for(j = 0; j < precision; j++) {
     printf(" ");
   }
 }
      }
    }
    printf("\n");
  }
  for(index = 0; index < numberInApplicationFilenameList(&application, argv, application.verbose_state); index++) {
    if(show_linenumber)
      printf("%-4d ", index+1);
    printf("%s", datain[index].filename);
    for(j = 0; j < maxfilenamelength - strlen(datain[index].filename); j++)
      printf(" ");
    if(c_index) {
      for(i = 0; i < nrwords; i++) {
 sscanf(pickWordFromString(argv[c_index], i+1, &nrwords, 0, ' ', application.verbose_state), "%s", cmd);
 if(strcasecmp(cmd, "p0") == 0 || strcasecmp(cmd, "period") == 0) {
   double period;
   int ret;
   ret = get_period(datain[index], 0, &period, application.verbose_state);
   if(ret == 2) {
     printerror(application.verbose_state.debug, "ERROR pheader (%s): Cannot obtain period", datain[index].filename);
     return 0;
   }
   if(period > 0) {
     printf(" %*.*lf", 7+precision, 4+precision, period);
   }else {
     printf(" SEARCH?");
     footnote_search = 1;
   }
 }else if(strcasecmp(cmd, "freq") == 0 || strcasecmp(cmd, "cfreq") == 0) {
   printf(" %*.*lf", 9+precision, 4+precision, get_centre_frequency(datain[index], application.verbose_state));
 }else if(strcasecmp(cmd, "reffreq") == 0) {
   printf(" %*.*lf", 9+precision, 4+precision, datain[index].freq_ref);
 }else if(strcasecmp(cmd, "npulses") == 0 || strcasecmp(cmd, "nsub") == 0 || strcasecmp(cmd, "nsubint") == 0 || strcasecmp(cmd, "nsubints") == 0 || strcasecmp(cmd, "nrsubint") == 0 || strcasecmp(cmd, "nrsubints") == 0 || strcasecmp(cmd, "nrpulses") == 0) {
   printf(" %9ld", datain[index].NrSubints);
 }else if(strcasecmp(cmd, "nbin") == 0 || strcasecmp(cmd, "nbins") == 0) {
   printf(" %9ld", datain[index].NrBins);
 }else if(strcasecmp(cmd, "npol") == 0) {
   printf(" %4ld", datain[index].NrPols);
 }else if(strcasecmp(cmd, "nchan") == 0 || strcasecmp(cmd, "nfreq") == 0) {
   printf(" %5ld", datain[index].NrFreqChan);
 }else if(strcasecmp(cmd, "nbits") == 0) {
   printf(" %5d", datain[index].NrBits);
 }else if(strcasecmp(cmd, "poltype") == 0) {
   printf(" %7d", datain[index].poltype);
 }else if(strcasecmp(cmd, "dt") == 0 || strcasecmp(cmd, "tsamp") == 0) {
   printf(" %*.*lf", 9+precision, 6+precision, get_tsamp(datain[index], 0, application.verbose_state));
 }else if(strcasecmp(cmd, "bw") == 0) {
   printf(" %*.*lf", 7+precision, 1+precision, get_bandwidth(datain[index], application.verbose_state));
 }else if(strcasecmp(cmd, "chbw") == 0 || strcasecmp(cmd, "chanbw") == 0) {
   double chanbw;
   if(get_channelbandwidth(datain[index], &chanbw, application.verbose_state) == 0) {
     printerror(application.verbose_state.debug, "ERROR pheader (%s): Cannot obtain channel bandwidth.", datain[index].filename);
     return 0;
   }
   printf(" %*.*lf", 9+precision, 3+precision, chanbw);
 }else if(strcasecmp(cmd, "dm") == 0) {
   printf(" %*.*lf", 8+precision, 3+precision, datain[index].dm);
 }else if(strcasecmp(cmd, "rm") == 0) {
   printf(" %*.*lf", 8+precision, 3+precision, datain[index].rm);
 }else if(strcasecmp(cmd, "ra") == 0) {
   printf(" %*.*lf", 8+precision, 3+precision, datain[index].ra);
 }else if(strcasecmp(cmd, "dec") == 0) {
   printf(" %*.*lf", 7+precision, 3+precision, datain[index].dec);
 }else if(strcasecmp(cmd, "long") == 0) {
   printf(" %*.*lf", 7+precision, 3+precision, observatory_long_geodetic(datain[index])*180.0/M_PI);
 }else if(strcasecmp(cmd, "lat") == 0) {
   printf(" %*.*lf", 7+precision, 3+precision, observatory_lat_geodetic(datain[index])*180.0/M_PI);
 }else if(strcasecmp(cmd, "mjd") == 0) {
   printf(" %*.*Lf", 9+precision, 3+precision, datain[index].mjd_start);
 }else if(strcasecmp(cmd, "length") == 0 || strcasecmp(cmd, "dur") == 0 || strcasecmp(cmd, "tobs") == 0) {
   printf(" %*.*lf", 8+precision, 1+precision, get_tobs(datain[index], application.verbose_state));
   if(get_tobs(datain[index], application.verbose_state) < 0.001)
     footnote_length = 1;
 }else if(strcasecmp(cmd, "length2") == 0) {
   printf(" %*.*lf", 8+precision, 1+precision, get_tsamp(datain[index], 0, application.verbose_state)*datain[index].NrSubints*datain[index].NrBins);
   footnote_length2 = 1;
 }else if(strcasecmp(cmd, "tsub") == 0 || strcasecmp(cmd, "tsubint") == 0 || strcasecmp(cmd, "t_sub") == 0) {
   printf(" %*.*lf", 7+precision, 1+precision, get_tobs(datain[index], application.verbose_state)/(double)datain[index].NrSubints);
 }else if(strcasecmp(cmd, "format") == 0) {
   printf(" %s", returnFileFormat_str(datain[index].format));
   for(j = 0; j < maxfileformatlength - strlen(returnFileFormat_str(datain[index].format)); j++)
     printf(" ");
 }else if(strcasecmp(cmd, "name") == 0) {
   printf(" %14s", datain[index].psrname);
 }else if(strcasecmp(cmd, "observatory") == 0) {
   printf(" %s", datain[index].observatory);
   for(j = 0; j < maxobservatorylength - strlen(datain[index].observatory); j++)
     printf(" ");
 }else if(strcasecmp(cmd, "gentype") == 0) {
   printf(" %s", returnGenType_str(datain[index].gentype));
   for(j = 0; j < maxgentypelength - strlen(returnGenType_str(datain[index].gentype)); j++)
     printf(" ");
 }else if(strcasecmp(cmd, "backend") == 0) {
   printf(" %s", datain[index].instrument);
   for(j = 0; j < maxinstrumentlength - strlen(datain[index].instrument); j++)
     printf(" ");
 }else if(strcasecmp(cmd, "scanid") == 0) {
   printf(" %s", datain[index].scanID);
   for(j = 0; j < maxscanidlength - strlen(datain[index].scanID); j++)
     printf(" ");
 }else if(strcasecmp(cmd, "observer") == 0 || strcasecmp(cmd, "observers") == 0) {
   printf(" %s", datain[index].observer);
   for(j = 0; j < maxobserverlength - strlen(datain[index].observer); j++)
     printf(" ");
 }else if(strcasecmp(cmd, "project") == 0 || strcasecmp(cmd,"projectid") == 0 || strcasecmp(cmd,"projid") == 0) {
   printf(" %s", datain[index].projectID);
   for(j = 0; j < maxprojectidlength - strlen(datain[index].projectID); j++)
     printf(" ");
 }else if(strcasecmp(cmd, "parang") == 0) {
   double parang;
   if(data_parang(datain[index], -1, &parang, application.verbose_state)) {
     printf(" %*.*lf", 23+precision, 2+precision, parang*180.0/M_PI);
   }else
     printf("                *");
   footnote_parang = 1;
 }else if(strcasecmp(cmd, "parang1") == 0) {
   double parang;
   if(data_parang(datain[index], 0, &parang, application.verbose_state)) {
     printf(" %*.*lf", 25+precision, 2+precision, parang*180.0/M_PI);
   }else
     printf("                *");
   footnote_parang = 1;
 }else if(strcasecmp(cmd, "dedisp") == 0) {
   printf(" %6d", datain[index].isDeDisp);
 }else if(strcasecmp(cmd, "defarad") == 0) {
   printf(" %7d", datain[index].isDeFarad);
 }else if(strcasecmp(cmd, "depar") == 0) {
   printf(" %5d", datain[index].isDePar);
 }else if(strcasecmp(cmd, "hist") == 0) {
 }else if(strcasecmp(cmd, "weights") == 0) {
 }else if(strcasecmp(cmd, "freqlist") == 0) {
 }else {
   printerror(application.verbose_state.verbose, "\nERROR pheader: Unknown header parameter (%s), specify -H for a list", cmd);
   return 0;
 }
      }
    }
    printf("\n");
  }
  if((footnote_length || footnote_length2 || footnote_search || footnote_parang) && showfootnotes) {
    printf("\nFootnotes:\n");
    if(footnote_length) {
      printf("- length: Observation duration is set to zero if undefined. Parameter length2 will give a lower limit for the observation duration.\n");
    }
    if(footnote_length2) {
      printf("- length2 = nbin*nsub*tsamp <= duration of the observation for folded data\n");
    }
    if(footnote_search) {
      printf("- period: is not set for search-mode data\n");
    }
    if(footnote_parang) {
      printf("- parallactic angle: This is a derived quantity\n");
    }
  }
  for(i = 0; i < numberInApplicationFilenameList(&application, argv, application.verbose_state); i++) {
    closePSRData(&datain[i], 0, application.verbose_state);
  }
  free(datain);
  terminateApplication(&application);
  return 0;
}

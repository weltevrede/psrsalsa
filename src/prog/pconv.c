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
#include <stdio.h>
#include <string.h>
#include <psrsalsa.h>
int pumawrite(void *prptr, int size, int nelem,FILE *out);
void print_help();
int main(int argc, char **argv)
{
  datafile_definition fin, fout;
  int read_wholefile, indx;
  long i, n, p, f, n1, n2;
  float *pulseData, sample;
  char outputname[MaxFilenameLength], *dummy_ptr;
  psrsalsaApplication application;
  initApplication(&application, "pconv", "[options] inputfile(s)");
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.switch_ext = 1;
  application.switch_output = 1;
  application.switch_fchan = 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_formatlist = 1;
  application.switch_tscr = 1;
  application.switch_tscr_complete = 1;
  application.switch_TSCR = 1;
  application.switch_FSCR = 1;
  application.switch_fscr = 1;
  application.switch_rebin = 1;
  application.switch_scale = 1;
  application.switch_polselect = 1;
  application.switch_nskip = 1;
  application.switch_nread = 1;
  application.switch_noweights = 1;
  application.switch_useweights = 1;
  application.switch_uniformweights = 1;
  application.switch_norm = 1;
  application.switch_debase = 1;
  application.switch_onpulse = 1;
  application.switch_onpulsef = 1;
  application.switch_nocounters = 1;
  application.switch_stokes = 1;
  application.switch_deparang = 1;
  application.switch_history_cmd_only = 1;
  read_wholefile = 1;
  if(argc <= 1) {
    printApplicationHelp(&application);
    print_help();
    terminateApplication(&application);
    return 0;
  }
  for(indx = 1; indx < argc; indx++) {
    if(processCommandLine(&application, argc, argv, &indx)) {
    }else if(strcmp(argv[indx], "-memsave") == 0) {
      read_wholefile = 0;
    }else {
      if(argv[indx][0] == '-') {
 printerror(application.verbose_state.debug, "ERROR pconv: Unknown option (%s).\n\nRun pconv without command line arguments to show help", argv[indx]);
 terminateApplication(&application);
 return 0;
      }else {
 if(applicationAddFilename(indx, application.verbose_state) == 0)
   return 0;
      }
    }
  }
  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pconv: No files specified");
    return 0;
  }
  if(isValidPSRDATA_format(application.oformat) == 0) {
    printerror(application.verbose_state.debug, "ERROR pconv: Please specify a valid output format with the -oformat option.");
    terminateApplication(&application);
    return 0;
  }
  cleanPSRData(&fout, application.verbose_state);
  while((dummy_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
    if(application.iformat <= 0) {
      application.iformat = guessPSRData_format(dummy_ptr, 0, application.verbose_state);
      if(application.iformat == -2 || application.iformat == -3) {
 closePSRData(&fout, 0, application.verbose_state);
 terminateApplication(&application);
 return 0;
      }
    }
    if(isValidPSRDATA_format(application.iformat) == 0) {
      printerror(application.verbose_state.debug, "ERROR pconv: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", dummy_ptr);
      closePSRData(&fout, 0, application.verbose_state);
      terminateApplication(&application);
      return 0;
    }
  if(!openPSRData(&fin, dummy_ptr, application.iformat, 0, read_wholefile, 0, application.verbose_state))
    return 0;
  if(read_wholefile == 0) {
    if(application.fchan_select != -1) {
      printerror(application.verbose_state.debug, "ERROR pconv: -fchan option doesn't work with -memsave option.");
      return 0;
    }
    if(!readHeaderPSRData(&fin, 0, 0, application.verbose_state))
      return 0;
  }
  if(PSRDataHeader_parse_commandline(&fin, argc, argv, application.verbose_state) == 0)
    return 0;
  region_frac_to_int(&(application.onpulse), fin.NrBins, 0);
  if(read_wholefile != 0) {
    for(i = 1; i < argc; i++) {
      if(strcmp(argv[i], "-header") == 0) {
 printwarning(application.verbose_state.debug, "WARNING pconv: If using the -header option, be aware it applied BEFORE the preprocessing.");
 break;
      }
    }
    if(preprocessApplication(&application, &fin) == 0) {
      return 0;
    }
  }
  copy_params_PSRData(fin, &fout, application.verbose_state);
  double period;
  int ret;
  if(fin.isFolded) {
    ret = get_period(fin, 0, &period, application.verbose_state);
    if(ret == 2) {
      printerror(application.verbose_state.debug, "ERROR pconv (%s): Cannot obtain period", fin.filename);
      return 0;
    }
  }else {
    ret = 1;
    period = -1;
  }
  if(ret == 0 && period >= 0 && period < 0.001) {
    fout.fixedPeriod = 1;
    fout.foldMode = FOLDMODE_FIXEDPERIOD;
  }
  if(getOutputName(&application, dummy_ptr, outputname, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pconv: Changing filename failed");
    return 0;
  }
  n1 = application.nskip;
  if(application.nread > 0)
    n2 = application.nskip + application.nread;
  else
    n2 = fin.NrSubints;
  fout.NrSubints = n2 - n1;
  if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
    return 0;
  if(!writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state))
    return 0;
  if(read_wholefile == 1) {
    if(writePSRData(&fout, fin.data, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR pconv: Cannot write data");
      return 0;
    }
  }else {
      if(application.iformat == SIGPROC_format && application.oformat == PUMA_format) {
      long j;
      unsigned char *sample_b;
      float *sample_f;
      void *data_ptr;
      if(fin.NrSubints > 1) {
 printerror(application.verbose_state.debug, "ERROR pconv: Non-folded sigproc data can be converted in PuMa format with the -memsave option. This data has %ld subints.", fin.NrSubints);
 return 0;
      }
      if(fin.NrBits == 32) {
 data_ptr = malloc(4*fin.NrFreqChan);
 sample_f = data_ptr;
      }else if(fin.NrBits == 8) {
 data_ptr = malloc(fin.NrFreqChan);
 sample_b = data_ptr;
      }else {
 printerror(application.verbose_state.debug, "ERROR pconv: Only 32-bit or 8-bit sigproc data can be converted in PuMa format. This is %ld bit data.", fin.NrBits);
 return 0;
      }
      if(data_ptr == NULL) {
 printerror(application.verbose_state.debug, "ERROR pconv: Cannot allocate memory.");
 return 0;
      }
      long long filepos;
      long polarization = 0;
      long pulsenr = 0;
      for(i = 0; i < fin.NrBins; i++) {
 for(j = 0; j < fin.NrFreqChan; j++) {
   if(fin.NrBits == 32) {
     fread(&sample_f[j], sizeof(float), 1, fin.fptr);
   }else if(fin.NrBits == 8) {
     fread(&sample_b[j], sizeof(unsigned char), 1, fin.fptr);
   }
 }
 filepos = polarization*fin.NrFreqChan*fin.NrSubints*fin.NrBins;
 filepos += pulsenr*fin.NrBins;
 filepos += i;
 filepos *= sizeof(float);
 filepos += fout.datastart;
 for(j = 0; j < fin.NrFreqChan; j++) {
   fseeko(fout.fptr, filepos, SEEK_SET);
   if(fin.NrBits == 32) {
     sample = sample_f[j];
   }else if(fin.NrBits == 8) {
     sample = sample_b[j];
   }
   pumawrite(&sample, sizeof(float), 1, fout.fptr);
   filepos += fin.NrSubints*fin.NrBins*sizeof(float);
 }
      }
      free(data_ptr);
    }else if(application.iformat == PSRCHIVE_ASCII_format && application.oformat == PUMA_format) {
      char txt[100];
      for(n = 0; n < fin.NrSubints; n++) {
 if(application.verbose_state.verbose) printf("pulse %ld/%ld\r", n+1-n1, fout.NrSubints);
 for(f = 0; f < fin.NrFreqChan; f++) {
   for(i = 0; i < fin.NrBins; i++) {
     fscanf(fin.fptr, "%s", txt);
     if(i == 0 && f == 0 && n == 0) {
       if(strcmp(txt, "0") != 0) {
  printerror(application.verbose_state.debug, "ERROR pconv: PSRCHIVE ascii header format has probably been changed.");
  return 0;
       }
     }
     fscanf(fin.fptr, "%s", txt);
     if(i == 0 && f == 0 && n == 0) {
       if(strcmp(txt, "0") != 0) {
  printerror(application.verbose_state.debug, "ERROR pconv: PSRCHIVE ascii header format has probably been changed.");
  return 0;
       }
     }
     fscanf(fin.fptr, "%s", txt);
     if(i == 0 && f == 0 && n == 0) {
       if(strcmp(txt, "0") != 0) {
  printerror(application.verbose_state.debug, "ERROR pconv: PSRCHIVE ascii header format has probably been changed.");
  return 0;
       }
     }
     for(p = 0; p < fin.NrPols; p++) {
       fscanf(fin.fptr, "%f", &sample);
       if(n >= n1 && n < n2) {
  if(writePulsePSRData(&fout, n-n1, p, f, i, 1, &sample, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pconv: Cannot write data");
    return 0;
  }
       }
     }
   }
 }
      }
    }else {
      pulseData = (float *)malloc(fin.NrBins*sizeof(float));
      if(pulseData == NULL) {
 printerror(application.verbose_state.debug, "ERROR pconv: Cannot allocate memory.");
 return 0;
      }
      for(n = n1; n < n2; n++) {
 if(application.verbose_state.verbose) {
   printf("pulse %ld/%ld\r", n+1-n1, fout.NrSubints);
   fflush(stdout);
 }
 for(p = 0; p < fin.NrPols; p++) {
   for(f = 0; f < fin.NrFreqChan; f++) {
     if(readPulsePSRData(&fin, n, p, f, 0, fin.NrBins, pulseData, application.verbose_state) == 0) {
       printerror(application.verbose_state.debug, "ERROR pconv: Cannot read individual pulses from data.");
       return 0;
     }
     if(writePulsePSRData(&fout, n-n1, p, f, 0, fin.NrBins, pulseData, application.verbose_state) == 0) {
       printerror(application.verbose_state.debug, "ERROR pconv: Cannot write individual pulses.");
       return 0;
     }
   }
 }
      }
      free(pulseData);
    }
  }
  closePSRData(&fin, 0, application.verbose_state);
  closePSRData(&fout, 0, application.verbose_state);
  }
  terminateApplication(&application);
  return 0;
}
void print_help()
{
  fprintf(stdout, "Other options:\n");
  fprintf(stdout, "  -memsave     Try to use less memory. Not all conversions will work with\n               or without this switch.\n");
  fprintf(stdout, "\n");
  printf("\n");
  printCitationInfo();
}

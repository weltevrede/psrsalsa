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
#include <string.h>
#include <math.h>
#include "psrsalsa.h"







int readSigprocHeader_readParamID(FILE *fin, char **id, verbose_definition verbose)
{
  int idlength;
  if(fread(&idlength, sizeof(int), 1, fin) != 1) {
    printerror(verbose.debug,"ERROR readSigprocHeader_readParamID: Read error");
    return 0;
  }
  if(idlength < 1 || idlength > 80) {
    printerror(verbose.debug,"ERROR readSigprocHeader_readParamID: ID length (%d) makes no sense", idlength);
    return 0;
  }
  if(*id != NULL) {
    free(*id);
  }
  *id = malloc(idlength+1);
  if(*id == NULL) {
    printerror(verbose.debug,"ERROR readSigprocHeader_readParamID: Memory allocation error");
    return 0;
  }
  if(fread(*id, 1, idlength, fin) != idlength) {
    printerror(verbose.debug,"ERROR readSigprocHeader_readParamID: Read error");
    return 0;
  }
  (*id)[idlength] = 0;
  return 1;
}

int readSigprocHeader(datafile_definition *datafile, verbose_definition verbose)
{
  char *id;
  int dummy_int, nifs;
  double freq_chan1, freq_chanbw, dummy_double;
  long int dummy_long;


  datafile->gentype = GENTYPE_SEARCHMODE;
  datafile->isFolded = 0;
  datafile->foldMode = FOLDMODE_UNKNOWN;
  datafile->fixedPeriod = -1;
  datafile->tsubMode = TSUBMODE_FIXEDTSUB;
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    printerror(verbose.debug,"ERROR readSigprocHeader: Memory allocation error");
    exit(-1);
  }
  datafile->tsub_list[0] = 0;

  id = NULL;
  freq_chan1 = -1;
  freq_chanbw = 0;
  nifs = 1;

  if(readSigprocHeader_readParamID(datafile->fptr, &id, verbose) == 0) {
    printerror(verbose.debug, "ERROR readSigprocHeader: Cannot read ID");
    return 0;
  }
  if(verbose.debug) {
    printf("DEBUG: ID=%s\n", id);
  }
  if(strcmp(id, "HEADER_START") != 0) {
    printerror(verbose.debug, "ERROR readSigprocHeader: First ID is not the expected ID (%s != HEADER_START).", id);
    return 0;
  }

  do {
    if(readSigprocHeader_readParamID(datafile->fptr, &id, verbose) == 0) {
      printerror(verbose.debug, "ERROR readSigprocHeader: Cannot read ID");
      return 0;
    }
    if(verbose.debug) {
      printf("DEBUG: ID=%s\n", id);
    }
    if(strcmp(id, "machine_id") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      char txt[100];
      sprintf(txt, "%d", dummy_int);
      if(set_instrument_PSRData(datafile, txt, verbose) == 0) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Setting instrument failed");
 return 0;
      }
    }else if(strcmp(id, "telescope_id") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      char txt[100];
      sprintf(txt, "%d", dummy_int);
      if(set_observatory_PSRData(datafile, txt, verbose) == 0) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Setting observatory failed");
 return 0;
      }
    }else if(strcmp(id, "data_type") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      if(dummy_int == 3) {
 datafile->gentype = GENTYPE_SUBINTEGRATIONS;
 datafile->isFolded = 1;
      }
    }else if(strcmp(id, "fch1") == 0) {
      if(fread(&freq_chan1, sizeof(double), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%lf\n", id, freq_chan1);
      }
    }else if(strcmp(id, "foff") == 0) {
      if(fread(&freq_chanbw, sizeof(double), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%lf\n", id, freq_chanbw);
      }
    }else if(strcmp(id, "nchans") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      datafile->NrFreqChan = dummy_int;
    }else if(strcmp(id, "source_name") == 0) {
      char *txt;
      txt = NULL;
      if(readSigprocHeader_readParamID(datafile->fptr, &txt, verbose) == 0) {
 printerror(verbose.debug, "ERROR readSigprocHeader: Cannot read string");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%s\n", id, txt);
      }
      if(set_psrname_PSRData(datafile, txt, verbose) == 0) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Setting observatory failed");
 return 0;
      }
      free(txt);
    }else if(strcmp(id, "src_raj") == 0) {
      int hour, min;
      double sec;
      if(fread(&(datafile->ra), sizeof(double), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%lf\n", id, datafile->ra);
      }
      hour = datafile->ra/10000.0;
      datafile->ra -= hour*10000.0;
      min = datafile->ra/100.0;
      datafile->ra -= min*100.0;
      sec = datafile->ra;
      datafile->ra = 2.0*M_PI*((double)(hour+(double)(min+(double)sec/60.0)/60.0)/24.0);
    }else if(strcmp(id, "src_dej") == 0) {
      int deg, min, sign;
      double sec;
      if(fread(&(datafile->dec), sizeof(double), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%lf\n", id, datafile->dec);
      }
      if(datafile->dec < 0) {
 sign = -1;
 datafile->dec = -datafile->dec;
      }else {
 sign = 1;
      }
      deg = datafile->dec/10000.0;
      datafile->dec -= deg*10000.0;
      min = datafile->dec/100.0;
      datafile->dec -= min*100.0;
      sec = datafile->dec;
      datafile->dec = (double)sign*M_PI*((double)(deg+(double)(min+(double)sec/60.0)/60.0))/180.0;
    }else if(strcmp(id, "refdm") == 0) {
      if(fread(&dummy_double, sizeof(double), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%lf\n", id, dummy_double);
      }
      datafile-> dm = dummy_double;
    }else if(strcmp(id, "nbits") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      datafile->NrBits = dummy_int;
      if(datafile->NrBits != 32 && datafile->NrBits != 8) {
 printerror(verbose.debug, "ERROR readSigprocHeader: Can only handle 32-bit or 8-bit data. Got %d bit data.", datafile->NrBits);
 return 0;
      }
    }else if(strcmp(id, "nifs") == 0) {
      if(fread(&nifs, sizeof(int), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%d\n", id, nifs);
      }
    }else if(strcmp(id, "tstart") == 0) {
      if(fread(&dummy_double, sizeof(double), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%lf\n", id, dummy_double);
      }
      datafile->mjd_start = dummy_double;
    }else if(strcmp(id, "tsamp") == 0) {
      if(fread(&dummy_double, sizeof(double), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%lf\n", id, dummy_double);
      }
      datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
      datafile->fixedtsamp = dummy_double;
    }else if(strcmp(id, "npuls") == 0) {
      if(fread(&(dummy_long), sizeof(long int), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%ld\n", id, dummy_long);
      }



    }else if(strcmp(id, "nbins") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      datafile->NrBins = dummy_int;
    }else if(strcmp(id, "period") == 0) {
      if(fread(&dummy_double, sizeof(double), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%lf\n", id, dummy_double);
      }
      datafile->isFolded = 1;
      datafile->foldMode = FOLDMODE_FIXEDPERIOD;
      datafile->fixedPeriod = dummy_double;
    }else if(strcmp(id, "barycentric") == 0 || strcmp(id, "pulsarcentric") == 0 || strcmp(id, "nsamples") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%d\n", id, dummy_int);
      }
    }else if(strcmp(id, "az_start") == 0 || strcmp(id, "za_start") == 0) {
      if(fread(&dummy_double, sizeof(double), 1, datafile->fptr) != 1) {
 printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
 return 0;
      }
      if(verbose.debug) {
 printf("DEBUG: %s=%lf\n", id, dummy_double);
      }
    }else if(strcmp(id, "HEADER_END") == 0) {
    }else {
      printerror(verbose.debug,"ERROR readSigprocHeader: ID=%s is not recognized", id);
      return 0;
    }
  }while(strcmp(id, "HEADER_END") != 0);

  if(datafile->NrFreqChan <= 0) {

    printwarning(verbose.debug, "WARNING readSigprocHeader: The number of frequency channels do not appear to be defined. It is assumed only one frequency channel is present.");
    datafile->NrFreqChan = 1;
  }

  if(datafile->isFolded) {
    printwarning(verbose.debug, "WARNING readSigprocHeader: For folded data the full baseline is assumed to be stored.");
    datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
    if(datafile->fixedPeriod <= 0) {
      printwarning(verbose.debug, "WARNING readSigprocHeader: Period is not set in the header, assume it is 1 sec.");
      datafile->foldMode = FOLDMODE_FIXEDPERIOD;
      datafile->fixedPeriod = 1;
    }
    datafile->fixedtsamp = datafile->fixedPeriod/(double)datafile->NrBins;
  }

  if(freq_chan1 > 0) {
    datafile->freqMode = FREQMODE_UNIFORM;
    datafile->uniform_bw = datafile->NrFreqChan*freq_chanbw;
    datafile->uniform_freq_cent = freq_chan1;
    datafile->uniform_freq_cent += 0.5*datafile->uniform_bw;
  }else {
    printwarning(verbose.debug, "WARNING readSigprocHeader: The frequency of the observation does not appear to be defined.");
    if(datafile->freq_ref > 0) {

      printwarning(verbose.debug, "WARNING readSigprocHeader: It is assumed that the frequency of the observationdefines the reference frequency.");
      datafile->freqMode = FREQMODE_UNIFORM;
      datafile->uniform_freq_cent = datafile->freq_ref;
    }
  }


  off_t curpos, datasize;
  curpos = ftello(datafile->fptr);
  fseeko(datafile->fptr, 0, SEEK_END);
  datasize = ftello(datafile->fptr);
  datasize -= curpos;
  fseeko(datafile->fptr, curpos, SEEK_SET);
  if(verbose.debug) {
    printf("DEBUG: data size=%ld byte\n", datasize);
  }

  if(datafile->isFolded && datafile->NrBins <= 0) {
    printwarning(verbose.debug, "WARNING readSigprocHeader: For folded data the number of bins should be defined in header. Assume the data is not folded");
    datafile->isFolded = 0;
  }

  if(datafile->isFolded == 0) {
    long double tmp;
    tmp = datasize*8.0;
    tmp /= (long double)(datafile->NrBits);
    tmp /= (long double)nifs;
    tmp /= (long double)datafile->NrFreqChan;
    datafile->NrBins = tmp;
    if(verbose.debug) {
      printf("DEBUG: Implied number of bins = %ld\n", datafile->NrBins);
    }

    datafile->NrSubints = 1;
  }else {

    long double tmp;
    long double subintsize;
    subintsize = ((long double)(datafile->NrFreqChan)*(long double)(datafile->NrBins)*(long double)(datafile->NrBits)/8.0 + (long double)curpos);
    tmp = datasize + curpos;
    tmp /= subintsize;
    datafile->NrSubints = roundl(tmp);
  }


  datafile->tsub_list[0] = datafile->NrBins * datafile->fixedtsamp;
  printwarning(verbose.debug, "WARNING readSigprocHeader: Assuming there is only one polarization channel in the data");
  datafile->NrPols = 1;
  return 1;
}
int readSigprocfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long n, f, i, p;
  float *sample_f;
  unsigned char *sample_b;
  int ret;
  if(datafile.NrPols > 1) {
    printerror(verbose.debug, "readSigprocfile: Data should have just one polarization");
    return 0;
  }
  if(datafile.NrBits != 32 && datafile.NrBits != 8) {
    printerror(verbose.debug, "ERROR readSigprocfile: Can only handle 32-bit or 8-bit data. Got %d bit data.", datafile.NrBits);
    return 0;
  }
  if(datafile.NrBits == 32) {
    sample_f = malloc(datafile.NrFreqChan*sizeof(float));
    if(sample_f == NULL) {
      printerror(verbose.debug, "readSigprocfile: Memory allocation error");
      return 0;
    }
  }else if(datafile.NrBits == 8) {
    sample_b = malloc(datafile.NrFreqChan);
    if(sample_b == NULL) {
      printerror(verbose.debug, "readSigprocfile: Memory allocation error");
      return 0;
    }
  }
  for(n = 0; n < datafile.NrSubints; n++) {
    for(i = 0; i < datafile.NrBins; i++) {
      if(datafile.NrBits == 32) {
 ret = fread(sample_f, sizeof(float), datafile.NrFreqChan, datafile.fptr);
      }else if(datafile.NrBits == 8) {
 ret = fread(sample_b, 1, datafile.NrFreqChan, datafile.fptr);
      }
      if(ret != datafile.NrFreqChan) {
 printerror(verbose.debug, "ERROR readSigprocfile: Cannot read data (sample %ld of subint %ld).", i, n);
 return 0;
      }
      p = 0;
      for(f = 0; f < datafile.NrFreqChan; f++) {
 if(datafile.NrBits == 32) {
   data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] = sample_f[f];
 }else if(datafile.NrBits == 8) {
   data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] = sample_b[f];
 }
      }
    }
    if(n != datafile.NrSubints - 1) {
      fseeko(datafile.fptr, datafile.datastart, SEEK_CUR);
    }
  }
  if(datafile.NrBits == 32) {
    free(sample_f);
  }else if(datafile.NrBits == 8) {
    free(sample_b);
  }
  if(verbose.verbose) printf("Reading is done.                           \n");
  return 1;
}
int readSigprocASCIIHeader(datafile_definition *datafile, verbose_definition verbose)
{
  int j, c;
  double sec, freq;
  long double mjd;
  char tmp1[10000];
  char tmp2[1000];
  char tmp3[1000];
  char tmp4[1000];
  char *txtptr;
  long pos;
  j = fscanf(datafile->fptr_hdr, "%s %Lf %lf %lf %s %lf %lf %ld %s %s %s", tmp1, &mjd, &sec, &(datafile->fixedPeriod), tmp2, &freq, &(datafile->dm), &(datafile->NrBins), tmp4, tmp3, &(datafile->psrname[0]));
  if(j != 11) {
    printerror(verbose.debug,"ERROR readSigprocASCIIHeader: Error reading first line (%d != 11).", j);
    exit(-1);
  }
  datafile->isFolded = 1;
  datafile->foldMode = FOLDMODE_FIXEDPERIOD;
  if(set_observatory_PSRData(datafile, tmp4, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readSigprocASCIIHeader: Setting observatory name failed.");
    return 0;
  }
  if(strcmp(tmp1, "#") != 0) {
    printerror(verbose.debug, "readSigprocASCIIHeader: I do not understand header (%s != #).", tmp1);
    return 0;
  }
  if(strcmp(tmp2, "1") != 0) {
    printerror(verbose.debug, "readSigprocASCIIHeader: I do not understand header (%s != 1).", tmp2);
    return 0;
  }
  if(strcmp(tmp3, "1") != 0) {
    printerror(verbose.debug, "readSigprocASCIIHeader: I do not understand header (%s != 1).", tmp3);
    return 0;
  }
  datafile->mjd_start = mjd;
  datafile->mjd_start += sec/(double)(24.0*60.0*60.0);
  datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
  datafile->fixedtsamp = 0;
  datafile->tsubMode = TSUBMODE_FIXEDTSUB;
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    printerror(verbose.debug,"ERROR readSigprocASCIIHeader: Memory allocation error");
    exit(-1);
  }
  datafile->tsub_list[0] = 0;
  datafile->freqMode = FREQMODE_UNIFORM;
  datafile->uniform_freq_cent = freq;
  printwarning(verbose.debug, "WARNING readSigprocASCIIHeader: Assuming there is only one polarization channel in the data");
  datafile->NrPols = 1;
  pos = ftell(datafile->fptr_hdr);
  datafile->NrSubints = 0;
  do {
    c = fgetc(datafile->fptr_hdr);
    if(c == '\n')
      datafile->NrSubints++;
  }while(c != EOF);
  datafile->NrSubints /= (datafile->NrBins + 1);
  fseek(datafile->fptr_hdr, pos, SEEK_SET);
  fgets(tmp1, 9999, datafile->fptr_hdr);
  if(strlen(tmp1) < 3)
    fgets(tmp1, 9999, datafile->fptr_hdr);
  datafile->NrFreqChan = 0;
  do {
    if(datafile->NrFreqChan == 0) {
      txtptr = strtok(tmp1, " ");
      if(
  strcmp(txtptr, "1")) {
 printerror(verbose.debug, "ERROR readSigprocASCIIHeader: Expected different bin number ('%s' != %d).", txtptr, 1);
 return 0;
      }
      txtptr = strtok(NULL, " ");
    }else {
      txtptr = strtok(NULL, " ");
    }
    if(txtptr != NULL) {
      datafile->NrFreqChan ++;
    }
  }while(txtptr != NULL);
  fseek(datafile->fptr_hdr, pos, SEEK_SET);
  return 1;
}
int writeSigprocASCIIHeader(datafile_definition datafile, verbose_definition verbose)
{
  int int_mjd;
  double sec, freq;
  int_mjd = datafile.mjd_start;
  sec = datafile.mjd_start;
  sec -= int_mjd;
  sec *= (24.0*60.0*60.0);
  freq = get_centre_freq(datafile, verbose);
  fprintf(datafile.fptr_hdr, "# %d %lf %f %d %lf %f %ld %s %d %s", int_mjd, sec, datafile.fixedPeriod, 1, freq, datafile.dm, datafile.NrBins, datafile.observatory, 1, datafile.psrname);
    fprintf(datafile.fptr_hdr, "\n");
  return 1;
}
int readSigprocASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  char txt[10000];
  long n, f, i, p;
  float sample;
  int ret;
  if(datafile.NrPols > 1) {
    printerror(verbose.debug, "readSigprocASCIIfile: Data should have just one polarization");
    return 0;
  }
  for(n = 0; n < datafile.NrSubints; n++) {
    if(verbose.verbose && verbose.nocounters == 0) printf("readSigprocASCIIfile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
    for(i = 0; i < datafile.NrBins; i++) {
      for(f = 0; f < datafile.NrFreqChan; f++) {
 ret = fscanf(datafile.fptr, "%s", txt);
 if(ret != 1) {
   printerror(verbose.debug, "ERROR readSigprocASCIIfile: Cannot read data.");
   return 0;
 }
 if(f == 0) {
   if(
      atoi(txt) != i+1) {
     printerror(verbose.debug, "ERROR readSigprocASCIIfile: Expected different bin number ('%s' != %ld).", txt, i+1);
     return 0;
   }
   ret = fscanf(datafile.fptr, "%s", txt);
   if(ret != 1) {
     printerror(verbose.debug, "ERROR readSigprocASCIIfile: Cannot read data.");
     return 0;
   }
 }
 sscanf(txt, "%f", &sample);
 p = 0;
 data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] = sample;
      }
    }
    fgets(txt, 9999, datafile.fptr_hdr);
    if(strlen(txt) < 2)
      fgets(txt, 9999, datafile.fptr_hdr);
  }
  if(verbose.verbose && verbose.nocounters == 0) printf("Reading is done.                           \n");
  return 1;
}
int writeSigprocASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long n, f, i, p;
  float sample;
  if(datafile.NrPols > 1) {
    printerror(verbose.debug, "writeSigprocASCIIfile: Data should have just one polarization");
    return 0;
  }
  p = 0;
  for(n = 0; n < datafile.NrSubints; n++) {
    if(verbose.verbose && verbose.nocounters == 0) printf("writeSigprocASCIIfile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
    for(i = 0; i < datafile.NrBins; i++) {
      for(f = 0; f < datafile.NrFreqChan; f++) {
 if(f == 0) {
     fprintf(datafile.fptr, "%5ld", i+1);
 }
 sample = data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i];
 fprintf(datafile.fptr, " %f", sample);
      }
      fprintf(datafile.fptr, "\n");
    }
    if(n != datafile.NrSubints - 1) {
      if(writeSigprocASCIIHeader(datafile, verbose) == 0) {
 printerror(verbose.debug, "writeSigprocASCIIfile: Writing header line failed.");
 return 0;
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) printf("Writing is done.                           \n");
  return 1;
}

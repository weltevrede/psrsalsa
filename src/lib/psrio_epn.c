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

#include <string.h>
#include <math.h>
#include "psrsalsa.h"







void stripspaces(char *txt)
{
  while(strlen(txt) > 0) {
    if(txt[strlen(txt)-1] == ' ')
      txt[strlen(txt)-1] = 0;
    else
      break;
  }
}

void fillspaces(char *txt, int n)
{
  memset(txt, ' ', n);
  txt[n] = 0;
}

double parse_unit_freq(char *txt, verbose_definition verbose)
{
  int i;
  stripspaces(txt);
  for(i = 0; i < 8; i++) {
    if(txt[i] != ' ')
      break;
  }
  if(i >= 7)
    i = 0;
  if(strcasecmp(&txt[i], "MHz") == 0) {
    return 1;
  }else if(strcasecmp(&txt[i], "GHz") == 0) {
    return 1000;
  }else {
    printwarning(verbose.debug, "WARNING readEPNsubHeader: Unknown frequency unit '%s', assume it is in GHz.", txt);
    return 1000;
  }
}


void write_epn_longheader(datafile_definition datafile, int version, verbose_definition verbose)
{
  int counter, nrlines, i;
  double dbl;
  char txt[100];

  nrlines = datafile.NrBins/20;
  if(datafile.NrBins > nrlines*20)
    nrlines++;

  counter=6+datafile.NrPols*datafile.NrFreqChan*(nrlines+2);

  if(version == 600) {
    fillspaces(txt, 68);
    sprintf(txt, "  Written using libpsrsalsa");
    txt[strlen(txt)] = ' ';
    txt[68] = 0;
    fprintf(datafile.fptr, "EPN 6.00%4d%s", counter, txt);
  }else {
    fillspaces(txt, 66);
    sprintf(txt, "  Written using libpsrsalsa");
    txt[strlen(txt)] = ' ';
    txt[66] = 0;
    fprintf(datafile.fptr, "EPN 6.30%6d%s", counter, txt);
  }



  fillspaces(txt, 12);
  if(strncmp(datafile.psrname, "PSR ", 4) == 0) {
    strncpy(txt, &datafile.psrname[4], 10);
  }else if(strncmp(datafile.psrname, "PSR_", 4) == 0 ) {
    strncpy(txt, &datafile.psrname[4], 10);
  }else if(strncmp(datafile.psrname, "PSR", 3) == 0 ) {
    strncpy(txt, &datafile.psrname[3], 10);
  }else {
    strncpy(txt, &datafile.psrname[0], 10);
  }
  txt[12] = 0;
  fprintf(datafile.fptr, "%12s", txt);
  fprintf(datafile.fptr, "%12s", txt);

  i = get_period(datafile, 0, &dbl, verbose);
  if(i == 2) {
    printerror(verbose.debug, "ERROR write_epn_longheader (%s): Cannot obtain period", datafile.filename);
    exit(0);
  }else if (i == 1) {
    printwarning(verbose.debug, "WARNING write_epn_longheader: Error obtaining the fold period");
  }
  if(dbl < 0) {
    printwarning(verbose.debug, "WARNING write_epn_longheader: Error obtaining the fold period");
    dbl = 0;
  }
  sprintf(txt, "%16.12lf", dbl);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%8.3f", datafile.dm);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%10.3f", datafile.rm);
  fprintf(datafile.fptr, "%s", txt);
  fillspaces(txt, 6);
  fprintf(datafile.fptr, "%s", txt);
  fillspaces(txt, 8);
  fprintf(datafile.fptr, "%s", txt);
  fillspaces(txt, 8);
  fprintf(datafile.fptr, "%s", txt);


  converthms_string(txt, datafile.ra*12.0/M_PI, 3, 4);
  fprintf(datafile.fptr, "%s", txt);
  if(datafile.ra >= 0)
    fprintf(datafile.fptr, "+");
  converthms_string(txt, datafile.dec*180.0/M_PI, 3, 4);
  fprintf(datafile.fptr, "%s", txt);

  fillspaces(txt, 8);
  if(datafile.observatory != NULL) {
    strncpy(txt, datafile.observatory, 8);


  }
  fprintf(datafile.fptr, "%8s", txt);

  sprintf(txt, "%10.3Lf", datafile.mjd_start);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%8.3f", 0.0);
  fprintf(datafile.fptr, "%s", txt);
  fprintf(datafile.fptr, " ");
  fprintf(datafile.fptr, " ");
  fillspaces(txt, 31);
  fprintf(datafile.fptr, "%s", txt);


  sprintf(txt, "%17.5lf", datafile.telescope_X);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%17.5lf", datafile.telescope_Y);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%17.5lf", datafile.telescope_Z);
  fprintf(datafile.fptr, "%s", txt);
  fillspaces(txt, 29);
  fprintf(datafile.fptr, "%s", txt);


  sprintf(txt, "01121978");
  fprintf(datafile.fptr, "%s", txt);

  i = 0;
  sscanf(datafile.scanID, "%d", &i);
  sprintf(txt, "%04d", i);
  fprintf(datafile.fptr, "%s", txt);
  if(version == 600) {
    sprintf(txt, "%04d", 0);
    fprintf(datafile.fptr, "%s", txt);
  }else {
    sprintf(txt, "%05d", 0);
    fprintf(datafile.fptr, "%s", txt);
  }
  sprintf(txt, "%02ld", datafile.NrPols);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%04ld", datafile.NrFreqChan);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%04ld", datafile.NrBins);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%12.6lf", 1e6*get_tsamp(datafile, 0, verbose));
  fprintf(datafile.fptr, "%s", txt);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%06d", 1);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%04d", 0);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%04d", 0);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, " ");
  fprintf(datafile.fptr, "%s", txt);
  if(version == 600) {
    fillspaces(txt, 15);
    fprintf(datafile.fptr, "%s", txt);
  }else {
    fillspaces(txt, 14);
    fprintf(datafile.fptr, "%s", txt);
  }
  fprintf(datafile.fptr, "--------------------------------------------------------------------------------");
}




void write_epn_shortheader(datafile_definition datafile, int version, char *idfield, int nband, verbose_definition verbose)
{
  char txt[100];


  fillspaces(txt, 8);
  strncpy(txt, idfield, 8);
  txt[8] = 0;
  fprintf(datafile.fptr, "%s", txt);

  sprintf(txt, "%04d", nband);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%04d", 1);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%12.8f", get_centre_freq(datafile, verbose)*0.001);
  fprintf(datafile.fptr, "%s", txt);
  if(version >= 610) {
    sprintf(txt, "%8s", "GHz");
    fprintf(datafile.fptr, "%s", txt);
  }
  sprintf(txt, "%12.6f", get_bw(datafile, verbose));
  fprintf(datafile.fptr, "%s", txt);
  if(version >= 610) {
    sprintf(txt, "%8s", "MHz");
    fprintf(datafile.fptr, "%s", txt);
  }
  sprintf(txt, "%17.5f", 0.0);
  fprintf(datafile.fptr, "%s", txt);
  if(version >= 620) {
    fprintf(datafile.fptr, " ");
    sprintf(txt, "%6.1f", 0.0);
    fprintf(datafile.fptr, "%s", txt);
  }else if(version == 610) {
    fillspaces(txt, 7);
    fprintf(datafile.fptr, "%s", txt);
  }else if(version == 600) {
    fillspaces(txt, 23);
    fprintf(datafile.fptr, "%s", txt);
  }
}

void write_epn_data(datafile_definition datafile, float scale, float offset, float rms, unsigned int *iprofile, verbose_definition verbose)
{
  int i;
  char txt[100];
  double dbl;

  sprintf(txt, "%#12G", scale);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%#12G", offset);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%#12G", rms);
  fprintf(datafile.fptr, "%s", txt);

  i = get_period(datafile, 0, &dbl, verbose);
  if(i == 2) {
    printerror(verbose.debug, "ERROR write_epn_data (%s): Cannot obtain period", datafile.filename);
    exit(0);
  }else if(i == 1) {
    printwarning(verbose.debug, "WARNING write_epn_data: Error obtaining the fold period");
  }
  if(dbl < 0) {
    printwarning(verbose.debug, "WARNING write_epn_data: Error obtaining the fold period");
    dbl = 0;
  }
  sprintf(txt, "%16.12lf", dbl);
  fprintf(datafile.fptr, "%s", txt);

  fillspaces(txt, 28);
  fprintf(datafile.fptr, "%s", txt);

  for(i = 0; i < datafile.NrBins; i++)
    fprintf(datafile.fptr, "%04X", iprofile[i]);
  for(i = 0; i < 20-(datafile.NrBins-(datafile.NrBins/20)*20); i++)
    fprintf(datafile.fptr, "0000");
}
int readEPNHeader(datafile_definition *datafile, verbose_definition verbose)
{
  int nlines, ret, i;
  double dbl_version;
  char txt[1000];
    ret = fread(txt, 1, 3, datafile->fptr);
    txt[3] = 0;
    if(ret != 3) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    if(strcmp(txt, "EPN") != 0) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 5, datafile->fptr);
    txt[5] = 0;
    if(ret != 5) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%lf", &dbl_version);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    datafile->version = round(dbl_version*10);
    if(verbose.debug) {
      printf("DEBUG readEPNHeader: EPN version = %lf (%d)\n", dbl_version, datafile->version);
    }
  if((datafile->version < 60) || (datafile->version > 63)) {
    printerror(verbose.debug, "ERROR readEPNHeader: EPN versions 6.0 - 6.3 are supported. This is version %lf.", dbl_version);
    return 0;
  }
    int counter;
    if(datafile->version == 63) {
      ret = fread(txt, 1, 6, datafile->fptr);
      txt[6] = 0;
      if(ret != 6) {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
      ret = sscanf(txt, "%d", &counter);
      if(ret != 1) {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
      ret = fread(txt, 1, 66, datafile->fptr);
      txt[66] = 0;
      if(ret != 66) {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
    }else {
      ret = fread(txt, 1, 4, datafile->fptr);
      txt[4] = 0;
      if(ret != 4) {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
      ret = sscanf(txt, "%d", &counter);
      if(ret != 1) {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
      ret = fread(txt, 1, 68, datafile->fptr);
      txt[68] = 0;
      if(ret != 68) {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
    }
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    stripspaces(txt);
    for(i = 0; i < 12; i++) {
      if(txt[i] != ' ')
 break;
    }
    if(i >= 11)
      i = 0;
    if(set_psrname_PSRData(datafile, &txt[i], verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readEPNHeader: Setting pulsar name failed.");
      return 0;
    }
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 16, datafile->fptr);
    txt[16] = 0;
    if(ret != 16) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    datafile->isFolded = 1;
    datafile->foldMode = FOLDMODE_FIXEDPERIOD;
    ret = sscanf(txt, "%lf", &datafile->fixedPeriod);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 8, datafile->fptr);
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%lf", &datafile->dm);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 10, datafile->fptr);
    txt[10] = 0;
    if(ret != 10) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%lf", &datafile->rm);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 6, datafile->fptr);
    txt[6] = 0;
    if(ret != 6) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 8, datafile->fptr);
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 8, datafile->fptr);
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    int hour, min;
    double sec, sign;
    ret = fread(txt, 1, 10, datafile->fptr);
    txt[10] = 0;
    if(ret != 10) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    sscanf(txt, "%2d%2d%lf", &hour, &min, &sec);
    datafile->ra = (hour+(min+sec/60.0)/60.0)*M_PI/12.0;
    ret = fread(txt, 1, 11, datafile->fptr);
    txt[11] = 0;
    if(ret != 11) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    sscanf(txt, "%3d%2d%lf", &hour, &min, &sec);
    if(hour < 0) {
      sign = -1.0;
      hour = -hour;
    }else {
      sign = +1.0;
    }
    datafile->dec = sign*(hour+(min+sec/60.0)/60.0)*M_PI/180.0;
    ret = fread(txt, 1, 8, datafile->fptr);
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    for(i = 0; i < 8; i++) {
      if(txt[i] != ' ')
 break;
    }
    if(i >= 7)
      i = 0;
    if(set_observatory_PSRData(datafile, &txt[i], verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readEPNHeader: Setting observatory name failed.");
      return 0;
    }
    ret = fread(txt, 1, 10, datafile->fptr);
    txt[10] = 0;
    if(ret != 10) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%Lf", &datafile->mjd_start);
    if(ret != 1) {
      datafile->mjd_start = 0;
    }
    ret = fread(txt, 1, 8, datafile->fptr);
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 1, datafile->fptr);
    txt[1] = 0;
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 1, datafile->fptr);
    txt[1] = 0;
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 31, datafile->fptr);
    txt[31] = 0;
    if(ret != 31) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 17, datafile->fptr);
    txt[17] = 0;
    if(ret != 17) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 17, datafile->fptr);
    txt[17] = 0;
    if(ret != 17) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 17, datafile->fptr);
    txt[17] = 0;
    if(ret != 17) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 29, datafile->fptr);
    txt[29] = 0;
    if(ret != 29) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    if(datafile->version > 60) {
      ret = fread(txt, 1, 17, datafile->fptr);
      txt[17] = 0;
      if(ret != 17) {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
    }else {
      ret = fread(txt, 1, 16, datafile->fptr);
      txt[16] = 0;
      if(ret != 16) {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
    }
    if(set_scanID_PSRData(datafile, &txt[8], verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readEPNHeader: Setting scan ID failed.");
      return 0;
    }
    ret = fread(txt, 1, 2, datafile->fptr);
    txt[2] = 0;
    if(ret != 2) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%ld", &datafile->NrPols);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    if(datafile->NrPols <= 0) {
      printwarning(verbose.debug, "WARNING readEPNHeader: The number of polarizations (%ld) does not make sense. I assume it should be one???", datafile->NrPols);
      datafile->NrPols = 1;
    }
    datafile->NrBits = 16;
    ret = fread(txt, 1, 4, datafile->fptr);
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%ld", &datafile->NrFreqChan);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 4, datafile->fptr);
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%ld", &datafile->NrBins);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
    datafile->tsubMode = TSUBMODE_FIXEDTSUB;
    if(datafile->tsub_list != NULL)
      free(datafile->tsub_list);
    datafile->tsub_list = (double *)malloc(sizeof(double));
    if(datafile->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readEPNHeader: Memory allocation error");
      return 0;
    }
    ret = sscanf(txt, "%lf", &(datafile->fixedtsamp));
    datafile->tsub_list[0] = datafile->fixedPeriod;
    datafile->fixedtsamp *= 1e-6;
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 6, datafile->fptr);
    txt[6] = 0;
    if(ret != 6) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%ld", &datafile->NrSubints);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 4, datafile->fptr);
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 4, datafile->fptr);
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 1, datafile->fptr);
    txt[1] = 0;
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    if(datafile->version > 60) {
      ret = fread(txt, 1, 14, datafile->fptr);
      txt[14] = 0;
      if(ret != 14) {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
    }else {
      ret = fread(txt, 1, 15, datafile->fptr);
      txt[15] = 0;
      if(ret != 15) {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
    }
    for(i = 0; i < 80; i++) {
      ret = fread(txt, 1, 1, datafile->fptr);
      if(ret != 1) {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
      if(txt[0] != '-') {
 printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
 return 0;
      }
    }
    nlines=datafile->NrBins/20+1;
    if(counter != 6+datafile->NrPols*datafile->NrFreqChan*(nlines+2)) {
      printerror(verbose.debug, "ERROR readEPNHeader: Something wrong with EPN file, counter doesn't make sense (%d != %ld). npol=%ld nfreq=%ld nbin=%ld", counter, 6+datafile->NrPols*datafile->NrFreqChan*(nlines+2), datafile->NrPols, datafile->NrFreqChan, datafile->NrBins);
      return 0;
    }
    int length_of_block;
    long size;
    length_of_block = counter*80;
    datafile->datastart = ftell(datafile->fptr);
    fseek(datafile->fptr, 0, SEEK_END);
    size = ftell(datafile->fptr);
    fseek(datafile->fptr, datafile->datastart, SEEK_SET);
    if(size < 0) {
      printerror(verbose.debug, "ERROR readEPNHeader: Determining file size failed.");
      return 0;
    }
    datafile->NrSubints = size/length_of_block;
    if(size != (long)datafile->NrSubints*(long)length_of_block) {
      printerror(verbose.debug, "ERROR readEPNHeader: File size does not match an integer number of subints.");
      return 0;
    }
  return 1;
}
int readEPNsubHeader(datafile_definition *datafile, float *scale, float *offset, verbose_definition verbose)
{
  int ret, version;
  double dbl_version;
  char txt[1000];
  version = datafile->version;
  version = 63;
  if((version < 60) || (version > 63)) {
    dbl_version = (version)*0.1;
    printerror(verbose.debug, "ERROR readEPNsubHeader: EPN versions 6.0 - 6.3 are supported. This is version %lf.", dbl_version);
    return 0;
  }
    ret = fread(txt, 1, 8, datafile->fptr);
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 4, datafile->fptr);
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 4, datafile->fptr);
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    if(datafile->freqMode == FREQMODE_UNKNOWN) {
      ret = sscanf(txt, "%lf", &(datafile->uniform_freq_cent));
      if(ret != 1) {
 printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
 return 0;
      }
    }
    if(version > 60) {
      ret = fread(txt, 1, 8, datafile->fptr);
      txt[8] = 0;
      if(ret != 8) {
 printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
 return 0;
      }
      stripspaces(txt);
    }else {
      txt[0] = 0;
    }
    if(datafile->freqMode == FREQMODE_UNKNOWN) {
      datafile->uniform_freq_cent *= parse_unit_freq(txt, verbose);
    }
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    if(datafile->freqMode == FREQMODE_UNKNOWN) {
      ret = sscanf(txt, "%lf", &(datafile->uniform_bw));
      if(ret != 1) {
 printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
 if(verbose.debug) {
   printerror(verbose.debug, "ERROR readEPNsubHeader: Cannot parse '%s' as a number.", txt);
 }
 return 0;
      }
    }
    if(version > 60) {
      ret = fread(txt, 1, 8, datafile->fptr);
      txt[8] = 0;
      if(ret != 8) {
 printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
 return 0;
      }
    }else {
      txt[0] = 0;
    }
    if(datafile->freqMode == FREQMODE_UNKNOWN) {
      datafile->uniform_bw *= parse_unit_freq(txt, verbose);
      datafile->freqMode = FREQMODE_UNIFORM;
    }
    ret = fread(txt, 1, 17, datafile->fptr);
    txt[17] = 0;
    if(ret != 17) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    if(version == 60) {
      ret = fread(txt, 1, 23, datafile->fptr);
      txt[23] = 0;
      if(ret != 23) {
 printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
 return 0;
      }
    }else if(version == 61) {
      ret = fread(txt, 1, 7, datafile->fptr);
      txt[7] = 0;
      if(ret != 7) {
 printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
 return 0;
      }
    }else if(version >= 62) {
      ret = fread(txt, 1, 1, datafile->fptr);
      txt[1] = 0;
      if(ret != 1) {
 printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
 return 0;
      }
      ret = fread(txt, 1, 6, datafile->fptr);
      txt[6] = 0;
      if(ret != 6) {
 printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
 return 0;
      }
    }
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%f", scale);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%f", offset);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 16, datafile->fptr);
    txt[16] = 0;
    if(ret != 16) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 28, datafile->fptr);
    txt[28] = 0;
    if(ret != 28) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
  return 1;
}
int writeEPNfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  float dmin, dmax, scale, offset;
  long n, p, i, f;
  int maxint;
  unsigned int *iprofile;
  maxint = (1 << 16);
  if(datafile.NrFreqChan != 1) {
    printerror(verbose.debug, "ERROR writeEPNfile: Only one freq channel is supported in this data format");
    return 0;
  }
  iprofile = (unsigned int *)calloc(datafile.NrBins, sizeof(unsigned int));
  if(iprofile == NULL) {
    printerror(verbose.debug, "ERROR writeEPNfile: Cannot allocate memory.");
    return 0;
  }
  f = 0;
  for(n = 0; n < datafile.NrSubints; n++) {
    if(verbose.verbose && verbose.nocounters == 0) {
      printf("writeEPNfile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
    }
    for(p = 0; p < datafile.NrPols; p++) {
      for(i = 0; i < datafile.NrPols*datafile.NrBins; i++) {
 if((i == 0 && p == 0) || dmax < data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i]) {
   dmax = data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i];
 }
 if((i == 0 && p == 0) || dmin > data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i]) {
   dmin = data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i];
 }
      }
    }
    if(dmax == dmin)
      dmax = dmin+1;
    offset = dmin;
    scale = (float)(maxint - 1)/(dmax - dmin);
    for(i = 0; i < datafile.NrBins; i++)
      iprofile[i] = (int)((data[datafile.NrBins*(0+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] - offset) * scale);
    write_epn_longheader(datafile, 600, verbose);
    write_epn_shortheader(datafile, 630, "I       ", 1, verbose);
    write_epn_data(datafile, scale, offset, 0.0, iprofile, verbose);
    if(datafile.NrPols == 4) {
      for(i = 0; i < datafile.NrBins; i++)
 iprofile[i] = (int)((data[datafile.NrBins*(1+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] - offset) * scale);
      write_epn_shortheader(datafile, 630, "Q       ", 2, verbose);
      write_epn_data(datafile, scale, offset, 0.0, iprofile, verbose);
      for(i = 0; i < datafile.NrBins; i++)
 iprofile[i] = (int)((data[datafile.NrBins*(2+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] - offset) * scale);
      write_epn_shortheader(datafile, 630, "U       ", 3, verbose);
      write_epn_data(datafile, scale, offset, 0.0, iprofile, verbose);
      for(i = 0; i < datafile.NrBins; i++)
 iprofile[i] = (int)((data[datafile.NrBins*(3+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] - offset) * scale);
      write_epn_shortheader(datafile, 630, "V       ", 4, verbose);
      write_epn_data(datafile, scale, offset, 0.0, iprofile, verbose);
    }
  }
  free(iprofile);
  printf("  Done                       \n");
  return 1;
}
int readPulseEPNData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  float scale, offset;
  int n;
  long pos, ret, nrrecords;
  nrrecords = datafile->NrBins/20;
  if(datafile->NrBins > nrrecords*20)
    nrrecords += 1;
  if(datafile->NrFreqChan > 1) {
    printerror(verbose.debug, "ERROR readPulseEPNData: Only one freq channel is supported in this data format");
    return 0;
  }
  pos = 480+(160+nrrecords*80)*(datafile->NrPols*datafile->NrFreqChan);
  pos *= pulsenr;
  pos += 480;
  pos += (160+nrrecords*80)*polarization;
  fseek(datafile->fptr, pos, SEEK_SET);
  if(readEPNsubHeader(datafile, &scale, &offset, verbose) == 0) {
    printerror(verbose.debug, "ERROR readPulseEPNData: Reading subheader failed.");
    return 0;
  }
  fseek(datafile->fptr, 4*binnr, SEEK_CUR);
  unsigned int iprofile;
  for(n = 0; n < nrSamples; n++) {
    ret = fscanf(datafile->fptr,"%04X", &iprofile);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readPulseEPNData: Reading data failed.");
      return 0;
    }
    pulse[n] = (float)iprofile / scale + offset;
  }
  return 1;
}
int readEPNfile(datafile_definition *datafile, float *data, verbose_definition verbose, long request_only_one_pulse)
{
  long n, f, p;
  if(verbose.verbose) {
    printf("Start reading EPN file\n");
  }
  for(n = 0; n < datafile->NrSubints; n++) {
    for(p = 0; p < datafile->NrPols; p++) {
      for(f = 0; f < datafile->NrFreqChan; f++) {
        if(verbose.verbose && verbose.nocounters == 0)
          printf("  Progress reading EPN file (%.1f%%)\r", 100.0*(n+(f+p*datafile->NrFreqChan)*datafile->NrSubints)/(float)(datafile->NrSubints*datafile->NrFreqChan*datafile->NrPols));
 if(readPulseEPNData(datafile, n, p, f, 0, datafile->NrBins, &data[datafile->NrBins*(p+datafile->NrPols*(f+n*datafile->NrFreqChan))], verbose) == 0) {
   printerror(verbose.debug, "ERROR readEPNfile: Read error.");
   return 0;
 }
      }
    }
  }
  if(verbose.verbose) printf("  Reading is done.                           \n");
  return 1;
}

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
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include "psrsalsa.h"
int readWSRTHeader(datafile_definition *datafile, verbose_definition verbose);
int readPulseWSRTData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse);
int writePulseWSRTData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse);
int writeWSRTHeader(datafile_definition datafile, verbose_definition verbose);
int writePuMafile(datafile_definition datafile, float *data, verbose_definition verbose);
int readPuMafile(datafile_definition datafile, float *data, verbose_definition verbose);
int readPSRFITSHeader(datafile_definition *datafile, int readnoscales, int nowarnings, verbose_definition verbose);
int readFITSpulse(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
int readFITSfile(datafile_definition *datafile, float *data, verbose_definition verbose);
int writePSRFITSHeader(datafile_definition *datafile, verbose_definition verbose);
int writeFITSpulse(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
int writeFITSfile(datafile_definition datafile, float *data, verbose_definition verbose);
int readPSRCHIVE_ASCIIHeader(datafile_definition *datafile, verbose_definition verbose);
int readPSRCHIVE_ASCIIfilepulse(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
int readPSRCHIVE_ASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose);
int writePSRCHIVE_ASCIIHeader(datafile_definition datafile, verbose_definition verbose);
int writePSRCHIVE_ASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose);
int readEPNHeader(datafile_definition *datafile, int what, verbose_definition verbose);
int readPulseEPNData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
int writeEPNfile(datafile_definition datafile, float *data, verbose_definition verbose);
int readEPNfile(datafile_definition *datafile, float *data, verbose_definition verbose, long request_only_one_pulse);
int readEPNsubHeader(datafile_definition *datafile, float *scale, float *offset, verbose_definition verbose);
int readSigprocHeader(datafile_definition *datafile, verbose_definition verbose);
int readPulseSigprocData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
int readPPOLHeader(datafile_definition *datafile, int extended, verbose_definition verbose);
int writePPOLHeader(datafile_definition datafile, int argc, char **argv, verbose_definition verbose);
int readHistoryFITS(datafile_definition *datafile, verbose_definition verbose);
int writeHistoryFITS(datafile_definition datafile, verbose_definition verbose);
int readHistoryPSRData(datafile_definition *datafile, verbose_definition verbose);
int writeHistoryPSRData(datafile_definition *datafile, int argc, char **argv, int cmdOnly, char *notes, verbose_definition verbose);
int writeHistoryPuma(datafile_definition datafile, verbose_definition verbose);
int readHistoryPuma(datafile_definition *datafile, verbose_definition verbose);
int readSigprocfile(datafile_definition datafile, float *data, verbose_definition verbose);
int readSigprocASCIIHeader(datafile_definition *datafile, verbose_definition verbose);
int writeSigprocASCIIHeader(datafile_definition datafile, verbose_definition verbose);
int writeSigprocASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose);
int readSigprocASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose);
void closeHistoryPSRData(datafile_definition *datafile, int remove_last_entry_only);
int isValidPSRDATA_format(int format)
{
  if(format == PUMA_format)
    return 1;
  if(format == PSRCHIVE_ASCII_format)
    return 1;
  if(format == EPN_format)
    return 1;
  if(format == FITS_format)
    return 1;
  if(format == SIGPROC_format)
    return 1;
  if(format == PPOL_format)
    return 1;
  if(format == PPOL_SHORT_format)
    return 1;
  if(format == SIGPROC_ASCII_format)
    return 1;
  if(format == PSRSALSA_BINARY_format)
    return 1;
  if(format == MEMORY_format)
    return 1;
  printerror(0, "ERROR isValidPSRDATA_format: specified data format is not recognized.");
  return 0;
}
void printPSRDataFormats(FILE *printdevice, int nrspaces)
{
  int i, nrspaces2;
  nrspaces2 = nrspaces + 17;
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "(PUMA)         - WSRT PuMa format\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "(ASCII)        - PSRCHIVE ascii dump file (generated by for instance\n");
  for(i = 0; i < nrspaces2; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "\"pdv -t\"). File has limited header information.\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "(EPN)          - EPN format\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "(PSRFITS)      - PSRFITS format (generated by for instance\n");
  for(i = 0; i < nrspaces2; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "\"pam -a PSRFITS\"). Note that the data files written out\n");
  for(i = 0; i < nrspaces2; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "do not 100%% conform with the PSRFITS definition, so\n");
  for(i = 0; i < nrspaces2; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "proper behavior in other software cannot be guaranteed.\n");
  for(i = 0; i < nrspaces2; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "Especially timing experiments are not recommended.\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "(SIGPROC)      - SIGPROC binary format\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "(PPOL)         - PPOL format\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "(PPOLSHORT)    - PPOL SHORT format (longitude, pa, pa error)\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "(SIGPROCASCII) - Sigproc ascii format\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "(PSRSALSA)     - PSRSALSA binary format\n");
}
int parsePSRDataFormats(char *cmd)
{
  if(strcasecmp(cmd, "ASCII") == 0 || atoi(cmd) == PSRCHIVE_ASCII_format)
    return PSRCHIVE_ASCII_format;
  else if(strcasecmp(cmd, "PSRFITS") == 0 || strcasecmp(cmd, "FITS") == 0 || atoi(cmd) == FITS_format)
    return FITS_format;
  else if(strcasecmp(cmd, "PSRSALSA") == 0 || strcasecmp(cmd, "SALSA") == 0 || atoi(cmd) == PSRSALSA_BINARY_format)
    return PSRSALSA_BINARY_format;
  else if(strcasecmp(cmd, "PUMA") == 0 || atoi(cmd) == PUMA_format)
    return PUMA_format;
  else if(strcasecmp(cmd, "EPN") == 0 || atoi(cmd) == EPN_format)
    return EPN_format;
  else if(strcasecmp(cmd, "SIGPROC") == 0 || atoi(cmd) == SIGPROC_format)
    return SIGPROC_format;
  else if(strcasecmp(cmd, "PPOL") == 0 || strcasecmp(cmd, "PASWING") == 0 || atoi(cmd) == PPOL_format)
    return PPOL_format;
  else if(strcasecmp(cmd, "PPOLSHORT") == 0 || strcasecmp(cmd, "PPOL_SHORT") == 0 || strcasecmp(cmd, "PASWINGSHORT") == 0 || atoi(cmd) == PPOL_SHORT_format)
    return PPOL_SHORT_format;
  else if(strcasecmp(cmd, "SIGPROCASCII") == 0 || strcasecmp(cmd, "SIGPROC_ASCII") == 0 || atoi(cmd) == SIGPROC_ASCII_format)
    return SIGPROC_ASCII_format;
  else {
    fflush(stdout);
    printerror(0, "parsePSRDataFormats: Cannot parse '%s' as a valid data format", cmd);
    return 0;
  }
  return 0;
}
void cleanPSRData(datafile_definition *datafile, verbose_definition verbose)
{
  memset(datafile, 0, sizeof(datafile_definition));
  datafile->fptr = NULL;
  datafile->fptr_hdr = NULL;
  datafile->fits_fptr = NULL;
  datafile->scales = NULL;
  datafile->offsets = NULL;
  datafile->weights = NULL;
  datafile->weight_stats_set = 0;
  datafile->data = NULL;
  datafile->format = 0;
  datafile->version = 0;
  datafile->opened_flag = 0;
  datafile->enable_write_flag = 0;
  datafile->dumpOnClose = 0;
  datafile->filename = malloc(1);
  if(datafile->filename == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cleanPSRData: Memory allocation error.");
    exit(0);
  }
  datafile->filename[0] = 0;
  datafile->NrSubints = 0;
  datafile->NrBins = 0;
  datafile->NrBits = 0;
  datafile->NrPols = 0;
  datafile->NrFreqChan = 0;
  datafile->isFolded = -1;
  datafile->foldMode = FOLDMODE_UNKNOWN;
  datafile->fixedPeriod = 0;
  datafile->tsampMode = TSAMPMODE_UNKNOWN;
  datafile->fixedtsamp = 0;
  datafile->tsamp_list = NULL;
  datafile->tsubMode = TSUBMODE_UNKNOWN;
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cleanPSRData: Memory allocation error.");
    exit(0);
  }
  datafile->tsub_list[0] = 0;
  datafile->freq_ref = -2;
  datafile->freqMode = FREQMODE_UNKNOWN;
  datafile->freqlabel_list = NULL;
  datafile->bandwidth = 0;
  datafile->centrefreq = 0;
  datafile->ra = 0;
  datafile->dec = 0;
  datafile->dm = 0;
  datafile->rm = 0;
  datafile->mjd_start = 0;
  datafile->psrname = malloc(1);
  datafile->observatory = malloc(1);
  datafile->institute = malloc(1);
  datafile->instrument = malloc(1);
  datafile->scanID = malloc(1);
  datafile->observer = malloc(1);
  datafile->projectID = malloc(1);
  if(datafile->psrname == NULL || datafile->observatory == NULL || datafile->institute == NULL || datafile->instrument == NULL || datafile->scanID == NULL || datafile->observer == NULL || datafile->projectID == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cleanPSRData: Memory allocation error.");
    exit(0);
  }
  datafile->psrname[0] = 0;
  datafile->observatory[0] = 0;
  datafile->institute[0] = 0;
  datafile->instrument[0] = 0;
  datafile->scanID[0] = 0;
  datafile->observer[0] = 0;
  datafile->projectID[0] = 0;
  datafile->offpulse_rms = NULL;
  datafile->feedtype = FEEDTYPE_UNKNOWN;
  datafile->poltype = POLTYPE_UNKNOWN;
  datafile->datastart = 0;
  datafile->isTransposed = 0;
  datafile->gentype = GENTYPE_UNDEFINED;
  datafile->xrangeset = 0;
  datafile->xrange[0] = 0;
  datafile->xrange[1] = 0;
  datafile->yrangeset = 0;
  datafile->yrange[0] = 0;
  datafile->yrange[1] = 0;
  datafile->isDeDisp = -1;
  datafile->isDeFarad = -1;
  datafile->isDePar = -1;
  datafile->isDebase = -1;
  datafile->telescope_X = 0;
  datafile->telescope_Y = 0;
  datafile->telescope_Z = 0;
  datafile->cableSwap = -1;
  datafile->cableSwapcor = -1;
  datafile->history.timestamp = NULL;
  datafile->history.cmd = NULL;
  datafile->history.user = NULL;
  datafile->history.hostname = NULL;
  datafile->history.notes = NULL;
  datafile->history.nextEntry = NULL;
}
int copy_params_PSRData(datafile_definition datafile_source, datafile_definition *datafile_dest, verbose_definition verbose)
{
  datafile_dest->fptr = NULL;
  datafile_dest->fptr_hdr = NULL;
  datafile_dest->fits_fptr = NULL;
  datafile_dest->scales = NULL;
  datafile_dest->offsets = NULL;
  datafile_dest->weights = NULL;
  datafile_dest->weight_stats_set = datafile_source.weight_stats_set;
  datafile_dest->weight_stats_zeroweightfound = datafile_source.weight_stats_zeroweightfound;
  datafile_dest->weight_stats_differentweights = datafile_source.weight_stats_differentweights;
  datafile_dest->weight_stats_negativeweights = datafile_source.weight_stats_negativeweights;
  datafile_dest->weight_stats_weightvalue = datafile_source.weight_stats_weightvalue;
  datafile_dest->offpulse_rms = NULL;
  datafile_dest->format = datafile_source.format;
  datafile_dest->version = datafile_source.version;
  datafile_dest->NrSubints = datafile_source.NrSubints;
  datafile_dest->NrBins = datafile_source.NrBins;
  datafile_dest->NrBits = datafile_source.NrBits;
  datafile_dest->NrPols = datafile_source.NrPols;
  datafile_dest->NrFreqChan = datafile_source.NrFreqChan;
  datafile_dest->isFolded = datafile_source.isFolded;
  datafile_dest->foldMode = datafile_source.foldMode;
  datafile_dest->fixedPeriod = datafile_source.fixedPeriod;
  datafile_dest->tsampMode = datafile_source.tsampMode;
  if(datafile_source.tsampMode == TSAMPMODE_LONGITUDELIST && datafile_source.tsamp_list != NULL) {
    if(datafile_dest->tsamp_list != NULL)
      free(datafile_dest->tsamp_list);
    datafile_dest->tsamp_list = (double *)malloc(datafile_source.NrBins*sizeof(double));
    if(datafile_dest->tsamp_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
      return 0;
    }
    memcpy(datafile_dest->tsamp_list, datafile_source.tsamp_list, datafile_source.NrBins*sizeof(double));
  }
  datafile_dest->fixedtsamp = datafile_source.fixedtsamp;
  datafile_dest->tsubMode = datafile_source.tsubMode;
  if(datafile_dest->tsub_list != NULL) {
    free(datafile_dest->tsub_list);
  }
  if(datafile_source.tsubMode == TSUBMODE_TSUBLIST) {
    datafile_dest->tsub_list = (double *)malloc(datafile_source.NrSubints*sizeof(double));
    if(datafile_dest->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
      return 0;
    }
    memcpy(datafile_dest->tsub_list, datafile_source.tsub_list, datafile_source.NrSubints*sizeof(double));
  }else {
    datafile_dest->tsub_list = (double *)malloc(sizeof(double));
    if(datafile_dest->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
      return 0;
    }
    memcpy(datafile_dest->tsub_list, datafile_source.tsub_list, sizeof(double));
  }
  datafile_dest->freq_ref = datafile_source.freq_ref;
  datafile_dest->freqMode = datafile_source.freqMode;
  if(datafile_dest->freqlabel_list != NULL) {
    free(datafile_dest->freqlabel_list);
    datafile_dest->freqlabel_list = NULL;
  }
  datafile_dest->bandwidth = datafile_source.bandwidth;
  datafile_dest->centrefreq = datafile_source.centrefreq;
  if(datafile_source.freqMode == FREQMODE_FREQTABLE) {
    datafile_dest->freqlabel_list = (double *)malloc(datafile_source.NrFreqChan*datafile_source.NrSubints*sizeof(double));
    if(datafile_dest->freqlabel_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
      return 0;
    }
    memcpy(datafile_dest->freqlabel_list, datafile_source.freqlabel_list, datafile_source.NrFreqChan*datafile_source.NrSubints*sizeof(double));
  }else {
    datafile_dest->freqlabel_list = NULL;
  }
  datafile_dest->ra = datafile_source.ra;
  datafile_dest->dec = datafile_source.dec;
  datafile_dest->dm = datafile_source.dm;
  datafile_dest->rm = datafile_source.rm;
  datafile_dest->mjd_start = datafile_source.mjd_start;
  if(set_filename_PSRData(datafile_dest, datafile_source.filename, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting file name failed.");
    return 0;
  }
  if(set_psrname_PSRData(datafile_dest, datafile_source.psrname, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting pulsar name failed.");
    return 0;
  }
  if(set_observatory_PSRData(datafile_dest, datafile_source.observatory, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting observatory name failed.");
    return 0;
  }
  if(set_institute_PSRData(datafile_dest, datafile_source.institute, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting institute name failed.");
    return 0;
  }
  if(set_instrument_PSRData(datafile_dest, datafile_source.instrument, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting instrument name failed.");
    return 0;
  }
  if(set_scanID_PSRData(datafile_dest, datafile_source.scanID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting scan ID name failed.");
    return 0;
  }
  if(set_observer_PSRData(datafile_dest, datafile_source.observer, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting name of observer failed.");
    return 0;
  }
  if(set_projectID_PSRData(datafile_dest, datafile_source.projectID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting project ID failed.");
    return 0;
  }
  datafile_dest->telescope_X = datafile_source.telescope_X;
  datafile_dest->telescope_Y = datafile_source.telescope_Y;
  datafile_dest->telescope_Z = datafile_source.telescope_Z;
  datafile_dest->feedtype = datafile_source.feedtype;
  datafile_dest->poltype = datafile_source.poltype;
  datafile_dest->isTransposed = datafile_source.isTransposed;
  datafile_dest->gentype = datafile_source.gentype;
  datafile_dest->xrangeset = datafile_source.xrangeset;
  datafile_dest->xrange[0] = datafile_source.xrange[0];
  datafile_dest->xrange[1] = datafile_source.xrange[1];
  datafile_dest->yrangeset = datafile_source.yrangeset;
  datafile_dest->yrange[0] = datafile_source.yrange[0];
  datafile_dest->yrange[1] = datafile_source.yrange[1];
  datafile_dest->isDeDisp = datafile_source.isDeDisp;
  datafile_dest->isDeFarad = datafile_source.isDeFarad;
  datafile_dest->isDePar = datafile_source.isDePar;
  datafile_dest->isDebase = datafile_source.isDebase;
  datafile_dest->cableSwap = datafile_source.cableSwap;
  datafile_dest->cableSwapcor = datafile_source.cableSwapcor;
  closeHistoryPSRData(datafile_dest, 0);
  datafile_history_entry_definition *dest_hist, *source_hist;
  dest_hist = &(datafile_dest->history);
  source_hist = &(datafile_source.history);
  int ok = 1;
  do {
    dest_hist->timestamp = NULL;
    dest_hist->cmd = NULL;
    dest_hist->user = NULL;
    dest_hist->hostname = NULL;
    dest_hist->notes = NULL;
    dest_hist->nextEntry = NULL;
    if(source_hist->timestamp != NULL) {
      dest_hist->timestamp = malloc(strlen(source_hist->timestamp)+1);
      if(dest_hist->timestamp == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
 return 0;
      }
      strcpy(dest_hist->timestamp, source_hist->timestamp);
    }
    if(source_hist->cmd != NULL) {
      dest_hist->cmd = malloc(strlen(source_hist->cmd)+1);
      if(dest_hist->cmd == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
 return 0;
      }
      strcpy(dest_hist->cmd, source_hist->cmd);
    }
    if(source_hist->user != NULL) {
      dest_hist->user = malloc(strlen(source_hist->user)+1);
      if(dest_hist->user == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
 return 0;
      }
      strcpy(dest_hist->user, source_hist->user);
    }
    if(source_hist->hostname != NULL) {
      dest_hist->hostname = malloc(strlen(source_hist->hostname)+1);
      if(dest_hist->hostname == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
 return 0;
      }
      strcpy(dest_hist->hostname, source_hist->hostname);
    }
    if(source_hist->notes != NULL) {
      dest_hist->notes = malloc(strlen(source_hist->notes)+1);
      if(dest_hist->notes == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
 return 0;
      }
      strcpy(dest_hist->notes, source_hist->notes);
    }
    if(source_hist->nextEntry != NULL) {
      dest_hist->nextEntry = malloc(sizeof(datafile_history_entry_definition));
      if(dest_hist->nextEntry == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
 return 0;
      }
      dest_hist = dest_hist->nextEntry;
      source_hist = source_hist->nextEntry;
    }else {
      ok = 0;
    }
  }while(ok);
  return 1;
}
int writePSRSALSAHeader(datafile_definition *datafile, verbose_definition verbose)
{
  int ret, dummyi;
  long dummyl;
  char identifier[] = "PSRSALSAdump";
  int version = 4;
  ret = fwrite(identifier, 12, 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&version, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  dummyi = strlen(datafile->psrname);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->psrname, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  dummyi = strlen(datafile->scanID);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->scanID, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  dummyi = strlen(datafile->observer);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->observer, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  dummyi = strlen(datafile->projectID);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->projectID, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  dummyi = strlen(datafile->observatory);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->observatory, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  dummyi = strlen(datafile->instrument);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->instrument, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  dummyi = strlen(datafile->institute);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->institute, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  ret = fwrite(&(datafile->telescope_X), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->telescope_Y), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->telescope_Z), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->NrBits), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->isDeDisp), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->isDeFarad), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->isDePar), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->isDebase), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->cableSwap), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->cableSwapcor), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->dm), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->rm), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->freq_ref), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->feedtype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->poltype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->mjd_start), sizeof(long double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->NrSubints), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->NrBins), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->NrPols), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->NrFreqChan), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->isFolded), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->foldMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->fixedPeriod), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->tsampMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->fixedtsamp), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(datafile->tsampMode == TSAMPMODE_LONGITUDELIST) {
    ret = fwrite(datafile->tsamp_list, sizeof(double), datafile->NrBins, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  ret = fwrite(&(datafile->tsubMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(datafile->tsubMode == TSUBMODE_TSUBLIST) {
    ret = fwrite(datafile->tsub_list, sizeof(double), datafile->NrSubints, datafile->fptr_hdr);
    if(ret != datafile->NrSubints) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s (return value is %ld)", datafile->filename, ret);
      return 0;
    }
  }else {
    ret = fwrite(datafile->tsub_list, sizeof(double), 1, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  ret = fwrite(&(datafile->ra), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->dec), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->freqMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->centrefreq), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->bandwidth), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(datafile->freqMode == FREQMODE_FREQTABLE) {
    ret = fwrite(datafile->freqlabel_list, sizeof(double), datafile->NrSubints*datafile->NrFreqChan, datafile->fptr_hdr);
    if(ret != datafile->NrSubints*datafile->NrFreqChan) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }else if(datafile->freqMode != FREQMODE_UNKNOWN && datafile->freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Frequency mode is not implemented");
    return 0;
  }
  ret = fwrite(&(datafile->gentype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->isTransposed), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->xrangeset), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(datafile->xrangeset) {
    ret = fwrite(datafile->xrange, sizeof(float), 2, datafile->fptr_hdr);
    if(ret != 2) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  ret = fwrite(&(datafile->yrangeset), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(datafile->yrangeset) {
    ret = fwrite(datafile->yrange, sizeof(float), 2, datafile->fptr_hdr);
    if(ret != 2) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  if(datafile->offpulse_rms == NULL) {
    dummyl = 0;
  }else {
    dummyl = datafile->NrSubints * datafile->NrFreqChan * datafile->NrPols;
  }
  ret = fwrite(&dummyl, sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyl != 0) {
    ret = fwrite(datafile->offpulse_rms, sizeof(float), dummyl, datafile->fptr_hdr);
    if(ret != dummyl) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s (return value is %ld)", datafile->filename, ret);
      return 0;
    }
  }
  datafile->datastart = ftell(datafile->fptr_hdr);
  return 1;
}
int writeHistoryPSRSALSA(datafile_definition *datafile, verbose_definition verbose)
{
  datafile_history_entry_definition *curHistoryEntry;
  int nrHistoryLines, dummyi;
  size_t ret;
  char *history_id = "HISTORY";
  curHistoryEntry = &(datafile->history);
  nrHistoryLines = 0;
  do {
    if(curHistoryEntry->timestamp != NULL || curHistoryEntry->cmd != NULL || curHistoryEntry->user != NULL || curHistoryEntry->hostname != NULL || curHistoryEntry->notes != NULL || curHistoryEntry->nextEntry != NULL) {
      nrHistoryLines++;
      curHistoryEntry = curHistoryEntry->nextEntry;
    }
  }while(curHistoryEntry != NULL);
  if(verbose.debug) {
    printf("Start writing history in PSRSALSA binary format (%d lines)\n", nrHistoryLines);
  }
  ret = fwrite(history_id, 1, strlen(history_id), datafile->fptr_hdr);
  if(ret != strlen(history_id)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(nrHistoryLines), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  curHistoryEntry = &(datafile->history);
  do {
    if(curHistoryEntry->timestamp != NULL || curHistoryEntry->cmd != NULL || curHistoryEntry->user != NULL || curHistoryEntry->hostname != NULL || curHistoryEntry->notes != NULL || curHistoryEntry->nextEntry != NULL) {
      if(curHistoryEntry->timestamp == NULL) {
 dummyi = 0;
      }else {
 dummyi = strlen(curHistoryEntry->timestamp);
      }
      ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
 return 0;
      }
      if(dummyi > 0) {
 ret = fwrite(curHistoryEntry->timestamp, sizeof(char), dummyi, datafile->fptr_hdr);
 if(ret != dummyi) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
   return 0;
 }
      }
      if(curHistoryEntry->user == NULL) {
 dummyi = 0;
      }else {
 dummyi = strlen(curHistoryEntry->user);
      }
      ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
 return 0;
      }
      if(dummyi > 0) {
 ret = fwrite(curHistoryEntry->user, sizeof(char), dummyi, datafile->fptr_hdr);
 if(ret != dummyi) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
   return 0;
 }
      }
      if(curHistoryEntry->hostname == NULL) {
 dummyi = 0;
      }else {
 dummyi = strlen(curHistoryEntry->hostname);
      }
      ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
 return 0;
      }
      if(dummyi > 0) {
 ret = fwrite(curHistoryEntry->hostname, sizeof(char), dummyi, datafile->fptr_hdr);
 if(ret != dummyi) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
   return 0;
 }
      }
      if(curHistoryEntry->cmd == NULL) {
 dummyi = 0;
      }else {
 dummyi = strlen(curHistoryEntry->cmd);
      }
      ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
 return 0;
      }
      if(dummyi > 0) {
 ret = fwrite(curHistoryEntry->cmd, sizeof(char), dummyi, datafile->fptr_hdr);
 if(ret != dummyi) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
   return 0;
 }
      }
      if(curHistoryEntry->notes == NULL) {
 dummyi = 0;
      }else {
 dummyi = strlen(curHistoryEntry->notes);
      }
      ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
 return 0;
      }
      if(dummyi > 0) {
 ret = fwrite(curHistoryEntry->notes, sizeof(char), dummyi, datafile->fptr_hdr);
 if(ret != dummyi) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
   return 0;
 }
      }
      curHistoryEntry = curHistoryEntry->nextEntry;
    }
  }while(curHistoryEntry != NULL);
  datafile->datastart = ftell(datafile->fptr_hdr);
  return 1;
}
int readPSRSALSAHeader(datafile_definition *datafile, int nohistory_expected, verbose_definition verbose)
{
  int ret, dummyi;
  long dummyl;
  char identifier[13], *txt;
  int version, maxversion_supported = 4;
  txt = malloc(10000);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error", datafile->filename);
    return 0;
  }
  ret = fread(identifier, 1, 12, datafile->fptr_hdr);
  if(ret != 12) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  identifier[12] = 0;
  if(strcmp(identifier, "PSRSALSAdump") != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: File does not appear to be in the expected format", datafile->filename);
    return 0;
  }
  ret = fread(&version, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(version < 1 || version > maxversion_supported) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: File %s is in an unsupported version number (%d). The maximum supported version number for reading in your installation is %d. An update of the software might be required to read this file.", datafile->filename, version, maxversion_supported);
    return 0;
  }
  if(verbose.debug) {
    printf("  PSRSALSA file is version %d\n", version);
  }
  ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(dummyi > 9999) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    txt[dummyi] = 0;
    if(set_psrname_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting pulsar name failed for %s", datafile->filename);
      return 0;
    }
  }
  ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(dummyi > 9999) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    txt[dummyi] = 0;
    if(set_scanID_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting scan ID failed for %s", datafile->filename);
      return 0;
    }
  }
  if(version > 1) {
    ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    if(dummyi > 9999) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
      return 0;
    }
    if(dummyi > 0) {
      ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
      if(ret != dummyi) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
 return 0;
      }
      txt[dummyi] = 0;
      if(set_observer_PSRData(datafile, txt, verbose) == 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting scan ID failed for %s", datafile->filename);
 return 0;
      }
    }
    ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    if(dummyi > 9999) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
      return 0;
    }
    if(dummyi > 0) {
      ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
      if(ret != dummyi) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
 return 0;
      }
      txt[dummyi] = 0;
      if(set_projectID_PSRData(datafile, txt, verbose) == 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting scan ID failed for %s", datafile->filename);
 return 0;
      }
    }
  }
  ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(dummyi > 9999) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    txt[dummyi] = 0;
    if(set_observatory_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting observatory name failed for %s", datafile->filename);
      return 0;
    }
  }
  ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(dummyi > 9999) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    txt[dummyi] = 0;
    if(set_instrument_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting instrument name failed for %s", datafile->filename);
      return 0;
    }
  }
  ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(dummyi > 9999) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    txt[dummyi] = 0;
    if(set_institute_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting institute name failed for %s", datafile->filename);
      return 0;
    }
  }
  ret = fread(&(datafile->telescope_X), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->telescope_Y), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->telescope_Z), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->NrBits), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->isDeDisp), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->isDeFarad), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->isDePar), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->isDebase), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->cableSwap), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->cableSwapcor), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->dm), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->rm), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->freq_ref), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->feedtype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->poltype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->mjd_start), sizeof(long double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->NrSubints), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->NrBins), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->NrPols), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->NrFreqChan), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->isFolded), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->foldMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->fixedPeriod), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->tsampMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->fixedtsamp), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  datafile->tsamp_list = NULL;
  if(datafile->tsampMode == TSAMPMODE_LONGITUDELIST) {
    datafile->tsamp_list = (double *)malloc(sizeof(double)*datafile->NrBins);
    if(datafile->tsamp_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading %s", datafile->filename);
      return 0;
    }
    ret = fread(datafile->tsamp_list, sizeof(double), datafile->NrBins, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  ret = fread(&(datafile->tsubMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(datafile->tsub_list != NULL)
    free(datafile->tsub_list);
  if(datafile->tsubMode == TSUBMODE_TSUBLIST) {
    datafile->tsub_list = (double *)malloc(sizeof(double)*datafile->NrSubints);
    if(datafile->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading %s", datafile->filename);
      return 0;
    }
    ret = fread(datafile->tsub_list, sizeof(double), datafile->NrSubints, datafile->fptr_hdr);
    if(ret != datafile->NrSubints) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
  }else {
    datafile->tsub_list = (double *)malloc(sizeof(double));
    if(datafile->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading %s", datafile->filename);
      return 0;
    }
    ret = fread(datafile->tsub_list, sizeof(double), 1, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
  }
  ret = fread(&(datafile->ra), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->dec), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->freqMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->centrefreq), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->bandwidth), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(datafile->freqMode == FREQMODE_FREQTABLE) {
    if(datafile->freqlabel_list != NULL) {
      free(datafile->freqlabel_list);
    }
    datafile->freqlabel_list = malloc(datafile->NrSubints*datafile->NrFreqChan*sizeof(double));
    if(datafile->freqlabel_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error");
      return 0;
    }
    ret = fread(datafile->freqlabel_list, sizeof(double), datafile->NrSubints*datafile->NrFreqChan, datafile->fptr_hdr);
    if(ret != datafile->NrSubints*datafile->NrFreqChan) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
  }
  ret = fread(&(datafile->gentype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->isTransposed), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->xrangeset), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(datafile->xrangeset) {
    ret = fread(datafile->xrange, sizeof(float), 2, datafile->fptr_hdr);
    if(ret != 2) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
  }
  ret = fread(&(datafile->yrangeset), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(datafile->yrangeset) {
    ret = fread(datafile->yrange, sizeof(float), 2, datafile->fptr_hdr);
    if(ret != 2) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
  }
  if(version > 2) {
    ret = fread(&(dummyl), sizeof(long), 1, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    if(datafile->offpulse_rms != NULL) {
      free(datafile->offpulse_rms);
      datafile->offpulse_rms = NULL;
    }
    if(dummyl > 0) {
      datafile->offpulse_rms = malloc(dummyl*sizeof(float));
      if(datafile->offpulse_rms == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error");
 return 0;
      }
      ret = fread(datafile->offpulse_rms, sizeof(float), dummyl, datafile->fptr_hdr);
      if(ret != dummyl) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
 return 0;
      }
    }
  }
  int readhistory;
  off_t filepos;
  filepos = ftello(datafile->fptr_hdr);
  readhistory = 1;
  ret = fread(identifier, 1, 7, datafile->fptr_hdr);
  if(ret != 7) {
    fflush(stdout);
    if(nohistory_expected == 0) {
      printwarning(verbose.debug, "WARNING readPSRSALSAHeader: No history table found in %s", datafile->filename);
    }
    readhistory = 0;
    fseeko(datafile->fptr_hdr, filepos, SEEK_SET);
  }else {
    identifier[7] = 0;
    if(strcmp(identifier, "HISTORY") != 0) {
      fflush(stdout);
      if(nohistory_expected == 0) {
 printwarning(verbose.debug, "WARNING readPSRSALSAHeader: No history table found in %s", datafile->filename);
      }
      if(verbose.debug) {
 printf("First bytes are '%c' '%c' and '%c'.\n", identifier[0], identifier[1], identifier[2]);
      }
      readhistory = 0;
      fseeko(datafile->fptr_hdr, filepos, SEEK_SET);
    }
  }
  if(nohistory_expected && readhistory) {
    printwarning(verbose.debug, "WARNING readPSRSALSAHeader: History table found in %s, while this was not expected.", datafile->filename);
  }
  if(readhistory) {
    datafile_history_entry_definition *curHistoryEntry;
    int curline, nrHistoryLines;
    ret = fread(&nrHistoryLines, sizeof(int), 1, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    curHistoryEntry = &(datafile->history);
    for(curline = 0; curline < nrHistoryLines; curline++) {
      ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
 return 0;
      }
      if(dummyi > 0) {
 curHistoryEntry->timestamp = malloc(dummyi+1);
 if(curHistoryEntry->timestamp == NULL) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading in %s", datafile->filename);
   return 0;
 }
 ret = fread(curHistoryEntry->timestamp, sizeof(char), dummyi, datafile->fptr_hdr);
 if(ret != dummyi) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
   return 0;
 }
 curHistoryEntry->timestamp[dummyi] = 0;
      }
      ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
 return 0;
      }
      if(dummyi > 0) {
 curHistoryEntry->user = malloc(dummyi+1);
 if(curHistoryEntry->user == NULL) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading in %s", datafile->filename);
   return 0;
 }
 ret = fread(curHistoryEntry->user, sizeof(char), dummyi, datafile->fptr_hdr);
 if(ret != dummyi) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
   return 0;
 }
 curHistoryEntry->user[dummyi] = 0;
      }
      ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
 return 0;
      }
      if(dummyi > 0) {
 curHistoryEntry->hostname = malloc(dummyi+1);
 if(curHistoryEntry->hostname == NULL) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading in %s", datafile->filename);
   return 0;
 }
 ret = fread(curHistoryEntry->hostname, sizeof(char), dummyi, datafile->fptr_hdr);
 if(ret != dummyi) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
   return 0;
 }
 curHistoryEntry->hostname[dummyi] = 0;
      }
      ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
 return 0;
      }
      if(dummyi > 0) {
 curHistoryEntry->cmd = malloc(dummyi+1);
 if(curHistoryEntry->cmd == NULL) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading in %s", datafile->filename);
   return 0;
 }
 ret = fread(curHistoryEntry->cmd, sizeof(char), dummyi, datafile->fptr_hdr);
 if(ret != dummyi) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
   return 0;
 }
 curHistoryEntry->cmd[dummyi] = 0;
      }
      if(version > 3) {
 ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
 if(ret != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
   return 0;
 }
 if(dummyi > 0) {
   curHistoryEntry->notes = malloc(dummyi+1);
   if(curHistoryEntry->notes == NULL) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading in %s", datafile->filename);
     return 0;
   }
   ret = fread(curHistoryEntry->notes, sizeof(char), dummyi, datafile->fptr_hdr);
   if(ret != dummyi) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
     return 0;
   }
   curHistoryEntry->notes[dummyi] = 0;
 }
      }
      if(curline < nrHistoryLines - 1) {
 curHistoryEntry->nextEntry = malloc(sizeof(datafile_history_entry_definition));
 if(curHistoryEntry->nextEntry == NULL) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading in %s", datafile->filename);
   return 0;
 }
 curHistoryEntry = curHistoryEntry->nextEntry;
 curHistoryEntry->timestamp = NULL;
 curHistoryEntry->cmd = NULL;
 curHistoryEntry->user = NULL;
 curHistoryEntry->hostname = NULL;
 curHistoryEntry->notes = NULL;
 curHistoryEntry->nextEntry = NULL;
      }
    }
  }
  datafile->datastart = ftell(datafile->fptr_hdr);
  free(txt);
  return 1;
}
int readPulsePSRSALSAData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  long long filepos;
  size_t ret;
  filepos = datafile.NrBins*(polarization+datafile.NrPols*(freq+pulsenr*datafile.NrFreqChan))+binnr;
  filepos *= sizeof(float);
  filepos += datafile.datastart;
  fseeko(datafile.fptr, filepos, SEEK_SET);
  ret = fread(pulse, sizeof(float), nrSamples, datafile.fptr);
  if(ret != nrSamples) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPulsePSRSALSAData: File read failed.");
    return 0;
  }
  return 1;
}
int writePulsePSRSALSAData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  long long filepos;
  size_t ret;
  filepos = datafile.NrBins*(polarization+datafile.NrPols*(freq+pulsenr*datafile.NrFreqChan))+binnr;
  filepos *= sizeof(float);
  filepos += datafile.datastart;
  fseeko(datafile.fptr, filepos, SEEK_SET);
  ret = fwrite(pulse, sizeof(float), nrSamples, datafile.fptr);
  if(ret != nrSamples) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePulsePSRSALSAData: File write failed.");
    return 0;
  }
  return 1;
}
int readPSRSALSAfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long n, f, p;
  size_t ret;
  if(verbose.verbose) {
    printf("Start reading PSRSALSA binary file\n");
  }
  fseeko(datafile.fptr, datafile.datastart, SEEK_SET);
  for(n = 0; n < datafile.NrSubints; n++) {
    for(f = 0; f < datafile.NrFreqChan; f++) {
      for(p = 0; p < datafile.NrPols; p++) {
 if(verbose.verbose && verbose.nocounters == 0)
   printf("  Progress reading PSRSALSA binary file (%.1f%%)\r", 100.0*(p+(f+n*datafile.NrFreqChan)*datafile.NrPols)/(float)(datafile.NrSubints*datafile.NrFreqChan*datafile.NrPols));
 ret = fread(&data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))], sizeof(float), datafile.NrBins, datafile.fptr);
 if(ret != datafile.NrBins) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRSALSAfile: File read failed.");
   return 0;
 }
      }
    }
  }
  if(verbose.verbose) printf("  Reading is done.                                \n");
  return 1;
}
int writePSRSALSAfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long n, f, p;
  size_t ret;
  fseeko(datafile.fptr, datafile.datastart, SEEK_SET);
  for(n = 0; n < datafile.NrSubints; n++) {
    for(f = 0; f < datafile.NrFreqChan; f++) {
      for(p = 0; p < datafile.NrPols; p++) {
 if(verbose.verbose && verbose.nocounters == 0)
   printf("  Progress writing PSRSALSA binary file (%.1f%%)\r", 100.0*(p+(f+n*datafile.NrFreqChan)*datafile.NrPols)/(float)(datafile.NrSubints*datafile.NrFreqChan*datafile.NrPols));
 ret = fwrite(&data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))], sizeof(float), datafile.NrBins, datafile.fptr);
 if(ret != datafile.NrBins) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writePSRSALSAfile: File write failed.");
   return 0;
 }
      }
    }
  }
  if(verbose.verbose) printf("  Writing is done.                                   \n");
  return 1;
}
int set_string_PSRData(char **dest, char *source, verbose_definition verbose)
{
  int i;
  if(*dest != NULL) {
    free(*dest);
    *dest = NULL;
  }
  if(source != NULL) {
    *dest = malloc(strlen(source)+1);
    if(*dest == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR set_string_PSRData: Memory allocation error.");
      return 0;
    }
    strcpy(*dest, source);
    for(i = 0; i < 2; i++) {
      if(strlen(*dest) > 0) {
 if((*dest)[strlen(*dest)-1] == '\n' || (*dest)[strlen(*dest)-1] == '\r') {
   (*dest)[strlen(*dest)-1] = 0;
 }
      }
    }
  }else {
    *dest = malloc(1);
    if(*dest == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR set_string_PSRData: Memory allocation error.");
      return 0;
    }
    (*dest)[0] = 0;
  }
  return 1;
}
int set_filename_PSRData(datafile_definition *datafile_dest, char *filename, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->filename), filename, verbose);
}
int set_psrname_PSRData(datafile_definition *datafile_dest, char *psrname, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->psrname), psrname, verbose);
}
int set_observatory_PSRData(datafile_definition *datafile_dest, char *observatory, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->observatory), observatory, verbose);
}
int set_institute_PSRData(datafile_definition *datafile_dest, char *institute, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->institute), institute, verbose);
}
int set_instrument_PSRData(datafile_definition *datafile_dest, char *instrument, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->instrument), instrument, verbose);
}
int set_scanID_PSRData(datafile_definition *datafile_dest, char *scanID, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->scanID), scanID, verbose);
}
int set_observer_PSRData(datafile_definition *datafile_dest, char *observer, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->observer), observer, verbose);
}
int set_projectID_PSRData(datafile_definition *datafile_dest, char *projectID, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->projectID), projectID, verbose);
}
int guessPSRData_format(char *filename, int noerror, verbose_definition verbose)
{
  FILE *fin;
  char *txt;
  int i, indent, ok;
  float fdummy;
  txt = malloc(10001);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR guessPSRData_format: Cannot allocate memory");
    return -3;
  }
  if(verbose.verbose) {
    for(indent = 0; indent < verbose.indent; indent++)
      printf(" ");
    printf("Trying to determine data format of file '%s'\n", filename);
    if(verbose.debug) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      printf("  Opening '%s'\n", filename);
    }
  }
  fin = fopen(filename, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR guessPSRData_format: opening file '%s' failed.", filename);
    free(txt);
    return -1;
  }
  i = fread(txt, 1, 3, fin);
  txt[3] = 0;
  if(i != 3) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR guessPSRData_format: reading file '%s' failed.", filename);
    free(txt);
    return -1;
  }
  if(verbose.debug) {
    printf("  First three bytes have values %d, %d and %d.\n", txt[0], txt[1], txt[2]);
    printf("  This corresponds to characters '%c', '%c' and '%c'.\n", txt[0], txt[1], txt[2]);
  }
  if(strcmp(txt, "PSR") == 0) {
    if(verbose.debug) {
      if(verbose.verbose) {
 for(indent = 0; indent < verbose.indent; indent++)
   printf(" ");
      }
      printf("  file might be in PSRSALSA binary format.\n");
    }
    rewind(fin);
    i = fread(txt, 1, 12, fin);
    txt[12] = 0;
    if(i != 12) {
      fflush(stdout);
      if(verbose.debug) {
 if(verbose.verbose) {
   for(indent = 0; indent < verbose.indent; indent++)
     printf(" ");
   printf("  file is not in PSRSALSA binary format.\n");
 }
      }
      rewind(fin);
      i = fread(txt, 1, 3, fin);
      txt[3] = 0;
    }else {
      if(strcmp(txt, "PSRSALSAdump") == 0) {
 if(verbose.verbose) {
   for(indent = 0; indent < verbose.indent; indent++)
     printf(" ");
   printf("  file is in PSRSALSA binary format.\n");
 }
 fclose(fin);
 free(txt);
 return PSRSALSA_BINARY_format;
      }else {
 fflush(stdout);
 if(verbose.debug) {
   for(indent = 0; indent < verbose.indent; indent++)
     printf(" ");
   printf("  file is not in PSRSALSA binary format.\n");
 }
 rewind(fin);
 i = fread(txt, 1, 3, fin);
 txt[3] = 0;
      }
    }
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    printf("  file is not in PSRSALSA binary format.\n");
  }
  if(strcmp(txt, "DPC") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      printf("  file is in WSRT PuMa 1 format.\n");
    }
    fclose(fin);
    free(txt);
    return PUMA_format;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    printf("  file is not in WSRT PuMa 1 format.\n");
  }
  if(strcmp(txt, "Fil") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    if(verbose.verbose)
      printf("  file is in PSRCHIVE ASCII format.\n");
    fclose(fin);
    free(txt);
    return PSRCHIVE_ASCII_format;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    printf("  file is not PSRCHIVE ASCII format.\n");
  }
  if(strcmp(txt, "EPN") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    if(verbose.verbose)
      printf("  file is in EPN format.\n");
    fclose(fin);
    free(txt);
    return EPN_format;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    printf("  file is not EPN format.\n");
  }
  if(txt[0] == 12 && txt[1] == 0 && txt[2] == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    if(verbose.verbose)
      printf("  file is in SIGPROC format.\n");
    fclose(fin);
    free(txt);
    return SIGPROC_format;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    printf("  file is not SIGPROC format.\n");
  }
  if(strcmp(txt, "SIM") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    if(verbose.verbose)
      printf("  file is in PSRFITS format.\n");
    fclose(fin);
    free(txt);
    return FITS_format;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    printf("  file is not PSRFITS format.\n");
  }
  if(strcmp(txt, "#pp") == 0
     ) {
    if(verbose.verbose && verbose.debug) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    if(verbose.debug)
      printf("  File might be in PPOL or PPOLSHORT format.\n");
    if(fgets(txt, 10000, fin) == NULL) {
      printf("ERROR guessPSRData_format: Cannot read first line\n");
      free(txt);
      return 0;
    }
    if(fgets(txt, 10000, fin) == NULL) {
      printf("ERROR guessPSRData_format: Cannot read second line\n");
      free(txt);
      return 0;
    }
    ok = 0;
    i = sscanf(txt, "%f %f %f", &fdummy, &fdummy, &fdummy);
    if(i == 3) {
      ok = 1;
      i = sscanf(txt, "%f %f %f %f %f %f %f %f %f", &fdummy, &fdummy, &fdummy, &fdummy, &fdummy, &fdummy, &fdummy, &fdummy, &fdummy);
      if(i == 9)
 ok = 2;
    }
    fclose(fin);
    if(ok == 2) {
      if(verbose.verbose) {
 for(indent = 0; indent < verbose.indent; indent++)
   printf(" ");
      }
      if(verbose.verbose) printf("  file is in PPOL format.\n");
      free(txt);
      return PPOL_format;
    }else if(ok == 1) {
      if(verbose.verbose) {
 for(indent = 0; indent < verbose.indent; indent++)
   printf(" ");
      }
      if(verbose.verbose) printf("  file is in PPOLSHORT format.\n");
      free(txt);
      return PPOL_SHORT_format;
    }else if(verbose.debug) {
      if(verbose.verbose) {
 for(indent = 0; indent < verbose.indent; indent++)
   printf(" ");
      }
      printf("  file is not in PPOL or PPOLSHORT format.\n");
    }
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    printf("  file is not in PPOL or PPOLSHORT format.\n");
  }
  int istimer = 1;
  if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    printf("  check if in TIMER format.\n");
  }
  if(fseek(fin, 32+32+8, SEEK_SET) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR guessPSRData_format: Seek in file '%s' failed.", filename);
    free(txt);
    return -1;
  }
  i = fread(&fdummy, sizeof(float), 1, fin);
  if(i != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR guessPSRData_format: reading file '%s' failed.", filename);
    free(txt);
    return -1;
  }
  if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
    }
    printf("    got version = %f\n", fdummy);
  }
  if(fabs(fdummy) < 1e-10) {
    i = fread(&fdummy, sizeof(float), 1, fin);
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR guessPSRData_format: reading file '%s' failed.", filename);
      free(txt);
      return -1;
    }
    if(verbose.debug) {
      if(verbose.verbose) {
 for(indent = 0; indent < verbose.indent; indent++)
   printf(" ");
      }
      printf("    got minor version = %f\n", fdummy);
    }
    if(fabs(fdummy) < 1e-10) {
      if(fseek(fin, 2*sizeof(int), SEEK_CUR) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR guessPSRData_format: Seek in file '%s' failed.", filename);
 free(txt);
 return -1;
      }
      i = fread(txt, 1, 16, fin);
      if(i != 16) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR guessPSRData_format: reading file '%s' failed.", filename);
 free(txt);
 return -1;
      }
      txt[16] = 0;
      if(verbose.debug) {
 if(verbose.verbose) {
   for(indent = 0; indent < verbose.indent; indent++)
     printf(" ");
 }
 printf("    got date = '%s'\n", txt);
      }
      if(txt[2] != '-' || txt[5] != '-') {
 istimer = 0;
      }
    }else {
      istimer = 0;
    }
  }else {
    istimer = 0;
  }
  fflush(stdout);
  fclose(fin);
  free(txt);
  if(istimer) {
    printerror(verbose.debug, "ERROR guessPSRData_format: file '%s' appears to be in TIMER format.", filename);
    printerror(verbose.debug, "ERROR guessPSRData_format: If you have psrchive, you can use 'pam -a psrfits -e fits %s' to convert to PSRFITS.", filename);
    return -2;
  }else {
    if(noerror == 0) {
      printerror(verbose.debug, "ERROR guessPSRData_format: determining file type of '%s' failed.", filename);
    }else {
      if(verbose.verbose) {
 for(indent = 0; indent < verbose.indent; indent++)
   printf(" ");
 printf("  file is in an undetermined format.\n");
      }
    }
  }
  return 0;
}
int openPSRData(datafile_definition *datafile, char *filename, int format, int enable_write, int read_in_memory, int nowarnings, verbose_definition verbose)
{
  int status = 0, iomode, i;
  char open_mode[100], filename2[1000];
  verbose_definition verbose2;
  if(enable_write == 0)
    cleanPSRData(datafile, verbose);
  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Opening file '%s'\n", filename);
  }
  if(set_filename_PSRData(datafile, filename, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR openPSRData: Making copy of data parameters failed.");
    return 0;
  }
  datafile->opened_flag = 0;
  datafile->format = format;
  if(enable_write) {
    datafile->enable_write_flag = 1;
  }
  if(format <= 0) {
    if(access(filename, F_OK) != 0) {
      fflush(stdout);
      if(enable_write == 0) {
 printerror(verbose.debug, "ERROR openPSRData: Cannot open file %s for reading.", filename);
      }else {
 printerror(verbose.debug, "ERROR openPSRData: File %s requested to be opened in write mode without file format being specified. File doesn't exist yet, so file format cannot be determined from existing file.", filename);
      }
      return 0;
    }
    format = guessPSRData_format(filename, 0, verbose2);
  }
  if(format == -1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR openPSRData:\n  Error reading file '%s' - file format couldn't be determined.", filename);
    return 0;
  }
  if(format == -2) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR openPSRData:\n  The file type of '%s' isn't supported for reading in.", filename);
    return 0;
  }
  if(format == 0 || format == -3) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR openPSRData:\n  Error determining file type of '%s'. It is either not supported, or you might need to specify the format on the command line.", filename);
    return 0;
  }
  datafile->format = format;
  open_mode[0] = 0;
  if(enable_write) {
    strcat(open_mode, "w+");
    if(verbose.debug) {
      if(verbose.verbose) {
 for(i = 0; i <= verbose2.indent; i++)
   printf(" ");
      }
      printf("Opening file '%s' for writing\n", filename);
    }
    if(format == PSRCHIVE_ASCII_format) {
      if(verbose.debug) {
 if(verbose.verbose) {
   for(i = 0; i <= verbose2.indent; i++)
     printf(" ");
 }
 printf("  The write commands will be buffered and executed once the file is being closed.\n");
      }
      datafile->dumpOnClose = 1;
    }else {
      datafile->dumpOnClose = 0;
    }
  }else {
    strcat(open_mode, "r");
    if(verbose.debug) {
      if(verbose.verbose) {
 for(i = 0; i <= verbose2.indent; i++)
   printf(" ");
      }
      printf("Opening file '%s' for reading\n", filename);
    }
  }
  if(format == PUMA_format || format == PSRSALSA_BINARY_format ||
     format == SIGPROC_format) {
    strcat(open_mode, "b");
    datafile->fptr = fopen(filename, open_mode);
    if(datafile->fptr == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR openPSRData: Cannot open file %s: %s", filename, strerror(errno));
      datafile->opened_flag = 0;
    }else {
      datafile->opened_flag = 1;
    }
  }else if(format == PSRCHIVE_ASCII_format || format == PPOL_format || format == PPOL_SHORT_format || format == SIGPROC_ASCII_format || format == EPN_format
       ) {
    datafile->fptr = fopen(filename, open_mode);
    if(datafile->fptr == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR openPSRData: Cannot open file %s: %s", filename, strerror(errno));
      datafile->opened_flag = 0;
    }else {
      datafile->opened_flag = 1;
    }
  }else if(format == FITS_format) {
    if(enable_write) {
      iomode = READWRITE;
      sprintf(filename2, "!%s", filename);
      if (!fits_create_file(&(datafile->fits_fptr), filename2, &status)) {
 if(verbose.debug) {
   if(verbose.verbose) {
     for(i = 0; i < verbose.indent; i++)
       printf(" ");
   }
   printf("  '%s' opened\n", filename);
 }
 datafile->opened_flag = 1;
      }else {
 datafile->opened_flag = 0;
      }
    }else {
      strcpy(filename2, filename);
      iomode = READONLY;
      if (!fits_open_file(&(datafile->fits_fptr), filename2, iomode, &status)) {
 if(verbose.debug) {
   if(verbose.verbose) {
     for(i = 0; i < verbose.indent; i++)
       printf(" ");
   }
   printf("  '%s' opened\n", filename);
 }
 datafile->opened_flag = 1;
      }else {
 datafile->opened_flag = 0;
      }
    }
    if(status) {
      fflush(stdout);
      fits_report_error(stderr, status);
    }
  }
    datafile->fptr_hdr = datafile->fptr;
  switch(format) {
  case PUMA_format: datafile->datastart = 4504; break;
  default:
    datafile->datastart = 0;
    break;
  }
  if(read_in_memory && datafile->opened_flag) {
    if(readHeaderPSRData(datafile, 0, nowarnings, verbose2)) {
      if(datafile->NrPols != 0) {
 long datasize = datafile->NrSubints*datafile->NrBins*datafile->NrPols*datafile->NrFreqChan*sizeof(float);
 datafile->data = (float *)malloc(datasize);
 if(datafile->data == NULL) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR openPSRData: Cannot allocate memory (data=%ld bytes=%.3fGB).", datasize, datasize/1073741824.0);
   closePSRData(datafile, 0, 0, verbose2);
   return 0;
 }
      }
      if(readPSRData(datafile, datafile->data, verbose2)) {
 closePSRData(datafile, 1, 1, verbose2);
 datafile->format = MEMORY_format;
 datafile->opened_flag = 1;
      }else {
 fflush(stdout);
 printerror(verbose.debug, "ERROR openPSRData: Cannot read data.");
 closePSRData(datafile, 0, 0, verbose2);
 return 0;
      }
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR openPSRData: Cannot read header.");
      closePSRData(datafile, 0, 0, verbose2);
      return 0;
    }
  }
  return datafile->opened_flag;
}
void closeHistoryPSRData(datafile_definition *datafile, int remove_last_entry_only)
{
  datafile_history_entry_definition *history_ptr, *history_ptr_next;
  int firsthistoryline;
  history_ptr = &(datafile->history);
  firsthistoryline = 1;
  do {
    history_ptr_next = history_ptr->nextEntry;
    if(history_ptr_next != NULL) {
      if(history_ptr_next->nextEntry == NULL) {
 history_ptr->nextEntry = NULL;
      }
    }
    if(remove_last_entry_only == 0 || history_ptr_next == NULL) {
      if(history_ptr->timestamp != NULL) {
 free(history_ptr->timestamp);
 history_ptr->timestamp = NULL;
      }
      if(history_ptr->cmd != NULL) {
 free(history_ptr->cmd);
 history_ptr->cmd = NULL;
      }
      if(history_ptr->user != NULL) {
 free(history_ptr->user);
 history_ptr->user = NULL;
      }
      if(history_ptr->hostname != NULL) {
 free(history_ptr->hostname);
 history_ptr->hostname = NULL;
      }
      if(history_ptr->notes != NULL) {
 free(history_ptr->notes);
 history_ptr->notes = NULL;
      }
      history_ptr->nextEntry = NULL;
      if(firsthistoryline == 0) {
 free(history_ptr);
      }
    }
    firsthistoryline = 0;
    history_ptr = history_ptr_next;
  }while(history_ptr_next != NULL);
}
int closePSRData(datafile_definition *datafile, int perserve_header, int perserve_data, verbose_definition verbose)
{
  int indent;
  int status = 0;
  if(verbose.debug) {
    printf("Closing file '%s'\n", datafile->filename);
  }
  if(datafile->opened_flag) {
    if(datafile->dumpOnClose) {
      if(verbose.debug) {
 printf("  - Dumping data to file before closing\n");
      }
      if(verbose.verbose) {
 for(indent = 0; indent < verbose.indent; indent++)
   printf(" ");
 printf("Writing buffer to file %s before it is closed\n", datafile->filename);
      }
      if(datafile->data == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR closePSRData: Although output is buffered, no memory is allocated?");
 return 1;
      }
      if(writePSRData(datafile, datafile->data, verbose) != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR closePSRData: Writing of buffered data failed.");
 return 1;
      }
    }
    if(datafile->format == FITS_format) {
      if(verbose.debug) {
 printf("  - Releasing FITS file pointer\n");
      }
      fits_close_file(datafile->fits_fptr, &status);
      if(status) {
 fflush(stdout);
 fits_report_error(stderr, status);
      }
      if(verbose.debug) {
 printf("  - Releasing memory related to scales/offsets/weights\n");
      }
      free(datafile->scales);
      free(datafile->offsets);
      free(datafile->weights);
      datafile->scales = NULL;
      datafile->offsets = NULL;
      datafile->weights = NULL;
    }else if(datafile->format != MEMORY_format){
      if(verbose.debug) {
 printf("  - Releasing file pointer\n");
      }
      fclose(datafile->fptr);
    }
    datafile->opened_flag = 0;
  }
  if(perserve_data == 0) {
    if(datafile->data != NULL) {
      if(verbose.debug) {
 printf("  - Releasing memory containing data\n");
      }
      free(datafile->data);
      datafile->data = NULL;
    }
  }
  if(perserve_header == 0) {
    if(verbose.debug) {
      printf("  - Releasing header related memory\n");
    }
    free(datafile->filename);
    datafile->filename = NULL;
    free(datafile->psrname);
    datafile->psrname = NULL;
    free(datafile->observatory);
    datafile->observatory = NULL;
    free(datafile->institute);
    datafile->institute = NULL;
    free(datafile->instrument);
    datafile->instrument = NULL;
    free(datafile->scanID);
    datafile->scanID = NULL;
    free(datafile->observer);
    datafile->observer = NULL;
    free(datafile->projectID);
    datafile->projectID = NULL;
    if(datafile->tsamp_list != NULL) {
      free(datafile->tsamp_list);
      datafile->tsamp_list = NULL;
    }
    if(datafile->tsub_list != NULL) {
      if(verbose.debug) {
 printf("  - Releasing tsub_list\n");
      }
      free(datafile->tsub_list);
      datafile->tsub_list = NULL;
    }
    if(datafile->freqlabel_list != NULL) {
      free(datafile->freqlabel_list);
      datafile->freqlabel_list = NULL;
    }
    if(datafile->offpulse_rms != NULL) {
      free(datafile->offpulse_rms);
      datafile->offpulse_rms = NULL;
    }
    closeHistoryPSRData(datafile, 0);
  }else if(perserve_header == 2) {
    closeHistoryPSRData(datafile, 1);
  }
  return status;
}
static char * internal_gentype_string_undefined = "Not set";
static char * internal_gentype_string_profile = "Profile";
static char * internal_gentype_string_pulsestack = "Pulsestack";
static char * internal_gentype_string_subints = "Subints";
static char * internal_gentype_string_searchmode = "Search mode";
static char * internal_gentype_string_bandpass = "Bandpass";
static char * internal_gentype_string_dynspec = "Dynamic spectrum";
static char * internal_gentype_string_polcalib = "Pol. calibrator";
static char * internal_gentype_string_lrfs = "LRFS";
static char * internal_gentype_string_2dfs = "2DFS";
static char * internal_gentype_string_s2dfsp3 = "S2DFS P3 map";
static char * internal_gentype_string_s2dfsp2 = "S2DFS P2 map";
static char * internal_gentype_string_p3fold = "P3 fold";
static char * internal_gentype_string_polarmap = "Polar map";
static char * internal_gentype_string_lrcc = "LRCC";
static char * internal_gentype_string_lrac = "LRAC";
static char * internal_gentype_string_gen_lr_dist = "Generic long. resolved distr.";
static char * internal_gentype_string_gen_lr_cumdist = "Generic long. resolved cum. distr.";
static char * internal_gentype_string_padist = "PA distr.";
static char * internal_gentype_string_elldist = "Ell distr.";
static char * internal_gentype_string_recmodel = "Receiver model";
static char * internal_gentype_string_recmodel2 = "Receiver model with chi^2";
static char * internal_gentype_string_rmmap = "RM map";
static char * internal_gentype_string_penergy = "penergy output";
static char * internal_gentype_string_hrfs_unfolded = "HRFS (unfolded)";
static char * internal_gentype_string_hrfs = "HRFS";
static char * internal_gentype_string_bug = "BUG, undefined????";
char *returnGenType_str(int gentype)
{
  switch(gentype) {
  case GENTYPE_UNDEFINED: return internal_gentype_string_undefined; break;
  case GENTYPE_PROFILE: return internal_gentype_string_profile; break;
  case GENTYPE_PULSESTACK: return internal_gentype_string_pulsestack; break;
  case GENTYPE_SUBINTEGRATIONS: return internal_gentype_string_subints; break;
  case GENTYPE_SEARCHMODE: return internal_gentype_string_searchmode; break;
  case GENTYPE_BANDPASS: return internal_gentype_string_bandpass; break;
  case GENTYPE_DYNAMICSPECTRUM: return internal_gentype_string_dynspec; break;
  case GENTYPE_POLNCAL: return internal_gentype_string_polcalib; break;
  case GENTYPE_LRFS: return internal_gentype_string_lrfs; break;
  case GENTYPE_2DFS: return internal_gentype_string_2dfs; break;
  case GENTYPE_S2DFSP3: return internal_gentype_string_s2dfsp3; break;
  case GENTYPE_S2DFSP2: return internal_gentype_string_s2dfsp2; break;
  case GENTYPE_P3FOLD: return internal_gentype_string_p3fold; break;
  case GENTYPE_POLARMAP: return internal_gentype_string_polarmap; break;
  case GENTYPE_LRCC: return internal_gentype_string_lrcc; break;
  case GENTYPE_LRAC: return internal_gentype_string_lrac; break;
  case GENTYPE_GEN_LR_DIST: return internal_gentype_string_gen_lr_dist; break;
  case GENTYPE_GEN_LR_CUMDIST: return internal_gentype_string_gen_lr_cumdist; break;
  case GENTYPE_PADIST: return internal_gentype_string_padist; break;
  case GENTYPE_ELLDIST: return internal_gentype_string_elldist; break;
  case GENTYPE_RECEIVERMODEL: return internal_gentype_string_recmodel; break;
  case GENTYPE_RECEIVERMODEL2: return internal_gentype_string_recmodel2; break;
  case GENTYPE_RMMAP: return internal_gentype_string_rmmap; break;
  case GENTYPE_PENERGY: return internal_gentype_string_penergy; break;
  case GENTYPE_HRFS_UNFOLDED: return internal_gentype_string_hrfs_unfolded; break;
  case GENTYPE_HRFS: return internal_gentype_string_hrfs; break;
  default: return internal_gentype_string_bug; break;
  }
}
void printGenType(int gentype, FILE *destination)
{
  fprintf(destination, "%s", returnGenType_str(gentype));
}
static char * internal_format_string_PuMa = "PuMa";
static char * internal_format_string_EPN = "EPN";
static char * internal_format_string_PPOL = "PAswing";
static char * internal_format_string_PPOL_short = "PAswing (short)";
static char * internal_format_string_PSRCHIVEAscii = "PSRCHIVE ascii";
static char * internal_format_string_PSRFITS = "PSRfits";
static char * internal_format_string_SIGPROC = "Sigproc";
static char * internal_format_string_SIGPROCAscii = "Sigproc (ascii)";
static char * internal_format_string_PSRSALSA = "PSRSALSA binary";
static char * internal_format_string_Memory = "Loaded in RAM";
static char * internal_format_string_bug = "BUG, undefined????";
char *returnFileFormat_str(int format)
{
  switch(format) {
  case PUMA_format: return internal_format_string_PuMa; break;
  case EPN_format: return internal_format_string_EPN; break;
  case SIGPROC_format: return internal_format_string_SIGPROC; break;
  case PPOL_format: return internal_format_string_PPOL; break;
  case PPOL_SHORT_format: return internal_format_string_PPOL_short; break;
  case SIGPROC_ASCII_format: return internal_format_string_SIGPROCAscii; break;
  case PSRCHIVE_ASCII_format: return internal_format_string_PSRCHIVEAscii; break;
  case FITS_format: return internal_format_string_PSRFITS; break;
  case PSRSALSA_BINARY_format: return internal_format_string_PSRSALSA; break;
  case MEMORY_format: return internal_format_string_Memory; break;
  default: return internal_format_string_bug; break;
  }
}
void printHeaderPSRData(datafile_definition datafile, int update, verbose_definition verbose)
{
  int i;
  char txt[1000];
  for(i = 0; i < verbose.indent; i++) printf(" ");
  if(update == 0)
    printf("========================= Dump of header parameters =========================\n");
  else
    printf("===================== Dump of updated header parameters =====================\n");
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("Pulsar=%s  ", datafile.psrname);
  printf("Ra=");
  converthms_string(txt, datafile.ra*12.0/M_PI, 2, 2);
  printf("%s ", txt);
  if(verbose.debug) {
    printf("= %f rad = %f deg = %f hours \n", datafile.ra, datafile.ra*180.0/M_PI, datafile.ra*12.0/M_PI);
    for(i = 0; i < verbose.indent; i++) printf(" ");
  }
  printf("Dec=");
  converthms_string(txt, datafile.dec*180.0/M_PI, 2, 3);
  printf("%s  ", txt);
  if(verbose.debug) {
    printf(" = %f rad = %f deg\n", datafile.dec, datafile.dec*180.0/M_PI);
    for(i = 0; i < verbose.indent; i++) printf(" ");
  }
  printf("\n");
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("tobs=%.1f sec", get_tobs(datafile, verbose));
  if(datafile.tsubMode == TSUBMODE_FIXEDTSUB) {
    printf(" (fixed subint length)\n");
  }else if(datafile.tsubMode == TSUBMODE_TSUBLIST) {
    printf(" (variable subint length)\n");
  }else {
    printf(" (unknown subint length)\n");
  }
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("DM=%f  RM=%f", datafile.dm, datafile.rm);
  printf(" RefFreq=");
  if((datafile.freq_ref < -0.9 && datafile.freq_ref >= -1.1) || (datafile.freq_ref > 0.99e10 && datafile.freq_ref < 1.01e10)){
    printf("infinity");
  }else if(datafile.freq_ref >= 0) {
    if(verbose.debug == 0)
      printf("%f MHz", datafile.freq_ref);
    else
      printf("%.9e MHz", datafile.freq_ref);
  }else {
    printf("Unknown");
  }
  printf("\n");
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("de-dispersed=");
  if(datafile.isDeDisp == 1)
    printf("yes");
  else if(datafile.isDeDisp == 0)
    printf("no");
  else
    printf("unknown");
  printf(" de-Faraday=");
  if(datafile.isDeFarad == 1)
    printf("yes");
  else if(datafile.isDeFarad == 0)
    printf("no");
  else
    printf("unknown");
  printf(" de-par. angle=");
  if(datafile.isDePar == 1)
    printf("yes");
  else if(datafile.isDePar == 0)
    printf("no");
  else
    printf("unknown");
  printf(" de-baselined=");
  if(datafile.isDebase == 1)
    printf("yes\n");
  else if(datafile.isDebase == 0)
    printf("no\n");
  else
    printf("unknown\n");
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("NrSubints=%ld NrBins=%ld NrPols=%ld NrFreqChan=%ld NrBits=%d\n",datafile.NrSubints,datafile.NrBins,datafile.NrPols, datafile.NrFreqChan,datafile.NrBits);
  for(i = 0; i < verbose.indent; i++) printf(" ");
  mjd2dateString(datafile.mjd_start, txt, 0, 1, " ");
  if(datafile.isFolded == 1) {
    double period;
    if(get_period(datafile, 0, &period, verbose) == 2) {
      printerror(verbose.debug, "ERROR printHeaderPSRData: Cannot obtain period");
      exit(0);
    }
    if(verbose.debug == 0)
      printf("period=%lf sec", period);
    else
      printf("period=%.9e sec", period);
  }else {
    printf("Not folded");
  }
  if(datafile.tsampMode == TSAMPMODE_LONGITUDELIST) {
    printf(" SampTime=longitude array");
  }else {
    if(verbose.debug == 0)
      printf(" SampTime=%f sec", get_tsamp(datafile, 0, verbose));
    else
      printf(" SampTime=%.9e sec", get_tsamp(datafile, 0, verbose));
  }
  if(verbose.debug == 0)
    printf(" mjd=%.9Lf (%s)\n", datafile.mjd_start, txt);
  else
    printf(" mjd=%.9Lf (%s)\n", datafile.mjd_start, txt);
  for(i = 0; i < verbose.indent; i++) printf(" ");
  if(datafile.freqMode == FREQMODE_UNKNOWN) {
    printwarning(verbose.debug, "OBSERVING FREQUENCY IS NOT SET");
  }else {
    double bw, chanbw, freq_cent;
    bw = get_bandwidth(datafile, verbose);
    if(get_channelbandwidth(datafile, &chanbw, verbose) == 0) {
      printerror(verbose.debug, "ERROR printHeaderPSRData: Cannot obtain channel bandwidth");
      exit(0);
    }
    freq_cent = get_centre_frequency(datafile, verbose);
    if(verbose.debug == 0)
      printf("freq_cent=%f MHz bw=%f MHz channelbw=%f MHz\n", freq_cent, bw, chanbw);
    else
      printf("freq_cent=%.9e MHz bw=%.9e MHz channelbw=%.9e MHz\n", freq_cent, bw, chanbw);
  }
  for(i = 0; i < verbose.indent; i++) printf(" ");
  if(verbose.debug) {
    printf("observatory=%s  ITRF (X,Y,Z)=(%lf,%lf,%lf) m\n", datafile.observatory, datafile.telescope_X, datafile.telescope_Y, datafile.telescope_Z);
    for(i = 0; i < verbose.indent; i++) printf(" ");
    printf("ITRF derived geocentric long=%lf deg lat=%lf deg\n", observatory_long_geocentric(datafile)*180.0/M_PI, observatory_lat_geocentric(datafile)*180.0/M_PI);
    for(i = 0; i < verbose.indent; i++) printf(" ");
    double longitude, latitude, height;
    tempo2_ITRF_to_GRS80(datafile.telescope_X, datafile.telescope_Y, datafile.telescope_Z, &longitude, &latitude, &height);
    printf("GRS80 derived geodetic  long=%lf deg lat=%lf deg height=%lf m\n", longitude*180.0/M_PI, latitude*180.0/M_PI, height);
  }else {
    printf("observatory=%s  long=%lf deg lat=%lf deg (geodetic derived)\n", datafile.observatory, observatory_long_geodetic(datafile)*180.0/M_PI, observatory_lat_geodetic(datafile)*180.0/M_PI);
  }
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("projectID=%s  observer=%s  institute=%s  instrument=%s\n", datafile.projectID, datafile.observer, datafile.institute, datafile.instrument);
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("feedtype=%d ", datafile.feedtype);
  if(datafile.feedtype == FEEDTYPE_LINEAR) {
    printf("(positive linear)");
  }else if(datafile.feedtype == FEEDTYPE_INV_LINEAR) {
    printf("(negative linear)");
  }else if(datafile.feedtype == FEEDTYPE_CIRCULAR) {
    printf("(positive circular)");
  }else if(datafile.feedtype == FEEDTYPE_INV_CIRCULAR) {
    printf("(negative circular)");
  }else {
    printf("(Unknown)");
  }
  printf("  poltype=%d ", datafile.poltype);
  if(datafile.poltype == POLTYPE_UNKNOWN) {
    printf("(undefined)\n");
  }else if(datafile.poltype == POLTYPE_STOKES) {
    printf("(Stokes)\n");
  }else if(datafile.poltype == POLTYPE_COHERENCY) {
    printf("(coherency)\n");
  }else if(datafile.poltype == POLTYPE_ILVPAdPA) {
    printf("(I,L,V,Pa and error)\n");
  }else if(datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
    printf("(I,L,V,Pa+error,tot pol,ell+error)\n");
  }else if(datafile.poltype == POLTYPE_PAdPA) {
    printf("(Pa and error)\n");
  }else {
    printf("BUG!!!!\n");
  }
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("FileName=");
  if(datafile.filename == NULL)
    printf(" ");
  else
    printf("%s", datafile.filename);
  printf("  ScanID=%s  Header length = %ld bytes\n", datafile.scanID, (long)datafile.datastart);
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("GenType=");
  printGenType(datafile.gentype, stdout);
  if(datafile.xrangeset)
    printf("  xrange=%f %f", datafile.xrange[0], datafile.xrange[1]);
  if(datafile.yrangeset)
    printf("  yrange=%f %f", datafile.yrange[0], datafile.yrange[1]);
  printf(" rms values=");
  if(datafile.offpulse_rms != NULL)
    printf("YES\n");
  else
    printf("NO\n");
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("====================== End of dump of header parameters ======================\n");
}
int setITRFlocation_by_name(datafile_definition *datafile, char *observatory, verbose_definition verbose)
{
  int indent;
  if(strcasecmp(observatory, "PARKES") == 0 || strcasecmp(observatory, "PKS") == 0 || strcasecmp(observatory, "7") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this Parkes telescope data");
    }
    datafile->telescope_X = -4554231.5;
    datafile->telescope_Y = 2816759.1;
    datafile->telescope_Z = -3454036.3;
  }else if(strcasecmp(observatory, "JODRELL") == 0 || strcasecmp(observatory, "JB") == 0 || strcasecmp(observatory, "JBO") == 0 || strcasecmp(observatory, "LOVELL") == 0 || strcasecmp(observatory, "JBDFB") == 0 || strcasecmp(observatory, "JBDFB") == 0 || strcasecmp(observatory, "JBODFB") == 0 || strcasecmp(observatory, "8") == 0 || strcasecmp(observatory, "q") == 0 || strcasecmp(observatory, "JB_MKII") == 0 || strcasecmp(observatory, "JBMK2") == 0 || strcasecmp(observatory, "h") == 0 || strcasecmp(observatory, "JB42") == 0 || strcasecmp(observatory, "JB_42ft") == 0) {
    if(strcasecmp(observatory, "JODRELL") == 0 || strcasecmp(observatory, "JB") == 0 || strcasecmp(observatory, "JBO") == 0 || strcasecmp(observatory, "LOVELL") == 0 || strcasecmp(observatory, "JBDFB") == 0 || strcasecmp(observatory, "JBDFB") == 0 || strcasecmp(observatory, "JBODFB") == 0 || strcasecmp(observatory, "8") == 0 || strcasecmp(observatory, "q") == 0) {
      if(verbose.verbose) {
 for(indent = 0; indent < verbose.indent; indent++)
   printf(" ");
 fflush(stdout);
 printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Jodrell Bank (Lovell) data");
      }
      datafile->telescope_X = 3822252.643;
      datafile->telescope_Y = -153995.683;
      datafile->telescope_Z = 5086051.443;
    }else if(strcasecmp(observatory, "JB_MKII") == 0 || strcasecmp(observatory, "JBMK2") == 0 || strcasecmp(observatory, "h") == 0) {
      if(verbose.verbose) {
 for(indent = 0; indent < verbose.indent; indent++)
   printf(" ");
 fflush(stdout);
 printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Jodrell Bank (MKII) data");
      }
      datafile->telescope_X = 3822473.365;
      datafile->telescope_Y = -153692.318;
      datafile->telescope_Z = 5085851.303;
    }else if(strcasecmp(observatory, "JB42") == 0 || strcasecmp(observatory, "JB_42ft") == 0) {
      if(verbose.verbose) {
 for(indent = 0; indent < verbose.indent; indent++)
   printf(" ");
 fflush(stdout);
 printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Jodrell Bank (42 ft telescope) data");
      }
      datafile->telescope_X = 3822294.825;
      datafile->telescope_Y = -153862.275;
      datafile->telescope_Z = 5085987.071;
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR setITRFlocation_by_name: BUG!!!!!!");
      return 0;
    }
  }else if(strcasecmp(observatory, "WSRT") == 0 || strcasecmp(observatory, "i") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is the WSRT data");
    }
    datafile->telescope_X = 3828445.659;
    datafile->telescope_Y = 445223.600000;
    datafile->telescope_Z = 5064921.5677;
  }else if(strcasecmp(observatory, "GBT") == 0 || strcasecmp(observatory, "1") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is GBT data");
    }
    datafile->telescope_X = 882589.65;
    datafile->telescope_Y = -4924872.32;
    datafile->telescope_Z = 3943729.348;
  }else if(strcasecmp(observatory, "ARECIBO") == 0 || strcasecmp(observatory, "AO") == 0 || strcasecmp(observatory, "3") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Arecibo data");
    }
    datafile->telescope_X = 2390490.0;
    datafile->telescope_Y = -5564764.0;
    datafile->telescope_Z = 1994727.0;
  }else if(strcasecmp(observatory, "EFFELSBERG") == 0 || strcasecmp(observatory, "eff") == 0 || strcasecmp(observatory, "g") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Effelsberg data");
    }
    datafile->telescope_X = 4033949.5;
    datafile->telescope_Y = 486989.4;
    datafile->telescope_Z = 4900430.8;
  }else if(strcasecmp(observatory, "GMRT") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is GMRT data");
    }
    datafile->telescope_X = 1656342.30;
    datafile->telescope_Y = 5797947.77;
    datafile->telescope_Z = 2073243.16;
  }else if(strcasecmp(observatory, "NANCAY") == 0 || strcasecmp(observatory, "ncy") == 0 || strcasecmp(observatory, "NUPPI") == 0 || strcasecmp(observatory, "ncyobs") == 0 || strcasecmp(observatory, "OP") == 0 || strcasecmp(observatory, "obspm") == 0 || strcasecmp(observatory, "f") == 0 || strcasecmp(observatory, "w") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Nancay data");
    }
    datafile->telescope_X = 4324165.81;
    datafile->telescope_Y = 165927.11;
    datafile->telescope_Z = 4670132.83;
  }else if(strcasecmp(observatory, "HARTRAO") == 0 || strcasecmp(observatory, "Hartebeesthoek") == 0 || strcasecmp(observatory, "hart") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is HARTRAO data");
    }
    datafile->telescope_X = 5085442.780;
    datafile->telescope_Y = 2668263.483;
    datafile->telescope_Z = -2768697.034;
  }else if(strcasecmp(observatory, "NANSHAN") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Nanshan data");
    }
    datafile->telescope_X = -228310.702;
    datafile->telescope_Y = 4631922.905;
    datafile->telescope_Z = 4367064.059;
  }else if(strcasecmp(observatory, "LOFAR") == 0 || strcasecmp(observatory, "t") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is central LOFAR data");
    }
    datafile->telescope_X = 3826577.462;
    datafile->telescope_Y = 461022.624;
    datafile->telescope_Z = 5064892.526;
  }else if(strcasecmp(observatory, "DE601LBA") == 0 || strcasecmp(observatory, "DE601LBH") == 0 || strcasecmp(observatory, "EFlfrlba") == 0 || strcasecmp(observatory, "EFlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Effelsberg international LOFAR station data (low band antenna's)");
    }
    datafile->telescope_X = 4034038.635;
    datafile->telescope_Y = 487026.223;
    datafile->telescope_Z = 4900280.057;
  }else if(strcasecmp(observatory, "DE601HBA") == 0 || strcasecmp(observatory, "DE601") == 0 || strcasecmp(observatory, "EFlfrhba") == 0 || strcasecmp(observatory, "EFlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Effelsberg international LOFAR station data (high band antenna's)");
    }
    datafile->telescope_X = 4034101.901;
    datafile->telescope_Y = 487012.401;
    datafile->telescope_Z = 4900230.210;
  }else if(strcasecmp(observatory, "DE602LBA") == 0 || strcasecmp(observatory, "DE602LBH") == 0 || strcasecmp(observatory, "UWlfrlba") == 0 || strcasecmp(observatory, "UWlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Garching/Unterweilenbach international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 4152561.068;
    datafile->telescope_Y = 828868.725;
    datafile->telescope_Z = 4754356.878;
  }else if(strcasecmp(observatory, "DE602HBA") == 0 || strcasecmp(observatory, "DE602") == 0 || strcasecmp(observatory, "UWlfrhba") == 0 || strcasecmp(observatory, "UWlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Garching/Unterweilenbach international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 4152568.416;
    datafile->telescope_Y = 828788.802;
    datafile->telescope_Z = 4754361.926;
  }else if(strcasecmp(observatory, "DE603LBA") == 0 || strcasecmp(observatory, "DE603LBH") == 0 || strcasecmp(observatory, "TBlfrlba") == 0 || strcasecmp(observatory, "TBlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Tautenburg international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 3940285.328;
    datafile->telescope_Y = 816802.001;
    datafile->telescope_Z = 4932392.757;
  }else if(strcasecmp(observatory, "DE603HBA") == 0 || strcasecmp(observatory, "DE603") == 0 || strcasecmp(observatory, "TBlfrhba") == 0 || strcasecmp(observatory, "TBlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Tautenburg international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 3940296.126;
    datafile->telescope_Y = 816722.532;
    datafile->telescope_Z = 4932394.152;
  }else if(strcasecmp(observatory, "DE604LBA") == 0 || strcasecmp(observatory, "DE604LBH") == 0 || strcasecmp(observatory, "POlfrlba") == 0 || strcasecmp(observatory, "POlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Potsdam international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 3796327.609;
    datafile->telescope_Y = 877591.315;
    datafile->telescope_Z = 5032757.252;
  }else if(strcasecmp(observatory, "DE604HBA") == 0 || strcasecmp(observatory, "DE604") == 0 || strcasecmp(observatory, "POlfrhba") == 0 || strcasecmp(observatory, "POlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Potsdam international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 3796380.254;
    datafile->telescope_Y = 877613.809;
    datafile->telescope_Z = 5032712.272;
  }else if(strcasecmp(observatory, "DE605LBA") == 0 || strcasecmp(observatory, "DE605LBH") == 0 || strcasecmp(observatory, "JUlfrlba") == 0 || strcasecmp(observatory, "JUlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Julich international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 4005681.742;
    datafile->telescope_Y = 450968.282;
    datafile->telescope_Z = 4926457.670;
  }else if(strcasecmp(observatory, "DE605HBA") == 0 || strcasecmp(observatory, "DE605") == 0 || strcasecmp(observatory, "JUlfrhba") == 0 || strcasecmp(observatory, "JUlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Julich international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 4005681.407;
    datafile->telescope_Y = 450968.304;
    datafile->telescope_Z = 4926457.940;
  }else if(strcasecmp(observatory, "FR606LBA") == 0 || strcasecmp(observatory, "FR606LBH") == 0 || strcasecmp(observatory, "FRlfrlba") == 0 || strcasecmp(observatory, "FRlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Nancey international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 4323980.155;
    datafile->telescope_Y = 165608.408;
    datafile->telescope_Z = 4670302.803;
  }else if(strcasecmp(observatory, "FR606HBA") == 0 || strcasecmp(observatory, "FR606") == 0 || strcasecmp(observatory, "FRlfrhba") == 0 || strcasecmp(observatory, "FRlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Nancey international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 4324017.054;
    datafile->telescope_Y = 165545.160;
    datafile->telescope_Z = 4670271.072;
  }else if(strcasecmp(observatory, "SE607LBA") == 0 || strcasecmp(observatory, "SE607LBH") == 0 || strcasecmp(observatory, "ONlfrlba") == 0 || strcasecmp(observatory, "ONlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Onsala international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 3370287.366;
    datafile->telescope_Y = 712053.586;
    datafile->telescope_Z = 5349991.228;
  }else if(strcasecmp(observatory, "SE607HBA") == 0 || strcasecmp(observatory, "SE607") == 0 || strcasecmp(observatory, "ONlfrhba") == 0 || strcasecmp(observatory, "ONlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Onsala international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 3370272.092;
    datafile->telescope_Y = 712125.596;
    datafile->telescope_Z = 5349990.934;
  }else if(strcasecmp(observatory, "UK608LBA") == 0 || strcasecmp(observatory, "UK608LBH") == 0 || strcasecmp(observatory, "UKlfrlba") == 0 || strcasecmp(observatory, "UKlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Chilbolton international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 4008438.796;
    datafile->telescope_Y = -100310.064;
    datafile->telescope_Z = 4943735.554;
  }else if(strcasecmp(observatory, "UK608HBA") == 0 || strcasecmp(observatory, "UK608") == 0 || strcasecmp(observatory, "UKlfrhba") == 0 || strcasecmp(observatory, "UKlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Chilbolton international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 4008462.280;
    datafile->telescope_Y = -100376.948;
    datafile->telescope_Z = 4943716.600;
  }else if(strcasecmp(observatory, "FI609LBA") == 0 || strcasecmp(observatory, "FI609LBH") == 0 || strcasecmp(observatory, "Filfrlba") == 0 || strcasecmp(observatory, "Filfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Kilpisjarvi international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 2136833.225;
    datafile->telescope_Y = 810088.740;
    datafile->telescope_Z = 5935285.279;
  }else if(strcasecmp(observatory, "FI609HBA") == 0 || strcasecmp(observatory, "FI609") == 0 || strcasecmp(observatory, "Filfrhba") == 0 || strcasecmp(observatory, "Filfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Kilpisjarvi international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 2136819.1940;
    datafile->telescope_Y = 810039.5757;
    datafile->telescope_Z = 5935299.0536;
  }else if(strcasecmp(observatory, "KAT7") == 0 || strcasecmp(observatory, "k7") == 0 || strcasecmp(observatory, "KAT-7") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is KAT-7 data");
    }
    datafile->telescope_X = 5109943.1050;
    datafile->telescope_Y = 2003650.7359;
    datafile->telescope_Z = -3239908.3195;
  }else if(strcasecmp(observatory, "MEERKAT") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is MeerKAT data");
    }
    datafile->telescope_X = 5109360.133;
    datafile->telescope_Y = 2006852.586;
    datafile->telescope_Z = -3238948.127;
  }else if(strcasecmp(observatory, "LWA") == 0 || strcasecmp(observatory, "LWA1") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is LWA data");
    }
    datafile->telescope_X = -1602196.60;
    datafile->telescope_Y = -5042313.47;
    datafile->telescope_Z = 3553971.51;
  }else if(strcasecmp(observatory, "FAST") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is FAST data");
    }
    datafile->telescope_X = -1668557.0;
    datafile->telescope_Y = 5506838.0;
    datafile->telescope_Z = 2744934.0;
  }else {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "WARNING setITRFlocation_by_name: Cannot guess location of telescope '%s'", observatory);
    }
    return 0;
  }
  return 1;
}
void determineWeightsStat(datafile_definition *datafile)
{
  int firstnonzeroweight = 1;
  long n;
  datafile->weight_stats_zeroweightfound = 0;
  datafile->weight_stats_differentweights = 0;
  datafile->weight_stats_negativeweights = 0;
  if(datafile->weights == NULL) {
    datafile->weight_stats_weightvalue = 1.0;
    return;
  }
  for(n = 0; n < datafile->NrSubints*datafile->NrFreqChan; n++) {
    if(datafile->weights[n] < 0) {
      datafile->weight_stats_negativeweights = 1;
    }
    if(datafile->weights[n] == 0.0) {
      datafile->weight_stats_zeroweightfound = 1;
    }else {
      if(firstnonzeroweight) {
 datafile->weight_stats_weightvalue = datafile->weights[n];
 firstnonzeroweight = 0;
      }else if(datafile->weights[n] != datafile->weight_stats_weightvalue) {
 datafile->weight_stats_differentweights = 1;
      }
    }
  }
  datafile->weight_stats_set = 1;
}
int readHeaderPSRData(datafile_definition *datafile, int readnoscales, int nowarnings, verbose_definition verbose)
{
  int i, nowarnings2, ret = 0;
  verbose_definition verbose2;
  if(nowarnings == 2) {
    nowarnings2 = 1;
  }else {
    nowarnings2 = 0;
  }
  if(verbose.debug) {
    printf("Entering readHeaderPSRData()\n");
  }
  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;
  if(datafile->format == PSRSALSA_BINARY_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Reading PSRSALSA binary header\n");
    }
    ret = readPSRSALSAHeader(datafile, 0, verbose);
  }else if(datafile->format == PUMA_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Reading WSRT header\n");
    }
    ret = readWSRTHeader(datafile, verbose);
  }else if(datafile->format == PSRCHIVE_ASCII_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Reading PSRCHIVE ascii header\n");
    }
    ret = readPSRCHIVE_ASCIIHeader(datafile, verbose);
  }else if(datafile->format == EPN_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Reading EPN header\n");
    }
    ret = readEPNHeader(datafile, 1, verbose);
    if(ret == 1) {
      float scale, offset;
      ret = 0;
      ret = readEPNsubHeader(datafile, &scale, &offset, verbose);
    }
  }else if(datafile->format == FITS_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Reading PSRFITS header\n");
    }
    ret = readPSRFITSHeader(datafile, readnoscales, nowarnings2, verbose);
    if(verbose.debug) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Reading PSRFITS header done\n");
    }
  }else if(datafile->format == SIGPROC_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Reading sigproc header\n");
    }
    ret = readSigprocHeader(datafile, verbose);
 }else if(datafile->format == PPOL_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Reading ppol file header\n");
    }
    ret = readPPOLHeader(datafile, 1, verbose);
  }else if(datafile->format == PPOL_SHORT_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Reading short ppol file header\n");
    }
    ret = readPPOLHeader(datafile, 0, verbose);
 }else if(datafile->format == SIGPROC_ASCII_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("Reading sigproc ascii header\n");
    }
    ret = readSigprocASCIIHeader(datafile, verbose);
  }else if(datafile->format == MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHeaderPSRData (%s): The file is already read in, cannot read header.", datafile->filename);
    return 0;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Data type (%s) not implemented in readHeaderPSRData", datafile->filename, returnFileFormat_str(datafile->format));
    return 0;
  }
  if(ret == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Error reading header.", datafile->filename);
    return 0;
  }
  if(datafile->format != FITS_format) {
    datafile->datastart = ftell(datafile->fptr);
  }else {
    datafile->datastart = 0;
  }
  if(convert_if_uniform_frequency_spacing(datafile, nowarnings2, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Attempt to convert dataset in an equally spaced frequency channel format failed.", datafile->filename);
    return 0;
  }
  if(fabs(datafile->telescope_X) < 1e-6 && fabs(datafile->telescope_Y) < 1e-6 && fabs(datafile->telescope_Z) < 1e-6) {
    fflush(stdout);
    if(strlen(datafile->observatory) == 0) {
      if(nowarnings == 0) {
 printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): Telescope location not set", datafile->filename);
      }
    }else {
      if(nowarnings == 0) {
 printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): Telescope location not set, try to determine location from name '%s'", datafile->filename, datafile->observatory);
      }
      if(setITRFlocation_by_name(datafile, datafile->observatory, verbose2) == 0) {
 fflush(stdout);
 if(nowarnings == 0) {
   printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): Telescope location not recognized by name", datafile->filename);
 }
      }
    }
  }
  double period;
  if(datafile->isFolded == 1) {
    if(get_period(*datafile, 0, &period, verbose) == 2) {
      printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Cannot obtain period", datafile->filename);
      return 0;
    }
  }else {
    period = -1;
  }
  if(period < 0.001 && datafile->isFolded != 0) {
    fflush(stdout);
    if(nowarnings == 0) {
      printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The period does not appear to be set in the header. Consider using the -header option.", datafile->filename);
      if(verbose.debug) {
 printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): isFolded=%d\n", datafile->filename, datafile->isFolded);
      }
    }
  }
  if(datafile->tsampMode != TSAMPMODE_LONGITUDELIST) {
    if(datafile->tsampMode == TSAMPMODE_UNKNOWN) {
      fflush(stdout);
      if(nowarnings == 0) {
 printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The sampling time is unknown.", datafile->filename);
      }
    }else {
      if(get_tsamp(*datafile, 0, verbose) < 0.0000001 || get_tsamp(*datafile, 0, verbose) >= 100) {
 fflush(stdout);
 if(nowarnings == 0) {
   printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The sampling time does not appear to be set correctly in the header. Consider using the -header option.", datafile->filename);
 }
      }
    }
  }
  if(datafile->isFolded == 0 && datafile->gentype == GENTYPE_UNDEFINED)
    datafile->gentype = GENTYPE_SEARCHMODE;
  if(datafile->tsampMode != TSAMPMODE_LONGITUDELIST && datafile->tsampMode != TSAMPMODE_UNKNOWN) {
    double diff_nbin;
    double period;
    if(datafile->isFolded == 1) {
      if(get_period(*datafile, 0, &period, verbose) == 2) {
 printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Cannot obtain period", datafile->filename);
 return 0;
      }
    }else {
      period = -1;
    }
    if(get_tsamp(*datafile, 0, verbose) != 0)
      diff_nbin = period/get_tsamp(*datafile, 0, verbose) - datafile->NrBins;
    else
      diff_nbin = 12345;
    if(datafile->isFolded && datafile->gentype != GENTYPE_2DFS && datafile->gentype != GENTYPE_RECEIVERMODEL && datafile->gentype != GENTYPE_RECEIVERMODEL2) {
      if(verbose.debug) {
 fflush(stdout);
 fprintf(stderr, "readHeaderPSRData (%s): Check if full period is stored - Tsamp=%lf, NBIN=%ld: Period=%lf suggests nr of bins is off by %lf.\n", datafile->filename, get_tsamp(*datafile, 0, verbose), datafile->NrBins, period, diff_nbin);
      }
      if(diff_nbin < -0.5 || diff_nbin > 0.5) {
 fflush(stdout);
 if(nowarnings == 0) {
   printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The sampling and period suggest that not the whole rotational phase range is stored. If not correct, consider using the -header option.", datafile->filename);
 }
      }
    }
  }
  if(datafile->gentype == GENTYPE_PULSESTACK) {
    double tobs_expected;
    double period;
    if(get_period(*datafile, 0, &period, verbose) == 2) {
      printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Cannot obtain period", datafile->filename);
      return 0;
    }
    tobs_expected = datafile->NrSubints * period;
    if(fabs(get_tsub(*datafile, 0, verbose)) < 0.0001) {
      datafile->tsubMode = TSUBMODE_FIXEDTSUB;
      if(datafile->tsub_list != NULL)
 free(datafile->tsub_list);
      datafile->tsub_list = (double *)malloc(sizeof(double));
      if(datafile->tsub_list == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Memory allocation error.", datafile->filename);
 return 0;
      }
      datafile->tsub_list[0] = period;
    }else if(get_tobs(*datafile, verbose)/tobs_expected > 1.02 || get_tobs(*datafile, verbose)/tobs_expected < 0.98) {
      fflush(stdout);
      if(nowarnings == 0) {
 printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The period and nr of pulses appear to be incompatable with tobs given that these are thought to be single pulses. Gentype is set to undefined. Consider using the -header option to fix problem.", datafile->filename);
      }
      datafile->gentype = GENTYPE_UNDEFINED;
    }
  }
  if(datafile->freq_ref < -1.1) {
    if(datafile->isDeDisp || datafile->isDeFarad) {
      if(nowarnings == 0) {
 printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The reference frequency used in dedispersion/de-Faraday rotation is unknown.", datafile->filename);
      }
    }else {
      datafile->freq_ref = -1;
    }
  }
  if(datafile->freq_ref < -0.9 && datafile->freq_ref >= -1.1) {
    datafile->freq_ref = 1e10;
  }
  if(readnoscales == 0) {
    determineWeightsStat(datafile);
    if(datafile->weight_stats_negativeweights) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): Found negative weights in input file, so probably the file is corrupted. You might want to use the -absweights option.", datafile->filename);
    }
    if(datafile->weight_stats_zeroweightfound && verbose.verbose) {
      printf("FITS file contains zero-weighted data. By default this is applied, unless the -noweights option is used.\n");
    }
    if(datafile->weight_stats_differentweights == 0 && (verbose.debug || datafile->weight_stats_weightvalue != 1.0)) {
      printf("FITS file contains uniform weights with a value %f", datafile->weight_stats_weightvalue);
      if(datafile->weight_stats_zeroweightfound)
 printf(" (appart from zero weights).\n");
      else
 printf("\n");
    }
    if(datafile->weight_stats_differentweights) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): FITS file contains data with different weights. Make sure it is used as desired.", datafile->filename);
    }
  }else {
    datafile->weight_stats_set = 1;
    datafile->weight_stats_zeroweightfound = 0;
    datafile->weight_stats_differentweights = 0;
    datafile->weight_stats_negativeweights = 0;
    datafile->weight_stats_weightvalue = 1.0;
  }
  if(readHistoryPSRData(datafile, verbose2) == 0) {
    printwarning(verbose.debug, "WARNING: Reading history failed.");
  }
  if(verbose.verbose) {
    printHeaderPSRData(*datafile, 0, verbose2);
  }
  if(verbose.debug) {
    printf("Exiting readHeaderPSRData()\n");
  }
  return 1;
}
int writeHeaderPSRData(datafile_definition *datafile, int argc, char **argv, int cmdOnly, char *notes, verbose_definition verbose)
{
  int i, ret;
  if(datafile->format == PSRSALSA_BINARY_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      if(verbose.verbose) printf("Write PSRSALSA header.\n");
    }
    ret = writePSRSALSAHeader(datafile, verbose);
    if(verbose.debug) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      if(verbose.verbose) printf("Writing PSRSALSA header done.\n");
    }
  }else if(datafile->format == PUMA_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      if(verbose.verbose) printf("Write PuMa header.\n");
    }
    ret = writeWSRTHeader(*datafile, verbose);
    if(verbose.debug) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      if(verbose.verbose) printf("Writing PuMa header done.\n");
    }
  }else if(datafile->format == PSRCHIVE_ASCII_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
    }
    if(verbose.verbose) printf("Write PSRCHIVE_ASCII header.\n");
    ret = writePSRCHIVE_ASCIIHeader(*datafile, verbose);
  }else if(datafile->format == FITS_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
    }
    if(verbose.verbose) printf("Write PSRFITS header.\n");
    ret = writePSRFITSHeader(datafile, verbose);
    if(verbose.debug) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      if(verbose.verbose) printf("Writing PSRFITS header done.\n");
    }
  }else if(datafile->format == EPN_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
    }
    if(verbose.verbose) {
      printwarning(verbose.debug, "WARNING writeHeaderPSRData: Write EPN header not necessary, will when writing data.");
    }
    ret = 1;
  }else if(datafile->format == SIGPROC_ASCII_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
    }
    if(verbose.verbose) printf("Write Sigproc ascii header.\n");
    ret = writeSigprocASCIIHeader(*datafile, verbose);
  }else if(datafile->format == PPOL_format || datafile->format == PPOL_SHORT_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
    }
    if(verbose.verbose) printf("Write ppol file header.\n");
    ret = writePPOLHeader(*datafile, argc, argv, verbose);
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeHeaderPSRData: Writing of this type of header is not implemented");
    ret = 0;
  }
  if(ret != 0) {
    if(datafile->format == PSRSALSA_BINARY_format || datafile->format == FITS_format || datafile->format == PUMA_format) {
      writeHistoryPSRData(datafile, argc, argv, cmdOnly, notes, verbose);
    }
  }
  if(ret == 0) {
    printerror(verbose.debug, "ERROR writeHeaderPSRData: Writing header failed.");
  }
  return ret;
}
int readPulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  if(datafile->format == PSRSALSA_BINARY_format)
    return readPulsePSRSALSAData(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  else if(datafile->format == PUMA_format)
    return readPulseWSRTData(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse);
  else if(datafile->format == FITS_format)
    return readFITSpulse(datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  else if(datafile->format == PSRCHIVE_ASCII_format)
    return readPSRCHIVE_ASCIIfilepulse(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  else if(datafile->format == EPN_format)
    return readPulseEPNData(datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  else if(datafile->format == SIGPROC_format)
    return readPulseSigprocData(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  else if(datafile->format == MEMORY_format) {
    memcpy(pulse, &datafile->data[datafile->NrBins*(polarization+datafile->NrPols*(freq+pulsenr*datafile->NrFreqChan))+binnr], sizeof(float)*nrSamples);
    return 1;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPulsePSRData: Reading of this format is not implemented. Maybe converting the data in a different format will solve this issue.");
  }
  return 0;
}
int writePulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  if(pulsenr < 0 || binnr < 0 || freq < 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePulsePSRData: Parameters outside boundaries.");
    return 0;
  }
  if(datafile->format == MEMORY_format || datafile->dumpOnClose) {
    if(datafile->dumpOnClose && datafile->data == NULL) {
      long datasize = datafile->NrSubints*datafile->NrBins*datafile->NrPols*datafile->NrFreqChan*sizeof(float);
      datafile->data = (float *)malloc(datasize);
      if(datafile->data == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePulsePSRData: Cannot allocate memory (data=%ld bytes=%.3fGB).", datasize, datasize/1073741824.0);
 datafile->dumpOnClose = 0;
 verbose_definition verbose2;
 copyVerboseState(verbose, &verbose2);
 if(verbose.debug == 0)
   verbose2.verbose = 0;
 verbose2.nocounters = 1;
 closePSRData(datafile, 0, 0, verbose2);
 return 0;
      }else if(verbose.debug) {
 printf("DEBUG: Allocated %ld bytes of memory for memory buffering.\n", datasize);
      }
    }
    memcpy(&datafile->data[datafile->NrBins*(polarization+datafile->NrPols*(freq+pulsenr*datafile->NrFreqChan))+binnr], pulse, sizeof(float)*nrSamples);
  }else if(datafile->format == PSRSALSA_BINARY_format) {
    return writePulsePSRSALSAData(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  }else if(datafile->format == PUMA_format) {
    return writePulseWSRTData(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse);
  }else if(datafile->format == FITS_format) {
    return writeFITSpulse(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  }else if(datafile->format == PSRCHIVE_ASCII_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePulsePSRData: Writing out individual subintegrations is not implemented for ASCII formats.");
    return 0;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePulsePSRData: Writing in this format (%s) is not implemented", returnFileFormat_str(datafile->format));
    return 0;
  }
  return 1;
}
int readPSRData(datafile_definition *datafile, float *data, verbose_definition verbose)
{
  if(datafile->format == PSRSALSA_BINARY_format)
    return readPSRSALSAfile(*datafile, data, verbose);
  else if(datafile->format == PUMA_format)
    return readPuMafile(*datafile, data, verbose);
  else if(datafile->format == PSRCHIVE_ASCII_format)
    return readPSRCHIVE_ASCIIfile(*datafile, data, verbose);
  else if(datafile->format == EPN_format)
    return readEPNfile(datafile, data, verbose, -1);
  else if(datafile->format == FITS_format) {
    int ret;
    ret = readFITSfile(datafile, data, verbose);
    return ret;
  }else if(datafile->format == PPOL_format)
    return readPPOLfile(datafile, data, 1, 0, verbose);
  else if(datafile->format == PPOL_SHORT_format)
    return readPPOLfile(datafile, data, 0, 0, verbose);
  else if(datafile->format == SIGPROC_ASCII_format)
    return readSigprocASCIIfile(*datafile, data, verbose);
  else if(datafile->format == SIGPROC_format)
    return readSigprocfile(*datafile, data, verbose);
  else
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRData: Reading whole dataset is not supported for this type of data.");
  return 0;
}
int writePSRData(datafile_definition *datafile, float *data, verbose_definition verbose)
{
  if(verbose.verbose) printf("Writing %ld x %ld x %ld x %ld samples\n", datafile->NrSubints, datafile->NrFreqChan, datafile->NrBins, datafile->NrPols);
  datafile->dumpOnClose = 0;
  if(datafile->format == PSRSALSA_BINARY_format) {
    return writePSRSALSAfile(*datafile, data, verbose);
  }else if(datafile->format == PUMA_format) {
    return writePuMafile(*datafile, data, verbose);
  }else if(datafile->format == EPN_format) {
    return writeEPNfile(*datafile, data, verbose);
  }else if(datafile->format == PSRCHIVE_ASCII_format) {
    return writePSRCHIVE_ASCIIfile(*datafile, data, verbose);
  }else if(datafile->format == FITS_format) {
    return writeFITSfile(*datafile, data, verbose);
  }else if(datafile->format == PPOL_format) {
    return writePPOLfile(*datafile, data, 1, 0, 0, 0, verbose);
  }else if(datafile->format == PPOL_SHORT_format) {
    return writePPOLfile(*datafile, data, 0, 0, 0, 0, verbose);
  }else if(datafile->format == SIGPROC_ASCII_format) {
    return writeSigprocASCIIfile(*datafile, data, verbose);
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRData: Writing whole dataset is not supported for this type of data (%s).", returnFileFormat_str(datafile->format));
  }
  return 0;
}
int read_profilePSRData(datafile_definition datafile, float *profileI, int *zapMask, int polchan, verbose_definition verbose)
{
  return read_partprofilePSRData(datafile, profileI, zapMask, polchan, 0, datafile.NrSubints, verbose);
}
int read_partprofilePSRData(datafile_definition datafile, float *profileI, int *zapMask, int polchan, long nskip, long nread, verbose_definition verbose)
{
  int zap;
  long i, j, k;
  float *data;
  if(verbose.verbose) printf("Generating average pulse profile (polarization channel %d)\n", polchan);
  data = (float *)malloc(datafile.NrBins*sizeof(float));
  if(data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR read_partprofilePSRData: Cannot allocate memory.");
    return 0;
  }
  for(j = 0; j < datafile.NrBins; j++) {
    profileI[j] = 0;
  }
  for(k = 0; k < datafile.NrFreqChan; k++) {
    for(i = nskip; i < nskip+nread; i++) {
      if(readPulsePSRData(&datafile, i, polchan, k, 0, datafile.NrBins, data, verbose) == 0)
 return 0;
      zap = 0;
      if(zapMask != NULL) {
 if(zapMask[i] != 0)
   zap = 1;
      }
      if(zap == 0) {
 for(j = 0; j < datafile.NrBins; j++) {
   profileI[j] += data[j];
 }
      }
    }
  }
  free(data);
  return 1;
}
int read_rmsPSRData(datafile_definition datafile, float *rms, float *avrg, int *zapMask, pulselongitude_regions_definition *regions, int invert, int polchan, int freqchan, verbose_definition verbose)
{
  long i, j, k, freq0, freq1;
  float *data;
  double *rms_double, *avrg_double;
  int nrOffpulseBins, zap, zap2;
  data = (float *)malloc(datafile.NrBins*sizeof(float));
  rms_double = (double *)malloc(datafile.NrSubints*sizeof(double));
  avrg_double = (double *)malloc(datafile.NrSubints*sizeof(double));
  if(data == NULL || rms_double == NULL || avrg_double == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR read_rmsPSRData: Cannot allocate memory.");
    return 0;
  }
  for(i = 0; i < datafile.NrSubints; i++) {
    rms_double[i] = 0;
    avrg_double[i] = 0;
  }
  nrOffpulseBins = datafile.NrBins;
  if(regions != NULL) {
    for(j = 0; j < datafile.NrBins; j++) {
      if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0))
 nrOffpulseBins--;
    }
  }
  if(nrOffpulseBins == 0) {
    printerror(verbose.debug, "ERROR read_rmsPSRData: An off-pulse rms was requested, but everything is defined as onpulse.");
    return 0;
  }
  if(freqchan < 0) {
    freq0 = 0;
    freq1 = datafile.NrFreqChan;
    freqchan = 0;
  }else {
    freq0 = freqchan;
    freq1 = freqchan+1;
    freqchan = 1;
  }
  for(k = freq0; k < freq1; k++) {
    for(i = 0; i < datafile.NrSubints; i++) {
      if(readPulsePSRData(&datafile, i, polchan, k, 0, datafile.NrBins, data, verbose) == 0)
 return 0;
      zap = 0;
      if(zapMask != NULL) {
 if(zapMask[i] != 0)
   zap = 1;
      }
      for(j = 0; j < datafile.NrBins; j++) {
 zap2 = zap;
 if(regions != NULL) {
   if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0))
     zap2 = 1;
 }
 if(zap2 == 0) {
   rms_double[i] += data[j]*data[j];
   avrg_double[i] += data[j];
 }
      }
    }
  }
  if(freqchan == 0) {
    freqchan = datafile.NrFreqChan;
  }else {
    freqchan = 1;
  }
  for(i = 0; i < datafile.NrSubints; i++) {
    avrg_double[i] /= (double)(nrOffpulseBins*freqchan);
    rms_double[i] -= nrOffpulseBins*freqchan*avrg_double[i]*avrg_double[i];
    rms_double[i] /= (double)(nrOffpulseBins*freqchan);
    if(rms != NULL)
      rms[i] = sqrt(rms_double[i]);
    if(avrg != NULL)
      avrg[i] = avrg_double[i];
  }
  free(data);
  free(rms_double);
  free(avrg_double);
  return 1;
}
int PSRDataHeader_parse_commandline(datafile_definition *psrdata, int argc, char **argv, verbose_definition verbose)
{
  int i, j, ok;
  char identifier[100], value[100], txt[100];
  ok = 0;
  for(i = 1; i < argc - 1; i++) {
    if(strcmp(argv[i], "-header") == 0 || strcmp(argv[i], "-headerUFL") == 0) {
      ok = 1;
    }
  }
  if(ok && verbose.verbose) {
    printf("Changing header parameters:\n");
  }
  for(i = 1; i < argc - 1; i++) {
    if(strcmp(argv[i], "-header") == 0) {
      i++;
      j = sscanf(argv[i], "%s %s", identifier, value);
      if(j != 2) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline:  Cannot parse -header option.");
 printHeaderCommandlineOptions(stderr);
 return 0;
      }else {
 if(strcasecmp(identifier,"name") == 0 || strcasecmp(identifier, "psrname") == 0 || strcasecmp(identifier, "pulsar") == 0 || strcasecmp(identifier, "psr") == 0) {
   if(set_psrname_PSRData(psrdata, value, verbose) == 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Setting pulsar name failed.");
     return 0;
   }
   if(verbose.verbose) printf("  hdr.psrname = %s\n", psrdata->psrname);
 }else if(strcasecmp(identifier, "freq") == 0 || strcasecmp(identifier, "cfreq") == 0 || strcasecmp(identifier, "freq_cent") == 0) {
   if(psrdata->freqMode != FREQMODE_UNIFORM) {
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Centre frequency can only be changed in the header if the data has uniformly distributed channels.");
       return 0;
   }else {
     double freq;
     sscanf(value, "%lf", &freq);
     set_centre_frequency(psrdata, freq, verbose);
     if(verbose.verbose) printf("  hdr.freq_cent = %lf MHz\n", freq);
   }
 }else if(strcasecmp(identifier, "bw") == 0) {
   if(psrdata->freqMode != FREQMODE_UNIFORM) {
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Bandwidth can only be changed in the header if the data has uniformly distributed channels.");
       return 0;
   }else {
     double bw;
     sscanf(value, "%lf", &bw);
     if(set_bandwidth(psrdata, bw, verbose) == 0) {
       printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Bandwidth changing failed.");
       return 0;
     }
     if(verbose.verbose) printf("  hdr.bw = %f MHz\n", bw);
     double chanbw;
     chanbw = bw/(double)psrdata->NrFreqChan;
     if(verbose.verbose) printf("  hdr.channelbw = %f MHz\n", chanbw);
   }
 }else if(strcasecmp(identifier,"chbw") == 0 || strcasecmp(identifier,"chanbw") == 0 || strcasecmp(identifier,"channelbw") == 0) {
   if(psrdata->freqMode != FREQMODE_UNIFORM) {
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Bandwidth can only be changed in the header if the data has uniformly distributed samples.");
       return 0;
   }else {
     double chanbw;
     sscanf(value, "%lf", &chanbw);
     double bw;
     bw = chanbw*psrdata->NrFreqChan;
     if(set_bandwidth(psrdata, bw, verbose) == 0) {
       printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Bandwidth changing failed.");
       return 0;
     }
     if(verbose.verbose) printf("  hdr.bw = %f MHz\n", bw);
     if(verbose.verbose) printf("  hdr.channelbw = %f MHz\n", chanbw);
   }
 }else if(strcasecmp(identifier, "freqref") == 0 || strcasecmp(identifier, "freq_ref") == 0 || strcasecmp(identifier, "ref_freq") == 0 || strcasecmp(identifier, "reffreq") == 0) {
   psrdata->freq_ref = atof(value);
   if(verbose.verbose) printf("  hdr.freq_ref = %f MHz\n", psrdata->freq_ref);
   if(psrdata->freq_ref > -1.01 && psrdata->freq_ref < -0.99) {
     psrdata->freq_ref = 1e10;
   }
 }else if(strcasecmp(identifier, "p0") == 0 || strcasecmp(identifier, "P0") == 0 || strcasecmp(identifier, "period") == 0) {
   if(psrdata->isFolded != 1 || psrdata->foldMode != FOLDMODE_FIXEDPERIOD) {
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Period can only be changed in the header if there is a fixed period throughout the whole dataset.");
       return 0;
   }else {
     psrdata->fixedPeriod = atof(value);
     double period;
     if(get_period(*psrdata, 0, &period, verbose) == 2) {
       printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline (%s): Cannot obtain period", psrdata->filename);
       return 0;
     }
     if(verbose.verbose) printf("  hdr.period = %lf s\n", period);
   }
 }else if(strcasecmp(identifier, "dt") == 0 || strcasecmp(identifier, "tsamp") == 0 || strcasecmp(identifier, "samptime") == 0) {
   if(psrdata->tsampMode != TSAMPMODE_FIXEDTSAMP) {
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Sampling time can only be changed in the header if there is a fixed sampling time throughout the whole dataset.");
       return 0;
   }else {
     psrdata->fixedtsamp = atof(value);
     if(verbose.verbose) printf("  hdr.SampTime = %lf s\n", get_tsamp(*psrdata, 0, verbose));
   }
 }else if(strcasecmp(identifier, "tsub") == 0 || strcasecmp(identifier, "tsubint") == 0 || strcasecmp(identifier, "t_sub") == 0) {
   char *substring;
   int nrwords, ret, nsub;
   substring = pickWordFromString(argv[i], 2, &nrwords, 1, ' ', verbose);
   if(nrwords == 2) {
     ret = sscanf(substring, "%lf", &(psrdata->tsub_list[0]));
     if(ret != 1) {
       printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Cannot parse %s as a value.", substring);
     }
     psrdata->tsubMode = TSUBMODE_FIXEDTSUB;
     if(verbose.verbose) printf("  hdr.tsub = %lf s (fixed for each subint)\n", get_tsub(*psrdata, 0, verbose));
   }else {
     if(nrwords != 1+psrdata->NrSubints) {
       printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: When setting individual subint durations, provide a list of space separated durations which defines a duration for each subint (got %d durations, need %d values). Example: -header 'tsub 30.0 30.0 30.0' when there are three subints present.", nrwords-1, psrdata->NrSubints);
       return 0;
     }
     if(psrdata->tsub_list != NULL)
       free(psrdata->tsub_list);
     psrdata->tsub_list = (double *)malloc(psrdata->NrSubints*sizeof(double));
     if(psrdata->tsub_list == NULL) {
       printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Memory allocation error.");
       return 0;
     }
     psrdata->tsubMode = TSUBMODE_TSUBLIST;
     if(verbose.verbose) printf("  hdr.tsub = ");
     for(nsub = 0; nsub < psrdata->NrSubints; nsub++) {
       substring = pickWordFromString(argv[i], 2+nsub, &nrwords, 1, ' ', verbose);
       ret = sscanf(substring, "%lf", &(psrdata->tsub_list[nsub]));
       if(ret != 1) {
  printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Cannot parse %s as a value.", substring);
       }
       if(verbose.verbose) {
  if(nsub != 0)
    printf(",");
  printf("%lf", get_tsub(*psrdata, nsub, verbose));
       }
     }
     if(verbose.verbose) {
       printf(" sec\n");
     }
   }
 }else if(strcasecmp(identifier,"mjd") == 0) {
   psrdata->mjd_start = atof(value);
   if(verbose.verbose) printf("  hdr.mjd = %Lf\n", psrdata->mjd_start);
 }else if(strcasecmp(identifier,"length") == 0 || strcasecmp(identifier,"tobs") == 0 || strcasecmp(identifier,"dur") == 0) {
   psrdata->tsubMode = TSUBMODE_FIXEDTSUB;
   if(psrdata->tsub_list != NULL)
     free(psrdata->tsub_list);
   psrdata->tsub_list = (double *)malloc(sizeof(double));
   if(psrdata->tsub_list == NULL) {
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Memory allocation error.");
     return 0;
   }
   psrdata->tsub_list[0] = atof(value)/(double)psrdata->NrSubints;
   if(verbose.verbose) printf("  hdr.tobs = %lf\n", get_tobs(*psrdata, verbose));
   printwarning(verbose.debug, "WARNING PSRDataHeader_parse_commandline: Assuming equal subint lengths.");
 }else if(strcasecmp(identifier,"loc") == 0 || strcasecmp(identifier,"location") == 0) {
   int ret;
   ret = sscanf(value, "%lf,%lf,%lf", &(psrdata->telescope_X), &(psrdata->telescope_Y), &(psrdata->telescope_Z));
   if(ret != 3) {
     if(verbose.verbose) {
       printf("  Looking up ITRF coordinates for site '%s'\n", value);
     }
     if(setITRFlocation_by_name(psrdata, value, verbose) == 0) {
       printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: '%s' not recognized as a telescope location. Consider entering location as XVALUE,YVALUE,ZVALUE (without spaces)", value);
       return 0;
     }
   }
   if(verbose.verbose) printf("  hdr.telescope_X = %lf m\n", psrdata->telescope_X);
   if(verbose.verbose) printf("  hdr.telescope_Y = %lf m\n", psrdata->telescope_Y);
   if(verbose.verbose) printf("  hdr.telescope_Z = %lf m\n", psrdata->telescope_Z);
 }else if(strcasecmp(identifier,"locationGEO") == 0 || strcasecmp(identifier,"locGEO") == 0) {
   int ret;
   double longitude, latitude, height;
   ret = sscanf(value, "%lf,%lf,%lf", &longitude, &latitude, &height);
   if(ret != 3) {
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: In option %s of -header, '%s' expected to be of the form LONGITUDE,LATITUDE,HEIGHT (without spaces).", identifier, value);
     return 0;
   }
   if(verbose.verbose) {
     printf("  Derive ITRF X,Y,Z coordinates from geodetic GRS80 longitude=%f deg, latitude=%f deg, height=%f m\n", longitude, latitude, height);
   }
   longitude *= M_PI/180.0;
   latitude *= M_PI/180.0;
   tempo2_GRS80_to_ITRF(longitude, latitude, height, &(psrdata->telescope_X), &(psrdata->telescope_Y), &(psrdata->telescope_Z));
   if(verbose.verbose) printf("  hdr.telescope_X = %lf m\n", psrdata->telescope_X);
   if(verbose.verbose) printf("  hdr.telescope_Y = %lf m\n", psrdata->telescope_Y);
   if(verbose.verbose) printf("  hdr.telescope_Z = %lf m\n", psrdata->telescope_Z);
 }else if(strcasecmp(identifier,"scan") == 0 || strcasecmp(identifier,"scanid") == 0) {
   if(set_scanID_PSRData(psrdata, value, verbose) == 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Setting scan ID failed.");
     return 0;
   }
   if(verbose.verbose) printf("  hdr.ScanID = %s\n", psrdata->scanID);
 }else if(strcasecmp(identifier,"project") == 0 || strcasecmp(identifier,"projectid") == 0 || strcasecmp(identifier,"projid") == 0) {
   if(set_projectID_PSRData(psrdata, value, verbose) == 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Setting project ID failed.");
     return 0;
   }
   if(verbose.verbose) printf("  hdr.projectID = %s\n", psrdata->projectID);
 }else if(strcasecmp(identifier, "observer") == 0 || strcasecmp(identifier, "observers") == 0) {
   if(set_observer_PSRData(psrdata, value, verbose) == 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Setting name of observer(s) failed.");
     return 0;
   }
   if(verbose.verbose) printf("  hdr.observer = %s\n", psrdata->observer);
 }else if(strcasecmp(identifier,"observatory") == 0 || strcasecmp(identifier,"telescope") == 0) {
   if(set_observatory_PSRData(psrdata, value, verbose) == 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Setting observatory name failed.");
     return 0;
   }
   if(verbose.verbose) printf("  hdr.observatory = %s\n", psrdata->observatory);
   int found;
   found = 0;
   for(j = 1; j < argc - 1; j++) {
     sscanf(argv[j], "%s", txt);
     if(strcasecmp(txt,"loc") == 0 || strcasecmp(txt,"location") == 0 || strcasecmp(txt,"locationGEO") == 0 || strcasecmp(txt,"locGEO") == 0) {
       found = 1;
       break;
     }
   }
   if(found == 0) {
     printwarning(verbose.debug, "WARNING PSRDataHeader_parse_commandline: When changing the name of the telescope, you might want to consider to change the location as well.");
   }
 }else if(strcasecmp(identifier,"nrpulses") == 0 || strcasecmp(identifier,"npulses") == 0 || strcasecmp(identifier,"pulses") == 0 || strcasecmp(identifier,"nrsub") == 0 || strcasecmp(identifier,"nsub") == 0 || strcasecmp(identifier, "nsubint") == 0 || strcasecmp(identifier, "subints") == 0) {
   psrdata->NrSubints = atol(value);
   if(verbose.verbose) printf("  hdr.NrSubints = %ld\n", psrdata->NrSubints);
 }else if(strcasecmp(identifier,"nrbin") == 0 || strcasecmp(identifier,"nbin") == 0) {
   psrdata->NrBins = atol(value);
   if(verbose.verbose) printf("  hdr.NrBins = %ld\n", psrdata->NrBins);
 }else if(strcasecmp(identifier,"nrbits") == 0 || strcasecmp(identifier,"nbits") == 0) {
   psrdata->NrBits = atoi(value);
   if(verbose.verbose) printf("  hdr.NrBits = %d\n", psrdata->NrBits);
 }else if(strcasecmp(identifier,"nrchan") == 0 || strcasecmp(identifier,"nchan") == 0 || strcasecmp(identifier,"nrfreq") == 0 || strcasecmp(identifier,"nfreq") == 0 || strcasecmp(identifier,"nrfreqchan") == 0 || strcasecmp(identifier,"nfreqchan") == 0) {
   psrdata->NrFreqChan = atol(value);
   if(verbose.verbose) printf("  hdr.NrFreqChan = %ld\n", psrdata->NrFreqChan);
 }else if(strcasecmp(identifier,"nrpol") == 0 || strcasecmp(identifier,"npol") == 0 || strcasecmp(identifier,"nrpols") == 0 || strcasecmp(identifier,"npols") == 0) {
   psrdata->NrPols = atol(value);
   if(verbose.verbose) printf("  hdr.NrPols = %ld\n", psrdata->NrPols);
 }else if(strcasecmp(identifier,"fdtype") == 0 || strcasecmp(identifier,"fd_type") == 0 || strcasecmp(identifier,"feedtype") == 0) {
   psrdata->feedtype = atoi(value);
   if(verbose.verbose) printf("  hdr.feedtype = %d\n", psrdata->feedtype);
 }else if(strcasecmp(identifier,"ra") == 0) {
   if(strstr(value, ":") != NULL) {
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: In option %s of -header, a single value in degrees is expected.", identifier);
     return 0;
   }
   psrdata->ra = atof(value)*M_PI/180.0;
   if(verbose.verbose) printf("  hdr.ra = %f rad = %f deg = %f hours = ", psrdata->ra, psrdata->ra*180.0/M_PI, psrdata->ra*12.0/M_PI);
   converthms_string(txt, psrdata->ra*12.0/M_PI, 2, 2);
   printf("%s\n", txt);
 }else if(strcasecmp(identifier, "dec") == 0) {
   if(strstr(value, ":") != NULL) {
     printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: In option %s of -header, a single value in degrees is expected.", identifier);
     return 0;
   }
   psrdata->dec = atof(value)*M_PI/180.0;
   if(verbose.verbose) printf("  hdr.dec = %f rad = %f deg = ", psrdata->dec, psrdata->dec*180.0/M_PI);
   converthms_string(txt, psrdata->dec*180.0/M_PI, 2, 3);
   printf("%s\n", txt);
 }else if(strcasecmp(identifier,"poltype") == 0 || strcasecmp(identifier,"pol_type") == 0) {
   psrdata->poltype = atoi(value);
   if(verbose.verbose) printf("  hdr.poltype = %d\n", psrdata->poltype);
 }else if(strcasecmp(identifier,"dedisp") == 0 || strcasecmp(identifier,"dedispersed") == 0 || strcasecmp(identifier,"isdedisp") == 0 || strcasecmp(identifier,"isdedispersed") == 0) {
   psrdata->isDeDisp = atoi(value);
   if(verbose.verbose) printf("  hdr.isDeDisp = %d\n", psrdata->isDeDisp);
 }else if(strcasecmp(identifier,"defarad") == 0 || strcasecmp(identifier,"isdefarad") == 0) {
   psrdata->isDeFarad = atoi(value);
   if(verbose.verbose) printf("  hdr.isDeFarad = %d\n", psrdata->isDeFarad);
   if(fabs(psrdata->rm) < 1e-6) {
     fflush(stdout);
     printwarning(verbose.debug, "WARNING PSRDataHeader_parse_commandline: Be cautious with changing the de-Faraday rotation state while the RM is not defined.");
   }
 }else if(strcasecmp(identifier,"depar") == 0 || strcasecmp(identifier,"isdepar") == 0) {
   psrdata->isDePar = atoi(value);
   if(verbose.verbose) printf("  hdr.isDePar = %d\n", psrdata->isDePar);
 }else if(strcasecmp(identifier,"debase") == 0 || strcasecmp(identifier,"isdebase") == 0) {
   psrdata->isDebase = atoi(value);
   if(verbose.verbose) printf("  hdr.isDebase = %d\n", psrdata->isDebase);
 }else if(strcasecmp(identifier,"dm") == 0) {
   psrdata->dm = atof(value);
   if(verbose.verbose) printf("  hdr.dm = %f\n", psrdata->dm);
 }else if(strcasecmp(identifier,"rm") == 0) {
   psrdata->rm = atof(value);
   if(verbose.verbose) printf("  hdr.rm = %f\n", psrdata->rm);
 }else if(strcasecmp(identifier,"cableswap") == 0) {
   psrdata->cableSwap = atoi(value);
   if(verbose.verbose) printf("  hdr.cableswap = %d\n", psrdata->cableSwap);
 }else if(strcasecmp(identifier,"cableswapcor") == 0) {
   psrdata->cableSwapcor = atoi(value);
   if(verbose.verbose) printf("  hdr.cableswapcor = %d\n", psrdata->cableSwapcor);
 }else if(strcasecmp(identifier,"gentype") == 0 || strcasecmp(identifier,"type") == 0) {
   int orig_gentype = psrdata->gentype;
   sscanf(value, "%d", &(psrdata->gentype));
   if(verbose.verbose) printf("  hdr.gentype = %d\n", psrdata->gentype);
   if(orig_gentype == GENTYPE_SEARCHMODE) {
     psrdata->isFolded = 1;
     psrdata->foldMode = FOLDMODE_FIXEDPERIOD;
     if(psrdata->fixedPeriod <= 0)
       psrdata->fixedPeriod = 1.0;
     fflush(stdout);
     printwarning(verbose.debug, "WARNING PSRDataHeader_parse_commandline:  Period is not set. Use -header \"period value\" to set it to something appropriate");
   }
   if(psrdata->gentype == GENTYPE_SEARCHMODE) {
     psrdata->fixedPeriod = 0;
     psrdata->foldMode = FOLDMODE_UNKNOWN;
     psrdata->isFolded = 0;
   }
 }else if(strcasecmp(identifier,"yrange") == 0) {
   int ret;
   double value1, value2;
   ret = sscanf(value, "%lf,%lf", &value1, &value2);
   if(ret == 2) {
     psrdata->yrangeset = 1;
     psrdata->yrange[0] = value1;
     psrdata->yrange[1] = value2;
     if(verbose.verbose) printf("  hdr.yrange1 = %lf\n", psrdata->yrange[0]);
     if(verbose.verbose) printf("  hdr.yrange2 = %lf\n", psrdata->yrange[1]);
   }else {
     if(strcasecmp(value, "x") == 0 || strcasecmp(value, "undefined") == 0 || strcasecmp(value, "empty") == 0 || strcasecmp(value, "nothing") == 0) {
       psrdata->yrangeset = 0;
       if(verbose.verbose) printf("  hdr.yrange = undefined\n");
     }else {
       printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: In option %s of -header, '%s' expected to be of the form VALUE1,VALUE2 (without spaces) or the word UNDEFINED.", identifier, value);
       return 0;
     }
   }
 }else {
   fflush(stdout);
   printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline:  '%s' not recognized as a header parameter.", identifier);
   printHeaderCommandlineOptions(stderr);
   return 0;
 }
      }
    }else if(strcmp(argv[i], "-headerUFL") == 0) {
      force_uniform_frequency_spacing(psrdata, verbose);
    }
  }
  if(ok && verbose.verbose) {
    verbose_definition verbose2;
    copyVerboseState(verbose, &verbose2);
    verbose2.verbose = 0;
    printHeaderPSRData(*psrdata, 1, verbose2);
  }
  return 1;
}
void printHeaderCommandlineOptions(FILE *printdevice)
{
  fprintf(printdevice, "Valid options for the -header option are:\n");
  fprintf(printdevice, "  bw           band width.\n");
  fprintf(printdevice, "  chbw         Channel band width.\n");
  fprintf(printdevice, "  dec          Declination (single number, in degrees).\n");
  fprintf(printdevice, "  dm           Dispersion Measure.\n");
  fprintf(printdevice, "  dt           Sampling time.\n");
  fprintf(printdevice, "  fdtype       1=lin 2=circ, sign is handiness.\n");
  fprintf(printdevice, "  freq         Centre frequency (i.e. defines frequency labeling).\n");
  fprintf(printdevice, "  debase       Data is baseline-subtracted (0 or 1).\n");
  fprintf(printdevice, "  dedisp       Data is de-dispersed (0 or 1).\n");
  fprintf(printdevice, "  defarad      Data is de-faraday rotated (0 or 1).\n");
  fprintf(printdevice, "  depar        Data is parallactic angle corrected (0 or 1).\n");
  fprintf(printdevice, "  gentype      Identifies the type of data, use -gentypelist to see options.\n");
  fprintf(printdevice, "  length       The duration of the observation.\n");
  fprintf(printdevice, "  location     The ITRF location of telescope in meters, e.g.\n");
  fprintf(printdevice, "               -header \"loc 3822252.643,-153995.683,5086051.443\" for the Lovell.\n");
  fprintf(printdevice, "               or -header \"loc Lovell\" for recognized site names.\n");
  fprintf(printdevice, "  locationGEO  Derive ITRF location of telescope via geodetic GRS80 coordinates:\n");
  fprintf(printdevice, "               -header \"locGEO long,lat,height\" where longitude/latitude are.\n");
  fprintf(printdevice, "               in degrees, and height in meters.\n");
  fprintf(printdevice, "  mjd          Start MJD of the observation.\n");
  fprintf(printdevice, "  name         Name of the pulsar.\n");
  fprintf(printdevice, "  nbin         Nr of time-bins (use at own risk).\n");
  fprintf(printdevice, "  nbits        Nr of bits (use at own risk).\n");
  fprintf(printdevice, "  nchan        Nr of frequency channels (use at own risk).\n");
  fprintf(printdevice, "  npol         Nr of polarization channels (use at own risk).\n");
  fprintf(printdevice, "  nsub         Nr of sub integrations in observation (use at own risk).\n");
  fprintf(printdevice, "  observatory  Change the telescope (name only, not location).\n");
  fprintf(printdevice, "  p0           Period.\n");
  fprintf(printdevice, "  poltype      Polarization type (%d = undefined, %d=Stokes, %d=Coherency,\n", POLTYPE_UNKNOWN, POLTYPE_STOKES, POLTYPE_COHERENCY);
  fprintf(printdevice, "               %d=I,L,V,Pa+error, %d=Pa+error, %d=I,L,V,Pa+error,tot pol,ell+error).\n", POLTYPE_ILVPAdPA, POLTYPE_PAdPA, POLTYPE_ILVPAdPATEldEl);
  fprintf(printdevice, "  ra           Right ascension (single number, in degrees).\n");
  fprintf(printdevice, "  reffreq      Reference frequency (for dedispersion/de Faraday rotation).\n");
  fprintf(printdevice, "               -1=Infinite freq -2=Unknown\n");
  fprintf(printdevice, "  rm           Rotation Measure.\n");
  fprintf(printdevice, "  scan         Scan ID of the observation.\n");
  fprintf(printdevice, "  tsub         Followed by single number: set subint durations to this number\n");
  fprintf(printdevice, "               in seconds.\n");
  fprintf(printdevice, "               If followed by space separated numbers: set individual subint\n");
  fprintf(printdevice, "               durations\n");
  fprintf(printdevice, "  yrange       Set the yrange to the range provided. Specify this as\n");
  fprintf(printdevice, "               -header \"yrange VALUE1,VALUE2\". To unset the use of\n");
  fprintf(printdevice, "               this parameter use -header \"yrange undefined\".\n");
}
void printHeaderGentypeOptions(FILE *printdevice)
{
  fprintf(printdevice, "Valid options for the -header \"gentype number\" option are:\n");
  fprintf(printdevice, "  %3d           Undefined.\n", GENTYPE_UNDEFINED);
  fprintf(printdevice, "  %3d           Profile.\n", GENTYPE_PROFILE);
  fprintf(printdevice, "  %3d           Pulse stack.\n", GENTYPE_PULSESTACK);
  fprintf(printdevice, "  %3d           Subintegration data.\n", GENTYPE_SUBINTEGRATIONS);
  fprintf(printdevice, "  %3d           Seach mode data (not folded).\n", GENTYPE_SEARCHMODE);
  fprintf(printdevice, "  %3d           Bandpass.\n", GENTYPE_BANDPASS);
  fprintf(printdevice, "  %3d           Dynamic spectrum.\n", GENTYPE_DYNAMICSPECTRUM);
  fprintf(printdevice, "  %3d           Polarization calibration signal (noise diode).\n", GENTYPE_POLNCAL);
  fprintf(printdevice, "  %3d           LRFS.\n", GENTYPE_LRFS);
  fprintf(printdevice, "  %3d           2DFS.\n", GENTYPE_2DFS);
  fprintf(printdevice, "  %3d           Sliding 2DFS (P3).\n", GENTYPE_S2DFSP3);
  fprintf(printdevice, "  %3d           Sliding 2DFS (P2).\n", GENTYPE_S2DFSP2);
  fprintf(printdevice, "  %3d           P3 fold.\n", GENTYPE_P3FOLD);
  fprintf(printdevice, "  %3d           Receiver model.\n", GENTYPE_RECEIVERMODEL);
  fprintf(printdevice, "  %3d           Receiver model with chi^2 and nfree.\n", GENTYPE_RECEIVERMODEL2);
}
int make_clone(datafile_definition original, datafile_definition *clone, verbose_definition verbose)
{
  int debug;
  debug = 0;
  cleanPSRData(clone, verbose);
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(debug, "ERROR make_clone: Only works if data is already loaded into memory.");
    return 0;
  }
  copy_params_PSRData(original, clone, verbose);
  clone->data = (float *)malloc((original.NrBins)*(original.NrPols)*(original.NrFreqChan)*(original.NrSubints)*sizeof(float));
  if(clone->data == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR make_clone: Memory allocation error.");
    return 0;
  }
  memcpy(clone->data, original.data, (original.NrBins)*(original.NrPols)*(original.NrFreqChan)*(original.NrSubints)*sizeof(float));
  if(original.offpulse_rms != NULL) {
    clone->offpulse_rms = (float *)malloc(original.NrPols*original.NrFreqChan*original.NrSubints*sizeof(float));
    if(clone->offpulse_rms == NULL) {
      fflush(stdout);
      printerror(debug, "ERROR make_clone: Memory allocation error.");
      return 0;
    }
    memcpy(clone->offpulse_rms, original.offpulse_rms, original.NrPols*original.NrFreqChan*original.NrSubints*sizeof(float));
  }
  clone->opened_flag = 1;
  return 1;
}
void swap_orig_clone(datafile_definition *original, datafile_definition *clone, verbose_definition verbose)
{
  verbose_definition verbose2;
  copyVerboseState(verbose, &verbose2);
  verbose2.verbose = 0;
  verbose2.nocounters = 1;
  closePSRData(original, 0, 0, verbose);
  memmove(original, clone, sizeof(datafile_definition));
}
int writeHistoryPSRData(datafile_definition *datafile, int argc, char **argv, int cmdOnly, char *notes, verbose_definition verbose)
{
  int ret;
  char txt[10000], txt2[1000], *username_ptr, hostname[1000];
  time_t curtime;
  datafile_history_entry_definition *curHistoryEntry;
  if(verbose.verbose)
    fprintf(stdout, "Writing history\n");
  if(argc > 0) {
    curHistoryEntry = &(datafile->history);
    if(curHistoryEntry->timestamp != NULL || curHistoryEntry->cmd != NULL || curHistoryEntry->user != NULL || curHistoryEntry->hostname != NULL || curHistoryEntry->notes != NULL || curHistoryEntry->nextEntry != NULL) {
      while(curHistoryEntry->nextEntry != NULL) {
 curHistoryEntry = curHistoryEntry->nextEntry;
      }
      curHistoryEntry->nextEntry = malloc(sizeof(datafile_history_entry_definition));
      if(curHistoryEntry->nextEntry == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error");
 return 0;
      }
      curHistoryEntry = curHistoryEntry->nextEntry;
      curHistoryEntry->timestamp = NULL;
      curHistoryEntry->cmd = NULL;
      curHistoryEntry->user = NULL;
      curHistoryEntry->hostname = NULL;
      curHistoryEntry->notes = NULL;
      curHistoryEntry->nextEntry = NULL;
    }
    constructCommandLineString(txt, 10000, argc, argv, verbose);
    curHistoryEntry->cmd = malloc(strlen(txt)+1);
    if(curHistoryEntry->cmd == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error");
      return 0;
    }
    strcpy(curHistoryEntry->cmd, txt);
    if(cmdOnly == 0) {
      curtime = time(NULL);
      strcpy(txt2, asctime(gmtime(&curtime)));
      if(txt2[strlen(txt2)-1] == '\n')
 txt2[strlen(txt2)-1] = 0;
      if(txt2[strlen(txt2)-1] == '\r')
 txt2[strlen(txt2)-1] = 0;
      if(txt2[strlen(txt2)-1] == '\n')
 txt2[strlen(txt2)-1] = 0;
      curHistoryEntry->timestamp = malloc(strlen(txt2)+1);
      if(curHistoryEntry->timestamp == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error");
 return 0;
      }
      strcpy(curHistoryEntry->timestamp, txt2);
      if(getUsername(&username_ptr, verbose) == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "writeHistoryPSRData: Cannot identify user.");
 username_ptr = malloc(8);
 if(username_ptr == NULL) {
   printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error");
   return 0;
 }
 sprintf(username_ptr, "Unknown");
      }
      curHistoryEntry->user = malloc(strlen(username_ptr)+1);
      if(curHistoryEntry->user == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error");
 return 0;
      }
      strcpy(curHistoryEntry->user, username_ptr);
      free(username_ptr);
      getMachinename(hostname, 1000, verbose);
      curHistoryEntry->hostname = malloc(strlen(hostname)+1);
      if(curHistoryEntry->hostname == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error");
 return 0;
      }
      strcpy(curHistoryEntry->hostname, hostname);
      if(notes != NULL) {
 curHistoryEntry->notes = malloc(strlen(notes)+1);
 if(curHistoryEntry->notes == NULL) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error");
   return 0;
 }
 strcpy(curHistoryEntry->notes, notes);
      }
    }
  }
  if(datafile->format == PSRSALSA_BINARY_format) {
    ret = writeHistoryPSRSALSA(datafile, verbose);
  }else if(datafile->format == FITS_format) {
    ret = writeHistoryFITS(*datafile, verbose);
  }else if(datafile->format == PUMA_format) {
    ret = writeHistoryPuma(*datafile, verbose);
  }else {
    fflush(stdout);
    printwarning(verbose.debug, "writeHistoryPSRData: Writing a history is not supported in this file format.");
    ret = 0;
  }
  if(verbose.verbose)
    fprintf(stdout, "  done\n");
  return ret;
}
int readHistoryPSRData(datafile_definition *datafile, verbose_definition verbose)
{
  int ret, indent, doread;
  if(verbose.verbose) {
    for(indent = 0; indent < verbose.indent; indent++)
      printf(" ");
    fprintf(stdout, "Reading history\n");
  }
  doread = 0;
  if(datafile->format == FITS_format) {
    ret = readHistoryFITS(datafile, verbose);
    doread = 1;
  }else if(datafile->format == PUMA_format) {
    ret = readHistoryPuma(datafile, verbose);
    doread = 1;
  }else {
    ret = -1;
  }
  if(doread) {
    if(ret == 0) {
      printwarning(verbose.debug, "WARNING: Reading history failed.");
    }
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fprintf(stdout, "  done\n");
    }
  }else {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      fprintf(stdout, "  skipped\n");
    }
  }
  return ret;
}
int showHistory(datafile_definition datafile, verbose_definition verbose)
{
  datafile_history_entry_definition *curHistoryEntry;
  long rownr;
  int indent;
  curHistoryEntry = &(datafile.history);
  rownr = 0;
  do {
    int foundsomething;
    foundsomething = 0;
    if(curHistoryEntry->timestamp != NULL || curHistoryEntry->cmd != NULL || curHistoryEntry->user != NULL || curHistoryEntry->hostname != NULL || curHistoryEntry->nextEntry != NULL) {
      for(indent = 0; indent < verbose.indent; indent++)
 printf(" ");
      printf("History line %ld: ", rownr+1);
      if(curHistoryEntry->timestamp != NULL)
 printf("%s ", curHistoryEntry->timestamp);
      if(curHistoryEntry->user != NULL)
 printf("%s ", curHistoryEntry->user);
      if(curHistoryEntry->hostname != NULL)
 printf("%s ", curHistoryEntry->hostname);
      if(curHistoryEntry->cmd != NULL)
 printf("%s\n", curHistoryEntry->cmd);
      foundsomething = 1;
    }
    if(curHistoryEntry->notes != NULL) {
      int oktoprint;
      oktoprint = 0;
      if(strlen(curHistoryEntry->notes) > 1) {
 oktoprint = 1;
      }else if(strlen(curHistoryEntry->notes) == 1) {
 if(curHistoryEntry->notes[0] != ' ') {
   oktoprint = 1;
 }
      }
      if(oktoprint) {
 for(indent = 0; indent < verbose.indent; indent++)
   printf(" ");
 printf("Notes line %ld:   ", rownr+1);
 if(curHistoryEntry->notes != NULL)
   printf("'%s'\n", curHistoryEntry->notes);
      }
      foundsomething = 1;
    }
    if(foundsomething) {
      curHistoryEntry = curHistoryEntry->nextEntry;
      rownr++;
    }
  }while(curHistoryEntry != NULL && rownr != 0);
  if(rownr > 0)
    return 1;
  else
    return 0;
}
int skipallhashedlines(datafile_definition *datafile)
{
  long pos;
  char c;
  int ret, debug;
  debug = 0;
  rewind(datafile->fptr);
  do {
    pos = ftell(datafile->fptr);
    ret = fscanf(datafile->fptr, "%c", &c);
    if(ret != 1) {
      fflush(stdout);
      printerror(debug, "ERROR skipallhashedlines: Cannot read data.");
      return 0;
    }
    if(c == '#') {
      do {
 fscanf(datafile->fptr, "%c", &c);
      }while(c != '\n');
      c = '#';
    }
  }while(c == '#');
  datafile->datastart = (long long)pos;
  fseek(datafile->fptr, datafile->datastart, SEEK_SET);
  return 1;
}
char *str_replace_header_params(datafile_definition data, char *text, verbose_definition verbose)
{
  int i, debug;
  char *newtext, *newtext2, headerparam[1000];
  if(verbose.debug) {
    fflush(stdout);
    printf("Entering str_replace_header_params()\n");
  }
  debug = 0;
  sprintf(headerparam, "%Lf", data.mjd_start);
  newtext = str_replace(text, "%MJD", headerparam, verbose);
  if(newtext == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  sprintf(headerparam, "%ld", data.NrSubints);
  newtext2 = str_replace(newtext, "%NRSUBINTS", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
  sprintf(headerparam, "%ld", data.NrFreqChan);
  newtext2 = str_replace(newtext, "%NRFREQCHAN", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
  double period;
  int ret;
  if(data.isFolded) {
    ret = get_period(data, 0, &period, verbose);
    if(ret == 2) {
      printerror(verbose.debug, "ERROR str_replace_header_params (%s): Cannot obtain period", data.filename);
      return 0;
    }
  }else {
    ret = 1;
    period = 0;
  }
  sprintf(headerparam, "%lf", period);
  newtext2 = str_replace(newtext, "%PERIOD", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
  if(data.freqMode != FREQMODE_UNKNOWN) {
    sprintf(headerparam, "%lf", get_centre_frequency(data, verbose));
  }
  newtext2 = str_replace(newtext, "%FREQ", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
  if((data.freq_ref > -1.1 && data.freq_ref < -0.9) || (data.freq_ref > 0.99e10 && data.freq_ref < 1.01e10)) {
    sprintf(headerparam, "Inf");
  }else if(data.freq_ref >= 0) {
    sprintf(headerparam, "%lf", data.freq_ref);
  }else {
    sprintf(headerparam, "Unknown");
  }
  newtext2 = str_replace(newtext, "%REFFREQ", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
  if(data.freqMode != FREQMODE_UNKNOWN) {
    sprintf(headerparam, "%lf", get_bandwidth(data, verbose));
  }
  newtext2 = str_replace(newtext, "%BW", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
  sprintf(headerparam, "%lf", data.dm);
  newtext2 = str_replace(newtext, "%DM", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
  sprintf(headerparam, "%lf", data.rm);
  newtext2 = str_replace(newtext, "%RM", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
  strncpy(headerparam, data.filename, 999);
  for(i = 0; i < strlen(headerparam); i++) {
    if(headerparam[i] == '.') {
      headerparam[i] = 0;
      break;
    }
  }
  newtext2 = str_replace(newtext, "%FILEBASE", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
  strncpy(headerparam, data.filename, 999);
  newtext2 = str_replace(newtext, "%FILE", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
  strncpy(headerparam, data.psrname, 999);
  newtext2 = str_replace(newtext, "%PSR", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
  if(verbose.debug) {
    fflush(stdout);
    printf("Exiting str_replace_header_params()\n");
  }
  return newtext;
}
void str_list_replace_keys(int nrspaces)
{
  int i;
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%BW             - Bandwidth in MHz\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%DM             - DM\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%FILE           - File name\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%FILEBASE       - File name (everything before first .)\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%FREQ           - Centre frequency in MHz\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%MJD            - MJD\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%NRFREQCHAN     - Nr of frequency channels\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%NRSUBINTS      - Nr of subintegrations\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%PERIOD         - Period in seconds\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%PSR            - Name of the pulsar\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%RM             - RM\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%REFFREQ        - Reference frequency for dedispersion/de-Faraday rotation in MHz\n");
}
void cleanVerboseState(verbose_definition *verbose_state)
{
  verbose_state->verbose = 0;
  verbose_state->debug = 0;
  verbose_state->nocounters = 0;
  verbose_state->indent = 0;
}
void copyVerboseState(verbose_definition verbose_state_src, verbose_definition *verbose_state_dst)
{
  verbose_state_dst->verbose = verbose_state_src.verbose;
  verbose_state_dst->debug = verbose_state_src.debug;
  verbose_state_dst->nocounters = verbose_state_src.nocounters;
  verbose_state_dst->indent = verbose_state_src.indent;
}
int convert_if_uniform_frequency_spacing(datafile_definition *datafile, int nowarnings, verbose_definition verbose)
{
  long subint, freqchan;
  double freq, freq_first, freq_prev, df_expected, actual_centre_freq;
  if(verbose.debug) {
    printf("Check if data has uniformly distributed channels\n");
  }
  if(datafile->freqMode != FREQMODE_FREQTABLE) {
    if(verbose.debug) {
      printf("  Data does not have frequency defined for separate channels/subintegrations, so nothing is done\n");
    }
    return 1;
  }
  freq_first = df_expected = actual_centre_freq = 0.0;
  for(subint = 0; subint < datafile->NrSubints; subint++) {
    freq_prev = 0.0;
    for(freqchan = 0; freqchan < datafile->NrFreqChan; freqchan++) {
      freq = get_weighted_channel_freq(*datafile, subint, freqchan, verbose);
      if(subint == 0)
 actual_centre_freq += freq;
      if(subint == 0 && freqchan == 0) {
 freq_first = freq;
 df_expected = -freq;
      }else if(subint == 0 && freqchan == 1) {
 df_expected += freq;
 freq_prev = freq;
      }else if(freqchan == 0 && subint != 0) {
 if(fabs(freq - freq_first) > 1e-6) {
   if(nowarnings == 0) {
     printwarning(verbose.debug, "WARNING (%s): Frequency channels do not appear to be equally spaced (%lf != %lf for first channels of subint 0 and %ld). Some operations might might not work with this data. Warnings for other channels are suppressed.", datafile->filename, freq, freq_first, subint);
   }
   return 2;
 }
 freq_prev = freq;
      }else {
 if(fabs(freq-freq_prev - df_expected) > 1e-6) {
   if(nowarnings == 0) {
     printwarning(verbose.debug, "WARNING (%s): Frequency channels do not appear to be equally spaced (%lf - %lf != %lf for channel %ld and subint %ld). Some operations might might not work with this data. Warnings for other channels are suppressed.", datafile->filename, freq, freq_prev, df_expected, freqchan, subint);
   }
   return 2;
 }
 freq_prev = freq;
      }
    }
  }
  if(verbose.debug) {
    printf("  The frequency channels appear to be uniformly separated, so convert in data-set described by only a bandwidth and a centre frequency.\n");
  }
  actual_centre_freq /= (double)datafile->NrFreqChan;
  double hdr_freq_cent;
  hdr_freq_cent = get_centre_frequency(*datafile, verbose);
  if(fabs(actual_centre_freq-hdr_freq_cent) > 1e-6) {
    if(hdr_freq_cent > 0 && (hdr_freq_cent < 0.99e10 || hdr_freq_cent > 1.01e10)) {
      double chanbw;
      if(get_channelbandwidth(*datafile, &chanbw, verbose) == 0) {
 printerror(verbose.debug, "ERROR (%s): Cannot obtain channel bandwidth.", datafile->filename);
 return 0;
      }
      if(fabs(hdr_freq_cent-actual_centre_freq-0.5*fabs(chanbw)) < 1e-6) {
 if(verbose.debug) {
   printf("  (%s): Updating centre frequency from %lf to %lf MHz as suggested by the frequency list in the fits table, corresponding to an offset expected from a dropped DC channel.\n", datafile->filename, hdr_freq_cent, actual_centre_freq);
 }
      }else {
 if(nowarnings == 0) {
   printwarning(verbose.debug, "WARNING (%s): Updating centre frequency from %lf to %lf MHz as suggested by the frequency list in the data.", datafile->filename, hdr_freq_cent, actual_centre_freq);
 }
      }
    }else {
      if(verbose.debug)
 printf("  %s: Use subint frequency table frequency as centre frequency = %f MHz.\n", datafile->filename, actual_centre_freq);
    }
    set_centre_frequency(datafile, actual_centre_freq, verbose);
  }
  datafile->freqMode = FREQMODE_UNIFORM;
  free(datafile->freqlabel_list);
  datafile->freqlabel_list = NULL;
  if(datafile->NrFreqChan > 1) {
    df_expected = (freq_prev-freq_first)/(double)(datafile->NrFreqChan-1);
    double chanbw;
    if(get_channelbandwidth(*datafile, &chanbw, verbose) == 0) {
      printerror(verbose.debug, "ERROR (%s): Cannot obtain channel bandwidth.", datafile->filename);
      return 0;
    }
    if(fabs(df_expected-chanbw) > 1e-6) {
      fflush(stdout);
      if(nowarnings == 0) {
 printwarning(verbose.debug, "WARNING (%s): Updating channel bandwidth from %lf to %lf MHz suggested by the frequency list in the subint table.", datafile->filename, chanbw, df_expected);
      }
      if(set_bandwidth(datafile, df_expected*datafile->NrFreqChan, verbose) == 0) {
 printerror(verbose.debug, "ERROR (%s): Bandwidth changing failed.", datafile->filename);
 return 0;
      }
    }
  }
  return 3;
}
int force_uniform_frequency_spacing(datafile_definition *datafile, verbose_definition verbose)
{
  if(verbose.debug) {
    printf("Force channels to be uniformly distributed in frequency\n");
  }
  if(datafile->freqMode != FREQMODE_FREQTABLE) {
    if(verbose.debug) {
      printf("  Data does not have frequency defined for separate channels/subintegrations, so nothing is done\n");
    }
    return 1;
  }
  free(datafile->freqlabel_list);
  datafile->freqlabel_list = NULL;
  datafile->freqMode = FREQMODE_UNIFORM;
  printwarning(verbose.debug, "WARNING (%s): Re-labelling frequency channels to be uniformly distributed and forcing them to be the same for each subintegration. Expect all frequency dependent effects to be wrong.", datafile->filename);
  return 3;
}

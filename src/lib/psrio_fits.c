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
#include "fitsio.h"
#include "psrsalsa.h"
static int psrfits_weightmode = 0;
static int psrfits_use_weighted_freq = 0;
static int psrfits_absweights = 0;
void print_fitsio_version_used(FILE *stream)
{
  float version;
  fprintf(stream, "%f (library) %f (header)", fits_get_version(&version), CFITSIO_VERSION);
}
void psrfits_set_noweights(int val)
{
  psrfits_weightmode = val;
}
void psrfits_set_absweights(int val)
{
  psrfits_absweights = val;
}
void psrfits_set_use_weighted_freq(int val)
{
  psrfits_use_weighted_freq = val;
}
void internalFITSscalePulse(float *pulse, long nrSamples, float *offset, float *scale, float maxvalue)
{
  long i;
  *offset = pulse[0];
  for(i = 0; i < nrSamples; i++) {
    if(pulse[i] < *offset)
      *offset = pulse[i];
  }
  *scale = (pulse[0]-*offset)/maxvalue;
  for(i = 0; i < nrSamples; i++) {
    if((pulse[i]-*offset)/(*scale) > maxvalue)
      *scale = (pulse[i]-*offset)/maxvalue;
  }
  if(*scale == 0.0)
    *scale = 1;
}
int lookupSubintTable(datafile_definition datafile, verbose_definition verbose)
{
  int ncols, status = 0;
  long nrows;
  if(fits_movnam_hdu(datafile.fits_fptr, BINARY_TBL, "SUBINT", 0, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR lookupSubintTable: Cannot move to subint HDU.");
    return 0;
  }
  fits_get_num_rows(datafile.fits_fptr, &nrows, &status);
  fits_get_num_cols(datafile.fits_fptr, &ncols, &status);
  if(ncols < 19) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR lookupSubintTable: Table has not the correct size, was writePSRFITSHeader not called?");
    return 0;
  }
  return 1;
}
int writeFITSpulse(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  int status = 0;
  long i, n;
  int *data, dummy_int;
  float scale, offset, weight;
  if(datafile.gentype == GENTYPE_RECEIVERMODEL || datafile.gentype == GENTYPE_RECEIVERMODEL2) {
    if(fits_movnam_hdu(datafile.fits_fptr, BINARY_TBL, "FEEDPAR", 0, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeFITSpulse: Cannot move to FEEDPAR HDU.");
      return 0;
    }
    for(n = binnr; n < binnr+nrSamples; n++) {
      if(datafile.gentype == GENTYPE_RECEIVERMODEL) {
 i = 1+freq*datafile.NrPols+polarization;
      }else if(datafile.gentype == GENTYPE_RECEIVERMODEL2) {
 i = 1+freq*(datafile.NrPols-2)+polarization;
      }
      if(n == 0) {
 if(datafile.gentype == GENTYPE_RECEIVERMODEL || polarization < datafile.NrPols-2) {
   if(fits_write_col(datafile.fits_fptr, TFLOAT, 3, 1+pulsenr, i, 1, &pulse[n-binnr], &status) != 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR writeFITSpulse: Error writing data (value pulse=%ld freq=%d pol=%d).", pulsenr, freq, polarization);
     fits_report_error(stderr, status);
     return 0;
   }
 }
 if(datafile.gentype == GENTYPE_RECEIVERMODEL2 && polarization == datafile.NrPols-2 && n == 0) {
   if(fits_write_col(datafile.fits_fptr, TFLOAT, 5, 1+pulsenr, 1+freq, 1, &pulse[n-binnr], &status) != 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR writeFITSpulse: Error writing data (value pulse=%ld freq=%d pol=%d).", pulsenr, freq, polarization);
     fits_report_error(stderr, status);
     return 0;
   }
 }
 if(datafile.gentype == GENTYPE_RECEIVERMODEL2 && polarization == datafile.NrPols-1 && n == 0) {
   dummy_int = pulse[n-binnr];
   if(fits_write_col(datafile.fits_fptr, TINT, 6, 1+pulsenr, 1+freq, 1, &dummy_int, &status) != 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR writeFITSpulse: Error writing data (value pulse=%ld freq=%d pol=%d).", pulsenr, freq, polarization);
     fits_report_error(stderr, status);
     return 0;
   }
 }
      }else if(n == 1) {
 if(datafile.gentype == GENTYPE_RECEIVERMODEL || polarization < datafile.NrPols-2) {
   if(fits_write_col(datafile.fits_fptr, TFLOAT, 4, 1+pulsenr, i, 1, &pulse[n-binnr], &status) != 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR writeFITSpulse: Error writing data (error pulse=%ld freq=%d pol=%d).", pulsenr, freq, polarization);
     fits_report_error(stderr, status);
     return 0;
   }
 }
      }else {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writeFITSpulse: Dont expect receiver model to contain more than two bins.");
 return 0;
      }
    }
  }else {
    data = (int *)malloc(datafile.NrBins*sizeof(int));
    if(data == NULL)
      return 0;
    if(lookupSubintTable(datafile, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeFITSpulse: Cannot move to subint table.");
      return 0;
    }
    double period;
    int ret;
    if(datafile.isFolded) {
      ret = get_period(datafile, 0, &period, verbose);
      if(ret == 2) {
 printerror(verbose.debug, "ERROR writeFITSpulse (%s): Cannot obtain period", datafile.filename);
 return 0;
      }
    }else {
      ret = 1;
      period = -1;
    }
    if(datafile.gentype == GENTYPE_SEARCHMODE || datafile.isFolded != 1 || ret == 1 || period < 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeFITSpulse: This function is not working properly for search mode data. Write out whole subints instead of per channel data.");
      return 0;
    }
    if(binnr != 0 || nrSamples != datafile.NrBins) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeFITSpulse: Can not write out less than a full period at a time in PSRFITS output (trying to write %ld samples out of %ld bins starting at position %ld).", nrSamples, datafile.NrBins, binnr);
      return 0;
    }
    internalFITSscalePulse(pulse, nrSamples, &offset, &scale, 32767);
    for(i = 0; i < nrSamples; i++) {
      data[i] = (pulse[i]-offset)/scale;
    }
    if(datafile.gentype != GENTYPE_SEARCHMODE && datafile.isFolded == 1 && ret == 0 && period >= 0) {
      i = 1+polarization*datafile.NrBins*datafile.NrFreqChan+freq*datafile.NrBins+binnr;
      if(fits_write_col(datafile.fits_fptr, TINT, 18, 1+pulsenr, i, nrSamples, data, &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writeFITSpulse: Error writing data: subint=%ld sample=%ld (pol=%ld freq=%ld bin=%d) nrSamples=%ld", pulsenr, i, polarization, freq, binnr, nrSamples);
 if(verbose.debug) {
   long sampnr;
   printerror(verbose.debug, "ERROR writeFITSpulse: Samples are: ");
   for(sampnr = 0; sampnr < nrSamples; sampnr++) {
     fprintf(stderr, "%d ", data[sampnr]);
   }
 }
 fits_report_error(stderr, status);
 for(i = 0; i < nrSamples; i++) {
   if(isnan(pulse[i]) || isinf(pulse[i])) {
     printerror(verbose.debug, "ERROR writeFITSpulse: Sample %ld: %f", i+1, pulse[i]);
     break;
   }
 }
 return 0;
      }
    }else {
    }
    free(data);
    if(fits_write_col(datafile.fits_fptr, TFLOAT, 17, 1+pulsenr, 1+polarization*datafile.NrFreqChan+freq, 1, &scale, &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeFITSpulse: Error writing scales.");
      fits_report_error(stderr, status);
      return 0;
    }
    if(fits_write_col(datafile.fits_fptr, TFLOAT, 16, 1+pulsenr, 1+polarization*datafile.NrFreqChan+freq, 1, &offset, &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeFITSpulse: Error writing offsets.");
      fits_report_error(stderr, status);
      return 0;
    }
    weight = 1;
    if(fits_write_col(datafile.fits_fptr, TFLOAT, 15, 1+pulsenr, 1+freq, 1, &weight, &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeFITSpulse: Error writing weights (subint=%d, freq=%d).", 1+pulsenr, 1+freq);
      fits_report_error(stderr, status);
      return 0;
    }
  }
  return 1;
}
int readFITSpulse_receivermodel(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  int status = 0;
  int anynul, colnum, datapoint_int, ncpar;
  float datapoint, datapoint_err;
  if(datafile->gentype == GENTYPE_RECEIVERMODEL2 && polarization == datafile->NrPols-2) {
    if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "CHISQ", &colnum, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: No chi^2 values defined in fits file?");
      return 0;
    }
    if(binnr == 0) {
      if(fits_read_col(datafile->fits_fptr, TFLOAT, colnum, 1+pulsenr, 1+freq, 1, NULL, &datapoint, &anynul, &status)) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: Cannot read CHISQ data");
 fits_report_error(stderr, status);
 return 0;
      }
      if(nrSamples == 1) {
 pulse[0] = datapoint;
      }else if(nrSamples == 2) {
 pulse[0] = datapoint;
 pulse[1] = 0;
      }else {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: Nr of bins is expected to be 1 or 2");
 return 0;
      }
    }else if(binnr == 1) {
      pulse[0] = 0;
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: Nr of bins is expected to be 1 or 2");
      return 0;
    }
    return 1;
  }
  if(datafile->gentype == GENTYPE_RECEIVERMODEL2 && polarization == datafile->NrPols-1) {
    if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "NFREE", &colnum, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: No NFREE values defined in fits file?");
      return 0;
    }
    if(binnr == 0) {
      if(fits_read_col(datafile->fits_fptr, TINT, colnum, 1+pulsenr, 1+freq, 1, NULL, &datapoint_int, &anynul, &status)) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: Cannot read NFREE data");
 fits_report_error(stderr, status);
 return 0;
      }
      if(nrSamples == 1) {
 pulse[0] = datapoint_int;
      }else if(nrSamples == 2) {
 pulse[0] = datapoint_int;
 pulse[1] = 0;
      }else {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: Nr of bins is expected to be 1 or 2");
 return 0;
      }
    }else if(binnr == 1) {
      pulse[0] = 0;
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: Nr of bins is expected to be 1 or 2");
      return 0;
    }
    return 1;
  }
  ncpar = datafile->NrPols;
  if(datafile->gentype == GENTYPE_RECEIVERMODEL2)
    ncpar -= 2;
  if(nrSamples == 2 || binnr == 1) {
    if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "DATAERR", &colnum, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: No receiver parameter error values defined in fits file?");
      return 0;
    }
    if(fits_read_col(datafile->fits_fptr, TFLOAT, colnum, 1+pulsenr, 1+ncpar*freq+polarization, 1, NULL, &datapoint_err, &anynul, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: Cannot read DATAERR data");
      fits_report_error(stderr, status);
      return 0;
    }
  }
  if(binnr == 0) {
    if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "DATA", &colnum, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: No receiver parameter values defined in fits file?");
      return 0;
    }
    if(fits_read_col(datafile->fits_fptr, TFLOAT, colnum, 1+pulsenr, 1+ncpar*freq+polarization, 1, NULL, &datapoint, &anynul, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: Cannot read DATA data");
      fits_report_error(stderr, status);
      return 0;
    }
    if(nrSamples == 1) {
      pulse[0] = datapoint;
    }else if(nrSamples == 2) {
      pulse[0] = datapoint;
      pulse[1] = datapoint_err;
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: Nr of bins is expected to be 1 or 2");
      return 0;
    }
  }else if(binnr == 1) {
    pulse[0] = datapoint_err;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readFITSpulse_receivermodel: Nr of bins is expected to be 1 or 2");
    return 0;
  }
  return 1;
}
int readFITSpulse(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  int status = 0;
  int anynul, ret, colnum;
  long i, istart, bstart, bsample;
  int *data;
  unsigned char singlebyte;
  if(datafile->gentype == GENTYPE_RECEIVERMODEL || datafile->gentype == GENTYPE_RECEIVERMODEL2) {
    return readFITSpulse_receivermodel(datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  }
  int weightmode;
  weightmode = 0;
  if(pulsenr == 0 && polarization == 0 && freq == 0 && binnr == 0) {
    if(psrfits_weightmode == 1) {
      printwarning(verbose.debug, "WARNING: The data is NOT multiplied with the weight, even when set to zero. This might undo any zapping done and avoids introducing artificial intensity fluctuations IF the data is written as the weighted average rather than the sum. This will not be ideal when summing data at a later stage.");
    }else if(psrfits_weightmode == 2) {
      printwarning(verbose.debug, "WARNING: The data is NOT multiplied with the weight, unless when zero. This might avoid introducing artificial intensity fluctuations IF the data is written as the weighted average rather than the sum. This will not be ideal when summing data at a later stage.");
    }else if(psrfits_weightmode == 3) {
      printwarning(verbose.debug, "WARNING: Data will be multiplied with the weight. This might be benificial when summing data at a later stage IF the data is written as the weighted average rather than the sum. However, it might introduce artificial intensity fluctuations as well.");
    }
  }
  if(datafile->NrFreqChan == 1) {
    if(psrfits_weightmode == 0) {
      if(datafile->weight_stats_differentweights && pulsenr == 0 && polarization == 0 && freq == 0 && binnr == 0) {
 printwarning(verbose.debug, "WARNING: FITS file contains data with different weights. For data with only one frequency channel the data is NOT multiplied with the weight if it is nonzero. This avoids introducing artificial intensity fluctuations IF the data is written as the weighted average rather than the sum. This will not be ideal when summing data at a later stage. Use the -useweights option to multiply the data with the weights.");
      }
      weightmode = 2;
    }else {
      weightmode = psrfits_weightmode;
    }
  }else {
    if(psrfits_weightmode == 0) {
      if(datafile->weight_stats_differentweights && pulsenr == 0 && polarization == 0 && freq == 0 && binnr == 0) {
 printwarning(verbose.debug, "WARNING: FITS file contains data with different weights. For data with multiple frequency channels the data is multiplied with the weight. This might be benificial when summing data at a later stage IF the data is written as the weighted average rather than the sum. However, it might introduce artificial intensity fluctuations as well. Use the -uniformweights to take the weights equal.");
      }
      weightmode = 3;
    }else {
      weightmode = psrfits_weightmode;
    }
  }
  data = (int *)malloc(datafile->NrBins*sizeof(int));
  if(data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readFITSpulse: Cannot allocate memory");
    return 0;
  }
  if(psrfits_use_weighted_freq) {
    if(fits_get_colnum(datafile->fits_fptr, CASEINSEN, "DAT_FREQ", &colnum, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readFITSpulse: No frequency column in fits file?");
      return 0;
    }
    datafile->freqMode = FREQMODE_UNIFORM;
    if(datafile->freqlabel_list != NULL) {
      free(datafile->freqlabel_list);
      datafile->freqlabel_list = NULL;
    }
    double freq;
    if(!fits_read_col(datafile->fits_fptr, TFLOAT, colnum, 1+pulsenr, 1, 1, NULL, &freq, &anynul, &status)) {
    }
    set_centre_frequency(datafile, freq, verbose);
    if(set_bandwidth(datafile, 0.0, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readFITSpulse: Bandwidth changing failed");
      return 0;
    }
  }
  if(fits_get_colnum(datafile->fits_fptr, CASEINSEN, "DATA", &colnum, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readFITSpulse: No data in fits file?");
    return 0;
  }
  ret = 0;
  if(datafile->NrBits == 16) {
    istart = polarization*datafile->NrBins*datafile->NrFreqChan+freq*datafile->NrBins+binnr;
    if(!fits_read_col(datafile->fits_fptr, TINT, colnum, 1+pulsenr, 1+istart, nrSamples, NULL, data, &anynul, &status)) {
      ret = 1;
    }
  }else if(datafile->NrBits == 2 || datafile->NrBits == 4 || datafile->NrBits == 8) {
    ret = 1;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readFITSpulse: Unsupported number of bits.");
  }
  if(ret == 1) {
    float scale, offset;
    long index = pulsenr*datafile->NrPols*datafile->NrFreqChan+polarization*datafile->NrFreqChan+freq;
    scale = datafile->scales[index];
    offset = datafile->offsets[index];
    for(i = 0; i < nrSamples; i++) {
      if(datafile->NrBits == 16) {
 pulse[i] = scale*(data[i]) + offset;
      }else if(datafile->NrBits == 2 || datafile->NrBits == 4 || datafile->NrBits == 8) {
 istart = (binnr+i)*datafile->NrPols*datafile->NrFreqChan+polarization*datafile->NrFreqChan+freq;
 bstart = istart*datafile->NrBits;
 istart = bstart/8;
 if(!fits_read_col(datafile->fits_fptr, TBYTE, colnum, 1+pulsenr, 1+istart, 1, NULL, &singlebyte, &anynul, &status)) {
   ret = 1;
 }
 bsample = bstart % 8;
 if(datafile->NrBits == 4) {
   if(bsample == 0) {
     singlebyte = singlebyte/16;
   }else if(bsample == 4) {
     singlebyte = singlebyte & 15;
   }else {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readFITSpulse: bin shift error.");
     return 0;
   }
 }else if(datafile->NrBits == 8) {
 }else if(datafile->NrBits == 2) {
   if(bsample == 0) {
     singlebyte = singlebyte/64;
   }else if(bsample == 2) {
     singlebyte = singlebyte & 48;
     singlebyte = singlebyte/16;
   }else if(bsample == 4) {
     singlebyte = singlebyte & 12;
     singlebyte = singlebyte/4;
   }else if(bsample == 6) {
     singlebyte = singlebyte & 3;
   }else {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readFITSpulse: bin shift error.");
     return 0;
   }
 }
 pulse[i] = scale*(singlebyte) + offset;
      }
      if(weightmode == 3 && psrfits_absweights == 0) {
 pulse[i] *= datafile->weights[pulsenr*datafile->NrFreqChan+freq];
      }else if(weightmode == 3 && psrfits_absweights != 0) {
 pulse[i] *= fabs(datafile->weights[pulsenr*datafile->NrFreqChan+freq]);
      }else if(weightmode == 2) {
 if(datafile->weights[pulsenr*datafile->NrFreqChan+freq] == 0.0)
 pulse[i] = 0;
      }
    }
  }
  free(data);
  if (status) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readFITSpulse: ");
    fits_report_error(stderr, status);
  }
  return ret;
}
int readPSRCHIVE_ASCIIHeader(datafile_definition *datafile, verbose_definition verbose)
{
  char *filename, *tmp, *psrname;
  int j;
  filename = malloc(1000);
  tmp = malloc(1000);
  psrname = malloc(1000);
  if(filename == NULL || tmp == NULL || psrname == NULL) {
    fflush(stdout);
    printerror(verbose.debug,"ERROR readPSRCHIVE_ASCIIHeader: Memory allocation error.");
    exit(-1);
  }
  j = fscanf(datafile->fptr_hdr, "File: %s Src: %s Nsub: %ld Nch: %ld Npol: %ld Nbin: %ld RMS: %s", filename, psrname, &(datafile->NrSubints), &(datafile->NrFreqChan), &(datafile->NrPols), &(datafile->NrBins), tmp);
  if(j != 7) {
    fflush(stdout);
    printerror(verbose.debug,"ERROR readPSRCHIVE_ASCIIHeader: Error reading first line.");
    exit(-1);
  }
  if(datafile->NrPols == 7 && strcmp(filename, "penergy") == 0) {
    datafile->gentype = GENTYPE_PENERGY;
  }
  if(set_psrname_PSRData(datafile, psrname, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIHeader: Setting pulsar name failed.");
    return 0;
  }
  datafile->isFolded = 1;
  datafile->foldMode = FOLDMODE_FIXEDPERIOD;
  datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
  datafile->tsubMode = TSUBMODE_FIXEDTSUB;
  if(datafile->tsub_list != NULL)
    free(datafile->tsub_list);
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIHeader: Memory allocation error");
    return 0;
  }
  datafile->tsub_list[0] = 0;
  datafile->freqMode = FREQMODE_UNIFORM;
  if(datafile->freqlabel_list != NULL) {
    free(datafile->freqlabel_list);
    datafile->freqlabel_list = NULL;
  }
  set_centre_frequency(datafile, 0.0, verbose);
  if(set_bandwidth(datafile, 0.0, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIHeader: Setting bandwidth failed");
    return 0;
  }
  free(filename);
  free(tmp);
  free(psrname);
  return 1;
}
int writePSRCHIVE_ASCIIHeader(datafile_definition datafile, verbose_definition verbose)
{
  int i;
  char *psrname;
  psrname = malloc(1000);
  if(psrname == NULL) {
    fflush(stdout);
    printerror(verbose.debug,"ERROR writePSRCHIVE_ASCIIHeader: Memory allocation error.");
  }
  strncpy(psrname, datafile.psrname, 999);
  if(strlen(psrname) > 0) {
    for(i = 0; i < strlen(psrname); i++)
      if(psrname[i] == ' ')
 psrname[i] = '_';
  }else {
    strcpy(psrname, "?");
  }
  fprintf(datafile.fptr_hdr, "File: %s Src: %s Nsub: %ld Nch: %ld Npol: %ld Nbin: %ld RMS: 0.0\n", datafile.filename, psrname, (datafile.NrSubints), (datafile.NrFreqChan), (datafile.NrPols), (datafile.NrBins));
  free(psrname);
  return 1;
}
int readPSRFITSscales(datafile_definition *datafile, verbose_definition verbose)
{
  char card[FLEN_CARD], comment[FLEN_COMMENT], value[FLEN_VALUE];
  int status = 0;
  int i, ncols, anynul, ret, fullscales, colnum_s, colnum_o, colnum_w;
  long nrows, n, p, datalength;
  char txt[100];
  ret = 1;
  datafile->scales = (float *)malloc(datafile->NrSubints*datafile->NrPols*datafile->NrFreqChan*sizeof(float));
  datafile->offsets = (float *)malloc(datafile->NrSubints*datafile->NrPols*datafile->NrFreqChan*sizeof(float));
  datafile->weights = (float *)malloc(datafile->NrSubints*datafile->NrFreqChan*sizeof(float));
  if(datafile->scales == NULL || datafile->offsets == NULL || datafile->weights == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSscales: Memory allocation failed");
    return 0;
  }
  fits_get_num_rows(datafile->fits_fptr, &nrows, &status);
  fits_get_num_cols(datafile->fits_fptr, &ncols, &status);
  if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "DAT_SCL", &colnum_s, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSscales: No scales in fits file?");
    return 0;
  }
  if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "DAT_OFFS", &colnum_o, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSscales: No offsets in fits file?");
    return 0;
  }
  if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "DAT_WTS", &colnum_w, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSscales: No weights in fits file?");
    return 0;
  }
  sprintf(txt, "TFORM%d", colnum_s);
  if (fits_read_card(datafile->fits_fptr,txt, card, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSscales: keyword does not exist");
    return 0;
  }
  fits_parse_value(card, value, comment, &status);
  if(value[1] == 'E') {
    datalength = 1;
  }else {
    sscanf(value, "'%ld'", &datalength);
  }
  if(verbose.verbose && verbose.debug) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("readPSRFITSscales: Datalength = %ld\n", datalength);
  }
  if(datalength == datafile->NrPols*datafile->NrFreqChan) {
    fullscales = 1;
  }else if(datalength == datafile->NrFreqChan) {
    fullscales = 0;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSscales: %ld is a weird number for the number of scales.", datalength);
    return 0;
  }
  if(verbose.verbose && verbose.debug) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("readPSRFITSscales: fullscales = %d\n", fullscales);
  }
  for(n = 0; n < datafile->NrSubints; n++) {
    for(p = 0; p < datafile->NrPols; p++) {
      if(!fits_read_col(datafile->fits_fptr, TFLOAT, colnum_s, 1+n, 1+fullscales*p*datafile->NrFreqChan, datafile->NrFreqChan, NULL, &datafile->scales[n*datafile->NrPols*datafile->NrFreqChan+p*datafile->NrFreqChan], &anynul, &status)) {
      }
      if(!fits_read_col(datafile->fits_fptr, TFLOAT, colnum_o, 1+n, 1+fullscales*p*datafile->NrFreqChan, datafile->NrFreqChan, NULL, &datafile->offsets[n*datafile->NrPols*datafile->NrFreqChan+p*datafile->NrFreqChan], &anynul, &status)) {
      }
      if(p == 0) {
 if(!fits_read_col(datafile->fits_fptr, TFLOAT, colnum_w, 1+n, 1+fullscales*p*datafile->NrFreqChan, datafile->NrFreqChan, NULL, &datafile->weights[n*datafile->NrFreqChan], &anynul, &status)) {
 }
      }
      if (status) {
 fflush(stdout);
 fits_report_error(stderr, status);
 ret = 0;
 break;
      }
    }
  }
  return ret;
}
int writePSRFITSHeader(datafile_definition *datafile, verbose_definition verbose)
{
  int status = 0;
  char dummy_txt[2000], dummy_txt2[2000], dummy_txt3[2000];
  int dummy_int, dummy_int2;
  float dummy_float;
  double dummy_double;
  long i, j, k, nrcolumns;
  if(fits_create_img(datafile->fits_fptr, 8, 0, NULL, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot create primary table.");
    return 0;
  }
  sprintf(dummy_txt, "4.1");
  if(fits_write_key(datafile->fits_fptr, TSTRING, "HDRVER", dummy_txt, "Header version", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  sprintf(dummy_txt, "PSRFITS");
  if(fits_write_key(datafile->fits_fptr, TSTRING, "FITSTYPE", dummy_txt, "FITS definition for pulsar data files", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(fits_write_key(datafile->fits_fptr, TSTRING, "TELESCOP", datafile->observatory, "Telescope name", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(fits_write_key(datafile->fits_fptr, TDOUBLE, "ANT_X", &datafile->telescope_X, "[m] Antenna ITRF X-coordinate", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(fits_write_key(datafile->fits_fptr, TDOUBLE, "ANT_Y", &datafile->telescope_Y, "[m] Antenna ITRF Y-coordinate", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(fits_write_key(datafile->fits_fptr, TDOUBLE, "ANT_Z", &datafile->telescope_Z, "[m] Antenna ITRF Z-coordinate", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  double period;
  int ret;
  if(datafile->isFolded) {
    ret = get_period(*datafile, 0, &period, verbose);
    if(ret == 2) {
      printerror(verbose.debug, "ERROR writePSRFITSHeader (%s): Cannot obtain period", datafile->filename);
      return 0;
    }
  }else {
    ret = 1;
    period = -1;
  }
  if(datafile->gentype == GENTYPE_RECEIVERMODEL || datafile->gentype == GENTYPE_RECEIVERMODEL2) {
    sprintf(dummy_txt, "PCM");
  }else if(datafile->gentype != GENTYPE_SEARCHMODE && datafile->isFolded == 1 && ret == 0 && period >= 0) {
    if(datafile->gentype == GENTYPE_POLNCAL) {
      sprintf(dummy_txt, "CAL");
    }else {
      sprintf(dummy_txt, "PSR");
    }
  }else {
    sprintf(dummy_txt, "SEARCH");
  }
  if(fits_write_key(datafile->fits_fptr, TSTRING, "OBS_MODE", dummy_txt, "(PSR, CAL, SEARCH)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(fits_write_key(datafile->fits_fptr, TSTRING, "SRC_NAME", datafile->psrname, "Source or scan ID", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword SRC_NAME (='%s').", datafile->psrname);
    return 0;
  }
  if(fits_write_key(datafile->fits_fptr, TSTRING, "BACKEND", datafile->instrument, "Backend ID", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(fits_write_key(datafile->fits_fptr, TSTRING, "OBSERVER", datafile->observer, "Observer name(s)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(fits_write_key(datafile->fits_fptr, TSTRING, "PROJID", datafile->projectID, "Project name", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  mjd2dateString(datafile->mjd_start, dummy_txt, 0, 1, "T");
  if(fits_write_key(datafile->fits_fptr, TSTRING, "DATE-OBS", dummy_txt, "Date of observation (YYYY-MM-DDThh:mm:ss UTC)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  converthms_string(dummy_txt, (12.0/M_PI)*datafile->ra, 4, 1);
  if(fits_write_key(datafile->fits_fptr, TSTRING, "RA", dummy_txt, "Right ascension (hh:mm:ss.ssss)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  converthms_string(dummy_txt, (180.0/M_PI)*datafile->dec, 3, 1);
  if(fits_write_key(datafile->fits_fptr, TSTRING, "DEC", dummy_txt, "Declination (-dd:mm:ss.sss)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_float = 2000;
  if(fits_write_key(datafile->fits_fptr, TFLOAT, "EQUINOX", &dummy_float, "Equinox of coords (e.g. 2000.0)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  sprintf(dummy_txt, "J2000");
  if(fits_write_key(datafile->fits_fptr, TSTRING, "COORD_MD", dummy_txt, "Coordinate mode (J2000, GALACTIC, ECLIPTIC)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_txt[0] = 0;
  if(fits_write_key(datafile->fits_fptr, TSTRING, "OBSERVER", dummy_txt, "Observer name(s)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(fits_write_key(datafile->fits_fptr, TSTRING, "FRONTEND", dummy_txt, "Receiver ID", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(fits_write_key(datafile->fits_fptr, TSTRING, "PROJID", dummy_txt, "Project name", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_int = datafile->NrFreqChan;
  if(fits_write_key(datafile->fits_fptr, TINT, "OBSNCHAN", &dummy_int, "Number of frequency channels (original)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_float = get_bandwidth(*datafile, verbose);
  if(fits_write_key(datafile->fits_fptr, TFLOAT, "OBSBW", &dummy_float, "[MHz] Bandwidth for observation", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_double = 0;
  if(fits_write_key(datafile->fits_fptr, TDOUBLE, "CHAN_DM", &dummy_double, "[cm-3 pc] DM used for on-line dedispersion", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  sprintf(dummy_txt, "TRACK");
  if(fits_write_key(datafile->fits_fptr, TSTRING, "TRK_MODE", dummy_txt, "Track mode (TRACK, SCANGC, SCANLAT)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_float = 0;
  if(fits_write_key(datafile->fits_fptr, TFLOAT, "BMIN", &dummy_float, "[deg] Beam minor axis length (needs to be set to something for PRESTO)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(datafile->feedtype == FEEDTYPE_LINEAR || datafile->feedtype == FEEDTYPE_INV_LINEAR) {
    sprintf(dummy_txt, "LIN");
  }else if(datafile->feedtype == FEEDTYPE_CIRCULAR || datafile->feedtype == FEEDTYPE_INV_CIRCULAR) {
    sprintf(dummy_txt, "CIRC");
  }
  if(datafile->feedtype != FEEDTYPE_UNKNOWN) {
    if(fits_write_key(datafile->fits_fptr, TSTRING, "FD_POLN", dummy_txt, "LIN or CIRC", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    if(datafile->feedtype == FEEDTYPE_LINEAR || datafile->feedtype == FEEDTYPE_CIRCULAR) {
      dummy_int = 1;
    }else {
      dummy_int = -1;
    }
    if(fits_write_key(datafile->fits_fptr, TINT, "FD_HAND", &dummy_int, "+/- 1. +1 is LIN:A=X,B=Y, CIRC:A=L,B=R (I)", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
  }else {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING writePSRFITSHeader: feedtype is not set.");
  }
  if(fits_write_key(datafile->fits_fptr, TDOUBLE, "OBSFREQ", &datafile->freq_ref, "[MHz] Centre frequency for observation", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_int = datafile->mjd_start;
  if(fits_write_key(datafile->fits_fptr, TINT, "STT_IMJD", &dummy_int, "Start MJD (UTC days) (J - long integer)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_int2 = (int)((datafile->mjd_start - (double)dummy_int)*86400.0 + 0.5);
  if(fits_write_key(datafile->fits_fptr, TINT, "STT_SMJD", &dummy_int2, "[s] Start time (sec past UTC 00h) (J)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_double = (datafile->mjd_start - (double)dummy_int)*(24.0*3600.0);
  dummy_double -= dummy_int2;
  if(fits_write_key(datafile->fits_fptr, TDOUBLE, "STT_OFFS", &dummy_double, "[s] Start time offset (D)", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(fits_write_key(datafile->fits_fptr, TINT, "PS_GTYPE", &datafile->gentype, "Describes the type of data in file", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(datafile->xrangeset) {
    if(fits_write_key(datafile->fits_fptr, TFLOAT, "PS_XMIN", &datafile->xrange[0], "Centre of leftmost bin", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    if(fits_write_key(datafile->fits_fptr, TFLOAT, "PS_XMAX", &datafile->xrange[1], "Centre of rightmost bin", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
  }
  if(datafile->yrangeset) {
    if(fits_write_key(datafile->fits_fptr, TFLOAT, "PS_YMIN", &datafile->yrange[0], "Centre of bottom bin", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    if(fits_write_key(datafile->fits_fptr, TFLOAT, "PS_YMAX", &datafile->yrange[1], "Centre of top bin", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
  }
  if(datafile->cableSwap != -1) {
    dummy_int = datafile->cableSwap;
    if(fits_write_key(datafile->fits_fptr, TINT, "PS_CSW", &dummy_int, "Flag is set if cable swap during observation", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
  }
  if(datafile->cableSwapcor != -1) {
    dummy_int = datafile->cableSwapcor;
    if(fits_write_key(datafile->fits_fptr, TINT, "PS_CSWC", &dummy_int, "Flag is set if cable swap is corrected", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
  }
  dummy_int = datafile->isDebase;
  if(fits_write_key(datafile->fits_fptr, TINT, "PS_DBASE", &dummy_int, "Describes if baseline is subtracted", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  char *ttypes_history[] = {"POL_TYPE", "NPOL", "NBIN", "TBIN", "CTR_FREQ", "NCHAN", "CHAN_BW", "RM_CORR", "DEDISP", "PR_CORR"};
  char *tform_history[] = {"8A", "1I", "1J", "1D", "1D", "1I", "1D", "1I", "1I", "1I"};
  char *tunit_history[] = {"", "", "", "s", "MHz", "", "MHz", "", "", ""};
  if(fits_create_tbl(datafile->fits_fptr, BINARY_TBL, 0, 10, ttypes_history, tform_history, tunit_history, "HISTORY", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot create table.");
    return 0;
  }
  if(datafile->poltype == POLTYPE_STOKES) {
    sprintf(dummy_txt, "STOKES");
  }else if(datafile->poltype == POLTYPE_COHERENCY) {
    sprintf(dummy_txt, "AABBCRCI");
  }else if(datafile->poltype == POLTYPE_ILVPAdPA) {
    sprintf(dummy_txt, "ILVPAdPA");
  }else if(datafile->poltype == POLTYPE_ILVPAdPATEldEl) {
    sprintf(dummy_txt, "ILVPATEL");
  }else if(datafile->poltype == POLTYPE_PAdPA) {
    sprintf(dummy_txt, "PAdPA");
  }else {
    sprintf(dummy_txt, "UNKNOWN");
  }
  char *tst = dummy_txt;
  if(fits_write_col(datafile->fits_fptr, TSTRING, 1, 1, 1, 1, &tst, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_int = datafile->NrPols;
  if(fits_write_col(datafile->fits_fptr, TINT, 2, 1, 1, 1, &dummy_int, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_int = datafile->NrBins;
  if(fits_write_col(datafile->fits_fptr, TINT, 3, 1, 1, 1, &dummy_int, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword: nrbins=%ld", dummy_int);
    return 0;
  }
  if(fits_write_col(datafile->fits_fptr, TDOUBLE, 4, 1, 1, 1, &datafile->fixedtsamp, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  double freq_ref_use;
  freq_ref_use = datafile->freq_ref;
  if((datafile->freq_ref < -0.9 && datafile->freq_ref > -1.1) || (datafile->freq_ref > 0.99e10 && datafile->freq_ref < 1.01e10)) {
    freq_ref_use = 1e10;
    fflush(stdout);
    printwarning(verbose.debug, "WARNING: The reference frequency of dedispersion is infinity. This is encoded in the header as a value of %.1lf MHz.", freq_ref_use);
 }else if(datafile->freq_ref < 0) {
    fflush(stdout);
    if(datafile->freq_ref < -1.9 && datafile->freq_ref > -2.1) {
      printwarning(verbose.debug, "WARNING: The reference frequency of dedispersion is unknown, and this is encoded in the header as a value of %.1lf MHz.", datafile->freq_ref);
    }else {
      printwarning(verbose.debug, "WARNING: The reference frequency of dedispersion appears to negative (%lf MHz), so it is unclear what this means.", freq_ref_use);
    }
  }
  if(fits_write_col(datafile->fits_fptr, TDOUBLE, 5, 1, 1, 1, &freq_ref_use, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_int = datafile->NrFreqChan;
  if(fits_write_col(datafile->fits_fptr, TINT, 6, 1, 1, 1, &dummy_int, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  if(get_channelbandwidth(*datafile, &dummy_double, verbose) == 0) {
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot obtain channel bandwidth.");
    return 0;
  }
  if(fits_write_col(datafile->fits_fptr, TDOUBLE, 7, 1, 1, 1, &dummy_double, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_int = datafile->isDeFarad;
  if(fits_write_col(datafile->fits_fptr, TINT, 8, 1, 1, 1, &dummy_int, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_int = datafile->isDeDisp;
  if(fits_write_col(datafile->fits_fptr, TINT, 9, 1, 1, 1, &dummy_int, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  dummy_int = datafile->isDePar;
  if(fits_write_col(datafile->fits_fptr, TINT, 10, 1, 1, 1, &dummy_int, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
    return 0;
  }
  char *ttypes_history2[] = {"DATE_PRO", "USER", "HOSTNAME", "PROC_CMD"};
  char *tform_history2[] = {"24A", "24A", "32A", "1024A"};
  char *tunit_history2[] = {"", "", "", ""};
  if(fits_create_tbl(datafile->fits_fptr, BINARY_TBL, 0, 4, ttypes_history2, tform_history2, tunit_history2, "HISTORY_NOT_PSRFITS", &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot create table.");
    return 0;
  }
  if(datafile->gentype != GENTYPE_RECEIVERMODEL && datafile->gentype != GENTYPE_RECEIVERMODEL2) {
    sprintf(dummy_txt, "%ldD", datafile->NrFreqChan);
    sprintf(dummy_txt2, "%ldE", datafile->NrFreqChan*datafile->NrPols);
    double period;
    int ret;
    if(datafile->isFolded) {
      ret = get_period(*datafile, 0, &period, verbose);
      if(ret == 2) {
 printerror(verbose.debug, "ERROR writePSRFITSHeader (%s): Cannot obtain period", datafile->filename);
 return 0;
      }
    }else {
      ret = 1;
      period = -1;
    }
    if(datafile->gentype != GENTYPE_SEARCHMODE && datafile->isFolded == 1 && ret == 0 && period >= 0) {
      sprintf(dummy_txt3, "%ldI", datafile->NrFreqChan*datafile->NrPols*datafile->NrBins);
    }else {
      long nrbytes = ceil((datafile->NrFreqChan*datafile->NrPols*datafile->NrBins*datafile->NrBits)/8.0);
      sprintf(dummy_txt3, "%ldb", nrbytes);
    }
    char *ttypes_subint[] = {"INDEXVAL", "TSUBINT", "OFFS_SUB", "LST_SUB", "RA_SUB", "DEC_SUB", "GLON_SUB", "GLAT_SUB", "FD_ANG", "POS_ANG", "PAR_ANG", "TEL_AZ", "TEL_ZEN", "DAT_FREQ", "DAT_WTS", "DAT_OFFS", "DAT_SCL", "DATA", "PERIOD"};
    char *tform_subint[] = {"1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1E", "1E", "1E", "1E", "1E", dummy_txt, dummy_txt, dummy_txt2, dummy_txt2, dummy_txt3, "1D"};
    char *tunit_subint[] = {"", "s", "s", "s", "deg", "deg", "deg", "deg", "deg", "deg", "deg", "deg", "deg", "MHz", "", "", "", "Jy", "s"};
    if(fits_create_tbl(datafile->fits_fptr, BINARY_TBL, datafile->NrSubints, 19, ttypes_subint, tform_subint, tunit_subint, "SUBINT", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot create table.");
      return 0;
    }
    sprintf(dummy_txt, "TIME");
    if(fits_write_key(datafile->fits_fptr, TSTRING, "INT_TYPE", dummy_txt, "Time axis (TIME, BINPHSPERI, BINLNGASC, etc)", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    sprintf(dummy_txt, "SEC");
    if(fits_write_key(datafile->fits_fptr, TSTRING, "INT_UNIT", dummy_txt, "Unit of time axis (SEC, PHS (0-1), DEG)", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    if(datafile->poltype == POLTYPE_STOKES) {
      sprintf(dummy_txt, "STOKES");
    }else if(datafile->poltype == POLTYPE_COHERENCY) {
      sprintf(dummy_txt, "AABBCRCI");
    }else if(datafile->poltype == POLTYPE_ILVPAdPA) {
      sprintf(dummy_txt, "ILVPAdPA");
    }else if(datafile->poltype == POLTYPE_ILVPAdPATEldEl) {
      sprintf(dummy_txt, "ILVPATEL");
    }else if(datafile->poltype == POLTYPE_PAdPA) {
      sprintf(dummy_txt, "PAdPA");
    }else {
      sprintf(dummy_txt, "UNKNOWN");
    }
    if(fits_write_key(datafile->fits_fptr, TSTRING, "POL_TYPE", dummy_txt, "Polarisation identifier (e.g., AABBCRCI, AA+BB)", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    if(datafile->poltype != POLTYPE_STOKES && datafile->poltype != POLTYPE_COHERENCY && datafile->poltype != POLTYPE_ILVPAdPA && datafile->poltype != POLTYPE_PAdPA && datafile->poltype != POLTYPE_ILVPAdPATEldEl) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING writePSRFITSHeader: poltype is not set.");
    }
    dummy_int = datafile->NrPols;
    if(fits_write_key(datafile->fits_fptr, TINT, "NPOL", &dummy_int, "Nr of polarisations", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    dummy_int = 0;
    if(fits_write_key(datafile->fits_fptr, TINT, "NCHNOFFS", &dummy_int, "Channel/sub-band offset for split files", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    if(fits_write_key(datafile->fits_fptr, TINT, "NSUBOFFS", &dummy_int, "Subint offset (Contiguous SEARCH-mode files)", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    dummy_double = datafile->fixedtsamp;
    if(fits_write_key(datafile->fits_fptr, TDOUBLE, "TBIN", &dummy_double, "[s] Time per bin or sample", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword TBIN=%lf.", dummy_double);
      return 0;
    }
    double period_set;
    int ret_prd;
    if(datafile->isFolded) {
      ret_prd = get_period(*datafile, 0, &period_set, verbose);
      if(ret_prd == 2) {
 printerror(verbose.debug, "ERROR writePSRFITSHeader (%s): Cannot obtain period", datafile->filename);
 return 0;
      }
    }else {
      ret_prd = 1;
      period = -1;
    }
    if(datafile->gentype != GENTYPE_SEARCHMODE && datafile->isFolded == 1 && ret_prd == 0 && period_set >= 0) {
      dummy_int = datafile->NrBins;
      if(fits_write_key(datafile->fits_fptr, TINT, "NBIN", &dummy_int, "Nr of bins (PSR/CAL mode; else 1)", &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
 return 0;
      }
    }else {
      dummy_int = 1;
      if(fits_write_key(datafile->fits_fptr, TINT, "NBIN", &dummy_int, "Nr of bins (PSR/CAL mode; else 1)", &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
 return 0;
      }
      dummy_int = datafile->NrBins;
      if(fits_write_key(datafile->fits_fptr, TINT, "NSBLK", &dummy_int, "Samples/row (SEARCH mode, else 1)", &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
 return 0;
      }
      if(datafile->NrBits <= 0 || datafile->NrBits > 8) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING writePSRFITSHeader: Assume you want to write out 8 bit search mode data.");
 datafile->NrBits = 8;
      }
      dummy_int = datafile->NrBits;
      if(fits_write_key(datafile->fits_fptr, TINT, "NBITS", &dummy_int, "Nr of bits/datum (SEARCH mode data, else 1)", &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
 return 0;
      }
    }
    dummy_int = datafile->NrFreqChan;
    if(fits_write_key(datafile->fits_fptr, TINT, "NCHAN", &dummy_int, "Number of channels/sub-bands in this file", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    if(get_channelbandwidth(*datafile, &dummy_double, verbose) == 0) {
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot obtain channel bandwidth.");
      return 0;
    }
    dummy_float = dummy_double;
    if(fits_write_key(datafile->fits_fptr, TFLOAT, "CHAN_BW", &dummy_float, "[MHz] Channel/sub-band width", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    dummy_double = (double)datafile->dm;
    if(fits_write_key(datafile->fits_fptr, TDOUBLE, "DM", &dummy_double, "[cm-3 pc] DM for post-detection dedisperion", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    dummy_double = (double)datafile->rm;
    if(fits_write_key(datafile->fits_fptr, TDOUBLE, "RM", &dummy_double, "[rad m-2] RM for post-detection deFaraday", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    for(i = 0; i < datafile->NrSubints; i++) {
      int ret_prd;
      if(datafile->isFolded) {
 ret_prd = get_period(*datafile, 0, &dummy_double, verbose);
 if(ret_prd == 2) {
   printerror(verbose.debug, "ERROR writePSRFITSHeader (%s): Cannot obtain period", datafile->filename);
   return 0;
 }
      }else {
 ret_prd = -1;
 dummy_double = -1;
      }
      if(ret_prd == 1) {
 dummy_double = -1;
      }
      if(fits_write_col(datafile->fits_fptr, TDOUBLE, 19, i+1, 1, 1, &dummy_double, &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
 return 0;
      }
      if(datafile->tsubMode != TSUBMODE_UNKNOWN) {
 dummy_double = get_tsub(*datafile, i, verbose);
      }else {
 dummy_double = 1.0/datafile->NrSubints;
 if(i == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING writePSRFITSHeader: Observation length seems to be undefined, set the length of the observation to 1 second.");
 }
      }
      if(fits_write_col(datafile->fits_fptr, TDOUBLE, 2, i+1, 1, 1, &dummy_double, &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
 return 0;
      }
      for(j = 0; j < datafile->NrFreqChan; j++) {
 dummy_double = get_weighted_channel_freq(*datafile, i, j, verbose);
 if(fits_write_col(datafile->fits_fptr, TDOUBLE, 14, i+1, 1+j, 1, &dummy_double, &status) != 0) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writePSRFITSHeader: Error writing frequency.");
   fits_report_error(stderr, status);
   return 0;
 }
 dummy_float = 1;
 if(fits_write_col(datafile->fits_fptr, TFLOAT, 15, i+1, 1+j, 1, &dummy_float, &status) != 0) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writePSRFITSHeader: Error writing weights.");
   fits_report_error(stderr, status);
   return 0;
 }
      }
    }
  }else if(datafile->gentype == GENTYPE_RECEIVERMODEL || datafile->gentype == GENTYPE_RECEIVERMODEL2) {
    sprintf(dummy_txt, "%ldE", datafile->NrFreqChan);
    sprintf(dummy_txt2, "%ldE", datafile->NrFreqChan*datafile->NrPols);
    sprintf(dummy_txt3, "%ldJ", datafile->NrFreqChan);
    char *ttypes_feedpar[] = {"DAT_FREQ", "DAT_WTS", "DATA", "DATAERR", "CHISQ", "NFREE"};
    char *tform_feedpar[] = {dummy_txt, dummy_txt, dummy_txt2, dummy_txt2, dummy_txt, dummy_txt3};
    char *tunit_feedpar[] = {"MHz", "", "", "", "", ""};
    if(datafile->gentype == GENTYPE_RECEIVERMODEL2) {
      nrcolumns = 6;
    }else {
      nrcolumns = 4;
    }
    if(fits_create_tbl(datafile->fits_fptr, BINARY_TBL, datafile->NrSubints, nrcolumns, ttypes_feedpar, tform_feedpar, tunit_feedpar, "FEEDPAR", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot create FEEDPAR table.");
      return 0;
    }
    if((datafile->NrPols == 3 && datafile->gentype == GENTYPE_RECEIVERMODEL) || (datafile->NrPols == 5 && datafile->gentype == GENTYPE_RECEIVERMODEL2)) {
      sprintf(dummy_txt, "single");
    }else if((datafile->NrPols == 7 && datafile->gentype == GENTYPE_RECEIVERMODEL) || (datafile->NrPols == 9 && datafile->gentype == GENTYPE_RECEIVERMODEL2)) {
      sprintf(dummy_txt, "van04e18");
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Expect 3 or 7 polarization channels for a receiver model.");
      return 0;
    }
    if(fits_write_key(datafile->fits_fptr, TSTRING, "CAL_MTHD", dummy_txt, "Cross-coupling method", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    dummy_int = datafile->NrPols;
    if(datafile->gentype == GENTYPE_RECEIVERMODEL2)
      dummy_int -= 2;
    if(fits_write_key(datafile->fits_fptr, TINT, "NCPAR", &dummy_int, "Number of coupling parameters", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    dummy_int = datafile->NrFreqChan;
    if(fits_write_key(datafile->fits_fptr, TINT, "NCHAN", &dummy_int, "Nr of channels in Feed coupling data", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    sprintf(dummy_txt, "G");
    if(fits_write_key(datafile->fits_fptr, TSTRING, "PAR_0000", dummy_txt, "scalar gain", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    sprintf(dummy_txt, "gamma");
    if(fits_write_key(datafile->fits_fptr, TSTRING, "PAR_0001", dummy_txt, "differential gain (hyperbolic radians)", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    sprintf(dummy_txt, "phi");
    if(fits_write_key(datafile->fits_fptr, TSTRING, "PAR_0002", dummy_txt, "differential phase (radians)", &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
      return 0;
    }
    if(datafile->NrPols >= 4) {
      sprintf(dummy_txt, "el0");
      if(fits_write_key(datafile->fits_fptr, TSTRING, "PAR_0003", dummy_txt, "ellipticity of receptor 0 (radians)", &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
 return 0;
      }
    }
    if(datafile->NrPols >= 5) {
      sprintf(dummy_txt, "or0");
      if(fits_write_key(datafile->fits_fptr, TSTRING, "PAR_0004", dummy_txt, "orientation of receptor 0 (radians)", &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
 return 0;
      }
    }
    if(datafile->NrPols >= 6) {
      sprintf(dummy_txt, "el1");
      if(fits_write_key(datafile->fits_fptr, TSTRING, "PAR_0005", dummy_txt, "ellipticity of receptor 1 (radians)", &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
 return 0;
      }
    }
    if(datafile->NrPols >= 7) {
      sprintf(dummy_txt, "or1");
      if(fits_write_key(datafile->fits_fptr, TSTRING, "PAR_0006", dummy_txt, "orientation of receptor 1 (radians)", &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writePSRFITSHeader: Cannot write keyword.");
 return 0;
      }
    }
    for(i = 0; i < datafile->NrSubints; i++) {
      for(j = 0; j < datafile->NrFreqChan; j++) {
 dummy_float = get_weighted_channel_freq(*datafile, i, j, verbose);
 if(fits_write_col(datafile->fits_fptr, TFLOAT, 1, i+1, 1+j, 1, &dummy_float, &status) != 0) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writePSRFITSHeader: Error writing frequency.");
   fits_report_error(stderr, status);
   return 0;
 }
 dummy_float = 1;
 if(fits_write_col(datafile->fits_fptr, TFLOAT, 2, i+1, 1+j, 1, &dummy_float, &status) != 0) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writePSRFITSHeader: Error writing weights.");
   fits_report_error(stderr, status);
   return 0;
 }
 dummy_float = 1e-5;
 for(k = 0; k < datafile->NrPols; k++) {
   if(fits_write_col(datafile->fits_fptr, TFLOAT, 4, i+1, datafile->NrPols*(j)+k+1, 1, &dummy_float, &status) != 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR writePSRFITSHeader: Error writing empty errors.");
     fits_report_error(stderr, status);
     return 0;
   }
 }
      }
    }
  }
  return 1;
}
void fits_strip_quotes(char *src, char *dst)
{
  int i;
  dst[0] = 0;
  if(strlen(src) == 0)
    return;
  for(i = 0; i < strlen(src); i++) {
    if(src[i] != '\'' && src[i] != ' ')
      break;
  }
  if(strlen(src) - i > 0) {
    strcpy(dst, src+i);
  }else {
    return;
  }
  for(i = strlen(dst)-1; i >= 0; i--) {
    if(dst[i] == '\'' || dst[i] == ' ') {
      dst[i] = 0;
    }else {
      return;
    }
  }
}
int parse_poltype_psrfits(char *poltypstr, int nowarnings, datafile_definition *datafile, verbose_definition verbose)
{
  long i;
  for(i = strlen(poltypstr)-1; i >= 0; i--) {
    if(i >= 0) {
      if(poltypstr[i] == '\'')
 poltypstr[i] = 0;
    }
  }
  if(strcmp(poltypstr, "STOKES") == 0 || strcmp(poltypstr, "STOKE") == 0 || strcmp(poltypstr, "INTEN") == 0) {
    return POLTYPE_STOKES;
  }else if(strcmp(poltypstr, "IQUV") == 0) {
    return POLTYPE_STOKES;
  }else if(strcmp(poltypstr, "AABBCRCI") == 0) {
    return POLTYPE_COHERENCY;
  }else if(strcmp(poltypstr, "ILVPAdPA") == 0 || strcmp(poltypstr, "ILVPADPA") == 0 || strcmp(poltypstr, "ILVPA") == 0) {
    return POLTYPE_ILVPAdPA;
  }else if(strcmp(poltypstr, "ILVPATEL") == 0 || strcmp(poltypstr, "PAELL+") == 0) {
    return POLTYPE_ILVPAdPATEldEl;
  }else if(strcmp(poltypstr, "PAdPA") == 0 || strcmp(poltypstr, "PADPA") == 0 || strcmp(poltypstr, "PA") == 0) {
    return POLTYPE_PAdPA;
  }else if(strcmp(poltypstr, "UNKNOWN") == 0) {
    return POLTYPE_UNKNOWN;
  }else {
    if(nowarnings == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): POL_TYPE '%s' not recognized", datafile->filename, poltypstr);
    }
    return POLTYPE_UNKNOWN;
  }
}
void update_freqs_using_DAT_FREQ_column(datafile_definition *datafile, verbose_definition verbose)
{
  int colnum, anynul;
  int status = 0;
  if(datafile->freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): The frequency channels should be assumed to be uniform at this point.", datafile->filename);
    exit(0);
  }
  if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "DAT_FREQ", &colnum, &status)) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find DAT_FREQ in table to determine observing frequency.", datafile->filename);
    status = 0;
  }else {
    datafile->freqMode = FREQMODE_FREQTABLE;
    if(datafile->freqlabel_list != NULL) {
      free(datafile->freqlabel_list);
    }
    datafile->freqlabel_list = malloc(datafile->NrSubints*datafile->NrFreqChan*sizeof(double));
    if(datafile->freqlabel_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Memory allocation error.", datafile->filename);
      exit(0);
    }
    long subint, freqchan;
    double freq;
    for(subint = 0; subint < datafile->NrSubints; subint++) {
      for(freqchan = 0; freqchan < datafile->NrFreqChan; freqchan++) {
 if(fits_read_col(datafile->fits_fptr, TDOUBLE, colnum, 1+subint, 1+freqchan, 1, NULL, &(freq), &anynul, &status)) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot read DAT_FREQ in subint table to determine observing frequency.", datafile->filename);
   status = 0;
   datafile->freqMode = FREQMODE_UNIFORM;
   free(datafile->freqlabel_list);
   datafile->freqlabel_list = NULL;
   return;
 }else {
   if(set_weighted_channel_freq(datafile, subint, freqchan, freq, verbose) == 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Error setting (weighted) frequency for channel/subint.", datafile->filename);
     exit(0);
   }
 }
      }
    }
  }
}
int readPSRFITSHeader(datafile_definition *datafile, int readnoscales, int nowarnings, verbose_definition verbose)
{
  char card[FLEN_CARD], comment[FLEN_COMMENT], value[FLEN_VALUE], value_tmp[FLEN_VALUE];
  int status = 0;
  int hdutype, nkeys, i, colnum, anynul, dmnotset, rmnotset, periodnotset, gentypenotset, tsubnotset, tsubsettotobs;
  int nodata, issearch, samptimenotset, dummy_int, ret;
  int version_major, version_minor;
  long numrows, dummy_long;
  double f, dummy_double, f1, f2, bandwidth;
  char command[4000], tmpfilename[2000];
  char dummy_txt[2000];
  FILE *fin;
  char *char_ptrptr[1];
  nodata = 0;
  dmnotset = 0;
  rmnotset = 0;
  periodnotset = 1;
  samptimenotset = 1;
  gentypenotset = 1;
  tsubnotset = 1;
  issearch = 0;
  tsubsettotobs = 0;
  version_major = 0;
  version_minor = 0;
  datafile->NrPols = -1;
  datafile->NrBins = -1;
  datafile->NrFreqChan = -1;
  bandwidth = 0.0;
  if(fits_movabs_hdu(datafile->fits_fptr, 1, &hdutype, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot move to primary HDU", datafile->filename);
    return 0;
  }
  fits_get_hdrspace(datafile->fits_fptr, &nkeys, NULL, &status);
  if(fits_read_card(datafile->fits_fptr, "HDRVER", card, &status)) {
    fflush(stdout);
    if(nowarnings == 0) {
      printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): HDRVER keyword does not exist. This could mean the version of the psrfits file is < 2.13, or this data is not in PSRFITS format.", datafile->filename);
    }
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    ret = sscanf(value, "'%d.%d", &version_major, &version_minor);
    if(ret != 2) {
      printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): HDRVER = '%s' cannot be interpreted as a major.minor version number.", datafile->filename, value);
      version_major = 0;
      version_minor = 0;
    }else if (verbose.debug) {
      printf("    Version is %d.%d\n", version_major, version_minor);
    }
  }
  if(fits_read_card(datafile->fits_fptr, "FITSTYPE", card, &status)) {
    fflush(stdout);
    if(nowarnings == 0) {
      if((version_major > 2) || (version_major == 2 && version_minor >= 13)) {
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): FITSTYPE keyword does not exist, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data.", datafile->filename);
      }else if(version_major == 0 && version_minor == 0) {
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): FITSTYPE keyword does not exist. Not sure this is a PSRFITS file.", datafile->filename);
      }
    }
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    while(value[strlen(value)-2] == ' ' && strlen(value) > 2) {
      value[strlen(value)-2] = '\'';
      value[strlen(value)-1] = 0;
    }
    if(strcmp(value, "'PSRFITS'") != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): This is not a PSRFITS file (%s)", datafile->filename, value);
      return 0;
    }
  }
  if(fits_read_card(datafile->fits_fptr, "TELESCOP", card, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): TELESCOP keyword does not exist", datafile->filename);
    return 0;
  }
  fits_parse_value(card, value, comment, &status);
  fits_strip_quotes(value, value_tmp);
  if(set_observatory_PSRData(datafile, value_tmp, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSHeader: Setting observatory name failed.");
    return 0;
  }
  if(fits_read_card(datafile->fits_fptr,"ANT_X", card, &status)) {
    if(nowarnings == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): ANT_X keyword does not exist", datafile->filename);
    }
    datafile->telescope_X = datafile->telescope_Y = datafile->telescope_Z = 0;
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%lf", &(datafile->telescope_X));
    if(fits_read_card(datafile->fits_fptr,"ANT_Y", card, &status)) {
      if(nowarnings == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): ANT_Y keyword does not exist", datafile->filename);
      }
      datafile->telescope_X = datafile->telescope_Y = datafile->telescope_Z = 0;
      status = 0;
    }else {
      fits_parse_value(card, value, comment, &status);
      sscanf(value, "%lf", &(datafile->telescope_Y));
      if(fits_read_card(datafile->fits_fptr,"ANT_Z", card, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): ANT_Z keyword does not exist", datafile->filename);
 }
 datafile->telescope_X = datafile->telescope_Y = datafile->telescope_Z = 0;
 status = 0;
      }else {
 fits_parse_value(card, value, comment, &status);
 sscanf(value, "%lf", &(datafile->telescope_Z));
      }
    }
  }
  if (fits_read_card(datafile->fits_fptr,"SRC_NAME", card, &status)) {
    fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): SRC_NAME keyword does not exist", datafile->filename);
    return 0;
  }
  fits_parse_value(card, value, comment, &status);
  fits_strip_quotes(value, value_tmp);
  if(set_psrname_PSRData(datafile, value_tmp, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSHeader: Setting pulsar name failed.");
    return 0;
  }
  if (fits_read_card(datafile->fits_fptr,"BACKEND", card, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): BACKEND keyword does not exist", datafile->filename);
    return 0;
  }
  fits_parse_value(card, value, comment, &status);
  fits_strip_quotes(value, value_tmp);
  if(set_instrument_PSRData(datafile, value_tmp, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readFITSHeader: Setting instrument name failed.");
    return 0;
  }
  if(fits_read_card(datafile->fits_fptr,"OBSERVER", card, &status)) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): OBSERVER keyword does not exist", datafile->filename);
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    fits_strip_quotes(value, value_tmp);
    if(set_observer_PSRData(datafile, value_tmp, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRFITSHeader: Setting name observer failed.");
      return 0;
    }
  }
  if(fits_read_card(datafile->fits_fptr,"PROJID", card, &status)) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): PROJID keyword does not exist", datafile->filename);
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    fits_strip_quotes(value, value_tmp);
    if(set_projectID_PSRData(datafile, value_tmp, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRFITSHeader: Setting project ID failed.");
      return 0;
    }
  }
  if(fits_read_card(datafile->fits_fptr, "RA", card, &status)) {
    fflush(stdout);
    if(nowarnings == 0) {
      if((version_major > 2) || (version_major == 2 && version_minor >= 8)) {
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): RA keyword does not exist, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data.", datafile->filename);
      }else if(version_major == 0 && version_minor == 0) {
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): RA keyword does not exist. Not sure this is a PSRFITS file.", datafile->filename);
      }
    }
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    converthms(value+1, &(datafile->ra));
    datafile->ra *= M_PI/12.0;
  }
  if(fits_read_card(datafile->fits_fptr, "DEC", card, &status)) {
    fflush(stdout);
    if(nowarnings == 0) {
      if((version_major > 2) || (version_major == 2 && version_minor >= 8)) {
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): DEC keyword does not exist, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data.", datafile->filename);
      }else if(version_major == 0 && version_minor == 0) {
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): DEC keyword does not exist. Not sure this is a PSRFITS file.", datafile->filename);
      }
    }
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    converthms(value+1, &(datafile->dec));
    datafile->dec *= M_PI/180.0;
  }
  if(fits_read_card(datafile->fits_fptr,"OBSFREQ", card, &status)) {
    fflush(stdout);
    if(nowarnings == 0) {
      if((version_major > 2) || (version_major == 2 && version_minor >= 8)) {
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): OBSFREQ keyword does not exist, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data.", datafile->filename);
      }else if(version_major == 0 && version_minor == 0) {
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): OBSFREQ keyword does not exist. Not sure this is a PSRFITS file.", datafile->filename);
      }
    }
    status = 0;
    dummy_double = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    ret = sscanf(value, "%lf", &dummy_double);
    if(ret != 1) {
      dummy_double = 0;
    }
  }
  datafile->freqMode = FREQMODE_UNIFORM;
  if(datafile->freqlabel_list != NULL) {
    free(datafile->freqlabel_list);
    datafile->freqlabel_list = NULL;
  }
  set_centre_frequency(datafile, dummy_double, verbose);
  if(dummy_double < 1 && (dummy_double < -1.1 || dummy_double > -0.9)) {
    if(nowarnings == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): OBSFREQ keyword is not defined or zero", datafile->filename);
    }
  }else {
    datafile->freq_ref = dummy_double;
  }
  if(verbose.debug)
    printf("  readPSRFITSHeader (%s): OBSFREQ = %lf (-1 or 1e10 is interpreted as infinity)\n", datafile->filename, dummy_double);
  if(fits_read_card(datafile->fits_fptr,"FD_POLN", card, &status)) {
    if(nowarnings == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): FD_POLN keyword does not exist", datafile->filename);
    }
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "'%s'", dummy_txt);
    if(strcmp(dummy_txt, "CIRC") == 0) {
      datafile->feedtype = FEEDTYPE_CIRCULAR;
    }else if(strcmp(dummy_txt, "LIN") == 0) {
      datafile->feedtype = FEEDTYPE_LINEAR;
    }else {
      if(nowarnings == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): FD_POLN '%s' not recognized", datafile->filename, dummy_txt);
      }
    }
  }
  if (fits_read_card(datafile->fits_fptr,"FD_HAND", card, &status)) {
    if(nowarnings == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): FD_HAND keyword does not exist", datafile->filename);
    }
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &i);
    if(i == 1 || i == -1)
      datafile->feedtype *= i;
    else {
      if(nowarnings == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): FD_HAND '%d' not recognized", datafile->filename, i);
      }
    }
  }
  if (fits_read_card(datafile->fits_fptr,"STT_IMJD", card, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): STT_IMJD keyword does not exist", datafile->filename);
    return 0;
  }
  fits_parse_value(card, value, comment, &status);
  sscanf(value, "%Lf", &(datafile->mjd_start));
  if (fits_read_card(datafile->fits_fptr,"STT_SMJD", card, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): STT_SMJD keyword does not exist", datafile->filename);
    return 0;
  }
  fits_parse_value(card, value, comment, &status);
  sscanf(value, "%lf", &f);
  datafile->mjd_start += f/(24.0*3600.0);
  if (fits_read_card(datafile->fits_fptr,"STT_OFFS", card, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): STT_OFFS keyword does not exist", datafile->filename);
    return 0;
  }
  fits_parse_value(card, value, comment, &status);
  sscanf(value, "%lf", &f);
  datafile->mjd_start += f/(24.0*3600.0);
  if(fits_read_card(datafile->fits_fptr,"OBS_MODE", card, &status)) {
    if(nowarnings == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): OBS_MODE keyword does not exist", datafile->filename);
    }
    status = 0;
  }
  fits_parse_value(card, value, comment, &status);
  if(strstr(value, "SEARCH") != 0) {
    fflush(stdout);
    if(verbose.debug)
      printf("  readPSRFITSHeader (%s): File is in search mode\n", datafile->filename);
    periodnotset = 0;
    issearch = 1;
    datafile->isFolded = 0;
    datafile->foldMode = FOLDMODE_UNKNOWN;
    datafile->fixedPeriod = -1;
    datafile->gentype = GENTYPE_SEARCHMODE;
    gentypenotset = 0;
  }else if(strstr(value, "PSR") != 0) {
    fflush(stdout);
    if(verbose.debug)
      printf("  readPSRFITSHeader (%s): File is in pulsar fold mode\n", datafile->filename);
    periodnotset = 1;
    datafile->isFolded = 1;
    datafile->foldMode = FOLDMODE_FIXEDPERIOD;
    datafile->fixedPeriod = 0;
  }else if(strstr(value, "CAL") != 0) {
    fflush(stdout);
    if(verbose.debug)
      printf("  readPSRFITSHeader (%s): File is in calibration mode\n", datafile->filename);
    periodnotset = 1;
    datafile->isFolded = 1;
    datafile->foldMode = FOLDMODE_FIXEDPERIOD;
    datafile->fixedPeriod = 0;
  }else if(strstr(value, "PCM") != 0) {
    fflush(stdout);
    if(verbose.debug)
      printf("  readPSRFITSHeader (%s): File appears to be a calibration solution\n", datafile->filename);
    periodnotset = 1;
    datafile->isFolded = 1;
    datafile->foldMode = FOLDMODE_FIXEDPERIOD;
    datafile->fixedPeriod = 0;
  }else {
    if(nowarnings == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): OBS_MODE keyword does not appear to be valid. It is set to %s.", datafile->filename, value);
    }
    periodnotset = 1;
    datafile->isFolded = 1;
    datafile->foldMode = FOLDMODE_FIXEDPERIOD;
    datafile->fixedPeriod = 0;
  }
  if(fits_read_card(datafile->fits_fptr,"PS_GTYPE", card, &status)) {
    if(issearch == 0) {
      datafile->gentype = GENTYPE_UNDEFINED;
      status = 0;
      if(verbose.debug) {
 fflush(stdout);
 printf("  readPSRFITSHeader (%s): gentype not defined in header\n", datafile->filename);
      }
    }
  }else {
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &datafile->gentype);
    gentypenotset = 0;
    if(verbose.debug) {
      fflush(stdout);
      printf("  readPSRFITSHeader (%s): gentype read in as: %s\n", datafile->filename, returnGenType_str(datafile->gentype));
    }
  }
  datafile->tsubMode = TSUBMODE_UNKNOWN;
  if(datafile->tsub_list != NULL) {
    free(datafile->tsub_list);
  }
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRFITSHeader: Memory allocation error");
    return 0;
  }
  if(verbose.debug)
    printf("  readPSRFITSHeader (%s): Initialise tsub_list\n", datafile->filename);
  datafile->tsub_list[0] = 0;
  if(fits_read_card(datafile->fits_fptr,"PS_XMIN", card, &status)) {
    datafile->xrangeset = 0;
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    datafile->xrangeset = 1;
    sscanf(value, "%f", &datafile->xrange[0]);
  }
  if (fits_read_card(datafile->fits_fptr,"PS_XMAX", card, &status)) {
    datafile->xrangeset = 0;
    status = 0;
  }else if (datafile->xrangeset) {
    fits_parse_value(card, value, comment, &status);
    datafile->xrangeset = 1;
    sscanf(value, "%f", &datafile->xrange[1]);
  }
  if (fits_read_card(datafile->fits_fptr,"PS_YMIN", card, &status)) {
    datafile->yrangeset = 0;
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    datafile->yrangeset = 1;
    sscanf(value, "%f", &datafile->yrange[0]);
  }
  if (fits_read_card(datafile->fits_fptr,"PS_YMAX", card, &status)) {
    datafile->yrangeset = 0;
    status = 0;
  }else if (datafile->yrangeset) {
    fits_parse_value(card, value, comment, &status);
    datafile->yrangeset = 1;
    sscanf(value, "%f", &datafile->yrange[1]);
  }
  if (fits_read_card(datafile->fits_fptr,"PS_CSW", card, &status)) {
    datafile->cableSwap = -1;
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &dummy_int);
    datafile->cableSwap = dummy_int;
  }
  if (fits_read_card(datafile->fits_fptr,"PS_CSWC", card, &status)) {
    datafile->cableSwapcor = -1;
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &dummy_int);
    datafile->cableSwapcor = dummy_int;
  }
  if(fits_read_card(datafile->fits_fptr,"PS_DBASE", card, &status)) {
    datafile->isDebase = -1;
    status = 0;
  }else {
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &dummy_int);
    datafile->isDebase = dummy_int;
  }
  if(fits_movnam_hdu(datafile->fits_fptr, BINARY_TBL, "HISTORY", 0, &status)) {
    if(nowarnings == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find dedispersion and de-Faraday rotation state in file.", datafile->filename);
    }
    datafile->isDeDisp = -1;
    datafile->isDeFarad = -1;
    datafile->isDePar = -1;
    if(periodnotset) {
      if(verbose.debug) {
 fflush(stdout);
 printf("  readPSRFITSHeader (%s): No HISTORY table - so cannot determine period from HISTORY:TBIN/NBIN_PRD.\n", datafile->filename);
      }
    }
    status = 0;
  }else {
    fits_get_num_rows (datafile->fits_fptr, &numrows, &status);
    if(numrows == 0 || fits_get_colnum (datafile->fits_fptr, CASEINSEN, "DEDISP", &colnum, &status)) {
      if(nowarnings == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find dedispersion state in history table.", datafile->filename);
      }
      datafile->isDeDisp = -1;
      status = 0;
    }else {
      if(fits_read_col(datafile->fits_fptr, TINT, colnum, numrows, 1, 1, NULL, &(dummy_int), &anynul, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find dedispersion state in history table.", datafile->filename);
 }
 datafile->isDeDisp = -1;
 status = 0;
      }else {
 datafile->isDeDisp = dummy_int;
      }
    }
    if(numrows == 0 || fits_get_colnum (datafile->fits_fptr, CASEINSEN, "RM_CORR", &colnum, &status)) {
      if(nowarnings == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find de-Faraday rotation state in history table.", datafile->filename);
      }
      datafile->isDeFarad = -1;
      status = 0;
    }else {
      if(fits_read_col(datafile->fits_fptr, TINT, colnum, numrows, 1, 1, NULL, &(dummy_int), &anynul, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find de-Faraday rotation state in history table.", datafile->filename);
 }
 datafile->isDeFarad = -1;
 status = 0;
      }else {
 datafile->isDeFarad = dummy_int;
      }
    }
    int ok = 0;
    if(numrows != 0) {
      ret = fits_get_colnum(datafile->fits_fptr, CASEINSEN, "PR_CORR", &colnum, &status);
      status = 0;
      if(ret == 0) {
 ok = 1;
      }else {
 ret = fits_get_colnum(datafile->fits_fptr, CASEINSEN, "PAR_CORR", &colnum, &status);
 if(ret == 0) {
   ok = 1;
 }
      }
    }
    if(ok == 0) {
      if(nowarnings == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find parallactic angle state in history table.", datafile->filename);
      }
      datafile->isDePar = -1;
      status = 0;
    }else {
      if(fits_read_col(datafile->fits_fptr, TINT, colnum, numrows, 1, 1, NULL, &(dummy_int), &anynul, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find parallactic angle state in history table.", datafile->filename);
 }
 datafile->isDePar = -1;
 status = 0;
      }else {
 datafile->isDePar = dummy_int;
      }
    }
    if(numrows == 0 || fits_get_colnum (datafile->fits_fptr, CASEINSEN, "NPOL", &colnum, &status)) {
      if(nowarnings == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find number polarizations in the history table.", datafile->filename);
      }
      datafile->NrPols = -1;
      status = 0;
    }else {
      if(fits_read_col(datafile->fits_fptr, TINT, colnum, numrows, 1, 1, NULL, &(dummy_int), &anynul, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find NPOL om the history table.", datafile->filename);
 }
 datafile->NrPols = -1;
 status = 0;
      }else {
 datafile->NrPols = dummy_int;
      }
    }
    if(numrows == 0 || fits_get_colnum (datafile->fits_fptr, CASEINSEN, "NBIN", &colnum, &status)) {
      if(nowarnings == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find number of bins in the history table.", datafile->filename);
      }
      datafile->NrBins = -1;
      status = 0;
    }else {
      if(fits_read_col(datafile->fits_fptr, TINT, colnum, numrows, 1, 1, NULL, &(dummy_int), &anynul, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find NBIN in the history table.", datafile->filename);
 }
 datafile->NrBins = -1;
 status = 0;
      }else {
 datafile->NrBins = dummy_int;
      }
    }
    if(numrows == 0 || fits_get_colnum (datafile->fits_fptr, CASEINSEN, "NCHAN", &colnum, &status)) {
      if(nowarnings == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find number of frequency channels in the history table.", datafile->filename);
      }
      datafile->NrFreqChan = -1;
      status = 0;
    }else {
      if(fits_read_col(datafile->fits_fptr, TINT, colnum, numrows, 1, 1, NULL, &(dummy_int), &anynul, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find NCHAN in the history table.", datafile->filename);
 }
 datafile->NrFreqChan = -1;
 status = 0;
      }else {
 datafile->NrFreqChan = dummy_int;
      }
    }
    if(datafile->NrFreqChan != -1) {
      if(numrows == 0 || fits_get_colnum (datafile->fits_fptr, CASEINSEN, "CHAN_BW", &colnum, &status)) {
 fflush(stdout);
 if(nowarnings == 0) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): CHAN_BW keyword does not exist in history table. Not sure this is a PSRFITS file.", datafile->filename);
 }
 bandwidth = 0.0;
 status = 0;
      }else {
 if(fits_read_col(datafile->fits_fptr, TDOUBLE, colnum, numrows, 1, 1, NULL, &(bandwidth), &anynul, &status)) {
   if(nowarnings == 0) {
     fflush(stdout);
     printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find CHAN_BW in the history table.", datafile->filename);
   }
   bandwidth = 0.0;
   status = 0;
 }else {
   bandwidth *= datafile->NrFreqChan;
 }
      }
    }
    if(numrows == 0 || fits_get_colnum (datafile->fits_fptr, CASEINSEN, "POL_TYPE", &colnum, &status)) {
      fflush(stdout);
      if(nowarnings == 0) {
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): POL_TYPE keyword does not exist in history table. Not sure this is a PSRFITS file.", datafile->filename);
      }
      datafile->poltype = POLTYPE_UNKNOWN;
      status = 0;
    }else {
      char_ptrptr[0] = dummy_txt;
      if(fits_read_col(datafile->fits_fptr, TSTRING, colnum, numrows, 1, 1, NULL, char_ptrptr, &anynul, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot read POL_TYPE from the history table.", datafile->filename);
 }
 datafile->poltype = POLTYPE_UNKNOWN;
 status = 0;
      }else {
 datafile->poltype = parse_poltype_psrfits(dummy_txt, nowarnings, datafile, verbose);
      }
    }
    if(periodnotset) {
      if(numrows == 0 || fits_get_colnum (datafile->fits_fptr, CASEINSEN, "TBIN", &colnum, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find TBIN in history table to obtain period.", datafile->filename);
 }
 periodnotset = 1;
 status = 0;
      }else {
 if(fits_read_col(datafile->fits_fptr, TDOUBLE, colnum, numrows, 1, 1, NULL, &(datafile->fixedPeriod), &anynul, &status)) {
   if(nowarnings == 0) {
     fflush(stdout);
     printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot read TBIN in history table to obtain period.", datafile->filename);
   }
   periodnotset = 1;
   status = 0;
 }else {
   double period;
   if(get_period(*datafile, 0, &period, verbose) == 2) {
     printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot obtain period", datafile->filename);
     return 0;
   }
   if(verbose.debug) {
     fflush(stdout);
     printf("  readPSRFITSHeader (%s): Period determination - Set sampling time to HISTORY:TBIN = %lf.\n", datafile->filename, period);
   }
   datafile->isFolded = 1;
   datafile->foldMode = FOLDMODE_FIXEDPERIOD;
   periodnotset = 0;
   datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
   datafile->fixedtsamp = period;
   samptimenotset = 0;
 }
      }
      if(periodnotset == 0) {
 int nbin_prd;
 if(numrows == 0 || fits_get_colnum (datafile->fits_fptr, CASEINSEN, "NBIN_PRD", &colnum, &status)) {
   if(verbose.debug) {
     fflush(stdout);
     printf("  readPSRFITSHeader (%s): Period determination - No HISTORY:NBIN_PRD to obtain period.\n", datafile->filename);
   }
   periodnotset = 1;
   status = 0;
 }else {
   if(fits_read_col(datafile->fits_fptr, TINT, colnum, numrows, 1, 1, NULL, &nbin_prd, &anynul, &status)) {
     if(verbose.debug) {
       fflush(stdout);
       printf("  readPSRFITSHeader (%s): Period determination - No HISTORY:NBIN_PRD to obtain period.\n", datafile->filename);
     }
     periodnotset = 1;
     status = 0;
   }else {
     if(nbin_prd == 0) {
       if(verbose.debug) {
  fflush(stdout);
  printf("  readPSRFITSHeader (%s): Period determination - HISTORY:NBIN_PRD = 0, suggest search-mode data?\n", datafile->filename);
       }
       periodnotset = 1;
     }else {
       datafile->fixedPeriod *= nbin_prd;
       datafile->isFolded = 1;
       datafile->foldMode = FOLDMODE_FIXEDPERIOD;
       if(verbose.debug) {
  fflush(stdout);
  double period;
  if(get_period(*datafile, 0, &period, verbose) == 2) {
    printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot obtain period", datafile->filename);
    return 0;
  }
  printf("  readPSRFITSHeader (%s): Period determination - multiplying HISTORY:TBIN with HISTORY:NBIN_PRD suggests %lf.\n", datafile->filename, period);
       }
     }
   }
 }
      }
    }
    if(numrows == 0 || fits_get_colnum (datafile->fits_fptr, CASEINSEN, "CTR_FREQ", &colnum, &status)) {
      if(nowarnings == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find centre frequency in history table.", datafile->filename);
      }
      status = 0;
    }else {
      double freq_history;
      if(fits_read_col(datafile->fits_fptr, TDOUBLE, colnum, numrows, 1, 1, NULL, &(freq_history), &anynul, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find centre frequency in history table.", datafile->filename);
 }
 set_centre_frequency(datafile, 0.0, verbose);
 status = 0;
      }
      if((freq_history > -1.01 && freq_history < -0.99)) {
 if(verbose.debug)
   printf("  readPSRFITSHeader (%s): Frequency defined in history table = %lf is interpreted as infinity\n", datafile->filename, freq_history);
 freq_history = 1e10;
      }
      if(get_centre_frequency(*datafile, verbose) < 1) {
 set_centre_frequency(datafile, freq_history, verbose);
 if(verbose.debug)
   printf("  readPSRFITSHeader (%s): Use history table frequency as centre frequency.\n", datafile->filename);
      }
      if(datafile->freq_ref < 1) {
 datafile->freq_ref = freq_history;
 if(verbose.debug)
   printf("  readPSRFITSHeader (%s): Use history table frequency as reference frequency.\n", datafile->filename);
      }else if(fabs(datafile->freq_ref - freq_history) > 0.001) {
 if(nowarnings == 0) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): OBSFREQ != Frequency defined in history table. Reference frequency set to %lf MHz", datafile->filename, freq_history);
 }
 datafile->freq_ref = freq_history;
      }
    }
  }
  if(get_centre_frequency(*datafile, verbose) < 1) {
    if(verbose.debug)
      printf("  readPSRFITSHeader (%s): Try to get centre frequency from polyco table\n", datafile->filename);
    if(fits_movnam_hdu(datafile->fits_fptr, BINARY_TBL, "POLYCO", 0, &status)) {
      if(verbose.debug)
 printf("  readPSRFITSHeader (%s): POLYCO table does not exist\n", datafile->filename);
      status = 0;
    }else {
      fits_get_hdrspace(datafile->fits_fptr, &nkeys, NULL, &status);
      fflush(stdout);
      if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "REF_FREQ", &colnum, &status)) {
 if(verbose.debug)
   printf("  readPSRFITSHeader (%s): No REF_FRQ in polyco table\n", datafile->filename);
      }else {
 if(fits_read_col(datafile->fits_fptr, TDOUBLE, colnum, 1, 1, 1, NULL, &dummy_double, &anynul, &status)) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot read column", datafile->filename);
   status = 0;
 }else {
   set_centre_frequency(datafile, dummy_double, verbose);
   if(nowarnings == 0) {
     printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): FREQUENCY=%f according to polyco", datafile->filename, dummy_double);
   }
 }
      }
    }
  }
  if(get_centre_frequency(*datafile, verbose) < 1) {
    if(verbose.debug)
      printf("  readPSRFITSHeader (%s): Try to get centre frequency from tempo2 predictor table\n", datafile->filename);
    if(fits_movnam_hdu(datafile->fits_fptr, BINARY_TBL, "T2PREDICT", 0, &status)) {
      if(verbose.debug)
 printf("  readPSRFITSHeader (%s): T2PREDICT table does not exist\n", datafile->filename);
      status = 0;
    }else {
      fits_get_hdrspace(datafile->fits_fptr, &nkeys, NULL, &status);
      fflush(stdout);
      fits_get_num_rows (datafile->fits_fptr, &numrows, &status);
      if(numrows == 0) {
 if(nowarnings == 0) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): No rows in T2predictor", datafile->filename);
 }
 status = 0;
      }else {
 if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "PREDICT", &colnum, &status)) {
   if(nowarnings == 0) {
     printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): PREDICT keyword does not exist", datafile->filename);
   }
   status = 0;
 }
 sprintf(command, "XXX");
 char_ptrptr[0] = command;
 for(i = 1; i <= numrows; i++) {
   fits_read_col (datafile->fits_fptr, TSTRING, colnum, i, 1, 1, 0, char_ptrptr, &anynul, &status);
   sscanf(command, "%s", dummy_txt);
   if(strcmp(dummy_txt, "FREQ_RANGE") == 0) {
     sscanf(command, "%s %lf %lf", dummy_txt, &f1, &f2);
     set_centre_frequency(datafile, 0.5*(f1+f2), verbose);
     if(nowarnings == 0) {
       printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): FREQUENCY=%lf according to tempo2 predictor", datafile->filename, get_centre_frequency(*datafile, verbose));
     }
   }
 }
      }
    }
  }
  status = 0;
  if(fits_movnam_hdu(datafile->fits_fptr, BINARY_TBL, "SUBINT", 0, &status)) {
    if(verbose.verbose)
      printf("  readPSRFITSHeader (%s): SUBINT table does not exist!\n", datafile->filename);
    status = 0;
    nodata = 1;
    if(periodnotset) {
      if(verbose.debug) {
 fflush(stdout);
 printf("  readPSRFITSHeader (%s): No SUBINT table - so cannot determine period from SUBINT:PERIOD.\n", datafile->filename);
      }
    }
  }else {
    int ncols;
    long nrows;
    if(periodnotset) {
      fits_get_num_rows(datafile->fits_fptr, &nrows, &status);
      fits_get_num_cols(datafile->fits_fptr, &ncols, &status);
      if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "PERIOD", &colnum, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find PERIOD in subint table, will use PSRCHIVE to obtain period.", datafile->filename);
 }
 status = 0;
 periodnotset = 1;
      }else {
 if(fits_read_col(datafile->fits_fptr, TDOUBLE, colnum, nrows, 1, 1, NULL, &(datafile->fixedPeriod), &anynul, &status)) {
   if(nowarnings == 0) {
     fflush(stdout);
     printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot read PERIOD in subint table, will use PSRCHIVE to obtain period.", datafile->filename);
   }
   periodnotset = 1;
   status = 0;
 }else {
   periodnotset = 0;
   datafile->isFolded = 1;
   datafile->foldMode = FOLDMODE_FIXEDPERIOD;
   if(verbose.debug) {
     fflush(stdout);
     double period;
     if(get_period(*datafile, 0, &period, verbose) == 2) {
       printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot obtain period", datafile->filename);
       return 0;
     }
     printf("  readPSRFITSHeader (%s): Period determination - SUBINT:PERIOD suggests %lf.\n", datafile->filename, period);
   }
 }
      }
    }
    if(periodnotset) {
      tmpfilename[0] = 0;
      char *username;
      if(getUsername(&username, verbose) == 0) {
 username = malloc(8);
 if(username == NULL) {
   printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Memory allocation error", datafile->filename);
   return 0;
 }
 sprintf(username, "Unknown");
      }
      sprintf(tmpfilename, "%sjunk_mjd_%s", PSRSALSA_TEMPDIRECTORY, username);
      free(username);
      remove(tmpfilename);
      sprintf(command, "vap -n -c period %s > %s", datafile->filename, tmpfilename);
      if(verbose.verbose) printf("  Running: %s\n", command);
      system(command);
      fin = fopen(tmpfilename, "r");
      if(fin == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot open %s", datafile->filename, tmpfilename);
 return 0;
      }
      i = fscanf(fin, "%s %lf", command, &(datafile->fixedPeriod));
      if(verbose.debug) {
 fflush(stdout);
 double period;
 if(get_period(*datafile, 0, &period, verbose) == 2) {
   printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot obtain period", datafile->filename);
   return 0;
 }
 printf("  readPSRFITSHeader (%s): System call to vap suggests period = %lf.\n", datafile->filename, period);
      }
      if(i != 2) {
 datafile->isFolded = 0;
 datafile->foldMode = FOLDMODE_UNKNOWN;
 datafile->fixedPeriod = -1;
 datafile->gentype = GENTYPE_SEARCHMODE;
 if(verbose.verbose) printf("readPSRFITSHeader (%s): Cannot parse %s\nAssume PSRFITS file is in search mode (or it is a receiver model).\n", datafile->filename, tmpfilename);
      }else {
 datafile->isFolded = 1;
 datafile->foldMode = FOLDMODE_FIXEDPERIOD;
      }
      fclose(fin);
      remove(tmpfilename);
    }
    fits_get_hdrspace(datafile->fits_fptr, &nkeys, NULL, &status);
    if(fits_read_card(datafile->fits_fptr, "NPOL", card, &status)) {
      fflush(stdout);
      if(datafile->NrPols == -1) {
 printerror(verbose.debug, "WARNING readPSRFITSHeader (%s): NPOL keyword does not exist in the subint table, and it was not found in a HISTORY table either. This file cannot be interpreted.", datafile->filename);
 return 0;
      }
      if(nowarnings == 0) {
 if((version_major > 2) || (version_major == 2 && version_minor >= 7)) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): NPOL keyword does not exist in the subint table, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data.", datafile->filename);
 }else if(version_major == 0 && version_minor == 0) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): NPOL keyword does not exist in the subint table. Not sure this is a PSRFITS file.", datafile->filename);
 }
      }
      status = 0;
    }else {
      fits_parse_value(card, value, comment, &status);
      long nrpols;
      sscanf(value, "%ld", &(nrpols));
      if(datafile->NrPols == -1) {
 datafile->NrPols = nrpols;
      }else if(datafile->NrPols != nrpols) {
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Updating the %ld polsarization channels as reported in the history table to %ld as suggested by the SUBINT table.", datafile->filename, datafile->NrPols, nrpols);
 datafile->NrPols = nrpols;
      }
    }
    if(fits_read_card(datafile->fits_fptr, "POL_TYPE", card, &status)) {
      if(nowarnings == 0) {
 fflush(stdout);
 if((version_major > 3) || (version_major == 3 && version_minor >= 0)) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): POL_TYPE keyword does not exist, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data.", datafile->filename);
 }else if(version_major == 0 && version_minor == 0) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): POL_TYPE keyword does not exist. Not sure this is a PSRFITS file.", datafile->filename);
 }
      }
      status = 0;
    }else {
      fits_parse_value(card, value, comment, &status);
      sscanf(value, "'%s'", dummy_txt);
    }
    int newpoltype = parse_poltype_psrfits(dummy_txt, nowarnings, datafile, verbose);
    if(datafile->poltype == POLTYPE_UNKNOWN) {
      datafile->poltype = newpoltype;
    }else {
      if(nowarnings == 0) {
 if(datafile->poltype != newpoltype) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Changing POL_TYPE from %d to %d, as suggested by the subint table.", datafile->filename, datafile->poltype, newpoltype);
 }
      }
      datafile->poltype = newpoltype;
    }
    if(fits_read_card(datafile->fits_fptr, "NBIN", card, &status)) {
      fflush(stdout);
      if(datafile->NrBins == -1) {
 printerror(verbose.debug, "WARNING readPSRFITSHeader (%s): NBIN keyword does not exist in the subint table, and it was not found in a HISTORY table either. This file cannot be interpreted.", datafile->filename);
 return 0;
      }
      if(nowarnings == 0) {
 if((version_major > 2) || (version_major == 2 && version_minor >= 7)) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): NBIN keyword does not exist in the subint table, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data.", datafile->filename);
 }else if(version_major == 0 && version_minor == 0) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): NBIN keyword does not exist in the subint table. Not sure this is a PSRFITS file.", datafile->filename);
 }
      }
      status = 0;
    }else {
      fits_parse_value(card, value, comment, &status);
      long nrbins;
      sscanf(value, "%ld", &(nrbins));
      if(datafile->NrBins == -1) {
 datafile->NrBins = nrbins;
      }else if(datafile->NrBins != nrbins) {
 printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Updating the %ld bins as reported in the history table to %ld as suggested by the SUBINT table.", datafile->filename, datafile->NrBins, nrbins);
 datafile->NrBins = nrbins;
      }
    }
    if(datafile->NrBins == 1) {
      if(fits_read_card(datafile->fits_fptr,"NSBLK", card, &status)) {
 if(nowarnings == 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): NSBLK keyword does not exist, while NrBins=1 suggests this is a search-mode file. Proceed with reading file as being in folded mode.", datafile->filename);
 }
 datafile->NrBits = 16;
      }else {
 fits_parse_value(card, value, comment, &status);
 sscanf(value, "%ld", &(datafile->NrBins));
 if (fits_read_card(datafile->fits_fptr,"NBITS", card, &status)) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): NBITS keyword does not exist, while NrBins=1 suggests this is a search-mode file.", datafile->filename);
   return 0;
 }
 fits_parse_value(card, value, comment, &status);
 sscanf(value, "%d", &(datafile->NrBits));
      }
    }else {
      datafile->NrBits = 16;
    }
    double period;
    int ret;
    if(datafile->isFolded == 1) {
      ret = get_period(*datafile, 0, &period, verbose);
      if(ret == 2) {
 printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot obtain period", datafile->filename);
 return 0;
      }
    }else {
      ret = 1;
      period = -1;
    }
    if(fits_read_card(datafile->fits_fptr,"TBIN", card, &status)) {
      if(samptimenotset == 0) {
 if(nowarnings == 0) {
   fflush(stdout);
   if((version_major > 3) || (version_major == 3 && version_minor >= 0)) {
     printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data.", datafile->filename);
   }else if(version_major == 0 && version_minor == 0) {
     printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table. Not sure this is a PSRFITS file.", datafile->filename);
   }
 }
      }else {
 if(periodnotset || period <= 0) {
   if(nowarnings == 0) {
     fflush(stdout);
     if((version_major > 3) || (version_major == 3 && version_minor >= 0)) {
       printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data. Period is not set, so cannot use period to get sampling time.", datafile->filename);
     }else if(version_major == 0 && version_minor == 0) {
       printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table. Not sure this is a PSRFITS file. Period is not set, so cannot use period to get sampling time.", datafile->filename);
     }
   }
 }else {
   if(nowarnings == 0) {
     fflush(stdout);
     if((version_major > 3) || (version_major == 3 && version_minor >= 0)) {
       printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data. Assuming whole period is stored and using period to get sampling time.", datafile->filename);
     }else if(version_major == 0 && version_minor == 0) {
       printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table. Not sure this is a PSRFITS file. Assuming whole period is stored and using period to get sampling time.", datafile->filename);
     }
   }
   datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
   datafile->fixedtsamp = period/(double)datafile->NrBins;
   samptimenotset = 0;
   if(verbose.debug) {
     fflush(stdout);
     printf("  readPSRFITSHeader (%s): Set sampling time to %lf.\n", datafile->filename, datafile->fixedtsamp);
   }
 }
      }
      status = 0;
    }else {
      fits_parse_value(card, value, comment, &status);
      if(value[1] == '*') {
 if(samptimenotset == 0) {
   if(nowarnings == 0) {
     fflush(stdout);
     if((version_major > 3) || (version_major == 3 && version_minor >= 0)) {
       printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data. It is assumed that stored in the history table is correct.", datafile->filename);
     }else if(version_major == 0 && version_minor == 0) {
       printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table. Not sure this is a PSRFITS file. It is assumed that stored in the history table is correct.", datafile->filename);
     }
   }
 }else {
   if(periodnotset || period <= 0) {
     if(nowarnings == 0) {
       fflush(stdout);
       if((version_major > 3) || (version_major == 3 && version_minor >= 0)) {
  printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data. Period is not set, so cannot use period to get sampling time.", datafile->filename);
       }else if(version_major == 0 && version_minor == 0) {
  printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table. Not sure this is a PSRFITS file. Period is not set, so cannot use period to get sampling time.", datafile->filename);
       }
     }
   }else {
     if(nowarnings == 0) {
       fflush(stdout);
       if((version_major > 3) || (version_major == 3 && version_minor >= 0)) {
  printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data. Assuming whole period is stored and using period to get sampling time.", datafile->filename);
       }else if(version_major == 0 && version_minor == 0) {
  printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): TBIN keyword does not exist in the subint table. Not sure this is a PSRFITS file. Assuming whole period is stored and using period to get sampling time.", datafile->filename);
       }
     }
     datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
     datafile->fixedtsamp = period/(double)datafile->NrBins;
     samptimenotset = 0;
     if(verbose.debug) {
       fflush(stdout);
       printf("  readPSRFITSHeader (%s): Set sampling time to %lf.\n", datafile->filename, datafile->fixedtsamp);
     }
   }
 }
      }else {
 if(samptimenotset) {
   datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
   sscanf(value, "%lf", &(datafile->fixedtsamp));
   samptimenotset = 0;
   if(verbose.debug) {
     fflush(stdout);
     printf("  readPSRFITSHeader (%s): Setting sampling time SUBINT:TBIN = %lf.\n", datafile->filename, get_tsamp(*datafile, 0, verbose));
   }
 }
      }
    }
    ret = -1;
    if(fits_read_card(datafile->fits_fptr, "NCHAN", card, &status)) {
      if(nowarnings == 0) {
 fflush(stdout);
 if((version_major > 3) || (version_major == 3 && version_minor >= 2)) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): NCHAN keyword does not exist in subint table, while this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data.", datafile->filename);
 }else if(version_major == 0 && version_minor == 0) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): NCHAN keyword does not exist in subint table, which might be true for old versions of PSRFITS.", datafile->filename);
 }
      }
      status = 0;
      if(fits_read_card(datafile->fits_fptr,"NCH_FILE", card, &status)) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): NCHAN and NCH_FILE keywords do not exist in subint table", datafile->filename);
 return 0;
      }else {
 fits_parse_value(card, value, comment, &status);
 ret = sscanf(value, "%ld", &(dummy_long));
      }
    }else {
      fits_parse_value(card, value, comment, &status);
      ret = sscanf(value, "%ld", &(dummy_long));
    }
    if(ret != -1) {
      if(ret != 1) {
 if(datafile->NrFreqChan > 0) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Parsing NCHAN or NCH_FILE keyword in subint table failed, assume %ld from history table is correct.", datafile->filename, datafile->NrFreqChan);
 }else {
   fflush(stdout);
   printerror(verbose.debug, "WARNING readPSRFITSHeader (%s): Parsing NCHAN or NCH_FILE keyword in subint table failed, and it was not defined in history table either.", datafile->filename);
   return 0;
 }
      }else {
 if(datafile->NrFreqChan == -1) {
   datafile->NrFreqChan = dummy_long;
 }else {
   if(datafile->NrFreqChan != dummy_long) {
     fflush(stdout);
     printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): The NCHAN or NCH_FILE keyword in subint table is different from that in the history table. Updating the number of frequency channels from %ld to %ld.", datafile->filename, datafile->NrFreqChan, dummy_long);
     datafile->NrFreqChan = dummy_long;
   }
 }
      }
    }
    if (fits_read_card(datafile->fits_fptr,"NAXIS2", card, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): NAXIS2 keyword does not exist", datafile->filename);
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%ld", &(datafile->NrSubints));
    if(tsubsettotobs) {
      datafile->tsub_list[0] /= (double)datafile->NrSubints;
      tsubsettotobs = 0;
    }
    if(fits_read_card(datafile->fits_fptr, "CHAN_BW", card, &status)) {
      fflush(stdout);
      if(nowarnings == 0) {
 if((version_major > 3) || (version_major == 3 && version_minor >= 2)) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): CHAN_BW keyword does not exist in the subint table, but this was expected to exist for claimed version of the PSRFITS file. Something might be wrong with the format data.", datafile->filename);
 }
      }
      if(bandwidth == 0.0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): CHAN_BW keyword does not exist in the subint table, and it was not set in the history table either. Header cannot be interpreted.", datafile->filename);
 return 0;
      }
      if(set_bandwidth(datafile, bandwidth, verbose) == 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Changing bandwidth failed.", datafile->filename);
 return 0;
      }
      status = 0;
    }else {
      fits_parse_value(card, value, comment, &status);
      double bw;
      sscanf(value, "%lf", &bw);
      bw *= datafile->NrFreqChan;
      if(bandwidth != 0.0) {
 if(fabs(bw/bandwidth -1.0) > 1e-5) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Bandwidth as reported by subint table (%lf MHz) is different from that in the history table (%lf MHz). Updating bandwidth to %lf MHz. Something might be wrong with the format data.", datafile->filename, bw, bandwidth);
 }
      }
      if(set_bandwidth(datafile, bw, verbose) == 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Changing bandwidth failed.", datafile->filename);
 return 0;
      }
    }
    if(fits_read_card(datafile->fits_fptr, "DM", card, &status)) {
      if(nowarnings == 0) {
 fflush(stdout);
 if((version_major > 4) || (version_major == 4 && version_minor >= 0)) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): DM not in SUBINT HDU, but this was expected to be set for the claimed PSRFITS version of this file.", datafile->filename);
 }else if(version_major == 0 && version_minor == 0) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): DM not in SUBINT HDU, which could be true for older versions of PSRFITS data type.", datafile->filename);
 }
      }
      status = 0;
      dmnotset = 1;
    }else {
      fits_parse_value(card, value, comment, &status);
      sscanf(value, "%lf", &(datafile->dm));
    }
    if(fits_read_card(datafile->fits_fptr, "RM", card, &status)) {
      if(nowarnings == 0) {
 fflush(stdout);
 if((version_major > 4) || (version_major == 4 && version_minor >= 0)) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): RM not in SUBINT HDU, but this was expected to be set for the claimed PSRFITS version of this file.", datafile->filename);
 }else if(version_major == 0 && version_minor == 0) {
   printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): RM not in SUBINT HDU, which could be true for older versions of PSRFITS data type.", datafile->filename);
 }
      }
      status = 0;
      rmnotset = 1;
    }else {
      fits_parse_value(card, value, comment, &status);
      sscanf(value, "%lf", &(datafile->rm));
    }
    if (status == END_OF_FILE) status = 0;
    if(verbose.debug) {
      fflush(stdout);
      printf("  readPSRFITSHeader (%s): Trying to use DAT_FREQ column in subint table to update centre frequency.\n", datafile->filename);
    }
    update_freqs_using_DAT_FREQ_column(datafile, verbose);
    if (status == END_OF_FILE) status = 0;
    if(issearch) {
      datafile->tsubMode = TSUBMODE_FIXEDTSUB;
      datafile->tsub_list[0] = datafile->NrBins * get_tsamp(*datafile, 0, verbose);
    }else {
      if(fits_get_colnum(datafile->fits_fptr, CASEINSEN, "TSUBINT", &colnum, &status)) {
 status = 0;
 fflush(stdout);
 if(gentypenotset) {
   if(nowarnings == 0) {
     printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find TSUBINT in subint table. Gentype cannot be determined.", datafile->filename);
   }
 }
 if(tsubnotset) {
   if(nowarnings == 0) {
     printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot find TSUBINT in subint table. Observation duration cannot be determined.", datafile->filename);
   }
 }
      }else {
 if(verbose.debug) fprintf(stderr, "  readPSRFITSHeader (%s): start determining total duration of observation\n", datafile->filename);
 double tot_duration, subint_duration;
 int ok;
 ok = 1;
 if(tsubnotset) {
   datafile->tsubMode = TSUBMODE_TSUBLIST;
   if(datafile->NrSubints == 1) {
     datafile->tsubMode = TSUBMODE_FIXEDTSUB;
   }
   if(datafile->tsub_list != NULL) {
     free(datafile->tsub_list);
   }
   if(verbose.debug)
     printf("  readPSRFITSHeader (%s): Re-initialise tsub_list as a list of length %ld\n", datafile->filename, datafile->NrSubints);
   datafile->tsub_list = (double *)malloc(datafile->NrSubints*sizeof(double));
   if(datafile->tsub_list == NULL) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPSRFITSHeader: Memory allocation error");
     return 0;
   }
 }
 tot_duration = subint_duration = 0;
 for(i = 0; i < datafile->NrSubints; i++) {
   if(fits_read_col(datafile->fits_fptr, TDOUBLE, colnum, 1+i, 1, 1, NULL, &(subint_duration), &anynul, &status)) {
     if(tsubnotset) {
       if(nowarnings == 0) {
  fflush(stdout);
  printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot read TSUBINT in subint table to determine observation duration.", datafile->filename);
       }
     }
     ok = 0;
     status = 0;
     if(gentypenotset) {
       if(nowarnings == 0) {
  fflush(stdout);
  printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Cannot read TSUBINT in subint table to determine observation duration. The gentype is not determined.", datafile->filename);
       }
     }
     break;
   }else {
     if(tsubnotset) {
       datafile->tsub_list[i] = subint_duration;
     }
     tot_duration += subint_duration;
   }
 }
 if(ok) {
   double expected_duration, testvalue;
   double period;
   if(get_period(*datafile, 0, &period, verbose) == 2) {
     printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot obtain period", datafile->filename);
     return 0;
   }
   expected_duration = period * datafile->NrSubints;
   testvalue = (tot_duration-expected_duration)/expected_duration;
   if(testvalue < 0.01) {
     if(verbose.debug) {
       fflush(stdout);
       printf("  readPSRFITSHeader (%s): Total duration of observation (%lf sec) suggests data are single pulses.\n", datafile->filename, tot_duration);
     }
     if(gentypenotset) {
       datafile->gentype = GENTYPE_PULSESTACK;
       gentypenotset = 0;
       if(verbose.debug) {
  fflush(stdout);
  printf("  readPSRFITSHeader (%s): gentype guessed to be: %s\n", datafile->filename, returnGenType_str(datafile->gentype));
       }
     }
     tsubnotset = 0;
   }else if(testvalue >= 0.01) {
     if(verbose.debug) {
       fflush(stdout);
       printf("  readPSRFITSHeader (%s): Total duration of observation (%lf sec) suggests multiple pulses are summed.\n", datafile->filename, tot_duration);
     }
     if(gentypenotset) {
       if(datafile->NrSubints == 1) {
  datafile->gentype = GENTYPE_PROFILE;
       }else {
  datafile->gentype = GENTYPE_SUBINTEGRATIONS;
       }
       if(verbose.debug) {
  fflush(stdout);
  printf("  readPSRFITSHeader (%s): gentype guessed to be: %s\n", datafile->filename, returnGenType_str(datafile->gentype));
       }
       gentypenotset = 0;
       tsubnotset = 0;
     }
   }else {
     fflush(stdout);
     double period;
     if(get_period(*datafile, 0, &period, verbose) == 2) {
       printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot obtain period", datafile->filename);
       return 0;
     }
     if(period != 0) {
       if(nowarnings == 0) {
  printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Total duration of observation (%lf sec) implies the nr of periods < nr of subintegrations. Something is wrong.", datafile->filename, tot_duration);
       }
       datafile->tsubMode = TSUBMODE_UNKNOWN;
       if(gentypenotset) {
  fflush(stdout);
  printf("readPSRFITSHeader (%s): Without an observation duration the gentype is not determined.\n", datafile->filename);
       }
     }else {
       if(nowarnings == 0) {
  printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Period does not appear to be set. Observation duration cannot be compared to number of subints.", datafile->filename);
       }
       if(gentypenotset) {
  fflush(stdout);
  printf("readPSRFITSHeader (%s): Without an observation duration the gentype is not determined.\n", datafile->filename);
       }
     }
   }
 }
 if(verbose.debug) fprintf(stderr, "  readPSRFITSHeader (%s): start determining total duration of observation done\n", datafile->filename);
      }
    }
    if (status == END_OF_FILE) status = 0;
    if((dmnotset || rmnotset) && issearch == 0) {
      if(fits_movnam_hdu(datafile->fits_fptr, BINARY_TBL, "PSREPHEM", 0, &status)) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): PSREPHEM table does not exist!", datafile->filename);
 return 0;
      }
      if(dmnotset) {
 if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "DM", &colnum, &status)) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): DM not in PSREPHEM as well.", datafile->filename);
   return 0;
 }
 if(fits_read_col(datafile->fits_fptr, TDOUBLE, colnum, 1, 1, 1, NULL, &dummy_double, &anynul, &status)) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot read DM.", datafile->filename);
   return 0;
 }
 datafile->dm = dummy_double;
      }
      if(rmnotset) {
 if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "RM", &colnum, &status)) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): RM not in PSREPHEM as well.", datafile->filename);
   return 0;
 }
 if(fits_read_col(datafile->fits_fptr, TDOUBLE, colnum, 1, 1, 1, NULL, &dummy_double, &anynul, &status)) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot read RM.", datafile->filename);
   return 0;
 }
 datafile->rm = dummy_double;
      }
    }
    if (status == END_OF_FILE) status = 0;
    if(fits_movnam_hdu(datafile->fits_fptr, BINARY_TBL, "SUBINT", 0, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): SUBINT table does not exist!", datafile->filename);
      return 0;
    }
    if (status) {
      fflush(stdout);
      fits_report_error(stderr, status);
    }
    if(readnoscales)
      return 1;
    else
      return readPSRFITSscales(datafile, verbose);
  }
  if(nodata) {
    if(fits_movnam_hdu(datafile->fits_fptr, BINARY_TBL, "FEEDPAR", 0, &status)) {
      if(verbose.verbose)
 printf("  readPSRFITSHeader (%s): Cannot move to FEEDPAR HDU.\n", datafile->filename);
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): No subint data and no receiver model, I don't know what to do with this file.", datafile->filename);
      return 0;
    }
    if(verbose.verbose) printf("  readPSRFITSHeader (%s): PSRFITS file contains a receiver model.\n", datafile->filename);
    if (fits_read_card(datafile->fits_fptr,"NCPAR", card, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): NCPAR keyword does not exist", datafile->filename);
      return 0;
    }
    if(status) {
      fflush(stdout);
      fits_report_error(stderr, status);
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%ld", &(datafile->NrPols));
    if (fits_read_card(datafile->fits_fptr,"NCHAN", card, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): NCHAN keyword does not exist", datafile->filename);
      return 0;
    }
    if(status) {
      fflush(stdout);
      fits_report_error(stderr, status);
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%ld", &(datafile->NrFreqChan));
    fits_get_num_rows(datafile->fits_fptr, &(datafile->NrSubints), &status);
    datafile->gentype = GENTYPE_RECEIVERMODEL2;
    if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "CHISQ", &colnum, &status)) {
      datafile->gentype = GENTYPE_RECEIVERMODEL;
      status = 0;
    }
    if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "NFREE", &colnum, &status)) {
      datafile->gentype = GENTYPE_RECEIVERMODEL;
      status = 0;
    }
    if(datafile->gentype == GENTYPE_RECEIVERMODEL2) {
      if(verbose.verbose) printf("  readPSRFITSHeader (%s): Found chi^2 and nr free parameters column.\n", datafile->filename);
      datafile->NrPols += 2;
    }else if(verbose.verbose) {
      printf("  readPSRFITSHeader (%s): Didn't find chi^2 and nr free parameters column.\n", datafile->filename);
    }
    datafile->NrBins = 2;
    if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "DATAERR", &colnum, &status)) {
      datafile->NrBins = 1;
      status = 0;
      if(verbose.verbose) printf("  readPSRFITSHeader (%s): Found no errorsbars on values.\n", datafile->filename);
    }else {
      if(verbose.verbose) printf("  readPSRFITSHeader (%s): Found errorsbars on values.\n", datafile->filename);
    }
    if(fits_get_colnum (datafile->fits_fptr, CASEINSEN, "DATA", &colnum, &status)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): No receiver solution data in fits file?", datafile->filename);
      return 0;
    }
    if(verbose.debug) {
      fflush(stdout);
      printf("  readPSRFITSHeader (%s): Trying to use DAT_FREQ column in FEEDPAR table to update centre frequency.\n", datafile->filename);
    }
    update_freqs_using_DAT_FREQ_column(datafile, verbose);
  }
  return 1;
}
int constructFITSsearchsubint(datafile_definition datafile, float *data, int subintnr, unsigned char **subintdata, float **scales, float **offsets, int alreadyscaled, int allocmem, int destroymem, verbose_definition verbose)
{
  long subintsize, f, p, b, samplenr, bitnr, bytenr;
  int bitoff, ivalue, bitofffac[8];
  float offset, scale, fvalue;
  if(datafile.NrBits == 4) {
    bitofffac[0] = 16;
    bitofffac[4] = 1;
  }else if(datafile.NrBits == 8) {
    bitofffac[0] = 1;
  }else if(datafile.NrBits == 2) {
    bitofffac[0] = 64;
    bitofffac[2] = 16;
    bitofffac[4] = 4;
    bitofffac[6] = 1;
  }else if(datafile.NrBits != 8) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR constructFITSsearchsubint: Writing of %d bits data is not supported", datafile.NrBits);
    return 0;
  }
  if(destroymem) {
    free(*subintdata);
    free(*scales);
    free(*offsets);
    return 1;
  }
  subintsize = datafile.NrPols*datafile.NrBins*datafile.NrFreqChan;
  subintsize *= datafile.NrBits;
  subintsize /= 8;
  if(allocmem) {
    *subintdata = (unsigned char *)malloc(subintsize);
    *scales = (float *)malloc(datafile.NrPols*datafile.NrFreqChan*sizeof(float));
    *offsets = (float *)malloc(datafile.NrPols*datafile.NrFreqChan*sizeof(float));
    if(*subintdata == NULL || *scales == NULL || *offsets == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR constructFITSsearchsubint: Cannot allocate %ld bytes", subintsize);
      return 0;
    }
    return 1;
  }
  memset(*subintdata, 0, subintsize);
  for(f = 0; f < datafile.NrFreqChan; f++) {
    for(p = 0; p < datafile.NrPols; p++) {
      if(alreadyscaled == 0) {
 internalFITSscalePulse(&data[datafile.NrBins*(p+datafile.NrPols*(f+subintnr*datafile.NrFreqChan))], datafile.NrBins, &offset, &scale, pow(2, datafile.NrBits)-1);
 (*offsets)[p*datafile.NrFreqChan+f] = offset;
 (*scales)[p*datafile.NrFreqChan+f] = scale;
      }
      for(b = 0; b < datafile.NrBins; b++) {
 if(alreadyscaled == 0)
   fvalue = (data[datafile.NrBins*(p+datafile.NrPols*(f+subintnr*datafile.NrFreqChan))+b]-(*offsets)[p*datafile.NrFreqChan+f])/(*scales)[p*datafile.NrFreqChan+f];
 else
   fvalue = data[datafile.NrBins*(p+datafile.NrPols*(f+subintnr*datafile.NrFreqChan))+b];
 ivalue = fvalue+0.5;
         samplenr = b*datafile.NrPols*datafile.NrFreqChan+p*datafile.NrFreqChan+f;
 bitnr = samplenr * datafile.NrBits;
 bytenr = bitnr/8;
 bitoff = bitnr % 8;
 if(datafile.NrBits != 16) {
   ivalue *= bitofffac[bitoff];
   if(ivalue > 255) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR constructFITSsearchsubint: Packing error (sample %ld: val=%d mult=%d nrbits=%d)", samplenr, ivalue, bitofffac[bitoff], datafile.NrBits);
     return 0;
   }
   (*subintdata)[bytenr] += ivalue;
 }else {
   int ivalue2;
   ivalue2 = ivalue/256;
   (*subintdata)[bytenr] += ivalue2;
   ivalue2 = ivalue % 256;
   (*subintdata)[bytenr+1] += ivalue2;
 }
      }
    }
  }
  return 1;
}
int writeFITSsubint(datafile_definition datafile, long subintnr, unsigned char *subintdata, float *scales, float *offsets, verbose_definition verbose)
{
  int status = 0;
  long subintsize, i;
  float weight;
  double offset;
  static int showErrorMessage = 1;
  if(lookupSubintTable(datafile, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeFITSsubint: Cannot mode to subint table.");
    return 0;
  }
  double period;
  int ret;
  if(datafile.isFolded) {
    ret = get_period(datafile, 0, &period, verbose);
    if(ret == 2) {
      printerror(verbose.debug, "ERROR writeFITSsubint (%s): Cannot obtain period", datafile.filename);
      return 0;
    }
  }else {
    ret = 1;
    period = -1;
  }
  if(datafile.gentype != GENTYPE_SEARCHMODE && datafile.isFolded && ret == 0 && period >= 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeFITSsubint: This function assumes file is search mode data, which it isn't.");
    return 0;
  }
  subintsize = datafile.NrPols*datafile.NrBins*datafile.NrFreqChan;
  subintsize *= datafile.NrBits;
  subintsize /= 8;
  if(fits_write_col(datafile.fits_fptr, TBYTE, 18, 1+subintnr, 1, subintsize, subintdata, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeFITSsubint: Error writing data.");
    fits_report_error(stderr, status);
    return 0;
  }
  for(i = 0; i < datafile.NrPols*datafile.NrFreqChan; i++) {
    if(isnan(offsets[i])) {
      offsets[i] = 0;
      scales[i] = 0;
      if(showErrorMessage) {
 fflush(stdout);
 printwarning(verbose.debug, "writeFITSsubint: Caught an NaN in an offset (duplicate messages are suppressed)");
      }
      showErrorMessage = 0;
    }
    if(isfinite(offsets[i]) == 0) {
      offsets[i] = 0;
      scales[i] = 0;
      if(showErrorMessage) {
 fflush(stdout);
 printwarning(verbose.debug, "writeFITSsubint: Caught an infinity in an offset (duplicate messages are suppressed)");
      }
      showErrorMessage = 0;
    }
    if(isnan(scales[i])) {
      offsets[i] = 0;
      scales[i] = 0;
      if(showErrorMessage) {
 fflush(stdout);
 printwarning(verbose.debug, "writeFITSsubint: Caught an NaN in a scale (duplicate messages are suppressed)");
      }
      showErrorMessage = 0;
    }
    if(isfinite(scales[i]) == 0) {
      offsets[i] = 0;
      scales[i] = 0;
      if(showErrorMessage) {
 fflush(stdout);
 printwarning(verbose.debug, "writeFITSsubint: Caught an infinity in a scale (duplicate messages are suppressed)");
      }
      showErrorMessage = 0;
    }
  }
  if(fits_write_col(datafile.fits_fptr, TFLOAT, 17, 1+subintnr, 1, datafile.NrPols*datafile.NrFreqChan, scales, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeFITSsubint: Error writing scales.");
    fits_report_error(stderr, status);
    return 0;
  }
  if(fits_write_col(datafile.fits_fptr, TFLOAT, 16, 1+subintnr, 1, datafile.NrPols*datafile.NrFreqChan, offsets, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeFITSsubint: Error writing offsets.");
    fits_report_error(stderr, status);
    return 0;
  }
  weight = 1;
  for(i = 0; i < datafile.NrFreqChan; i++) {
    if(fits_write_col(datafile.fits_fptr, TFLOAT, 15, 1+subintnr, 1+i, 1, &weight, &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeFITSsubint: Error writing weights.");
      fits_report_error(stderr, status);
      return 0;
    }
  }
  offset = ((double)((subintnr + 0.5)*datafile.NrBins))*get_tsamp(datafile, 0, verbose);
  if(fits_write_col(datafile.fits_fptr, TDOUBLE, 3, 1+subintnr, 1, 1, &offset, &status) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeFITSsubint: Error writing subint offset.");
    fits_report_error(stderr, status);
    return 0;
  }
  return 1;
}
int writeFITSfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  unsigned char *subintdata;
  float *scales, *offsets;
  long n, f, p;
  double period;
  int ret;
  if(datafile.isFolded) {
    ret = get_period(datafile, 0, &period, verbose);
    if(ret == 2) {
      printerror(verbose.debug, "ERROR writeFITSfile (%s): Cannot obtain period", datafile.filename);
      return 0;
    }
  }else {
    ret = 1;
    period = -1;
  }
  if((datafile.gentype != GENTYPE_SEARCHMODE && datafile.isFolded && ret == 0 && period >= 0) || datafile.gentype == GENTYPE_RECEIVERMODEL || datafile.gentype == GENTYPE_RECEIVERMODEL2) {
    for(n = 0; n < datafile.NrSubints; n++) {
      if(verbose.verbose && verbose.nocounters == 0) printf("writeFITSfile: pulse %ld/%ld          \r", n+1, datafile.NrSubints);
      for(f = 0; f < datafile.NrFreqChan; f++) {
 for(p = 0; p < datafile.NrPols; p++) {
   if(writeFITSpulse(datafile, n, p, f, 0, datafile.NrBins, &data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))], verbose) == 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR writeFITSfile: Cannot write data               ");
     return 0;
   }
 }
      }
    }
  }else {
    n = 0;
    if(constructFITSsearchsubint(datafile, data, n, &subintdata, &scales, &offsets, 0, 1, 0, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeFITSfile: Cannot allocate temporary memory               ");
      return 0;
    }
    for(n = 0; n < datafile.NrSubints; n++) {
      if(verbose.verbose && verbose.nocounters == 0) printf("writeFITSfile: pulse %ld/%ld             \r", n+1, datafile.NrSubints);
      if(constructFITSsearchsubint(datafile, data, n, &subintdata, &scales, &offsets, 0, 0, 0, verbose) != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writeFITSfile: Cannot construct subint data               ");
 return 0;
      }
      if(writeFITSsubint(datafile, n, subintdata, scales, offsets, verbose) != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writeFITSfile: Cannot write subint data               ");
 return 0;
      }
    }
    if(constructFITSsearchsubint(datafile, data, n, &subintdata, &scales, &offsets, 0, 0, 1, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeFITSfile: Cannot free temporary memory               ");
      return 0;
    }
  }
  if(verbose.verbose) printf("Writing is done.                                  \n");
  return 1;
}
int readFITSfile(datafile_definition *datafile, float *data, verbose_definition verbose)
{
  long n, f, p, i;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Reading in entire PSRFITS file\n");
  }
  for(n = 0; n < datafile->NrSubints; n++) {
    if(verbose.verbose && verbose.nocounters == 0) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  subint %ld/%ld (%.1f%%)\r", n+1, datafile->NrSubints, 100.0*(n+1)/(float)datafile->NrSubints);
      fflush(stdout);
    }
    for(f = 0; f < datafile->NrFreqChan; f++) {
      for(p = 0; p < datafile->NrPols; p++) {
 if(readFITSpulse(datafile, n, p, f, 0, datafile->NrBins, &data[datafile->NrBins*(p+datafile->NrPols*(f+n*datafile->NrFreqChan))], verbose) == 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readFITSfile: Cannot read data  (n=%ld f=%ld p=%ld)             ", n, f, p);
     return 0;
   }
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  Reading is done.                      \n");
  }
  return 1;
}
int readPSRCHIVE_ASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  char txt[100];
  int ret;
  long n, f, i, p;
  float sample;
  for(n = 0; n < datafile.NrSubints; n++) {
    if(verbose.verbose && verbose.nocounters == 0) printf("readPSRCHIVE_ASCIIfile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
    for(f = 0; f < datafile.NrFreqChan; f++) {
      for(i = 0; i < datafile.NrBins; i++) {
 ret = fscanf(datafile.fptr, "%s", txt);
 if(ret != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIfile: Read error. File not in right format (not generated with \"pdv -t\"?) (%s != 0).", txt);
   return 0;
 }
 if(i == 0 && f == 0 && n == 0) {
   if(strcmp(txt, "0") != 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIfile: File not in right format (not generated with \"pdv -t\"?) (%s != 0).", txt);
     return 0;
   }
 }
 ret = fscanf(datafile.fptr, "%s", txt);
 if(ret != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIfile: Read error. File not in right format (not generated with \"pdv -t\"?) (%s != 0).", txt);
   return 0;
 }
 if(i == 0 && f == 0 && n == 0) {
   if(strcmp(txt, "0") != 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIfile: File not in right format (not generated with \"pdv -t\"?) (%s != 0).",txt );
     return 0;
   }
 }
 ret = fscanf(datafile.fptr, "%s", txt);
 if(ret != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIfile: Read error. File not in right format (not generated with \"pdv -t\"?) (%s != 0).", txt);
   return 0;
 }
 if(i == 0 && f == 0 && n == 0) {
   if(strcmp(txt, "0") != 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIfile: File not in right format (not generated with \"pdv -t\"?) (%s != 0).", txt);
     return 0;
   }
 }
 for(p = 0; p < datafile.NrPols; p++) {
   ret = fscanf(datafile.fptr, "%f", &sample);
   if(ret != 1) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIfile: Read error. File not in right format (not generated with \"pdv -t\"?) (%s != 0).", txt);
     return 0;
   }
   data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] = sample;
 }
      }
    }
  }
  if(verbose.verbose) printf("Reading is done.                           \n");
  return 1;
}
int writePSRCHIVE_ASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long n, f, i, p;
  float sample;
  for(n = 0; n < datafile.NrSubints; n++) {
    if(verbose.verbose && verbose.nocounters == 0) printf("writePSRCHIVE_ASCIIfile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
    for(f = 0; f < datafile.NrFreqChan; f++) {
      for(i = 0; i < datafile.NrBins; i++) {
 fprintf(datafile.fptr, "%ld %ld %ld ", n, f, i);
 for(p = 0; p < datafile.NrPols; p++) {
   sample = data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i];
   fprintf(datafile.fptr, "%e ", sample);
 }
 fprintf(datafile.fptr, "\n");
      }
    }
  }
  if(verbose.verbose) printf("Writing is done.                           \n");
  return 1;
}
int readPSRCHIVE_ASCIIfilepulse(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  char txt[100];
  long n, f, i, p, pn, fn, bn;
  int ret;
  float sample;
  if(pulsenr != 0 || freq != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIfilepulse: I'm not sure which column is freq and pulse, so I can only read in first pulse and first frequency channel at this point.");
    return 0;
  }
  rewind(datafile.fptr);
  readPSRCHIVE_ASCIIHeader(&datafile, verbose);
  for(n = 0; n < datafile.NrSubints; n++) {
    for(f = 0; f < datafile.NrFreqChan; f++) {
      for(i = 0; i < datafile.NrBins; i++) {
 fscanf(datafile.fptr, "%s", txt);
 if(i == 0 && f == 0 && n == 0) {
   if(strcmp(txt, "0") != 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIfilepulse: File not in right format (not generated with \"pdv -t\"?) (%s).", txt);
     return 0;
   }
 }
 ret = sscanf(txt, "%ld", &pn);
 if(ret != 1)
   pn = -1;
 fscanf(datafile.fptr, "%s", txt);
 if(i == 0 && f == 0 && n == 0) {
   if(strcmp(txt, "0") != 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIfilepulse: File not in right format (not generated with \"pdv -t\"?) (%s).", txt);
     return 0;
   }
 }
 ret = sscanf(txt, "%ld", &fn);
 if(ret != 1)
   fn = -1;
 fscanf(datafile.fptr, "%s", txt);
 if(i == 0 && f == 0 && n == 0) {
   if(strcmp(txt, "0") != 0) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPSRCHIVE_ASCIIfilepulse: File not in right format (not generated with \"pdv -t\"?) (%s).", txt);
     return 0;
   }
 }
 ret = sscanf(txt, "%ld", &bn);
 if(ret != 1)
   bn = -1;
 for(p = 0; p < datafile.NrPols; p++) {
   fscanf(datafile.fptr, "%f", &sample);
   if(pn == pulsenr && fn == freq && (bn >= 0 || bn < datafile.NrBins) && p == polarization)
   pulse[bn] = sample;
 }
      }
    }
  }
  return 1;
}
int writeHistoryFITS(datafile_definition datafile, verbose_definition verbose)
{
  int status = 0;
  int colnum_date, colnum_cmd, colnum_user, colnum_hostname;
  datafile_history_entry_definition *curhistoryEntry;
  long nrows;
  int nstart, n;
  char questionmark[2], *str_ptr, *cmd_ptr, txt2[2000];
  questionmark[0] = '?';
  questionmark[1] = 0;
  if(fits_movnam_hdu(datafile.fits_fptr, BINARY_TBL, "HISTORY_NOT_PSRFITS", 0, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeHistoryFITS: Cannot move to history HDU.");
    return 0;
  }
  if(fits_get_colnum (datafile.fits_fptr, CASEINSEN, "DATE_PRO", &colnum_date, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeHistoryFITS: No DATE_PRO column is history table?");
    return 0;
  }
  if(fits_get_colnum (datafile.fits_fptr, CASEINSEN, "PROC_CMD", &colnum_cmd, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeHistoryFITS: No PROC_CMD column is history table?");
    return 0;
  }
  if(fits_get_colnum (datafile.fits_fptr, CASEINSEN, "USER", &colnum_user, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeHistoryFITS: No HOSTNAME column is history table?");
    return 0;
  }
  if(fits_get_colnum (datafile.fits_fptr, CASEINSEN, "HOSTNAME", &colnum_hostname, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeHistoryFITS: No HOSTNAME column is history table?");
    return 0;
  }
  curhistoryEntry = &(datafile.history);
  do {
    if(curhistoryEntry->cmd != NULL)
      cmd_ptr = curhistoryEntry->cmd;
    else
      cmd_ptr = questionmark;
    nstart = 0;
    long totallength = strlen(cmd_ptr);
    do {
      fits_get_num_rows(datafile.fits_fptr, &nrows, &status);
      if(fits_insert_rows(datafile.fits_fptr, nrows, 1, &status)) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writeHistoryFITS: Cannot add a row to history table");
 fits_report_error(stderr, status);
 return 0;
      }
      fits_get_num_rows(datafile.fits_fptr, &nrows, &status);
      if(nstart == 0) {
 if(curhistoryEntry->timestamp != NULL)
   str_ptr = curhistoryEntry->timestamp;
 else
   str_ptr = questionmark;
 if(fits_write_col(datafile.fits_fptr, TSTRING, colnum_date, nrows, 1, 1, &str_ptr, &status) != 0) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writeHistoryFITS: Error writing date to history table.");
   fits_report_error(stderr, status);
   return 0;
 }
 if(curhistoryEntry->user != NULL)
   str_ptr = curhistoryEntry->user;
 else
   str_ptr = questionmark;
 if(fits_write_col(datafile.fits_fptr, TSTRING, colnum_user, nrows, 1, 1, &str_ptr, &status) != 0) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writeHistoryFITS: Error writing user to history table.");
   fits_report_error(stderr, status);
   return 0;
 }
 if(curhistoryEntry->hostname != NULL)
   str_ptr = curhistoryEntry->hostname;
 else
   str_ptr = questionmark;
 if(fits_write_col(datafile.fits_fptr, TSTRING, colnum_hostname, nrows, 1, 1, &str_ptr, &status) != 0) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR writeHistoryFITS: Error writing hostname to history table.");
   fits_report_error(stderr, status);
   return 0;
 }
      }
      txt2[1024] = 0;
      strncpy(txt2, cmd_ptr, 1024);
      n = strlen(txt2);
      cmd_ptr += n;
      nstart += n;
      str_ptr = txt2;
      if(fits_write_col(datafile.fits_fptr, TSTRING, colnum_cmd, nrows, 1, 1, &str_ptr, &status) != 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR writeHistoryFITS: Error writing data to history table.");
 fits_report_error(stderr, status);
 return 0;
      }
    }while(nstart < totallength);
    curhistoryEntry = curhistoryEntry->nextEntry;
  }while(curhistoryEntry != NULL);
  return 1;
}
int readHistoryFITS_psrsalsa(datafile_definition *datafile, int surpressHDUerror, verbose_definition verbose)
{
  int status = 0;
  int colnum_date, colnum_cmd, colnum_user, colnum_hostname, anynul;
  long nrows;
  char txt2[2000];
  int rownr;
  char *char_ptrptr[1];
  datafile_history_entry_definition *curHistoryEntry;
  if(datafile->opened_flag == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryFITS: File already closed?");
    return 1;
  }
  if(fits_movnam_hdu(datafile->fits_fptr, BINARY_TBL, "HISTORY_NOT_PSRFITS", 0, &status)) {
    if(surpressHDUerror == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS: Cannot move to HISTORY_NOT_PSRFITS HDU, file not generated by psrsalsa code?");
    }
    return 2;
  }
  if(fits_get_colnum(datafile->fits_fptr, CASEINSEN, "DATE_PRO", &colnum_date, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryFITS: No DATE_PRO column is history table?");
    return 3;
  }
  if(fits_get_colnum(datafile->fits_fptr, CASEINSEN, "PROC_CMD", &colnum_cmd, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryFITS: No PROC_CMD column is history table?");
    return 3;
  }
  if(fits_get_colnum(datafile->fits_fptr, CASEINSEN, "USER", &colnum_user, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryFITS: No HOSTNAME column is history table?");
    return 3;
  }
  if(fits_get_colnum(datafile->fits_fptr, CASEINSEN, "HOSTNAME", &colnum_hostname, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryFITS: No HOSTNAME column is history table?");
    return 3;
  }
  fits_get_num_rows(datafile->fits_fptr, &nrows, &status);
  if(nrows <= 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryFITS: No rows present in history table?");
    return 4;
  }
  curHistoryEntry = &(datafile->history);
  for(rownr = 0; rownr < nrows; rownr++) {
    if(curHistoryEntry->timestamp != NULL || curHistoryEntry->cmd != NULL || curHistoryEntry->user != NULL || curHistoryEntry->hostname != NULL || curHistoryEntry->nextEntry != NULL) {
      curHistoryEntry->nextEntry = malloc(sizeof(datafile_history_entry_definition));
      if(curHistoryEntry->nextEntry == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readHistoryFITS: Memory allocation error");
 return 6;
      }
      curHistoryEntry = curHistoryEntry->nextEntry;
      curHistoryEntry->timestamp = NULL;
      curHistoryEntry->cmd = NULL;
      curHistoryEntry->user = NULL;
      curHistoryEntry->hostname = NULL;
      curHistoryEntry->nextEntry = NULL;
    }
    char_ptrptr[0] = txt2;
    if(fits_read_col(datafile->fits_fptr, TSTRING, colnum_date, rownr+1, 1, 1, NULL, char_ptrptr, &anynul, &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS: Error reading date from history table.");
      fits_report_error(stderr, status);
      return 5;
    }
    curHistoryEntry->timestamp = malloc(strlen(txt2)+1);
    if(curHistoryEntry->timestamp == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS: Memory allocation error");
      return 6;
    }
    strcpy(curHistoryEntry->timestamp, txt2);
    if(fits_read_col(datafile->fits_fptr, TSTRING, colnum_user, rownr+1, 1, 1, NULL, char_ptrptr, &anynul, &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS: Error reading user from history table.");
      fits_report_error(stderr, status);
      return 5;
    }
    curHistoryEntry->user = malloc(strlen(txt2)+1);
    if(curHistoryEntry->user == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS: Memory allocation error");
      return 6;
    }
    strcpy(curHistoryEntry->user, txt2);
    if(fits_read_col(datafile->fits_fptr, TSTRING, colnum_hostname, rownr+1, 1, 1, NULL, char_ptrptr, &anynul, &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS: Error reading hostname from history table.");
      fits_report_error(stderr, status);
      return 5;
    }
    curHistoryEntry->hostname = malloc(strlen(txt2)+1);
    if(curHistoryEntry->hostname == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS: Memory allocation error");
      return 6;
    }
    strcpy(curHistoryEntry->hostname, txt2);
    if(fits_read_col(datafile->fits_fptr, TSTRING, colnum_cmd, rownr+1, 1, 1, NULL, char_ptrptr, &anynul, &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS: Error reading data from history table.");
      fits_report_error(stderr, status);
      return 5;
    }
    curHistoryEntry->cmd = malloc(strlen(txt2)+1);
    if(curHistoryEntry->cmd == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS: Memory allocation error");
      return 6;
    }
    strcpy(curHistoryEntry->cmd, txt2);
  }
  return 0;
}
int readHistoryFITS_psrfits(datafile_definition *datafile, int surpressHDUerror, int surpressColumnNotPresent, verbose_definition verbose)
{
  int status = 0;
  int colnum_date, colnum_cmd, anynul;
  long nrows;
  char txt2[2000];
  int rownr;
  char *char_ptrptr[1];
  datafile_history_entry_definition *curHistoryEntry;
  if(datafile->opened_flag == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryFITS_psrfits: File already closed?");
    return 1;
  }
  if(fits_movnam_hdu(datafile->fits_fptr, BINARY_TBL, "HISTORY", 0, &status)) {
    if(surpressHDUerror == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readHistoryFITS_psrfits: Cannot move to HISTORY HDU, no history stored in file?");
    }
    return 2;
  }
  if(fits_get_colnum(datafile->fits_fptr, CASEINSEN, "DATE_PRO", &colnum_date, &status)) {
    if(surpressColumnNotPresent == 0 || verbose.debug) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readHistoryFITS_psrfits: No DATE_PRO column is history table?");
    }
    return 3;
  }
  if(fits_get_colnum(datafile->fits_fptr, CASEINSEN, "PROC_CMD", &colnum_cmd, &status)) {
    if(surpressColumnNotPresent == 0 || verbose.debug) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readHistoryFITS_psrfits: No PROC_CMD column is history table?");
    }
    return 3;
  }
  fits_get_num_rows(datafile->fits_fptr, &nrows, &status);
  if(nrows <= 0) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING readHistoryFITS_psrfits: No rows present in history table?");
    return 4;
  }
  curHistoryEntry = &(datafile->history);
  for(rownr = 0; rownr < nrows; rownr++) {
    if(curHistoryEntry->timestamp != NULL || curHistoryEntry->cmd != NULL || curHistoryEntry->user != NULL || curHistoryEntry->hostname != NULL || curHistoryEntry->nextEntry != NULL) {
      curHistoryEntry->nextEntry = malloc(sizeof(datafile_history_entry_definition));
      if(curHistoryEntry->nextEntry == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR readHistoryFITS: Memory allocation error");
 return 6;
      }
      curHistoryEntry = curHistoryEntry->nextEntry;
      curHistoryEntry->timestamp = NULL;
      curHistoryEntry->cmd = NULL;
      curHistoryEntry->user = NULL;
      curHistoryEntry->hostname = NULL;
      curHistoryEntry->nextEntry = NULL;
    }
    char_ptrptr[0] = txt2;
    if(fits_read_col(datafile->fits_fptr, TSTRING, colnum_date, rownr+1, 1, 1, NULL, char_ptrptr, &anynul, &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS_psrfits: Error reading date from history table.");
      fits_report_error(stderr, status);
      return 5;
    }
    curHistoryEntry->timestamp = malloc(strlen(txt2)+1);
    if(curHistoryEntry->timestamp == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS: Memory allocation error");
      return 6;
    }
    strcpy(curHistoryEntry->timestamp, txt2);
    if(fits_read_col(datafile->fits_fptr, TSTRING, colnum_cmd, rownr+1, 1, 1, NULL, char_ptrptr, &anynul, &status) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS_psrfits: Error reading data from history table.");
      fits_report_error(stderr, status);
      return 5;
    }
    curHistoryEntry->cmd = malloc(strlen(txt2)+1);
    if(curHistoryEntry->cmd == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readHistoryFITS: Memory allocation error");
      return 6;
    }
    strcpy(curHistoryEntry->cmd, txt2);
  }
  return 0;
}
int readHistoryFITS(datafile_definition *datafile, verbose_definition verbose)
{
  int ret, ret2, hdunum, hdutype;
  int status = 0;
  if(fits_get_hdu_num(datafile->fits_fptr, &hdunum)) {
  }
  ret = readHistoryFITS_psrfits(datafile, 1, 1, verbose);
  ret2 = readHistoryFITS_psrsalsa(datafile, 1, verbose);
  if(fits_movabs_hdu(datafile->fits_fptr, hdunum, &hdutype, &status)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryFITS (%s): Cannot goto original HDU", datafile->filename);
    return 0;
  }
  if(ret == 0 && ret2 == 0)
    return 1;
  if(ret == 0 && ret2 == 2)
    return 1;
  if(ret == 2 && ret2 == 0)
    return 1;
  if(ret == 2 && ret2 == 2) {
    fflush(stdout);
    printwarning(verbose.debug, "ERROR readHistoryFITS: File doesn't contain a history HDU, so reading history failed.");
    return 0;
  }
  if(ret == 3 && ret2 == 0)
    return 1;
  return 0;
}

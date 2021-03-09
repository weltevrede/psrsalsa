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
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sort_float.h>
#include <gsl/gsl_statistics_float.h>
#include <gsl/gsl_integration.h>
#include "psrsalsa.h"
int filterPApoints(datafile_definition *datafile, verbose_definition verbose)
{
  int dPa_polnr;
  long i, j, nrpoints;
  float *olddata;
  if(datafile->poltype != POLTYPE_ILVPAdPA && datafile->poltype != POLTYPE_PAdPA && datafile->poltype != POLTYPE_ILVPAdPATEldEl) {
    printerror(verbose.debug, "ERROR filterPApoints: Data doesn't appear to have poltype ILVPAdPA, PAdPA or ILVPAdPATEldEl.");
    return -1;
  }
  if(datafile->poltype == POLTYPE_ILVPAdPA && datafile->NrPols != 5) {
    printerror(verbose.debug, "ERROR filterPApoints: 5 polarization channels were expected, but there are only %ld.", datafile->NrPols);
    return -1;
  }else if(datafile->poltype == POLTYPE_ILVPAdPATEldEl && datafile->NrPols != 8) {
    printerror(verbose.debug, "ERROR filterPApoints: 8 polarization channels were expected, but there are only %ld.", datafile->NrPols);
    return -1;
  }else if(datafile->poltype == POLTYPE_PAdPA && datafile->NrPols != 2) {
    printerror(verbose.debug, "ERROR filterPApoints: 2 polarization channels were expected, but there are only %ld.", datafile->NrPols);
    return -1;
  }
  if(datafile->NrSubints > 1 || datafile->NrFreqChan > 1) {
    printerror(verbose.debug, "ERROR filterPApoints: Can only do this opperation if there is one subint and one frequency channel.");
    return -1;
  }
  if(datafile->tsampMode != TSAMPMODE_LONGITUDELIST) {
    printerror(verbose.debug, "ERROR filterPApoints: Expected pulse longitudes to be defined.");
    return -1;
  }
  if(datafile->poltype == POLTYPE_ILVPAdPA || datafile->poltype == POLTYPE_ILVPAdPATEldEl) {
    dPa_polnr = 4;
  }else if(datafile->poltype == POLTYPE_PAdPA) {
    dPa_polnr = 1;
  }
  nrpoints = 0;
  for(i = 0; i < datafile->NrBins; i++) {
    if(datafile->data[i+dPa_polnr*datafile->NrBins] > 0 && isfinite(datafile->data[i+dPa_polnr*datafile->NrBins])) {
      nrpoints++;
    }
  }
  if(verbose.verbose)
    printf("Keeping %ld significant PA points\n", nrpoints);
  olddata = datafile->data;
  datafile->data = (float *)malloc(nrpoints*datafile->NrPols*sizeof(float));
  if(datafile->data == NULL) {
    printerror(verbose.debug, "ERROR filterPApoints: Memory allocation error.");
    return -1;
  }
  j = 0;
  for(i = 0; i < datafile->NrBins; i++) {
    if(olddata[i+dPa_polnr*datafile->NrBins] > 0 && isfinite(olddata[i+dPa_polnr*datafile->NrBins])) {
      if(datafile->poltype == POLTYPE_ILVPAdPA) {
 datafile->data[j+0*nrpoints] = olddata[i+0*datafile->NrBins];
 datafile->data[j+1*nrpoints] = olddata[i+1*datafile->NrBins];
 datafile->data[j+2*nrpoints] = olddata[i+2*datafile->NrBins];
 datafile->data[j+3*nrpoints] = olddata[i+3*datafile->NrBins];
 datafile->data[j+4*nrpoints] = olddata[i+4*datafile->NrBins];
      }else if(datafile->poltype == POLTYPE_ILVPAdPATEldEl) {
 datafile->data[j+0*nrpoints] = olddata[i+0*datafile->NrBins];
 datafile->data[j+1*nrpoints] = olddata[i+1*datafile->NrBins];
 datafile->data[j+2*nrpoints] = olddata[i+2*datafile->NrBins];
 datafile->data[j+3*nrpoints] = olddata[i+3*datafile->NrBins];
 datafile->data[j+4*nrpoints] = olddata[i+4*datafile->NrBins];
 datafile->data[j+5*nrpoints] = olddata[i+5*datafile->NrBins];
 datafile->data[j+6*nrpoints] = olddata[i+6*datafile->NrBins];
 datafile->data[j+7*nrpoints] = olddata[i+7*datafile->NrBins];
      }else {
 datafile->data[j+0*nrpoints] = olddata[i+0*datafile->NrBins];
 datafile->data[j+1*nrpoints] = olddata[i+1*datafile->NrBins];
      }
      datafile->tsamp_list[j] = datafile->tsamp_list[i];
      j++;
    }
  }
  free(olddata);
  datafile->NrBins = nrpoints;
  return datafile->NrBins;
}
int make_paswing_fromIQUV_remove_lowS2N_points_sp_isLsignificant(float sigma_limit, int sigmaI, float dataL, float rmsI, float rmsL, verbose_definition verbose)
{
  if(sigmaI == 0) {
    if(dataL < sigma_limit*rmsL) {
      return 0;
    }
    return 1;
  }else {
    if(dataL < sigma_limit*rmsI) {
      return 0;
    }
    return 1;
  }
}
int make_paswing_fromIQUV_remove_lowS2N_points_sp_isPsignificant(float sigma_limit, int sigmaI, float dataP, float rmsI, float rmsP, verbose_definition verbose)
{
  if(sigmaI == 0) {
    if(dataP < sigma_limit*rmsP) {
      return 0;
    }
    return 1;
  }else {
    if(dataP < sigma_limit*rmsI) {
      return 0;
    }
    return 1;
  }
}
int make_paswing_fromIQUV_remove_lowS2N_points_sp_isVsignificant(float sigma_limit, int sigmaI, float dataV, float rmsI, float rmsV, verbose_definition verbose)
{
  if(sigmaI == 0) {
    if(fabs(dataV) < sigma_limit*rmsV) {
      return 0;
    }
    return 1;
  }else {
    if(fabs(dataV) < sigma_limit*rmsI) {
      return 0;
    }
    return 1;
  }
}
void make_paswing_fromIQUV_remove_lowS2N_points_sp(float sigma_limit, int sigmaI, int nrBins, float *dataL, float *dataP, float *dataPA, float *dataPaErr, float *dataEll, float *dataEllErr, float rmsI, float rmsL, float rmsP, pulselongitude_regions_definition *onpulse, verbose_definition verbose)
{
  long j;
  if(sigma_limit < 0) {
    return;
  }
  for(j = 0; j < nrBins; j++) {
    int issignificant;
    issignificant = 1;
    if(onpulse != NULL) {
      if(checkRegions(j, onpulse, 0, verbose) == 0) {
 issignificant = 0;
      }
    }
    if(dataL != NULL && (dataPA != NULL || dataPaErr != NULL)) {
      if(issignificant == 0 || make_paswing_fromIQUV_remove_lowS2N_points_sp_isLsignificant(sigma_limit, sigmaI, dataL[j], rmsI, rmsL, verbose) == 0) {
 if(dataPA != NULL) {
   dataPA[j] = 0;
 }
 if(dataPaErr != NULL) {
   dataPaErr[j] = -1;
 }
      }
    }
    if(dataP != NULL && (dataEll != NULL || dataEllErr != NULL)) {
      if(issignificant == 0 || make_paswing_fromIQUV_remove_lowS2N_points_sp_isPsignificant(sigma_limit, sigmaI, dataP[j], rmsI, rmsP, verbose) == 0) {
 if(dataEll != NULL) {
   dataEll[j] = 0;
 }
 if(dataEllErr != NULL) {
   dataEllErr[j] = -1;
 }
      }
    }
  }
}
int make_paswing_fromIQUV_sp(float *dataI, float *dataQ, float *dataU, float *dataV, int nrBins, float *dataL, float *dataP, float *dataPA, float *dataPaErr, float *dataEll, float *dataEllErr, float *baseline_intensity, float *rmsI, float *rmsQ, float *rmsU, float *rmsV, float *rmsL, float *rmsP, float *medianL, float *medianP, pulselongitude_regions_definition onpulse, int normalize, int correctLbias, int correctPbias, float correctQV, float correctV, float paoffset, int rms_file_nrBins, float *rms_file_I, float *rms_file_Q, float *rms_file_U, float *rms_file_V, float rebin_factor, verbose_definition verbose)
{
  int rms_file_specified;
  float *Loffpulse, *Poffpulse, median_L, median_P;
  double ymax, avrgI, RMSI, RMSQ, RMSU, RMSV, RMSL, RMSP;
  long i, nrOffpulseBins;
  rms_file_specified = 1;
  if(rms_file_I == NULL) {
    rms_file_I = dataI;
    rms_file_Q = dataQ;
    rms_file_U = dataU;
    rms_file_V = dataV;
    rms_file_nrBins = nrBins;
    rms_file_specified = 0;
  }
  Loffpulse = (float *)malloc(rms_file_nrBins*sizeof(float));
  Poffpulse = (float *)malloc(rms_file_nrBins*sizeof(float));
  if(Loffpulse == NULL || Poffpulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV_sp: Memory allocation error.");
    return 0;
  }
  if(normalize == 0) {
    ymax = 1;
  }else {
    ymax = dataI[0];
    for(i = 1; i < nrBins; i++) {
      if(dataI[i] > ymax) {
 ymax = dataI[i];
      }
    }
    if(ymax == 0) {
      ymax = 1;
    }
  }
  if(ymax != 1.0 || correctQV != 1.0 || correctV != 1.0) {
    for(i = 0; i < nrBins; i++) {
      dataI[i] /= ymax;
      dataQ[i] /= correctQV*ymax;
      dataU[i] /= ymax;
      dataV[i] /= correctV*correctQV*ymax;
    }
    if(rms_file_specified) {
      for(i = 0; i < rms_file_nrBins; i++) {
 rms_file_I[i] /= ymax;
 rms_file_Q[i] /= correctQV*ymax;
 rms_file_U[i] /= ymax;
 rms_file_V[i] /= correctV*correctQV*ymax;
      }
    }
  }
  nrOffpulseBins = 0;
  avrgI = 0;
  RMSI = 0;
  RMSQ = 0;
  RMSU = 0;
  RMSV = 0;
  RMSL = 0;
  RMSP = 0;
  for(i = 0; i < rms_file_nrBins; i++) {
    if(checkRegions(i, &onpulse, 0, verbose) == 0) {
      double sampleI2, sampleQ2, sampleU2, sampleV2;
      avrgI += rms_file_I[i];
      sampleI2 = rms_file_I[i]*rms_file_I[i];
      sampleQ2 = rms_file_Q[i]*rms_file_Q[i];
      sampleU2 = rms_file_U[i]*rms_file_U[i];
      sampleV2 = rms_file_V[i]*rms_file_V[i];
      RMSI += sampleI2;
      RMSQ += sampleQ2;
      RMSU += sampleU2;
      RMSV += sampleV2;
      RMSL += sampleQ2+sampleU2;
      RMSP += sampleQ2+sampleU2+sampleV2;
      Loffpulse[nrOffpulseBins] = sqrt(sampleQ2+sampleU2);
      Poffpulse[nrOffpulseBins] = sqrt(sampleQ2+sampleU2+sampleV2);
      nrOffpulseBins++;
    }
  }
  avrgI /= (double)nrOffpulseBins;
  RMSI = sqrt(RMSI/(double)nrOffpulseBins);
  RMSQ = sqrt(RMSQ/(double)nrOffpulseBins);
  RMSU = sqrt(RMSU/(double)nrOffpulseBins);
  RMSV = sqrt(RMSV/(double)nrOffpulseBins);
  RMSL = sqrt(RMSL/(double)nrOffpulseBins);
  RMSP = sqrt(RMSP/(double)nrOffpulseBins);
  if(rms_file_specified) {
    double scale = 1.0/sqrt(rebin_factor);
    RMSI *= scale;
    RMSQ *= scale;
    RMSU *= scale;
    RMSV *= scale;
    RMSL *= scale;
    RMSP *= scale;
  }
  gsl_sort_float(Loffpulse, 1, nrOffpulseBins);
  median_L = gsl_stats_float_median_from_sorted_data(Loffpulse, 1, nrOffpulseBins);
  gsl_sort_float(Poffpulse, 1, nrOffpulseBins);
  median_P = gsl_stats_float_median_from_sorted_data(Poffpulse, 1, nrOffpulseBins);
  if(baseline_intensity != NULL)
    *baseline_intensity = avrgI;
  if(rmsI != NULL)
    *rmsI = RMSI;
  if(rmsQ != NULL)
    *rmsQ = RMSQ;
  if(rmsU != NULL)
    *rmsU = RMSU;
  if(rmsV != NULL)
    *rmsV = RMSV;
  if(rmsL != NULL)
    *rmsL = RMSL;
  if(rmsP != NULL)
    *rmsP = RMSP;
  if(medianL != NULL)
    *medianL = median_L;
  if(medianP != NULL)
    *medianP = median_P;
  for(i = 0; i < nrBins; i++) {
    double sampleQ2, sampleU2, sampleV2, sampleL, sampleP;
    sampleQ2 = dataQ[i]*dataQ[i];
    sampleU2 = dataU[i]*dataU[i];
    sampleV2 = dataV[i]*dataV[i];
    sampleL = sqrt(sampleQ2+sampleU2);
    if(correctLbias == 0) {
      sampleL -= median_L;
    }else if(correctLbias == 1) {
      double junk = (sqrt(0.5*(RMSQ*RMSQ+RMSU*RMSU))/sampleL);
      if(junk < 1.0) {
 sampleL *= sqrt(1.0-junk*junk);
      }else {
 sampleL = 0.0;
      }
    }
    if(dataL != NULL) {
      dataL[i] = sampleL;
    }
    sampleP = sqrt(sampleQ2+sampleU2+sampleV2);
    if(correctPbias == 0) {
      sampleP -= median_P;
    }
    if(dataP != NULL) {
      dataP[i] = sampleP;
    }
    if(dataPA != NULL) {
      dataPA[i] = 90.0*atan2(dataU[i], dataQ[i])/M_PI;
      if(paoffset != 0.0) {
 dataPA[i] += paoffset;
 dataPA[i] = derotate_180_small_double(dataPA[i]);
      }
    }
    if(dataPaErr != NULL) {
      if(sampleQ2 == 0 && sampleU2 == 0) {
 dataPaErr[i] = 0;
      }else {
 dataPaErr[i] = sqrt(sampleQ2*RMSU*RMSU + sampleU2*RMSQ*RMSQ);
 dataPaErr[i] /= 2.0*(sampleQ2 + sampleU2);
 dataPaErr[i] *= 180.0/M_PI;
      }
    }
    if(dataEll != NULL) {
      if(sampleQ2 == 0 && sampleU2 == 0 && sampleV2 == 0) {
 dataEll[i] = 0;
      }
      double value;
      value = dataV[i]/sqrt(sampleQ2+sampleU2+sampleV2);
      if(value < -1.0) {
 value = -1.0;
      }else if(value > 1.0) {
 value = 1.0;
      }
      dataEll[i] = 90.0*asin(value)/M_PI;
    }
    if(dataEllErr != NULL) {
      if(dataQ[i] == 0 && dataU[i] == 0 && dataV[i] == 0) {
 dataEllErr[i] = -1;
      }
      dataEllErr[i] = sampleV2*(sampleQ2*RMSQ*RMSQ+sampleU2*RMSU*RMSU);
      dataEllErr[i] += (sampleQ2+sampleU2)*(sampleQ2+sampleU2)*RMSV*RMSV;
      dataEllErr[i] /= 4.0*(sampleQ2+sampleU2);
      dataEllErr[i] = sqrt(dataEllErr[i]);
      dataEllErr[i] /= (sampleQ2+sampleU2+sampleV2);
      dataEllErr[i] *= 180.0/M_PI;
    }
  }
  free(Loffpulse);
  free(Poffpulse);
  return 1;
}
void make_paswing_fromIQUV_reportRMS(long pulsenr, long freqnr, int extended, int spstat, float rmsI, float rmsQ, float rmsU, float rmsV, float rmsL, float rmsP, float medianL, float medianP, float baseline_intensity, verbose_definition verbose)
{
  int indent;
  if(spstat == 0) {
    for(indent = 0; indent < verbose.indent; indent++) printf(" ");
    fprintf(stdout, "  PA conversion output for subint %ld frequency channel %ld:\n", pulsenr, freqnr);
  }
  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
  fprintf(stdout, "    Avrg baseline Stokes I:   %f (only reported, not subtracted)\n", baseline_intensity);
  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
  fprintf(stdout, "    RMS I:                    %f\n", rmsI);
  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
  fprintf(stdout, "    RMS Q:                    %f\n", rmsQ);
  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
  fprintf(stdout, "    RMS U:                    %f\n", rmsU);
  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
  fprintf(stdout, "    RMS V:                    %f\n", rmsV);
  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
  fprintf(stdout, "    RMS L (before de-bias):   %f\n", rmsL);
  if(extended) {
    for(indent = 0; indent < verbose.indent; indent++) printf(" ");
    fprintf(stdout, "    RMS sqrt(Q^2+U^2+V^2):    %f\n", rmsP);
  }
  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
  fprintf(stdout, "    Median L:                 %f\n", medianL);
  if(extended) {
    for(indent = 0; indent < verbose.indent; indent++) printf(" ");
    fprintf(stdout, "    Median sqrt(Q^2+U^2+V^2): %f\n", medianP);
  }
}
int make_paswing_fromIQUV(datafile_definition *datafile, int extended, int spstat, float sigma_limit, int sigmaI, pulselongitude_regions_definition onpulse, int normalize, int correctLbias, int correctPbias, float correctQV, float correctV, int nolongitudes, float loffset, float paoffset, datafile_definition *rms_file, float rebin_factor, int onpulseonly, verbose_definition verbose)
{
  int indent, output_nr_pols;
  long i, pulsenr, freqnr;
  float *newdata, *newdata_current_pulse;
  if(verbose.verbose) {
    for(indent = 0; indent < verbose.indent; indent++) printf(" ");
    printf("Constructing PA and degree of linear polarization");
    if(extended)
      printf(", total polarization and ellipticity");
    if(rms_file != NULL)
      printf(" (using a seperate file to determine the off-pulse rms)");
    printf("\n");
    for(indent = 0; indent < verbose.indent; indent++) printf(" ");
    printf("  Reference frequency for PA is ");
    if(datafile->isDeFarad) {
      if((datafile->freq_ref > -1.1 && datafile->freq_ref < -0.9) || (datafile->freq_ref > 0.99e10 && datafile->freq_ref < 1.01e10))
 printf("infinity\n");
      else if(datafile->freq_ref < 0)
 printf("unknown\n");
      else
 printf("%f MHz\n", datafile->freq_ref);
    }else {
      if(datafile->NrFreqChan == 1)
 printf("%lf MHz\n", get_centre_frequency(*datafile, verbose));
      else
 printf("observing frequencies of individual frequency channels\n");
    }
    for(indent = 0; indent < verbose.indent; indent++) printf(" ");
    printf("  ");
    switch(correctLbias) {
    case -1: printf("No L de-bias applied"); break;
    case 0: printf("De-bias L using median noise subtraction"); break;
    case 1: printf("De-bias L using Wardle & Kronberg correction"); break;
    default: printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Undefined L de-bias method specified."); return 0;
    }
    if(extended) {
      switch(correctPbias) {
      case -1: printf(", no P de-bias applied"); break;
      case 0: printf(", de-bias P using median noise subtraction"); break;
      case 1: printf(", de-bias P using Wardle & Kronberg like correction"); break;
      default: printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Undefined P de-bias method specified."); return 0;
      }
    }
    if(correctQV != 1 || correctV != 1)
      printf(", Q correction factor %f, V correction factor %f", 1.0/correctQV, 1.0/(correctQV*correctV));
    if(normalize)
      printf(", output is normalised");
    if(loffset != 0)
      printf(", pulse longitude shifted by %f deg\n", loffset);
    if(paoffset != 0)
      printf(", PA shifted by %f deg\n", paoffset);
    printf("\n");
  }
  if(datafile->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Expected 4 input polarizations.");
    return 0;
  }
  if(rms_file != NULL) {
    if(rms_file->NrPols != 4) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Expected 4 input polarizations.");
      return 0;
    }
  }
  if(datafile->poltype != POLTYPE_STOKES) {
    if(datafile->poltype == POLTYPE_UNKNOWN) {
      printwarning(verbose.debug, "WARNING make_paswing_fromIQUV: Polarization state unknown, it is assumed the data are Stokes parameters.");
    }else {
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Please convert data into Stokes parameters first.");
      return 0;
    }
  }
  if(rms_file != NULL) {
    if(rms_file->poltype != POLTYPE_STOKES) {
      if(rms_file->poltype == POLTYPE_UNKNOWN) {
 printwarning(verbose.debug, "WARNING make_paswing_fromIQUV: Polarization state of the data to be used to determine the off-pulse rms is unknown. It is assumed the data are Stokes parameters.");
      }else {
 printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Please convert data to be used to determine the off-pulse rms into Stokes parameters first.");
 return 0;
      }
    }
  }
  if(datafile->tsampMode != TSAMPMODE_FIXEDTSAMP) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: It is expected that the input has a uniform time sampling.");
    return 0;
  }
  if(correctQV == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: correctQV is set to zero, you probably want this to be 1.");
    return 0;
  }
  if(correctV == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: correctV is set to zero, you probably want this to be 1.");
    return 0;
  }
  if(datafile->isDebase == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Please remove baseline first, i.e. use pmod -debase.");
    return 0;
  }else if(datafile->isDebase != 1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING make_paswing_fromIQUV: Unknown baseline removal state. It will be assumed the baseline has already removed from the data.");
  }
  if(rms_file != NULL) {
    if(rms_file->isDebase == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Please remove baseline first, i.e. use pmod -debase.");
      return 0;
    }else if(rms_file->isDebase != 1) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING make_paswing_fromIQUV: Unknown baseline removal state. It is assumed the baseline has already removed from the data.");
    }
    if(datafile->NrSubints != rms_file->NrSubints) {
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Number of subintegrations is different in the data to be used to determine the off-pulse rms compared to the data used to compute the polarization information.");
 return 0;
    }
    if(datafile->NrFreqChan != rms_file->NrFreqChan) {
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Number of frequency channels is different in the data to be used to determine the off-pulse rms compared to the data used to compute the polarization information (%ld != %ld).", rms_file->NrFreqChan, datafile->NrFreqChan);
 return 0;
    }
    if(correctLbias == 0) {
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Subtracting the median of L is not supported when a separate file is used for the off-pulse statistics.");
      return 0;
    }
    if(extended && correctPbias == 0) {
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Subtracting the median of P is not supported when a separate file is used for the off-pulse statistics.");
      return 0;
    }
  }
  if(normalize && (datafile->NrSubints > 1 || datafile->NrFreqChan > 1)) {
    if(spstat == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING make_paswing_fromIQUV: Normalization of the polarization information will cause all subintegrations/frequency channels to be normalised individually. This may not be desired.");
    }
  }
  if(extended) {
    output_nr_pols = 8;
  }else {
    output_nr_pols = 5;
  }
  if(spstat == 0) {
    newdata = (float *)malloc(datafile->NrBins*datafile->NrSubints*datafile->NrFreqChan*output_nr_pols*sizeof(float));
  }else {
    newdata = (float *)malloc(datafile->NrBins*output_nr_pols*sizeof(float));
    newdata_current_pulse = (float *)malloc(datafile->NrBins*output_nr_pols*sizeof(float));
  }
  if(datafile->offpulse_rms != NULL) {
    free(datafile->offpulse_rms);
  }
  if(spstat == 0) {
    datafile->offpulse_rms = (float *)malloc(datafile->NrSubints*datafile->NrFreqChan*output_nr_pols*sizeof(float));
  }else {
    datafile->offpulse_rms = (float *)malloc(output_nr_pols*sizeof(float));
  }
  if(newdata == NULL || datafile->offpulse_rms == NULL
     ) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Memory allocation error.");
    return 0;
  }
  if(nolongitudes == 0) {
    datafile->tsamp_list = (double *)malloc(datafile->NrBins*sizeof(double));
    if(datafile->tsamp_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Memory allocation error.");
      return 0;
    }
    for(i = 0; i < datafile->NrBins; i++) {
      datafile->tsamp_list[i] = get_pulse_longitude(*datafile, 0, i, verbose);
      datafile->tsamp_list[i] += loffset;
    }
  }
  float *dataI, *dataQ, *dataU, *dataV, *newdataL, *newdataP, *newdataPa, *newdataPaErr, *newdataEll, *newdataEllErr;
  float baseline_intensity, rmsI, rmsQ, rmsU, rmsV, rmsL, rmsP, medianL, medianP;
  for(pulsenr = 0; pulsenr < datafile->NrSubints; pulsenr++) {
    for(freqnr = 0; freqnr < datafile->NrFreqChan; freqnr++) {
      int normalize_sp;
      dataI = &(datafile->data[datafile->NrBins*(0+datafile->NrPols*(freqnr+pulsenr*datafile->NrFreqChan))]);
      dataQ = &(datafile->data[datafile->NrBins*(1+datafile->NrPols*(freqnr+pulsenr*datafile->NrFreqChan))]);
      dataU = &(datafile->data[datafile->NrBins*(2+datafile->NrPols*(freqnr+pulsenr*datafile->NrFreqChan))]);
      dataV = &(datafile->data[datafile->NrBins*(3+datafile->NrPols*(freqnr+pulsenr*datafile->NrFreqChan))]);
      if(spstat == 0) {
 newdataL = &(newdata[datafile->NrBins*(1+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr))]);
 newdataPa = &(newdata[datafile->NrBins*(3+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr))]);
 newdataPaErr = &(newdata[datafile->NrBins*(4+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr))]);
 normalize_sp = normalize;
      }else {
 newdataL = &(newdata_current_pulse[datafile->NrBins*1]);
 newdataPa = NULL;
 newdataPaErr = NULL;
 normalize_sp = 0;
      }
      if(extended == 0) {
 newdataP = NULL;
 newdataEll = NULL;
 newdataEllErr = NULL;
      }else {
 if(spstat == 0) {
   newdataP = &(newdata[datafile->NrBins*(5+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr))]);
   newdataEll = &(newdata[datafile->NrBins*(6+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr))]);
   newdataEllErr = &(newdata[datafile->NrBins*(7+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr))]);
 }else {
   newdataP = &(newdata_current_pulse[datafile->NrBins*5]);
   newdataEll = NULL;
   newdataEllErr = NULL;
 }
      }
      long rms_file_nrBins;
      float *rms_file_I, *rms_file_Q, *rms_file_U, *rms_file_V;
      if(rms_file != NULL) {
 rms_file_nrBins = rms_file->NrBins;
 rms_file_I = &(rms_file->data[rms_file->NrBins*(0+rms_file->NrPols*(freqnr+pulsenr*rms_file->NrFreqChan))]);
 rms_file_Q = &(rms_file->data[rms_file->NrBins*(1+rms_file->NrPols*(freqnr+pulsenr*rms_file->NrFreqChan))]);
 rms_file_U = &(rms_file->data[rms_file->NrBins*(2+rms_file->NrPols*(freqnr+pulsenr*rms_file->NrFreqChan))]);
 rms_file_V = &(rms_file->data[rms_file->NrBins*(3+rms_file->NrPols*(freqnr+pulsenr*rms_file->NrFreqChan))]);
      }else {
 rms_file_nrBins = 0;
 rms_file_I = NULL;
 rms_file_Q = NULL;
 rms_file_U = NULL;
 rms_file_V = NULL;
      }
      if(make_paswing_fromIQUV_sp(dataI, dataQ, dataU, dataV, datafile->NrBins, newdataL, newdataP, newdataPa, newdataPaErr, newdataEll, newdataEllErr, &baseline_intensity, &rmsI, &rmsQ, &rmsU, &rmsV, &rmsL, &rmsP, &medianL, &medianP, onpulse, normalize_sp, correctLbias, correctPbias, correctQV, correctV, paoffset, rms_file_nrBins, rms_file_I, rms_file_Q, rms_file_U, rms_file_V, rebin_factor, verbose) == 0) {
 printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Calculating polarization products failed.");
 return 0;
      }
      pulselongitude_regions_definition *onpulse_ptr;
      onpulse_ptr = NULL;
      if(onpulseonly) {
 onpulse_ptr = &onpulse;
      }
      make_paswing_fromIQUV_remove_lowS2N_points_sp(sigma_limit, sigmaI, datafile->NrBins, newdataL, newdataP, newdataPa, newdataPaErr, newdataEll, newdataEllErr, rmsI, rmsL, rmsP, onpulse_ptr, verbose);
      if(spstat == 0) {
 for(i = 0; i < datafile->NrBins; i++) {
   newdata[datafile->NrBins*(0+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr))+i] = dataI[i];
   newdata[datafile->NrBins*(2+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr))+i] = dataV[i];
        }
 datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = rmsI;
 datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = rmsL;
 datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = rmsV;
 datafile->offpulse_rms[3+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = -1;
 datafile->offpulse_rms[4+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = -1;
 if(extended) {
   datafile->offpulse_rms[5+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = rmsP;
   datafile->offpulse_rms[6+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = -1;
   datafile->offpulse_rms[7+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = -1;
 }
      }
      if(verbose.verbose) {
 if(spstat == 0 && ((freqnr == 0 && pulsenr == 0) || verbose.debug)) {
   if(datafile->NrFreqChan > 1 || datafile->NrSubints > 1) {
     for(indent = 0; indent < verbose.indent; indent++) printf(" ");
     fprintf(stdout, "    Statistics based on first processed pulse\n");
   }
   if(extended) {
     make_paswing_fromIQUV_reportRMS(pulsenr, freqnr, extended, spstat, datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)], rmsQ, rmsU, datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)], datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)], datafile->offpulse_rms[5+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)], medianL, medianP, baseline_intensity, verbose);
   }else {
     make_paswing_fromIQUV_reportRMS(pulsenr, freqnr, extended, spstat, datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)], rmsQ, rmsU, datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)], datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)], 0.0, medianL, 0.0, baseline_intensity, verbose);
   }
 }
      }
    }
  }
  free(datafile->data);
  datafile->data = newdata;
  if(nolongitudes == 0) {
    datafile->tsampMode = TSAMPMODE_LONGITUDELIST;
  }
  datafile->NrPols = output_nr_pols;
  if(extended) {
    datafile->poltype = POLTYPE_ILVPAdPATEldEl;
  }else {
    datafile->poltype = POLTYPE_ILVPAdPA;
  }
  return 1;
}
void paswing_remove_observed_PA_swing_sp(float *dataPA, float *dataPAerr, float *dataPAref, float *dataPArefErr, int nrBins, int add, verbose_definition verbose)
{
  int ok;
  long j;
  for(j = 0; j < nrBins; j++) {
    ok = 1;
    if(dataPAerr != NULL) {
      if(dataPAerr[j] < 0) {
 ok = 0;
      }
    }
    if(dataPArefErr != NULL) {
      if(dataPArefErr[j] < 0) {
 ok = 0;
 dataPA[j] = 0;
 if(dataPAerr != NULL) {
   dataPAerr[j] = -1;
 }
      }
    }
    if(ok) {
      if(add) {
 dataPA[j] += dataPAref[j];
      }else {
 dataPA[j] -= dataPAref[j];
      }
      dataPA[j] = derotate_180(dataPA[j]) - 90;
      if(dataPAerr != NULL && dataPArefErr != NULL) {
 dataPAerr[j] = sqrt(dataPAerr[j]*dataPAerr[j]+dataPArefErr[j]*dataPArefErr[j]);
      }
    }else {
      dataPA[j] = 0;
      if(dataPAerr != NULL) {
 dataPAerr[j] = -1;
      }
    }
  }
}
int paswing_remove_observed_PA_swing(datafile_definition *datafile, datafile_definition datafile_reference, int add, verbose_definition verbose)
{
  int pachannel_ref, pachannelerr_ref, pachannel, pachannelerr;
  if(datafile->poltype != POLTYPE_ILVPAdPA && datafile->poltype != POLTYPE_PAdPA && datafile->poltype != POLTYPE_ILVPAdPATEldEl) {
    printerror(verbose.debug, "ERROR paswing_remove_observed_PA_swing: Data doesn't appear to have poltype ILVPAdPA, PAdPA or ILVPAdPATEldEl.");
    return 0;
  }
  if(datafile_reference.poltype != POLTYPE_ILVPAdPA && datafile_reference.poltype != POLTYPE_PAdPA && datafile_reference.poltype != POLTYPE_ILVPAdPATEldEl) {
    printerror(verbose.debug, "ERROR paswing_remove_observed_PA_swing: Data in the reference doesn't appear to have poltype ILVPAdPA, PAdPA or ILVPAdPATEldEl.");
    return 0;
  }
  if(datafile->poltype == POLTYPE_ILVPAdPA && datafile->NrPols != 5) {
    printerror(verbose.debug, "ERROR paswing_remove_observed_PA_swing: 5 polarization channels were expected, but there are only %ld.", datafile->NrPols);
    return 0;
  }else if(datafile->poltype == POLTYPE_ILVPAdPATEldEl && datafile->NrPols != 8) {
    printerror(verbose.debug, "ERROR paswing_remove_observed_PA_swing: 8 polarization channels were expected, but there are only %ld.", datafile->NrPols);
    return 0;
  }else if(datafile->poltype == POLTYPE_PAdPA && datafile->NrPols != 2) {
    printerror(verbose.debug, "ERROR paswing_remove_observed_PA_swing: 2 polarization channels were expected, but there are only %ld.", datafile->NrPols);
    return 0;
  }
  if(datafile_reference.poltype == POLTYPE_ILVPAdPA && datafile_reference.NrPols != 5) {
    printerror(verbose.debug, "ERROR paswing_remove_observed_PA_swing: 5 polarization channels were expected in the reference, but there are only %ld.", datafile_reference.NrPols);
    return 0;
  }else if(datafile_reference.poltype == POLTYPE_ILVPAdPATEldEl && datafile_reference.NrPols != 8) {
    printerror(verbose.debug, "ERROR paswing_remove_observed_PA_swing: 8 polarization channels were expected in the reference, but there are only %ld.", datafile_reference.NrPols);
    return 0;
  }else if(datafile_reference.poltype == POLTYPE_PAdPA && datafile_reference.NrPols != 2) {
    printerror(verbose.debug, "ERROR paswing_remove_observed_PA_swing: 2 polarization channels were expected in the reference, but there are only %ld.", datafile_reference.NrPols);
    return 0;
  }
  if(datafile_reference.NrBins != datafile->NrBins) {
    printerror(verbose.debug, "ERROR paswing_remove_observed_PA_swing: Mismatch in the number of bins in the data and the reference (%ld != %ld).", datafile->NrBins, datafile_reference.NrBins);
    return 0;
  }
  if(datafile_reference.NrFreqChan != 1) {
    printerror(verbose.debug, "ERROR paswing_remove_observed_PA_swing: The number of frequency channels in the reference should be 1 (it is %ld).", datafile_reference.NrFreqChan);
    return 0;
  }
  if(datafile_reference.NrSubints != 1) {
    printerror(verbose.debug, "ERROR paswing_remove_observed_PA_swing: The number of subints in the reference should be 1 (it is %ld).", datafile_reference.NrSubints);
    return 0;
  }
  if(datafile->poltype == POLTYPE_ILVPAdPA || datafile->poltype == POLTYPE_ILVPAdPATEldEl) {
    pachannel = 3;
    pachannelerr = 4;
  }else if(datafile->poltype == POLTYPE_PAdPA) {
    pachannel = 0;
    pachannelerr = 1;
  }
  if(datafile_reference.poltype == POLTYPE_ILVPAdPA || datafile_reference.poltype == POLTYPE_ILVPAdPATEldEl) {
    pachannel_ref = 3;
    pachannelerr_ref = 4;
  }else if(datafile_reference.poltype == POLTYPE_PAdPA) {
    pachannel_ref = 0;
    pachannelerr_ref = 1;
  }
  long i, f;
  for(i = 0; i < datafile->NrSubints; i++) {
    for(f = 0; f < datafile->NrFreqChan; f++) {
      paswing_remove_observed_PA_swing_sp(&(datafile->data[datafile->NrBins*(pachannel + datafile->NrPols*(f+datafile->NrFreqChan*i))]), &(datafile->data[datafile->NrBins*(pachannelerr + datafile->NrPols*(f+datafile->NrFreqChan*i))]), &(datafile_reference.data[datafile->NrBins*pachannel_ref]), &(datafile_reference.data[datafile->NrBins*pachannelerr_ref]), datafile->NrBins, add, verbose);
    }
  }
  return 1;
}
int writePPOLHeader(datafile_definition datafile, int argc, char **argv, verbose_definition verbose)
{
  char *txt;
  txt = malloc(10000);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePPOLHeader: Memory allocation error.");
    return 0;
  }
  constructCommandLineString(txt, 10000, argc, argv, verbose);
  fprintf(datafile.fptr_hdr, "#ppol file: %s\n", txt);
  free(txt);
  return 1;
}
int readPPOLHeader(datafile_definition *datafile, int extended, verbose_definition verbose)
{
  float dummy_float;
  int ret, maxlinelength, nrwords;
  char *txt, *ret_ptr, *word_ptr;
  maxlinelength = 2000;
  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLHeader: Memory allocation error.");
    return 0;
  }
  datafile->isFolded = 1;
  datafile->foldMode = FOLDMODE_FIXEDPERIOD;
  datafile->fixedPeriod = 0;
  datafile->tsampMode = TSAMPMODE_LONGITUDELIST;
  datafile->fixedtsamp = 0;
  datafile->tsubMode = TSUBMODE_FIXEDTSUB;
  if(datafile->tsub_list != NULL)
    free(datafile->tsub_list);
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLHeader: Memory allocation error");
    return 0;
  }
  datafile->tsub_list[0] = 0;
  datafile->NrSubints = 1;
  datafile->NrFreqChan = 1;
  datafile->datastart = 0;
  rewind(datafile->fptr);
  ret = fread(txt, 1, 3, datafile->fptr);
  txt[3] = 0;
  if(ret != 3) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLHeader: cannot read from file.");
    free(txt);
    return 0;
  }
  if(strcmp(txt, "#pp") != 0
) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING readPPOLHeader: File does not appear to be in PPOL or PPOLSHORT format. I will try to load file, but this will probably fail. Did you run ppol first?");
  }
  skipallhashedlines(datafile);
  datafile->NrBins = 0;
  dummy_float = 0;
  do {
    ret_ptr = fgets(txt, maxlinelength, datafile->fptr);
    if(ret_ptr != NULL) {
      if(txt[0] != '#') {
 if(extended) {
   word_ptr = pickWordFromString(txt, 2, &nrwords, 1, ' ', verbose);
   if(nrwords != 10 && nrwords != 14) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPPOLHeader: Line should have 10 or 14 words, got %d", nrwords);
     if(nrwords == 3)
       printerror(verbose.debug, "                             Maybe file is in format %s?", returnFileFormat_str(PPOL_SHORT_format));
     printerror(verbose.debug, "                             Line: '%s'.", txt);
     free(txt);
     return 0;
   }
   if(nrwords == 10) {
     datafile->poltype = POLTYPE_ILVPAdPA;
     datafile->NrPols = 5;
   }else {
     datafile->poltype = POLTYPE_ILVPAdPATEldEl;
     datafile->NrPols = 8;
   }
 }else {
   word_ptr = pickWordFromString(txt, 1, &nrwords, 1, ' ', verbose);
   if(nrwords != 3) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR readPPOLHeader: Line should have 3 words, got %d", nrwords);
     if(nrwords == 10)
       printerror(verbose.debug, "                             Maybe file is in format %s?", returnFileFormat_str(PPOL_format));
     printerror(verbose.debug, "                             Line: '%s'.", txt);
     free(txt);
     return 0;
   }
 }
 ret = sscanf(word_ptr, "%f", &dummy_float);
 if(ret != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR readPPOLHeader: Cannot interpret as a float: '%s'.", txt);
   free(txt);
   return 0;
 }
 if(dummy_float >= 360) {
   fflush(stdout);
   printwarning(verbose.debug, "WARNING: IGNORING POINTS AT PULSE LONGITUDES > 360 deg.");
 }else {
   (datafile->NrBins)++;
 }
      }
    }
  }while(ret_ptr != NULL && dummy_float < 360);
  if(extended == 0) {
    datafile->poltype = POLTYPE_PAdPA;
    datafile->NrPols = 2;
  }
  fflush(stdout);
  if(verbose.verbose) fprintf(stdout, "Going to load %ld points from %s\n", datafile->NrBins, datafile->filename);
  if(datafile->NrBins == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLHeader: No data in %s", datafile->filename);
    free(txt);
    return 0;
  }
  fseek(datafile->fptr, datafile->datastart, SEEK_SET);
  free(txt);
  if(datafile->offpulse_rms != NULL) {
    free(datafile->offpulse_rms);
    datafile->offpulse_rms = NULL;
  }
  if(extended) {
    datafile->offpulse_rms = (float *)malloc(datafile->NrSubints*datafile->NrFreqChan*datafile->NrPols*sizeof(float));
    if(datafile->offpulse_rms == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPPOLHeader: Memory allocation error");
      return 0;
    }
  }
  datafile->tsamp_list = (double *)malloc(datafile->NrBins*sizeof(double));
  if(datafile->tsamp_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLHeader: Memory allocation error");
    return 0;
  }
  return 1;
}
int readPPOLfile(datafile_definition *datafile, float *data, int extended, float add_longitude_shift, verbose_definition verbose)
{
  int maxlinelength;
  long i, k, dummy_long;
  char *txt, *ret_ptr;
  if(datafile->NrBins == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLfile: No data in %s", datafile->filename);
    return 0;
  }
  maxlinelength = 2000;
  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLfile: Memory allocation error.");
    return 0;
  }
  fseek(datafile->fptr, datafile->datastart, SEEK_SET);
  k = 0;
  if(extended) {
    datafile->offpulse_rms[3] = -1;
    datafile->offpulse_rms[4] = -1;
    if(datafile->NrPols == 8) {
      datafile->offpulse_rms[6] = -1;
      datafile->offpulse_rms[7] = -1;
    }
  }
  for(i = 0; i < datafile->NrBins; i++) {
    ret_ptr = fgets(txt, maxlinelength, datafile->fptr);
    if(ret_ptr == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPPOLfile: Cannot read next line, should not happen after successfully reading in header");
      free(txt);
      return 0;
    }
    if(txt[0] != '#') {
      if(extended == 0) {
 sscanf(txt, "%lf %f %f", &(datafile->tsamp_list[k]), &(data[k]), &(data[k+datafile->NrBins]));
      }else {
 if(datafile->NrPols == 8) {
   sscanf(txt, "%ld %lf %f %f %f %f %f %f %f %f %f %f %f %f", &dummy_long, &(datafile->tsamp_list[k]), &(data[k]), &(datafile->offpulse_rms[0]), &(data[k+datafile->NrBins]), &(datafile->offpulse_rms[1]), &(data[k+2*datafile->NrBins]), &(datafile->offpulse_rms[2]), &(data[k+3*datafile->NrBins]), &(data[k+4*datafile->NrBins]), &(data[k+5*datafile->NrBins]), &(datafile->offpulse_rms[5]), &(data[k+6*datafile->NrBins]), &(data[k+7*datafile->NrBins]));
 }else {
   sscanf(txt, "%ld %lf %f %f %f %f %f %f %f %f", &dummy_long, &(datafile->tsamp_list[k]), &(data[k]), &(datafile->offpulse_rms[0]), &(data[k+datafile->NrBins]), &(datafile->offpulse_rms[1]), &(data[k+2*datafile->NrBins]), &(datafile->offpulse_rms[2]), &(data[k+3*datafile->NrBins]), &(data[k+4*datafile->NrBins]));
 }
      }
      datafile->tsamp_list[k] += add_longitude_shift;
      if(datafile->tsamp_list[k] >= 0 && datafile->tsamp_list[k] < 360) {
   k++;
      }else {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING readPPOLfile: IGNORING POINTS AT PULSE LONGITUDES outside range 0 ... 360 deg.");
      }
    }
  }
  if(k != datafile->NrBins) {
    fflush(stdout);
    printerror(verbose.debug, "WARNING readPPOLfile: The nr of bins read in is different as determined from header. Something is wrong.");
    return 0;
  }
  fflush(stdout);
  if(verbose.verbose) fprintf(stdout, "readPPOLfile: Accepted %ld points\n", datafile->NrBins);
     free(txt);
  return 1;
}
int writePPOLfile(datafile_definition datafile, float *data, int extended, int onlysignificantPA, int twoprofiles, float PAoffset, verbose_definition verbose)
{
  long j;
  if(datafile.poltype != POLTYPE_ILVPAdPA && datafile.poltype != POLTYPE_PAdPA && datafile.poltype != POLTYPE_ILVPAdPATEldEl) {
    printerror(verbose.debug, "ERROR writePPOLfile: Data doesn't appear to have poltype ILVPAdPA, PAdPA or ILVPAdPATEldEl (it is %d).", datafile.poltype);
    return 0;
  }
  if(datafile.poltype == POLTYPE_ILVPAdPA && datafile.NrPols != 5) {
    printerror(verbose.debug, "ERROR writePPOLfile: 5 polarization channels were expected, but there are %ld.", datafile.NrPols);
    return 0;
  }else if(datafile.poltype == POLTYPE_PAdPA && datafile.NrPols != 2) {
    printerror(verbose.debug, "ERROR writePPOLfile: 2 polarization channels were expected, but there are %ld.", datafile.NrPols);
    return 0;
  }else if(datafile.poltype == POLTYPE_ILVPAdPATEldEl && datafile.NrPols != 8) {
    printerror(verbose.debug, "ERROR writePPOLfile: 8 polarization channels were expected, but there are %ld.", datafile.NrPols);
    return 0;
  }
  if(datafile.NrSubints > 1 || datafile.NrFreqChan > 1) {
    printerror(verbose.debug, "ERROR writePPOLfile: Can only do this opperation if there is one subint and one frequency channel.");
    return 0;
  }
  if(datafile.tsampMode != TSAMPMODE_LONGITUDELIST) {
    printerror(verbose.debug, "ERROR writePPOLfile: Expected pulse longitudes to be defined.");
    return 0;
  }
  int pa_offset, dpa_offset;
  if(datafile.poltype == POLTYPE_ILVPAdPA || datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
    pa_offset = 3;
    dpa_offset = 4;
  }else if(datafile.poltype == POLTYPE_PAdPA) {
    pa_offset = 0;
    dpa_offset = 1;
  }
  for(j = 0; j < datafile.NrBins; j++) {
    if(data[j+dpa_offset*datafile.NrBins] > 0 || onlysignificantPA == 0) {
      if(extended) {
 fprintf(datafile.fptr, "%ld %e %e %e %e %e %e %e %e %e", j, datafile.tsamp_list[j], data[j], datafile.offpulse_rms[0], data[j+datafile.NrBins], datafile.offpulse_rms[1], data[j+2*datafile.NrBins], datafile.offpulse_rms[2], data[j+pa_offset*datafile.NrBins]+PAoffset, data[j+dpa_offset*datafile.NrBins]);
 if(datafile.poltype == POLTYPE_ILVPAdPATEldEl)
   fprintf(datafile.fptr, " %e %e %e %e", data[j+5*datafile.NrBins], datafile.offpulse_rms[5], data[j+6*datafile.NrBins], data[j+7*datafile.NrBins]);
 fprintf(datafile.fptr, "\n");
      }else {
 fprintf(datafile.fptr, "%e %e %e\n", datafile.tsamp_list[j], data[j+pa_offset*datafile.NrBins]+PAoffset, data[j+dpa_offset*datafile.NrBins]);
      }
    }
  }
  if(twoprofiles) {
    for(j = 0; j < datafile.NrBins; j++) {
      if(data[j+dpa_offset*datafile.NrBins] > 0 || onlysignificantPA == 0) {
 if(extended) {
   fprintf(datafile.fptr, "%ld %e %e %e %e %e %e %e %e %e", j, datafile.tsamp_list[j]+360, data[j], datafile.offpulse_rms[0], data[j+datafile.NrBins], datafile.offpulse_rms[1], data[j+2*datafile.NrBins], datafile.offpulse_rms[2], data[j+pa_offset*datafile.NrBins]+PAoffset, data[j+dpa_offset*datafile.NrBins]);
   if(datafile.poltype == POLTYPE_ILVPAdPATEldEl)
     fprintf(datafile.fptr, " %e %e %e %e", data[j+5*datafile.NrBins], datafile.offpulse_rms[5], data[j+6*datafile.NrBins], data[j+7*datafile.NrBins]);
   fprintf(datafile.fptr, "\n");
 }else {
   fprintf(datafile.fptr, "%e %e %e\n", datafile.tsamp_list[j]+360, data[j]+PAoffset, data[j+dpa_offset*datafile.NrBins]);
 }
      }
    }
  }
  return 1;
}
int make_pa_distribution(datafile_definition datain, datafile_definition *dataout, int nrbins, int normalise, int weighttype, datafile_definition *pamask, float pamask_value, int ellipticity, verbose_definition verbose)
{
  long i, j, f, nrpointsadded, nrpointsadded_max, binnr;
  float dpa;
  if(datain.NrSubints <= 1 && datain.NrFreqChan <= 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_pa_distribution: Need more than a single subints and frequency channel to make a PA distribution");
    return 0;
  }
  if(datain.poltype != POLTYPE_ILVPAdPA && datain.poltype != POLTYPE_PAdPA && datain.poltype != POLTYPE_ILVPAdPATEldEl) {
    printerror(verbose.debug, "ERROR make_pa_distribution: Data doesn't appear to have poltype ILVPAdPA, ILVPAdPATEldEl or PAdPA.");
    return 0;
  }
  if(ellipticity && datain.poltype != POLTYPE_ILVPAdPATEldEl) {
    printerror(verbose.debug, "ERROR make_pa_distribution: The data isn't of polarization type ILVPAdPATEldEl, while the ellipticity distribution was requested.");
    return 0;
  }
  if(datain.poltype == POLTYPE_ILVPAdPA && datain.NrPols != 5) {
    printerror(verbose.debug, "ERROR make_pa_distribution: 5 polarization channels were expected, but there are %ld.", datain.NrPols);
    return 0;
  }else if(datain.poltype == POLTYPE_ILVPAdPATEldEl && datain.NrPols != 8) {
    printerror(verbose.debug, "ERROR make_pa_distribution: 8 polarization channels were expected, but there are %ld.", datain.NrPols);
    return 0;
  }else if(datain.poltype == POLTYPE_PAdPA && datain.NrPols != 2) {
    printerror(verbose.debug, "ERROR make_pa_distribution: 2 polarization channels were expected, but there are %ld.", datain.NrPols);
    return 0;
  }
  if(pamask != NULL) {
    if(datain.NrBins != pamask->NrBins) {
      printerror(verbose.debug, "ERROR make_pa_distribution: Applying a PA mask only works when the input data has the same number of pulse longitude bins compared to that of the provided mask. (the input data has %ld pulse longitude bins, while the mask has %ld).", datain.NrBins, pamask->NrBins);
      return 0;
    }
    if(nrbins != pamask->NrSubints) {
      printerror(verbose.debug, "ERROR make_pa_distribution: Applying a PA mask only works when generating a PA-distribution with an equal number of PA bins compared to that of the provided mask. (now %ld pa bins are requested, while the mask has %ld pa-bins defined).", nrbins, pamask->NrSubints);
      return 0;
    }
    if(pamask->NrFreqChan > 1) {
      printerror(verbose.debug, "ERROR make_pa_distribution: Applying a PA mask only works when the mask has one frequency channel defined. There are currently %ld channels defined).", pamask->NrFreqChan);
      return 0;
    }
  }
  cleanPSRData(dataout, verbose);
  copy_params_PSRData(datain, dataout, verbose);
  dataout->format = MEMORY_format;
  dataout->NrSubints = nrbins;
  dataout->NrPols = 1;
  dataout->NrFreqChan = 1;
  if(ellipticity == 0)
    dataout->gentype = GENTYPE_PADIST;
  else
    dataout->gentype = GENTYPE_ELLDIST;
  dataout->tsubMode = TSUBMODE_FIXEDTSUB;
  if(dataout->tsub_list != NULL)
    free(dataout->tsub_list);
  dataout->tsub_list = (double *)malloc(sizeof(double));
  if(dataout->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_pa_distribution: Memory allocation error");
    return 0;
  }
  dataout->tsub_list[0] = get_tobs(datain, verbose);
  dataout->yrangeset = 1;
  if(ellipticity == 0) {
    dataout->yrange[0] = -90.0+0.5*180.0/(float)(nrbins);
    dataout->yrange[1] = 90.0-0.5*180.0/(float)(nrbins);
  }else {
    dataout->yrange[0] = -45.0+0.5*90.0/(float)(nrbins);
    dataout->yrange[1] = 45.0-0.5*90.0/(float)(nrbins);
  }
  dataout->data = (float *)calloc(dataout->NrSubints*dataout->NrBins*dataout->NrPols*dataout->NrFreqChan, sizeof(float));
  if(dataout->data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_pa_distribution: Cannot allocate memory for data.");
    return 0;
  }
  if(ellipticity == 0) {
    dpa = 180.0/(float)nrbins;
  }else {
    dpa = 90.0/(float)nrbins;
  }
  int pa_chan, dpa_chan, weight_chan;
  if(datain.poltype == POLTYPE_ILVPAdPA || datain.poltype == POLTYPE_ILVPAdPATEldEl) {
    if(ellipticity == 0) {
      pa_chan = 3;
      dpa_chan = 4;
    }else {
      pa_chan = 6;
      dpa_chan = 7;
    }
    if(weighttype == 0) {
      weight_chan = -1;
    }else if(weighttype == 1) {
      weight_chan = 1;
    }else if(weighttype == 2) {
      weight_chan = 2;
    }else if(weighttype == 3) {
      weight_chan = 0;
    }else if(weighttype == 4) {
      weight_chan = 5;
      if(datain.poltype != POLTYPE_ILVPAdPATEldEl) {
 printerror(verbose.debug, "ERROR make_pa_distribution: Weighting by the total polarization is requested, but that is not appears to be defined in the input data.");
 return 0;
      }
    }else {
      printerror(verbose.debug, "ERROR make_pa_distribution: Unsupported weighttype is specified.");
      return 0;
    }
  }else {
    pa_chan = 0;
    dpa_chan = 1;
    if(weighttype != 0) {
      printerror(verbose.debug, "ERROR make_pa_distribution: Data only has PA values defined, so weighting is not supported.");
      return 0;
    }
  }
  nrpointsadded_max = 0;
  for(j = 0; j < datain.NrBins; j++) {
    nrpointsadded = 0;
    for(i = 0; i < datain.NrSubints; i++) {
      for(f = 0; f < datain.NrFreqChan; f++) {
 float paerr;
 paerr = datain.data[j+datain.NrBins*(dpa_chan+datain.NrPols*(f+datain.NrFreqChan*i))];
 if(paerr > 0) {
   float pa = derotate_180_small_double(datain.data[j+datain.NrBins*(pa_chan+datain.NrPols*(f+datain.NrFreqChan*i))]);
   float weight = 1.0;
   if(weighttype != 0) {
     weight = datain.data[j+datain.NrBins*(weight_chan+datain.NrPols*(f+datain.NrFreqChan*i))];
     if(weighttype == 2) {
       weight = fabs(weight);
     }
   }
   if(ellipticity == 0) {
     if(pa == 90.0)
       binnr = 0;
     else
       binnr = (pa + 90.0)/dpa;
   }else {
     if(pa == 45.0)
       binnr = 0;
     else
       binnr = (pa + 45.0)/dpa;
   }
   if(binnr < 0 || binnr >= nrbins) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR make_pa_distribution: %ld %f %f BUG!!!!!!!!!!!!!!!", binnr, datain.data[j+datain.NrBins*(pa_chan+datain.NrPols*(f+datain.NrFreqChan*i))], dpa);
     return 0;
   }
   dataout->data[j+dataout->NrBins*binnr] += weight;
   nrpointsadded++;
 }
 if(pamask != NULL) {
   float value;
   if(paerr <= 0) {
     value = 0;
   }else {
     value = pamask->data[j+pamask->NrBins*(0+pamask->NrPols*(0+pamask->NrFreqChan*binnr))];
   }
   if(isnan(pamask_value)) {
     int curpol;
     for(curpol = 0; curpol < datain.NrPols; curpol++) {
       datain.data[j+datain.NrBins*(curpol+datain.NrPols*(f+datain.NrFreqChan*i))] = value;
     }
   }else {
     if(value < pamask_value-0.01 || value > pamask_value+0.01 || paerr <= 0) {
       int curpol;
       for(curpol = 0; curpol < datain.NrPols; curpol++) {
  datain.data[j+datain.NrBins*(curpol+datain.NrPols*(f+datain.NrFreqChan*i))] = 0;
       }
     }
   }
 }
      }
    }
    if(nrpointsadded > nrpointsadded_max)
      nrpointsadded_max = nrpointsadded;
  }
  if(normalise) {
    if(nrpointsadded_max > 0) {
      for(i = 0; i < nrbins; i++) {
 for(j = 0; j < datain.NrBins; j++) {
   dataout->data[j+datain.NrBins*i] /= (float)nrpointsadded_max;
 }
      }
    }
  }
  return 1;
}
int make_polarization_projection_map(datafile_definition datafile, float *map, int nrx, int nry, float background, int binnr, pulselongitude_regions_definition onpulse, int weighting, float threshold, int projection, float rot_long, float rot_lat, float conalselection, datafile_definition *subtract_pa_data, verbose_definition verbose)
{
  int ok;
  long i, xi, yi, pulsenr;
  float longitude, stokesI, L, P, latitude, x, y, weight;
  if(projection < 1 || projection > 3) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_projection_map_formIQUV: Projection type is not implemented.");
    return 0;
  }
  if(datafile.NrFreqChan != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_projection_map_formIQUV: Expected 1 frequency channel.");
    return 0;
  }
  if(datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
    if(datafile.NrPols != 8) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_projection_map_formIQUV: Expected 8 input polarizations when reading in data containing PA's and ellipticities.");
      return 0;
    }
  }else {
    if(datafile.NrPols != 4) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_projection_map_formIQUV: Expected 4 input polarizations.");
      return 0;
    }
  }
  if(subtract_pa_data != NULL) {
    if(subtract_pa_data->NrBins != datafile.NrBins) {
      printerror(verbose.debug, "ERROR make_projection_map_formIQUV: The reference PA-swing has a different number of pulse phase bins compared to the input data.");
      return 0;
    }
    if(subtract_pa_data->NrFreqChan != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_projection_map_formIQUV: The reference PA-swing should have a single frequency channel.");
      return 0;
    }
    if(subtract_pa_data->NrSubints != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_projection_map_formIQUV: The reference PA-swing should have a single subint.");
      return 0;
    }
  }
  rot_long *= M_PI/180.0;
  rot_lat *= M_PI/180.0;
  float *normmap;
  if(weighting == 2) {
    normmap = malloc(nrx*nry*sizeof(float));
    if(normmap == NULL) {
      printerror(verbose.debug, "ERROR make_projection_map_formIQUV: Memory allocation error.");
      return 0;
    }
  }
  for(xi = 0; xi < nrx; xi++) {
    for(yi = 0; yi < nry; yi++) {
      map[xi+nrx*yi] = 0;
      if(weighting == 2) {
 normmap[xi+nrx*yi] = 0;
      }
    }
  }
  for(pulsenr = 0; pulsenr < datafile.NrSubints; pulsenr++) {
    for(i = 0; i < (datafile.NrBins); i++) {
      ok = 1;
      if(binnr >= 0) {
 if(i != binnr)
   ok = 0;
      }else if(binnr == -1) {
 if(checkRegions(i, &onpulse, 0, verbose) == 0) {
   ok = 0;
 }
      }
      if(ok) {
 if(datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
   longitude = 2.0*datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+3*(datafile.NrBins)]*M_PI/180.0;
   L = datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+1*(datafile.NrBins)];
   latitude = 2.0*datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+6*(datafile.NrBins)]*M_PI/180.0;
   if(datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+4*(datafile.NrBins)] < 0 || datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+7*(datafile.NrBins)] < 0) {
     longitude = latitude = sqrt(-1.0);
   }
 }else {
   longitude = atan2(datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+2*(datafile.NrBins)],datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+(datafile.NrBins)]);
   L = sqrt(datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+2*(datafile.NrBins)]*datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+2*(datafile.NrBins)] + datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+(datafile.NrBins)]*datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+(datafile.NrBins)]);
   latitude = atan(datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+3*(datafile.NrBins)]/L);
 }
 if(subtract_pa_data != NULL) {
   int pachannel_subtract_fin;
   if(subtract_pa_data->poltype == POLTYPE_ILVPAdPATEldEl) {
     pachannel_subtract_fin = 3;
   }else {
     pachannel_subtract_fin = subtract_pa_data->NrPols-2;
   }
   longitude -= 2.0*subtract_pa_data->data[i+subtract_pa_data->NrBins*(pachannel_subtract_fin)]*M_PI/180.0;
 }
 if(!isnan(latitude)) {
   if(conalselection > 0) {
     double sphericaldistance;
     sphericaldistance = acos(cos(latitude-rot_lat)*cos(longitude+rot_long))*180.0/M_PI;
     if(sphericaldistance > conalselection) {
       latitude = sqrt(-1.0);
     }
   }
 }
 if(!isnan(latitude)) {
   if(weighting) {
     if(weighting != 3) {
       if(datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
  P = datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+5*(datafile.NrBins)];
       }else {
  P = sqrt(L*L+datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+3*(datafile.NrBins)]*datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+3*(datafile.NrBins)]);
       }
     }
     if(weighting == 2 || weighting == 3) {
       stokesI = datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr];
     }
   }
   if(projection == 1) {
     projectionHammerAitoff_xy(longitude, latitude, rot_long, rot_lat, &x, &y);
     weight = 1;
   }else if(projection == 2) {
     projection_sphere_xy(longitude, latitude, rot_long, rot_lat, &x, &y, &weight);
   }else if(projection == 3) {
     projection_longlat_xy(longitude, latitude, rot_long, rot_lat, &x, &y);
     weight = 1;
     x /= 80.0;
     y /= 80.0;
   }
   xi = 0.5*nrx + x*nrx/4.5;
   yi = 0.5*nry + y*nry/2.25;
   if(weighting == 0) {
     map[xi+nrx*yi] += 1.0*weight;
   }else {
     if(weighting == 3) {
       map[xi+nrx*yi] += stokesI*weight;
     }else {
       map[xi+nrx*yi] += P*weight;
     }
     if(weighting == 2) {
       normmap[xi+nrx*yi] += stokesI*weight;
     }
   }
 }
      }
    }
  }
  if((weighting == 1 || weighting == 2) && threshold > 0) {
    float maxvalue = map[0];
    for(xi = 0; xi < nrx; xi++) {
      for(yi = 0; yi < nry; yi++) {
 if(map[xi+nrx*yi] > maxvalue) {
   maxvalue = map[xi+nrx*yi];
 }
      }
    }
    for(xi = 0; xi < nrx; xi++) {
      for(yi = 0; yi < nry; yi++) {
 if(map[xi+nrx*yi] < threshold*maxvalue) {
   map[xi+nrx*yi] = 0;
 }
      }
    }
  }
  if(weighting == 2) {
    for(xi = 0; xi < nrx; xi++) {
      for(yi = 0; yi < nry; yi++) {
 if(normmap[xi+nrx*yi] != 0.0) {
   map[xi+nrx*yi] /= normmap[xi+nrx*yi];
   if(map[xi+nrx*yi] < 0) {
     map[xi+nrx*yi] = 0;
   }
   if(map[xi+nrx*yi] > 1) {
     map[xi+nrx*yi] = 1;
   }
 }else {
   map[xi+nrx*yi] = 0;
 }
      }
    }
    free(normmap);
  }
  return 1;
}

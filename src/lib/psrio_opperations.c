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
#include <math.h>
#include "psrsalsa.h"
int rebinPulse(float *Ipulse, long NrBins, float *Ipulse2, long NrBins2, int noDependencyWarning, verbose_definition verbose)
{
  long j, i1, i2;
  float x, x2;
  if(noDependencyWarning == 0) {
    if(NrBins % NrBins2 != 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING rebinPulse: Rebinning from %ld to %ld bins implies that separate bins are not entirely independent.", NrBins, NrBins2);
    }
  }
  for(j = 0; j < NrBins2; j++)
    Ipulse2[j] = 0;
  for(j = 0; j < NrBins; j++) {
    x = (j)/(float)NrBins;
    x *= NrBins2;
    x2 = (j+1)/(float)NrBins;
    x2 *= NrBins2;
    i1 = x;
    i2 = x2;
    if(i1 == i2) {
      Ipulse2[i1] += Ipulse[j]*(x2-x);
    }else if(i2-i1 == 1) {
      Ipulse2[i1] += Ipulse[j]*(i2-x);
      if(i2 < NrBins2)
 Ipulse2[i2] += Ipulse[j]*(x2-i2);
    }else if(i2-i1 == 2) {
      Ipulse2[i1] += Ipulse[j]*(i1+1-x);
      Ipulse2[i1+1] += Ipulse[j];
      if(i2 < NrBins2)
 Ipulse2[i2] += Ipulse[j]*(x2-i2);
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR rebinPulse: Error in rebinning function.");
      return 0;
    }
  }
  return 1;
}
int continuous_shift(datafile_definition fin, datafile_definition *fout, int shift, int circularShift, char *output_name, int oformat, int argc, char **argv, verbose_definition verbose, int verbose2)
{
  int i;
  long nout, p, f, n;
  float *Ipulse, *Ifirst, *Ilast;
  verbose_definition verbose_counters_verbose2;
  if(fin.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR continuous_shift: Frequency channels need to be uniformely separated.");
    return 0;
  }
  copyVerboseState(verbose, &verbose_counters_verbose2);
  if(verbose2)
    verbose_counters_verbose2.nocounters = 0;
  else
    verbose_counters_verbose2.nocounters = 1;
  cleanPSRData(fout, verbose);
  copy_params_PSRData(fin, fout, verbose);
  fout->NrSubints = fin.NrSubints-1;
  if(circularShift != 0)
    (fout->NrSubints)++;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Continuous shift: shift=%d bins circularShift=%d.\n", shift, circularShift);
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Output data will contain %ld pulses.\n", fout->NrSubints);
  }
  if(fout->NrFreqChan > 1 && fout->isDeDisp == 0) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING continuous_shift: You might want to dedisperse the data first.");
  }
  if(fout->NrSubints == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR continuous_shift: If there is only one pulse circular shifting must be used.");
    return 0;
  }
  if(oformat == MEMORY_format) {
    fout->format = oformat;
    fout->data = (float *)malloc(fout->NrPols*fout->NrBins*fout->NrSubints*fout->NrFreqChan*sizeof(float));
    if(fout->data == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR continuous_shift: Cannot allocate memory.");
      return 0;
    }
  }else {
    if(!openPSRData(fout, output_name, oformat, 1, 0, 0, verbose_counters_verbose2))
      return 0;
    int cmdOnly = 0;
    if(!writeHeaderPSRData(fout, argc, argv, cmdOnly, NULL, verbose))
      return 0;
  }
  Ipulse = (float *)malloc(fout->NrPols*fout->NrBins*sizeof(float));
  Ifirst = (float *)malloc(fout->NrPols*fout->NrBins*sizeof(float));
  Ilast = (float *)malloc(fout->NrPols*fout->NrBins*sizeof(float));
  if(Ipulse == NULL || Ifirst == NULL || Ilast == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR continuous_shift: Cannot allocate memory.");
    return 0;
  }
  if(shift < 0)
    shift += fout->NrBins;
  if(shift >= fout->NrBins || shift < 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR continuous_shift: Shift is not valid.");
    return 0;
  }
  for(p = 0; p < fin.NrPols; p++) {
    for(n = 0; n < fin.NrSubints; n++) {
      if(shift >= 0)
 nout = n;
      else
 nout = n-1;
      for(f = 0; f < fin.NrFreqChan; f++) {
 if(readPulsePSRData(&fin, n, p, f, 0, fin.NrBins, Ipulse, verbose) != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR continuous_shift: Read error.");
   return 0;
 }
 if(n == 0) {
   if(shift == 0 && circularShift != 0) {
     if(verbose2) printf("\ncontinuous_shift: Write out whole pulse %ld\n", n);
     if(shift >= 0) {
       if(writePulsePSRData(fout, nout, p, f, shift, fin.NrBins-shift, Ipulse, verbose) != 1) {
  fflush(stdout);
  printerror(verbose.debug, "ERROR continuous_shift: Write error.");
  return 0;
       }
       if(writePulsePSRData(fout, nout+1, p, f, 0, shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
  fflush(stdout);
  printerror(verbose.debug, "ERROR continuous_shift: Write error.");
  return 0;
       }
     }else {
       if(writePulsePSRData(fout, nout, p, f, fin.NrBins-shift, shift, Ipulse, verbose) != 1) {
  fflush(stdout);
  printerror(verbose.debug, "ERROR continuous_shift: Write error.");
  return 0;
       }
       if(writePulsePSRData(fout, nout+1, p, f, 0, fin.NrBins-shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
  fflush(stdout);
  printerror(verbose.debug, "ERROR continuous_shift: Write error.");
  return 0;
       }
     }
   }
   if(shift > 0) {
     if(circularShift != 0) {
       if(verbose2) printf("\ncontinuous_shift: Write %d bins of last pulse\n", shift);
       if(readPulsePSRData(&fin, fin.NrSubints-1, p, f, 0, fin.NrBins, Ilast, verbose) != 1) {
  fflush(stdout);
  printerror(verbose.debug, "ERROR continuous_shift: Read error.");
  return 0;
       }
       if(writePulsePSRData(fout, nout, p, f, 0, shift, Ilast+fout->NrBins-shift, verbose) != 1) {
  fflush(stdout);
  printerror(verbose.debug, "ERROR continuous_shift: Write error.");
  return 0;
       }
       if(n != fin.NrSubints-1) {
  if(verbose2) printf("\ncontinuous_shift: Write out whole pulse %ld\n", n);
  if(shift >= 0) {
    if(writePulsePSRData(fout, nout, p, f, shift, fin.NrBins-shift, Ipulse, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR continuous_shift: Write error.");
      return 0;
    }
    if(writePulsePSRData(fout, nout+1, p, f, 0, shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR continuous_shift: Write error.");
      return 0;
    }
  }else {
    if(writePulsePSRData(fout, nout, p, f, fin.NrBins-shift, shift, Ipulse, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR continuous_shift: Write error.");
      return 0;
    }
    if(writePulsePSRData(fout, nout+1, p, f, 0, fin.NrBins-shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR continuous_shift: Write error.");
      return 0;
    }
  }
       }
     }else {
       if(verbose2) printf("\ncontinuous_shift: Write %d bins of pulse %ld (freq = %ld)\n", shift, n, f);
       if(writePulsePSRData(fout, nout, p, f, 0, shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
  fflush(stdout);
  printerror(verbose.debug, "ERROR continuous_shift: Write error.");
  return 0;
       }
     }
   }else if(shift < 0) {
     if(verbose2) printf("\ncontinuous_shift: Write %ld bins of pulse %ld (freq = %ld)\n", fout->NrBins+shift, n, f);
     if(writePulsePSRData(fout, nout, p, f, 0, fout->NrBins+shift, Ipulse-shift, verbose) != 1) {
       fflush(stdout);
       printerror(verbose.debug, "ERROR continuous_shift: Write error.");
       return 0;
     }
   }
 }
 if(n == fin.NrSubints-1) {
   if(shift >= 0) {
     if(circularShift == 0 && f == 0)
       nout -= 1;
     if(verbose2) printf("\ncontinuous_shift: Write %ld bins of pulse %ld\n", fout->NrBins-shift, n);
     if(writePulsePSRData(fout, nout, p, f, shift, fout->NrBins-shift, Ipulse, verbose) != 1) {
       fflush(stdout);
       printerror(verbose.debug, "ERROR continuous_shift: Write error.");
       return 0;
     }
   }else if(shift < 0) {
     if(circularShift != 0) {
       if(n != 0) {
  if(verbose2) printf("\ncontinuous_shift: continuous_shift: Write out whole pulse %ld (freq = %ld)\n", n, f);
  if(shift >= 0) {
    if(writePulsePSRData(fout, nout, p, f, shift, fin.NrBins-shift, Ipulse, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR continuous_shift: Write error.");
      return 0;
    }
    if(writePulsePSRData(fout, nout+1, p, f, 0, shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR continuous_shift: Write error.");
      return 0;
    }
  }else {
    if(writePulsePSRData(fout, nout, p, f, fin.NrBins-shift, shift, Ipulse, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR continuous_shift: Write error.");
      return 0;
    }
    if(writePulsePSRData(fout, nout+1, p, f, 0, fin.NrBins-shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR continuous_shift: Write error.");
      return 0;
    }
  }
       }
       if(verbose2) printf("\ncontinuous_shift: Write %d bins of first pulse\n", -shift);
       if(readPulsePSRData(&fin, 0, p, f, 0, fin.NrBins, Ifirst, verbose) != 1) {
  fflush(stdout);
  printerror(verbose.debug, "ERROR continuous_shift: Read error.");
  return 0;
       }
       if(writePulsePSRData(fout, nout+1, p, f, fin.NrBins-shift, -shift, Ipulse, verbose) != 1) {
  fflush(stdout);
  printerror(verbose.debug, "ERROR continuous_shift: Write error.");
  return 0;
       }
     }else {
       if(verbose2) printf("\ncontinuous_shift: Write %d bins of pulse %ld\n", -shift, n);
       if(writePulsePSRData(fout, nout, p, f, fin.NrBins-shift, -shift, Ipulse, verbose) != 1) {
  fflush(stdout);
  printerror(verbose.debug, "ERROR continuous_shift: Write error.");
  return 0;
       }
     }
   }
 }
 if(n != fin.NrSubints-1 && n != 0) {
   if(n == 1)
     if(verbose2) printf("\ncontinuous_shift: Write %ld pulses\n", fout->NrSubints-2);
   if(shift >= 0) {
     if(circularShift == 0 && f == 0)
       nout -= 1;
     if(writePulsePSRData(fout, nout, p, f, shift, fin.NrBins-shift, Ipulse, verbose) != 1) {
       fflush(stdout);
       printerror(verbose.debug, "ERROR continuous_shift: Write error (pulse=%ld pol=%ld freq=%ld, start=%d, length=%ld).", nout, p, f, shift, fin.NrBins-shift);
       return 0;
     }
     if(writePulsePSRData(fout, nout+1, p, f, 0, shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
       fflush(stdout);
       printerror(verbose.debug, "ERROR continuous_shift: Write error (pulse=%ld pol=%ld freq=%ld, start=%d, length=%d).", nout, p, f, 0, shift);
       return 0;
     }
   }else {
     if(writePulsePSRData(fout, nout, p, f, fin.NrBins-shift, shift, Ipulse, verbose) != 1) {
       fflush(stdout);
       printerror(verbose.debug, "ERROR continuous_shift: Write error (pulse=%ld pol=%ld freq=%ld, start=%ld, length=%d).", nout, p, f, fin.NrBins-shift, shift);
       return 0;
     }
     if(writePulsePSRData(fout, nout+1, p, f, 0, fin.NrBins-shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
       fflush(stdout);
       printerror(verbose.debug, "ERROR continuous_shift: Write error (pulse=%ld pol=%ld freq=%ld, start=%d, length=%ld).", nout, p, f, 0, fin.NrBins-shift);
       return 0;
     }
   }
 }
 if(verbose2 && verbose.nocounters == 0) {
   printf("continuous_shift: Writing file: %.1f%%     \r", (100.0*(n+1))/(float)(fout->NrSubints-2));
   fflush(stdout);
 }
      }
    }
  }
  free(Ipulse);
  free(Ifirst);
  free(Ilast);
  if(verbose2 && verbose.nocounters == 0)
    printf("\n");
  return 1;
}
int get_period(datafile_definition datafile, long subint, double *period, verbose_definition verbose)
{
  *period = 0;
  if(datafile.isFolded != 1) {
    printwarning(verbose.debug, "WARNING get_period: Data does not appear to be folded");
    return 1;
  }
  if(datafile.foldMode != FOLDMODE_FIXEDPERIOD) {
    printerror(verbose.debug, "ERROR get_period: Unknown folding mode");
    return 2;
  }
  *period = datafile.fixedPeriod;
  return 0;
}
double get_tsamp(datafile_definition datafile, long subint, verbose_definition verbose)
{
  if(datafile.tsampMode == TSAMPMODE_LONGITUDELIST) {
    printerror(verbose.debug, "ERROR get_tsamp: Data is not defined to have a regular sampling interval.");
    exit(0);
  }
  if(datafile.tsampMode != TSAMPMODE_FIXEDTSAMP) {
    printerror(verbose.debug, "ERROR get_tsamp: Unknown sampling mode");
    exit(0);
  }
  return datafile.fixedtsamp;
}
double get_pulse_longitude(datafile_definition datafile, long subint, long binnr, verbose_definition verbose)
{
  if(datafile.tsampMode == TSAMPMODE_LONGITUDELIST) {
    if(datafile.tsamp_list == NULL) {
      printerror(verbose.debug, "ERROR get_pulse_longitude: Longitudes appear to be undifined in the data");
      exit(0);
    }
    return datafile.tsamp_list[binnr];
  }else if(datafile.tsampMode == TSAMPMODE_FIXEDTSAMP) {
    if(datafile.fixedtsamp == 0) {
      printerror(verbose.debug, "ERROR get_pulse_longitude (%s): Sampling time appears to be set to zero.", datafile.filename);
      exit(0);
    }
    double longitude;
    double period;
    int ret;
    ret = get_period(datafile, subint, &period, verbose);
    if(ret != 0) {
      printerror(verbose.debug, "ERROR get_pulse_longitude (%s): Cannot obtain period", datafile.filename);
      exit(0);
    }
    longitude = binnr;
    longitude *= datafile.fixedtsamp;
    longitude /= period;
    return 360.0*longitude;
  }else {
    printerror(verbose.debug, "ERROR get_pulse_longitude: Unknown sampling mode");
    exit(0);
  }
  return datafile.fixedtsamp;
}
int convert_to_fixed_tsamp(datafile_definition *datafile, verbose_definition verbose)
{
  int ret;
  long binnr;
  double period, delta, delta_exp, diff;
  if(verbose.debug)
    printf("Entering convert_to_fixed_tsamp()\n");
  if(datafile->tsampMode == TSAMPMODE_FIXEDTSAMP) {
    printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: Data already has a fixed sampling time.");
    return 0;
  }
  ret = get_period(*datafile, 0, &period, verbose);
  if(ret == 2) {
    printerror(verbose.debug, "ERROR convert_to_fixed_tsamp (%s): Cannot obtain period", datafile->filename);
    return 0;
  }
  if(period <= 0.0) {
    printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: Period is not set, command is ignored.");
    return 0;
  }
  if(datafile->tsamp_list == NULL) {
    printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: The tsamp_list does not appear to be initialised, command is ignored.");
    return 0;
  }
  if(datafile->NrBins <= 1) {
    printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: More than 1 bin is required, command is ignored.");
    return 0;
  }
  delta_exp = get_pulse_longitude(*datafile, 0, 1, verbose);
  delta_exp -= get_pulse_longitude(*datafile, 0, 0, verbose);
  if(delta_exp <= 0.0) {
    printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: Sampling time seems to be negative or zero, command is ignored.");
    return 0;
  }
  if(datafile->NrBins > 2) {
    for(binnr = 1; binnr < datafile->NrBins-1; binnr++) {
      delta = get_pulse_longitude(*datafile, 0, binnr+1, verbose);
      delta -= get_pulse_longitude(*datafile, 0, binnr, verbose);
      diff = fabs(delta-delta_exp)/delta_exp;
      if(diff > 1.0001 && diff < 0.9999) {
 printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: Sampling time does not appear to be regular.");
 return 0;
      }
    }
  }
  free(datafile->tsamp_list);
  datafile->tsamp_list = NULL;
  datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
  datafile->fixedtsamp = period*delta_exp/360.0;
  if(verbose.debug) {
    printf("  convert_to_fixed_tsamp: new sampling time is %lf sec\n", datafile->fixedtsamp);
    printf("Exiting convert_to_fixed_tsamp()\n");
  }
  return 1;
}
double get_tsub(datafile_definition datafile, long subint, verbose_definition verbose)
{
  if(datafile.tsubMode == TSUBMODE_FIXEDTSUB) {
    return datafile.tsub_list[0];
  }else if(datafile.tsubMode == TSUBMODE_TSUBLIST) {
    return datafile.tsub_list[subint];
  }else {
    printerror(verbose.debug, "ERROR get_tsub: Unknown subint duration mode", datafile.tsubMode);
    if(verbose.debug) {
      printerror(verbose.debug, "ERROR get_tsub: Subint duration mode is set to %d", datafile.tsubMode);
    }
    exit(0);
  }
}
double get_tobs(datafile_definition datafile, verbose_definition verbose)
{
  long i;
  double tobs;
  if(datafile.gentype == GENTYPE_LRFS || datafile.gentype == GENTYPE_2DFS || datafile.gentype == GENTYPE_S2DFSP3 || datafile.gentype == GENTYPE_S2DFSP2 || datafile.gentype == GENTYPE_P3FOLD || datafile.gentype == GENTYPE_LRCC || datafile.gentype == GENTYPE_RMMAP || datafile.gentype == GENTYPE_PADIST) {
    return get_tsub(datafile, 0, verbose);
  }else if(datafile.gentype == GENTYPE_RECEIVERMODEL || datafile.gentype == GENTYPE_RECEIVERMODEL2) {
    return 0;
  }else {
    tobs = 0;
    for(i = 0; i < datafile.NrSubints; i++)
      tobs += get_tsub(datafile, i, verbose);
    return tobs;
  }
}
long double get_mjd_subint(datafile_definition datafile, long subint, verbose_definition verbose)
{
  int i;
  long double mjd = datafile.mjd_start;
  long double offset;
  offset = 0;
  if(subint > 0) {
    for(i = 0; i < subint; i++)
      offset += get_tsub(datafile, i, verbose);
  }
  offset += 0.5*get_tsub(datafile, subint, verbose);
  mjd += offset/(3600.0*24.0);
  return mjd;
}
int get_channelbandwidth(datafile_definition datafile, double *channelbw, verbose_definition verbose)
{
  if(datafile.isTransposed == 0)
    *channelbw = get_bandwidth(datafile, verbose)/(double)datafile.NrFreqChan;
  else
    *channelbw = get_bandwidth(datafile, verbose)/(double)datafile.NrSubints;
  return 1;
}
void set_channelbandwidth(datafile_definition *datafile, double channelbw, verbose_definition verbose)
{
  datafile->bandwidth = channelbw*datafile->NrFreqChan;
}
double get_bandwidth(datafile_definition datafile, verbose_definition verbose)
{
  return datafile.bandwidth;
}
int set_bandwidth(datafile_definition *datafile, double bw, verbose_definition verbose)
{
  datafile->bandwidth = bw;
  return 1;
}
double get_centre_frequency(datafile_definition datafile, verbose_definition verbose)
{
  return datafile.centrefreq;
}
void set_centre_frequency(datafile_definition *datafile, double freq, verbose_definition verbose)
{
  datafile->centrefreq = freq;
}
double get_nonweighted_channel_freq(datafile_definition psrdata, long channel, verbose_definition verbose)
{
  double freq_bottom, cfreq;
  if(psrdata.isTransposed == 0) {
    if(channel < 0 || channel >= psrdata.NrFreqChan) {
      printerror(verbose.debug, "ERROR get_nonweighted_channel_freq: channel out of range (channel=%ld/%ld).", channel, psrdata.NrFreqChan);
      exit(0);
    }
  }else {
    if(channel < 0 || channel >= psrdata.NrSubints) {
      printerror(verbose.debug, "ERROR get_nonweighted_channel_freq: channel out of range (channel=%ld/%ld).", channel, psrdata.NrFreqChan);
      exit(0);
    }
  }
  double chanbw, cfreq_total;
  if(get_channelbandwidth(psrdata, &chanbw, verbose) == 0) {
    printerror(verbose.debug, "ERROR get_nonweighted_channel_freq (%s): Cannot obtain channel bandwidth.", psrdata.filename);
    exit(0);
  }
  cfreq_total = get_centre_frequency(psrdata, verbose);
  if(psrdata.isTransposed == 0) {
    freq_bottom = cfreq_total - 0.5*(double)psrdata.NrFreqChan*chanbw;
    cfreq = freq_bottom + chanbw*(channel+0.5);
  }else {
    freq_bottom = cfreq_total - 0.5*(double)psrdata.NrSubints*chanbw;
    cfreq = freq_bottom + chanbw*(channel+0.5);
  }
  return cfreq;
}
double get_weighted_channel_freq(datafile_definition psrdata, long subint, long channel, verbose_definition verbose)
{
  if(psrdata.freqMode != FREQMODE_FREQTABLE) {
    return get_nonweighted_channel_freq(psrdata, channel, verbose);
  }
  if(subint < 0 || channel < 0 || subint >= psrdata.NrSubints || channel >= psrdata.NrFreqChan) {
      printerror(verbose.debug, "ERROR get_weighted_channel_freq: subint/channel out of range (subint=%ld/%ld, channel=%ld/%ld).", subint, psrdata.NrSubints, channel, psrdata.NrFreqChan);
      exit(0);
  }
  return psrdata.freqlabel_list[subint*psrdata.NrFreqChan+channel];
}
int set_weighted_channel_freq(datafile_definition *psrdata, long subint, long channel, double freq, verbose_definition verbose)
{
  if(psrdata->freqMode != FREQMODE_FREQTABLE) {
    printerror(verbose.debug, "ERROR set_weighted_channel_freq: freqMode is set to %d, expected %d.\n", psrdata->freqMode, FREQMODE_FREQTABLE);
    return 0;
  }
  if(subint < 0 || channel < 0 || subint >= psrdata->NrSubints || channel >= psrdata->NrFreqChan) {
      printerror(verbose.debug, "ERROR set_weighted_channel_freq: subint/channel out of range (subint=%ld/%ld, channel=%ld/%ld).", subint, psrdata->NrSubints, channel, psrdata->NrFreqChan);
      return 0;
  }
  psrdata->freqlabel_list[subint*psrdata->NrFreqChan+channel] = freq;
  return 1;
}
char * get_history_notes_last(datafile_definition *psrdata)
{
  datafile_history_entry_definition *source_hist;
  char *ret;
  ret = NULL;
  source_hist = &(psrdata->history);
  while(source_hist->nextEntry != NULL) {
    source_hist = source_hist->nextEntry;
  }
  if(source_hist != NULL) {
    ret = source_hist->notes;
  }
  return ret;
}
int data_parang(datafile_definition data, long subintnr, double *parang, verbose_definition verbose)
{
  if(fabs(data.telescope_X) < 1e-6 && fabs(data.telescope_Y) < 1e-6 && fabs(data.telescope_Z) < 1e-6) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR data_parang: Position of telescope appears to be undefined.");
    return 0;
  }
  if(fabs(data.ra) < 1e-6 && fabs(data.dec) < 1e-6) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR data_parang: Position of source appears to be undefined.");
    return 0;
  }
  double telescope_long = observatory_long_geodetic(data);
  double telescope_lat = observatory_lat_geodetic(data);
  long double mjd;
  if(subintnr < 0) {
    mjd = data.mjd_start + 0.5*get_tobs(data, verbose)/(double)(3600.0*24.0);
  }else {
    mjd = get_mjd_subint(data, subintnr, verbose);
  }
  *parang = calc_parang(telescope_long, telescope_lat, data.ra, data.dec, mjd, 1);
  return 1;
}
int check_baseline_subtracted(datafile_definition data, verbose_definition verbose)
{
  long f, n, b;
  float sample, miny, maxy;
  miny = maxy = 0;
  if(data.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR check_baseline_subtracted: Data should be loaded into memory.");
    return 0;
  }
  for(f = 0; f < data.NrFreqChan; f++) {
    for(n = 0; n < data.NrSubints; n++) {
      for(b = 0; b < data.NrBins; b++) {
 if(readPulsePSRData(&data, n, 0, f, b, 1, &sample, verbose) != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR check_baseline_subtracted: Cannot read data.");
   exit(-1);
 }
 if(b == 0) {
   miny = sample;
   maxy = sample;
 }else {
   if(sample < miny)
     miny = sample;
   if(sample > maxy)
     maxy = sample;
 }
      }
      if(maxy < 0 || miny > 0) {
 return 0;
      }
    }
  }
  return 1;
}

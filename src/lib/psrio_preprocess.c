/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <math.h>
#include <string.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_randist.h"
#include "psrsalsa.h"




int preprocess_rebin(datafile_definition original, datafile_definition *clone, long NrBins, verbose_definition verbose)
{
  long p, f, n;
  int i;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Rebinning data to %ld bins\n", NrBins);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rebin: Rebinning only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rebin: Cannot handle PA data.");
    return 0;
  }
  if(NrBins > original.NrBins) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rebin: Cannot rebin to a larger amount of bins.");
    return 0;
  }
  if(original.NrBins % NrBins != 0) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_rebin: Rebinning from %ld to %ld bins implies that separate bins are not entirely independent.", original.NrBins, NrBins);
  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->fixedtsamp *= original.NrBins/(double)NrBins;
  clone->NrBins = NrBins;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  if(clone->data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rebin: Memory allocation error.");
    return 0;
  }


  for(p = 0; p < clone->NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(n = 0; n < clone->NrSubints; n++) {
 if(!rebinPulse(&(original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))]), original.NrBins, &(clone->data[clone->NrBins*(p+clone->NrPols*(f+n*clone->NrFreqChan))]), clone->NrBins, 1, verbose)) {
   return 0;
 }
        if(verbose.verbose && verbose.nocounters == 0) {
   long doprint;
   doprint = 1;
   if(original.NrFreqChan > 4 && n != 0)
     doprint = 0;
   if(doprint) {
     for(i = 0; i < verbose.indent; i++)
       printf(" ");
     printf("  %.1f%%     \r", (100.0*((p*clone->NrFreqChan+f)*clone->NrSubints+n))/(float)(clone->NrSubints*clone->NrFreqChan*clone->NrPols));
     fflush(stdout);
   }
 }
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done                               \n");
  }
  return 1;
}






int preprocess_fftshift(datafile_definition original, float shiftPhase, verbose_definition verbose)
{
  long p, f, n;
  int i;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Rotating data by %f phase\n", shiftPhase);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_fftshift: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_fftshift: Cannot handle PA data.");
    return 0;
  }


  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
 if(rotateSinglepulse(&(original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))]), original.NrBins, shiftPhase*original.NrBins, verbose) == 0)
   return 0;
        if(verbose.verbose && verbose.nocounters == 0) {
   long doprint;
   doprint = 1;
   if(original.NrFreqChan > 4 && n != 0)
     doprint = 0;
   if(doprint) {
     for(i = 0; i < verbose.indent; i++)
       printf(" ");
     printf("  %.1f%%     \r", (100.0*((p*original.NrFreqChan+f)*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan*original.NrPols));
     fflush(stdout);
   }
 }
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(verbose.verbose) printf("  done                            \n");
  }
  return 1;
}






int preprocess_channelselect(datafile_definition original, datafile_definition *clone, long chanelnr, verbose_definition verbose)
{
  long p, f, n;
  int i;
  float *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Selecting frequency channel %ld\n", chanelnr);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_chanelselect: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_chanelselect: Cannot handle PA data.");
    return 0;
  }
  if(chanelnr < 0 || chanelnr >= original.NrFreqChan) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_chanelselect: Invalid frequency chanel number.");
    return 0;
  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->NrFreqChan = 1;

  clone->uniform_bw = get_bw(original, verbose);
  clone->uniform_freq_cent = get_channel_freq(original, chanelnr, verbose);
  clone->freqMode = FREQMODE_UNIFORM;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_chanelselect: Memory allocation error.");
    return 0;
  }



  for(p = 0; p < clone->NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < clone->NrSubints; n++) {
 if(f == chanelnr) {
   if(readPulsePSRData(&original, n, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR preprocess_chanelselect: Error reading pulse.");
     return 0;
   }
   if(writePulsePSRData(clone, n, p, 0, 0, clone->NrBins, pulse, verbose) != 1) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR preprocess_chanelselect: Error writing pulse.");
     return 0;
   }
 }
      }
    }
  }
  free(pulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done\n");
  }
  return 1;
}
int preprocess_polselect(datafile_definition original, datafile_definition *clone, long polnr, verbose_definition verbose)
{
  long p, f, n;
  int i;
  float *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Selecting polarization channel %ld\n", polnr);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_polselect: only works if data is loaded into memory.");
    return 0;
  }
  if(polnr < 0 || polnr >= original.NrPols) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_polselect: Invalid polarization chanel number.");
    return 0;
  }
  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->NrPols = 1;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_polselect: Memory allocation error.");
    return 0;
  }
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(n = 0; n < clone->NrSubints; n++) {
 if(p == polnr) {
   if(readPulsePSRData(&original, n, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR preprocess_polselect: Error reading pulse.");
     return 0;
   }
   if(writePulsePSRData(clone, n, 0, f, 0, clone->NrBins, pulse, verbose) != 1) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR preprocess_polselect: Error writing pulse.");
     return 0;
   }
 }
      }
    }
  }
  free(pulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done              \n");
  }
  return 1;
}
int preprocess_pulsesselect(datafile_definition original, datafile_definition *clone, long nskip, long nread, verbose_definition verbose)
{
  long p, f, n, i;
  float *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Selecting pulses %ld-%ld\n", nskip, nread+nskip-1);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_pulsesselect: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_pulsesselect: Cannot handle PA data.");
    return 0;
  }
  if(nskip < 0 || nskip >= original.NrSubints) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_pulsesselect: Invalid nskip.");
    return 0;
  }
  if(nread < 0 || nskip+nread > original.NrSubints) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_pulsesselect: Invalid nread.");
    return 0;
  }
  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->NrSubints = nread;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_pulsesselect: Memory allocation error.");
    return 0;
  }
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(n = nskip; n < nskip+nread; n++) {
 if(readPulsePSRData(&original, n, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR preprocess_pulsesselect: Error reading pulse.");
   return 0;
 }
 if(writePulsePSRData(clone, n-nskip, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR preprocess_pulsesselect: Error writing pulse.");
   return 0;
 }
      }
    }
  }
  free(pulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done            \n");
  }
  return 1;
}
int preprocess_blocksize(datafile_definition original, datafile_definition *clone, int blocksize, verbose_definition verbose)
{
  int i;
  long nread;
  verbose_definition verbose2;
  nread = original.NrSubints/blocksize;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Selecting %ld blocks of %d subints\n", nread, blocksize);
  }
  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;
  return preprocess_pulsesselect(original, clone, 0, nread*blocksize, verbose2);
}
int preprocess_transposeRawFBdata(datafile_definition original, datafile_definition *clone, verbose_definition verbose)
{
  long p, f, n;
  float *pulse;
  int i;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Transposing data\n");
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_transposeRawFBdata: only works if data is loaded into memory.");
    return 0;
  }
  if(original.isTransposed) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_transposeRawFBdata: Data already appears to be transposed.");
    return 0;
  }
  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_transposeRawFBdata: Memory allocation error.");
    return 0;
  }
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(n = 0; n < clone->NrSubints; n++) {
 if(readPulsePSRData(&original, n, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR preprocess_transposeRawFBdata: Error reading pulse.");
   return 0;
 }
 memcpy(&(clone->data[original.NrBins*(n+original.NrSubints*(f+p*original.NrFreqChan))]), pulse, sizeof(float)*(clone->NrBins));
      }
    }
  }
  clone->isTransposed = 1;
  free(pulse);
  return 1;
}
int preprocess_addsuccessivepulses(datafile_definition original, datafile_definition *clone, long nrpulses, int complete, verbose_definition verbose)
{
  long p, f, n, n2, b;
  float *pulse, *addedpulse;
  int i, use_depar;
  datafile_definition clone_depar;
  verbose_definition verbose2;
  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;
  use_depar = 0;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(nrpulses >= 0)
      printf("Add each %ld succesive subints\n", nrpulses);
    else
      printf("Write out each subint %ld times\n", -nrpulses);
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Cannot handle position angle data when adding subints. Add subints using Stokes parameter data.");
    return 0;
  }
  if(nrpulses > original.NrSubints || nrpulses == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Invalid number of subints to add.");
    return 0;
  }
  if(original.NrPols == 4 && original.isDePar == -1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_addsuccessivepulses: Parallactic angle correction state is unknown, no correction will be done.");
  }else if(original.NrPols == 4 && original.isDePar == 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Correcting for parallactic angle changes first.\n");
    }
    if(preprocess_corrParAng(&original, &clone_depar, 0, verbose2) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Parallactic angle correction failed.");
      return 0;
    }
    use_depar = 1;
  }
  cleanPSRData(clone, verbose);
  if(use_depar)
    copy_params_PSRData(clone_depar, clone, verbose);
  else
    copy_params_PSRData(original, clone, verbose);
  if(nrpulses > 0) {
    clone->NrSubints = original.NrSubints/nrpulses;
    if(complete == 0) {
      if(clone->NrSubints * nrpulses != original.NrSubints) {
 printwarning(verbose.debug, "WARNING preprocess_addsuccessivepulses: Last subint has a different duration.");
 clone->NrSubints += 1;
      }
    }
  }else {
    clone->NrSubints = original.NrSubints*(-nrpulses);
  }
  clone->tsubMode = TSUBMODE_TSUBLIST;
  if(clone->tsub_list != NULL)
    free(clone->tsub_list);
  clone->tsub_list = (double *)malloc(clone->NrSubints * sizeof(double));
  if(clone->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Memory allocation error");
    return 0;
  }
  clone->format = MEMORY_format;
  if(clone->gentype == GENTYPE_PULSESTACK) {
    clone->gentype = GENTYPE_SUBINTEGRATIONS;
  }
  if(clone->gentype == GENTYPE_SUBINTEGRATIONS && clone->NrSubints == 1) {
    clone->gentype = GENTYPE_PROFILE;
  }
  if(clone->gentype != GENTYPE_PROFILE && clone->gentype != GENTYPE_SUBINTEGRATIONS && clone->gentype != GENTYPE_PULSESTACK && clone->gentype != GENTYPE_UNDEFINED && clone->gentype != GENTYPE_POLNCAL) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_addsuccessivepulses: Unsure about adding subints for a %s file. Setting gentype to undefined.", returnGenType_str(clone->gentype));
    clone->gentype = GENTYPE_UNDEFINED;
  }
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  addedpulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL || addedpulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Memory allocation error.");
    return 0;
  }
  double curtsub;
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      curtsub = 0;
      for(n = 0; n < clone->NrSubints; n++) {
 if(nrpulses > 0) {
   for(b = 0; b < original.NrBins; b++)
     addedpulse[b] = 0;
   for(n2 = 0; n2 < nrpulses; n2++) {
     if(n*nrpulses+n2 < original.NrSubints) {
       if(use_depar) {
  if(readPulsePSRData(&clone_depar, n*nrpulses+n2, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Error reading pulse.");
    return 0;
  }
  curtsub += get_tsub(clone_depar, n*nrpulses+n2, verbose);
       }else {
  if(readPulsePSRData(&original, n*nrpulses+n2, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Error reading pulse.");
    return 0;
  }
  curtsub += get_tsub(original, n*nrpulses+n2, verbose);
       }
       for(b = 0; b < original.NrBins; b++)
  addedpulse[b] += pulse[b];
     }
   }
 }else {
   if(use_depar) {
     if(readPulsePSRData(&clone_depar, n/(-nrpulses), p, f, 0, clone->NrBins, addedpulse, verbose) != 1) {
       fflush(stdout);
       printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Error reading pulse.");
       return 0;
     }
     curtsub += get_tsub(clone_depar, n/(-nrpulses), verbose)/(double)(-nrpulses);
   }else {
     if(readPulsePSRData(&original, n/(-nrpulses), p, f, 0, clone->NrBins, addedpulse, verbose) != 1) {
       fflush(stdout);
       printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Error reading pulse.");
       return 0;
     }
     curtsub += get_tsub(original, n/(-nrpulses), verbose)/(double)(-nrpulses);
   }
 }
 clone->tsub_list[n] = curtsub;
 if(writePulsePSRData(clone, n, p, f, 0, clone->NrBins, addedpulse, verbose) != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Error writing pulse.");
   return 0;
 }
 curtsub = 0;
        if(verbose.verbose && verbose.nocounters == 0) {
   long doprint;
   doprint = 1;
   if(clone->NrFreqChan > 4 && n != 0)
     doprint = 0;
   if(doprint) {
     for(i = 0; i < verbose.indent; i++)
       printf(" ");
     printf("  %.1f%%     \r", (100.0*((p*clone->NrFreqChan+f)*clone->NrSubints+n))/(float)(clone->NrSubints*clone->NrFreqChan*clone->NrPols));
     fflush(stdout);
   }
 }
      }
    }
  }
  free(pulse);
  free(addedpulse);
  if(use_depar) {
    closePSRData(&clone_depar, verbose);
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(verbose.verbose) printf("  done                            \n");
  }
  return 1;
}
int preprocess_addsuccessiveFreqChans(datafile_definition original, datafile_definition *clone, long nrfreq, int *fzapMask, verbose_definition verbose)
{
  long p, f, n, n2, b;
  float *pulse, *addedpulse;
  int i, ok;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(nrfreq >= 0)
      printf("Add %ld succesive frequency channels\n", nrfreq);
    else
      printf("Write out each frequency channel %ld times\n", -nrfreq);
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Cannot handle position angle data when adding frequency channels. Add subints using Stokes parameter data.");
    return 0;
  }
  if(nrfreq > original.NrFreqChan || nrfreq == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Invalid number of frequency channels to add.");
    return 0;
  }
  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  if(nrfreq > 0) {
    clone->NrFreqChan = original.NrFreqChan/nrfreq;
    if(clone->NrFreqChan * nrfreq != original.NrFreqChan) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_addsuccessiveFreqChans: Last channel will be the sum of a different number of channels compared to the others. The frequency labeling will no longer be correct.");
    }
  }else {
    clone->NrFreqChan = original.NrFreqChan*(-nrfreq);
  }
  clone->format = MEMORY_format;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  addedpulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL || addedpulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Memory allocation error.");
    return 0;
  }
  for(p = 0; p < clone->NrPols; p++) {
    for(n = 0; n < clone->NrSubints; n++) {
      for(f = 0; f < clone->NrFreqChan; f++) {
 if(nrfreq > 0) {
   for(b = 0; b < original.NrBins; b++)
     addedpulse[b] = 0;
   for(n2 = 0; n2 < nrfreq; n2++) {
     ok = 1;
     if(fzapMask != NULL) {
       if(fzapMask[f*nrfreq+n2] != 0) {
  ok = 0;
       }
     }
     if(ok) {
       if(readPulsePSRData(&original, n, p, f*nrfreq+n2, 0, clone->NrBins, pulse, verbose) != 1) {
  fflush(stdout);
  printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Error reading pulse.");
  return 0;
       }
       for(b = 0; b < original.NrBins; b++)
  addedpulse[b] += pulse[b];
     }
   }
 }else {
   ok = 1;
   if(fzapMask != NULL) {
     if(fzapMask[f/(-nrfreq)] != 0) {
       ok = 0;
     }
   }
   if(ok) {
     if(readPulsePSRData(&original, n, p, f/(-nrfreq), 0, clone->NrBins, addedpulse, verbose) != 1) {
       fflush(stdout);
       printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Error reading pulse.");
       return 0;
     }
   }else {
     for(b = 0; b < original.NrBins; b++)
       addedpulse[b] = 0;
   }
 }
 if(writePulsePSRData(clone, n, p, f, 0, clone->NrBins, addedpulse, verbose) != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Error writing pulse.");
   return 0;
 }
        if(verbose.verbose && verbose.nocounters == 0) {
   long doprint;
   doprint = 1;
   if(clone->NrFreqChan > 4 && n != 0)
     doprint = 0;
   if(doprint) {
     for(i = 0; i < verbose.indent; i++)
       printf(" ");
     printf("  %.1f%%     \r", (100.0*((p*clone->NrFreqChan+f)*clone->NrSubints+n))/(float)(clone->NrSubints*clone->NrFreqChan*clone->NrPols));
     fflush(stdout);
   }
 }
      }
    }
  }
  free(pulse);
  free(addedpulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done                            \n");
  }
  return 1;
}
int preprocess_debase(datafile_definition *original, regions_definition onpulse, verbose_definition verbose)
{
  long p, f, n, j, nrOffpulseBins;
  int i;
  float avrg, *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Subtracting baseline\n");
  }
  region_frac_to_int(&onpulse, original->NrBins, 0);
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_debase: only works if data is loaded into memory.");
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_debase: Cannot handle PA data.");
    return 0;
  }
  pulse = (float *)malloc((original->NrBins)*sizeof(float));
  if(pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_debase: Memory allocation error.");
    return 0;
  }
  for(p = 0; p < original->NrPols; p++) {
    for(f = 0; f < original->NrFreqChan; f++) {
      for(n = 0; n < original->NrSubints; n++) {
 if(readPulsePSRData(original, n, p, f, 0, original->NrBins, pulse, verbose) != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR preprocess_debase: Error reading data.");
   return 0;
 }
 nrOffpulseBins = 0;
 avrg = 0;
 for(j = 0; j < original->NrBins; j++) {
   if(checkRegions(j, &onpulse, 0, verbose) == 0 || onpulse.nrRegions == 0) {
     avrg += pulse[j];
     nrOffpulseBins++;
   }
 }
 if(nrOffpulseBins > 0) {
   avrg /= (float)nrOffpulseBins;
   for(j = 0; j < original->NrBins; j++) {
     pulse[j] -= avrg;
   }
   if(writePulsePSRData(original, n, p, f, 0, original->NrBins, pulse, verbose) != 1) {
     fflush(stdout);
     printerror(verbose.debug, "ERROR preprocess_debase: Error writing data.");
     return 0;
   }
 }
        if(verbose.verbose && verbose.nocounters == 0) {
   long doprint;
   doprint = 1;
   if(original->NrFreqChan > 4 && n != 0)
     doprint = 0;
   if(doprint) {
     for(i = 0; i < verbose.indent; i++)
       printf(" ");
     printf("  %.1f%%     \r", (100.0*((p*original->NrFreqChan+f)*original->NrSubints+n))/(float)(original->NrSubints*original->NrFreqChan*original->NrPols));
     fflush(stdout);
   }
 }
      }
    }
  }
  original->isDebase = 1;
  free(pulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(verbose.verbose) printf("  done                            \n");
  }
  return 1;
}
static int preprocess_addNoise_random_nr_generater_initialized = 0;
static gsl_rng *preprocess_addNoiserand_num_gen;
static const gsl_rng_type *rand_num_gen_type;
int preprocess_addNoise(datafile_definition original, datafile_definition *clone, float rms, verbose_definition verbose)
{
  long p, f, n, b;
  int i;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Adding noise to data with rms %f\n", rms);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addNoise: Adding noise only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA) {
    if(original.NrPols == 1) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_addNoise: The polarization state suggests this is not necessarily a Stokes parameter, hence the error distribution is not necessarily Gaussian, which is assumed to be the case. This might be a problem.");
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_addNoise: The polarization channels appear to include L and PA, hence the error distribution is not necessarily Gaussian.");
      return 0;
    }
  }
  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  if(clone->data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addNoise: Memory allocation error.");
    return 0;
  }
  if(preprocess_addNoise_random_nr_generater_initialized == 0) {
    gsl_rng_env_setup();
    rand_num_gen_type = gsl_rng_default;
    preprocess_addNoiserand_num_gen = gsl_rng_alloc(rand_num_gen_type);
    preprocess_addNoise_random_nr_generater_initialized = 1;
  }
  for(p = 0; p < clone->NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(n = 0; n < clone->NrSubints; n++) {
 if(writePulsePSRData(clone, n, p, f, 0, clone->NrBins, &(original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))]), verbose) != 1) {
   return 0;
 }
 for(b = 0; b < clone->NrBins; b++) {
   clone->data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] += gsl_ran_gaussian(preprocess_addNoiserand_num_gen, rms);
 }
 if(verbose.verbose && verbose.nocounters == 0) {
   long doprint;
   doprint = 1;
   if(clone->NrFreqChan > 4 && n != 0)
     doprint = 0;
   if(doprint) {
     for(i = 0; i < verbose.indent; i++)
       printf(" ");
     printf("  %.1f%%     \r", (100.0*((p*clone->NrFreqChan+f)*clone->NrSubints+n))/(float)(clone->NrSubints*clone->NrFreqChan*clone->NrPols));
     fflush(stdout);
   }
 }
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(verbose.verbose) printf("done                              \n");
  }
  return 1;
}
static int preprocess_shuffle_random_nr_generater_initialized = 0;
static gsl_rng *preprocess_shufflerand_num_gen;
static const gsl_rng_type *rand_num_gen_type;
int preprocess_shuffle(datafile_definition original, datafile_definition *clone, int fixseed, verbose_definition verbose)
{
  long p, f, n;
  int i;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Shuffeling order of subints\n");
  }
  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->format = MEMORY_format;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  if(clone->data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_shuffle: Memory allocation error.");
    return 0;
  }
  if(preprocess_shuffle_random_nr_generater_initialized == 0) {
    long idnum;
    gsl_rng_env_setup();
    rand_num_gen_type = gsl_rng_default;
    preprocess_shufflerand_num_gen = gsl_rng_alloc(rand_num_gen_type);
    if(fixseed)
      idnum = 1;
    else
      randomize_idnum(&idnum);
    gsl_rng_set(preprocess_shufflerand_num_gen, idnum);
    preprocess_shuffle_random_nr_generater_initialized = 1;
  }
  long *subintlist;
  subintlist = (long *)malloc(clone->NrSubints*sizeof(long));
  if(subintlist == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_shuffle: Memory allocation error.");
    return 0;
  }
  for (n = 0; n < clone->NrSubints; n++) {
    subintlist[n] = n;
  }
  gsl_ran_shuffle (preprocess_shufflerand_num_gen, subintlist, clone->NrSubints, sizeof(long));
  for(n = 0; n < clone->NrSubints; n++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(p = 0; p < clone->NrPols; p++) {
 if(readPulsePSRData(&original, n, p, f, 0, clone->NrBins, &(clone->data[clone->NrBins*(p+clone->NrPols*(f+subintlist[n]*clone->NrFreqChan))]), verbose) != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR preprocess_shuffle: Read error.");
   return 0;
 }
 if(verbose.verbose && verbose.nocounters == 0) {
   long doprint;
   doprint = 1;
   if(clone->NrFreqChan > 4 && n != 0)
     doprint = 0;
   if(doprint) {
     for(i = 0; i < verbose.indent; i++)
       printf(" ");
     printf("  %.1f%%     \r", (100.0*((p*clone->NrFreqChan+f)*clone->NrSubints+n))/(float)(clone->NrSubints*clone->NrFreqChan*clone->NrPols));
     fflush(stdout);
   }
 }
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(verbose.verbose) printf("done                              \n");
  }
  return 1;
}
int preprocess_rotateStokes(datafile_definition *original, datafile_definition *clone, int inplace, int subint, float angle, int stokes1, int stokes2, verbose_definition verbose)
{
  int i;
  long f, n, b, p;
  float x, y;
  if(verbose.verbose) {
    for(n = 0; n < verbose.indent; n++)
      printf(" ");
    printf("Rotating ");
    switch(stokes1) {
    case 0: printf("I"); break;
    case 1: printf("Q"); break;
    case 2: printf("U"); break;
    case 3: printf("V"); break;
    default:
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_rotateStokes: Non-existing stokes parameter was defined.");
      return 0;
    }
    printf(" and ");
    switch(stokes2) {
    case 0: printf("I"); break;
    case 1: printf("Q"); break;
    case 2: printf("U"); break;
    case 3: printf("V"); break;
    default:
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_rotateStokes: Non-existing stokes parameter was defined.");
      return 0;
    }
    printf(" by %f degrees\n", angle);
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rotateStokes: Rotating Stokes parameters only works if data is loaded into memory.");
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rotateStokes: Cannot handle PA data.");
    return 0;
  }
  if(original->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rotateStokes: There should be 4 polarization channels in the data.");
    return 0;
  }
  if(original->poltype == POLTYPE_COHERENCY) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rotateStokes: Date is written as coherency parameters. Please convert into Stokes parameters first.");
    return 0;
  }else if(original->poltype == POLTYPE_UNKNOWN) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_rotateStokes: It is assumed the data are Stokes parameters.");
  }else if(original->poltype != POLTYPE_STOKES) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_rotateStokes: Cannot process this polarization type of the data.");
    return 0;
  }
  if(verbose.debug) {
    for(n = 0; n < verbose.indent; n++)
      printf(" ");
    if(original->isDebase != 1) {
      printf("WARNING preprocess_rotateStokes: Baseline variations in time/freq will be introduced if baseline is not removed\n");
    }
  }
  if(inplace == 0) {
    cleanPSRData(clone, verbose);
    copy_params_PSRData(*original, clone, verbose);
    clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
    if(clone->data == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_rotateStokes: Memory allocation error.");
      return 0;
    }
  }
  angle *= M_PI/180.0;
  for(f = 0; f < original->NrFreqChan; f++) {
    for(n = 0; n < original->NrSubints; n++) {
      if(n == subint || subint < 0 || inplace == 0) {
 for(b = 0; b < original->NrBins; b++) {
   for(p = 0; p < 4; p++) {
     if((p != stokes1 && p != stokes2 && inplace == 0) || (inplace == 0 && subint != n)) {
       clone->data[original->NrBins*(p+original->NrPols*(f+n*original->NrFreqChan))+b] = original->data[original->NrBins*(p+original->NrPols*(f+n*original->NrFreqChan))+b];
     }else if(p == stokes1) {
       x = original->data[original->NrBins*(stokes1+original->NrPols*(f+n*original->NrFreqChan))+b]*cos(angle)-original->data[original->NrBins*(stokes2+original->NrPols*(f+n*original->NrFreqChan))+b]*sin(angle);
       y = original->data[original->NrBins*(stokes1+original->NrPols*(f+n*original->NrFreqChan))+b]*sin(angle)+original->data[original->NrBins*(stokes2+original->NrPols*(f+n*original->NrFreqChan))+b]*cos(angle);
       if(inplace == 0) {
  clone->data[original->NrBins*(stokes2+original->NrPols*(f+n*original->NrFreqChan))+b] = y;
  clone->data[original->NrBins*(stokes1+original->NrPols*(f+n*original->NrFreqChan))+b] = x;
       }else {
  original->data[original->NrBins*(stokes2+original->NrPols*(f+n*original->NrFreqChan))+b] = y;
  original->data[original->NrBins*(stokes1+original->NrPols*(f+n*original->NrFreqChan))+b] = x;
       }
     }
   }
 }
      }
    }
    if(verbose.verbose && verbose.nocounters == 0) {
      long doprint;
      doprint = 1;
      if(original->NrFreqChan > 4 && n != 0)
 doprint = 0;
      if(doprint) {
 for(i = 0; i < verbose.indent; i++)
   printf(" ");
 printf("  %.1f%%     \r", (100.0*(f*original->NrSubints+n))/(float)(original->NrSubints*original->NrFreqChan));
 fflush(stdout);
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(verbose.verbose) printf("  done                              \n");
  }
  return 1;
}
int preprocess_norm(datafile_definition original, float normvalue, regions_definition *onpulse, verbose_definition verbose)
{
  int ok, first;
  long p, f, n, b, i;
  float max, fac;
  regions_definition onpulse_converted;
  if(onpulse != NULL) {
    memcpy(&onpulse_converted, onpulse, sizeof(regions_definition));
    region_frac_to_int(&onpulse_converted, original.NrBins, 0);
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Normalize data to %f (individual channels/subints)\n", normvalue);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_norm: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_norm: Cannot handle PA data.");
    return 0;
  }
  for(f = 0; f < original.NrFreqChan; f++) {
    for(n = 0; n < original.NrSubints; n++) {
      p = 0;
      max = original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))];
      first = 1;
      for(b = 0; b < original.NrBins; b++) {
 ok = 0;
 if(onpulse == NULL) {
   ok = 1;
 }else {
   if(checkRegions(b, &onpulse_converted, 0, verbose) != 0 || onpulse_converted.nrRegions == 0) {
     ok = 1;
   }
 }
 if((original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] > max || first == 1) && ok == 1) {
   max = original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b];
   first = 0;
 }
      }
      if(max != 0)
 fac = normvalue/max;
      else
 fac = 1;
      for(p = 0; p < original.NrPols; p++) {
 for(b = 0; b < original.NrBins; b++) {
   original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] *= fac;
 }
      }
      if(verbose.verbose && verbose.nocounters == 0) {
 long doprint;
 doprint = 1;
 if(original.NrFreqChan > 4 && n != 0)
   doprint = 0;
 if(doprint) {
   for(i = 0; i < verbose.indent; i++)
     printf(" ");
   printf("  %.1f%%     \r", (100.0*(f*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan));
   fflush(stdout);
 }
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done                              \n");
  }
  return 1;
}
int preprocess_stokes(datafile_definition *original, verbose_definition verbose)
{
  long f, n, b, i;
  float I, Q, U, V, c1, c2, c3, c4;
  int basis;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Forming Stokes parameters\n");
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: only works if data is loaded into memory.");
    return 0;
  }
  if(original->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Expected 4 polarization channels, only %ld present.", original->NrPols);
    return 0;
  }
  if(original->poltype == POLTYPE_STOKES) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_stokes: Data already in Stokes parameters, nothing will be done.");
    return 1;
  }else {
    if(original->poltype != POLTYPE_COHERENCY) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_stokes: Expected poltype=%d (coherency parameters), got %d.", POLTYPE_COHERENCY, original->poltype);
      return 0;
    }
  }
  if(original->feedtype == FEEDTYPE_LINEAR || original->feedtype == FEEDTYPE_INV_LINEAR) {
    basis = 1;
  }else if(original->feedtype == FEEDTYPE_CIRCULAR || original->feedtype == FEEDTYPE_INV_CIRCULAR) {
    basis = 2;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Expected feedtype=%d, this basis is not implemented.", original->feedtype);
    return 0;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(basis == 1)
      printf("  Using a linear basis\n");
    else if(basis == 2)
      printf("  Using a circular basis\n");
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Cannot handle PA data.");
    return 0;
  }
  for(f = 0; f < original->NrFreqChan; f++) {
    for(n = 0; n < original->NrSubints; n++) {
      for(b = 0; b < original->NrBins; b++) {
 c1 = original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b];
 c2 = original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b];
 c3 = original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b];
 c4 = original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b];
 if(basis == 1) {
   I = c1+c2;
   Q = c1-c2;
   U = 2.0*c3;
   V = 2.0*c4;
 }else if(basis == 2) {
   I = c1+c2;
   Q = 2.0*c3;
   U = 2.0*c4;
   V = c1-c2;
 }
 original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b] = I;
 original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b] = Q;
 original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b] = U;
 original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b] = V;
      }
    }
  }
  original->poltype = POLTYPE_STOKES;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done                    \n");
  }
  return 1;
}
int preprocess_coherency(datafile_definition *original, verbose_definition verbose)
{
  long f, n, b, i;
  float I, Q, U, V, c1, c2, c3, c4;
  int basis;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Forming coherency parameters\n");
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: only works if data is loaded into memory.");
    return 0;
  }
  if(original->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Expected 4 polarization channels, only %ld present.", original->NrPols);
    return 0;
  }
  if(original->poltype != POLTYPE_STOKES) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Expected poltype=%d (Stokes parameters), got %d.", POLTYPE_STOKES, original->poltype);
    return 0;
  }
  if(original->feedtype == FEEDTYPE_LINEAR || original->feedtype == FEEDTYPE_INV_LINEAR) {
    basis = 1;
  }else if(original->feedtype == FEEDTYPE_CIRCULAR || original->feedtype == FEEDTYPE_INV_CIRCULAR) {
    basis = 2;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Expected feedtype=%d, this basis is not implemented.", original->feedtype);
    return 0;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(basis == 1)
      printf("  Using a linear basis\n");
    else if(basis == 2)
      printf("  Using a circular basis\n");
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Cannot handle PA data.");
    return 0;
  }
  for(f = 0; f < original->NrFreqChan; f++) {
    for(n = 0; n < original->NrSubints; n++) {
      for(b = 0; b < original->NrBins; b++) {
 I = original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b];
 Q = original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b];
 U = original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b];
 V = original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b];
 if(basis == 1) {
          c1 = 0.5*(I+Q);
   c2 = 0.5*(I-Q);
   c3 = 0.5*U;
   c4 = 0.5*V;
 }else if(basis == 2) {
          c1 = 0.5*(I+V);
   c2 = 0.5*(I-V);
   c3 = 0.5*Q;
   c4 = 0.5*U;
 }
 original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b] = c1;
 original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b] = c2;
 original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b] = c3;
 original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b] = c4;
      }
    }
  }
  original->poltype = POLTYPE_COHERENCY;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done                      \n");
  }
  return 1;
}
int preprocess_scale(datafile_definition original, float factor, float offset, verbose_definition verbose)
{
  long p, f, n, b, i;
  float sample;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Scale data with %f and offset %f\n", factor, offset);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_scale: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_scale: Cannot handle PA data.");
    return 0;
  }
  for(f = 0; f < original.NrFreqChan; f++) {
    for(n = 0; n < original.NrSubints; n++) {
      for(p = 0; p < original.NrPols; p++) {
 for(b = 0; b < original.NrBins; b++) {
   sample = original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b];
   original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] = (sample+offset)* factor;
 }
      }
      if(verbose.verbose && verbose.nocounters == 0) {
 long doprint;
 doprint = 1;
 if(original.NrFreqChan > 4 && n != 0)
   doprint = 0;
 if(doprint) {
   for(i = 0; i < verbose.indent; i++)
     printf(" ");
   printf("  %.1f%%     \r", (100.0*(f*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan));
   fflush(stdout);
 }
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done                              \n");
  }
  return 1;
}
int preprocess_dedisperse(datafile_definition *original, int update, double freq_ref, verbose_definition verbose)
{
  long p, f, n;
  int i, inffreq, inffreq_old;
  long double dt, dt_samples;
  double freq;
  if(original->freq_ref < -1.1) {
    printwarning(verbose.debug, "WARNING preprocess_dedisperse (%s): Reference frequency is unknown. The reference frequency is set to infinite frequency.", original->filename);
    original->freq_ref = -1;
  }
  if(original->freq_ref < 0)
    inffreq = 1;
  else
    inffreq = 0;
  if(freq_ref > -1.1 && freq_ref < -0.9)
    inffreq_old = 1;
  else if(freq_ref < 0) {
    printerror(verbose.debug, "ERROR preprocess_dedisperse (%s): Requested frequency is invalid (%f).", freq_ref);
    return 0;
  }else
    inffreq_old = 0;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("De-dispersing %s with DM=%f with reference frequency=", original->filename, original->dm);
    if(inffreq)
      printf("infinity\n");
    else
      printf("%f MHz\n", original->freq_ref);
  }
  if(original->isDeDisp == 1) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Data in %s is already dedispersed\n", original->filename);
    }
  }else if(original->isDeDisp == -1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_dedisperse (%s): unknown dedispersion state. Data is assumed to be already dedispersed.", original->filename);
  }
  if(original->isDeDisp == 1 || original->isDeDisp == -1) {
    if(update) {
      if(verbose.verbose) {
 for(i = 0; i < verbose.indent; i++)
   printf(" ");
 printf("  Updating reference frequency from ");
 if(inffreq_old)
   printf("infinity");
 else
   printf("%f MHz", freq_ref);
 if(inffreq)
   printf(" to infinity\n");
 else
   printf(" to %f MHz\n", original->freq_ref);
      }
      if((inffreq == 1 && inffreq_old == 1) || fabs(original->freq_ref - freq_ref) < 0.001) {
 if(verbose.verbose) {
   for(i = 0; i < verbose.indent; i++)
     printf(" ");
   printf("  Reference frequency is the same, nothing done\n");
 }
 return 1;
      }
    }else {
      return 1;
    }
  }else {
    update = 0;
  }
  dt_samples = 0;
  if(update) {
    freq = get_channel_freq(*original, 0, verbose);
    dt = -calcDMDelay(freq, freq_ref, inffreq_old, original->dm);
    dt += calcDMDelay(freq, original->freq_ref, inffreq, original->dm);
    dt_samples = dt/get_tsamp(*original, 0, verbose);
  }
  if(fabs(original->dm) < 1e-6) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_dedisperse (%s): DM appears to be not set", original->filename);
    return 1;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    f = original->NrFreqChan/2;
    freq = get_channel_freq(*original, f, verbose);
    if(update == 0)
      printf("  Rotating central frequency channel (%lf MHz) by %Lf phase\n", freq, calcDMDelay(freq, original->freq_ref, inffreq, original->dm)/(original->NrBins*get_tsamp(*original, 0, verbose)));
    else
      printf("  Rotating central frequency channel (%lf MHz) by %Lf phase\n", freq, dt_samples/(double)(original->NrBins));
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_dedisperse (%s): only works if data is loaded into memory.", original->filename);
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_dedisperse (%s): Cannot handle PA data.", original->filename);
    return 0;
  }
  for(p = 0; p < original->NrPols; p++) {
    for(f = 0; f < original->NrFreqChan; f++) {
      for(n = 0; n < original->NrSubints; n++) {
 if(update == 0) {
   dt = calcDMDelay(get_channel_freq(*original, f, verbose), original->freq_ref, inffreq, original->dm);
   dt /= get_tsamp(*original, 0, verbose);
 }else {
   dt = dt_samples;
 }
 if(rotateSinglepulse(&(original->data[original->NrBins*(p+original->NrPols*(f+n*original->NrFreqChan))]), original->NrBins, -dt, verbose) == 0)
   return 0;
      }
    }
  }
  original->isDeDisp = 1;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done              \n");
  }
  return 2;
}
int preprocess_deFaraday(datafile_definition *original, int undo, int update, double freq_ref, verbose_definition verbose)
{
  long f, n, b;
  int i, inffreq, inffreq_old;
  float dphi, phi, L, *pulseQ, *pulseU;
  verbose_definition verbose2;
  if(original->freq_ref < -1.1) {
    printwarning(verbose.debug, "WARNING preprocess_deFaraday (%s): Reference frequency is unknown. The reference frequency is set to infinite frequency.", original->filename);
    original->freq_ref = -1;
  }
  if(original->freq_ref < 0)
    inffreq = 1;
  else
    inffreq = 0;
  if(freq_ref > -1.1 && freq_ref < -0.9)
    inffreq_old = 1;
  else if(freq_ref < 0) {
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Requested original reference frequency is invalid (%f).", original->filename, freq_ref);
    return 0;
  }else
    inffreq_old = 0;
  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(undo == 0) {
      printf("De-Faraday rotate frequency channels of %s with RM=%lf with reference frequency=", original->filename, original->rm);
      if(inffreq)
 printf("infinity\n");
      else
 printf("%f MHz\n", original->freq_ref);
    }else {
      printf("Undo de-Faraday rotation of frequency channels of %s with RM=%lf\n", original->filename, original->rm);
    }
  }
  if(update && undo) {
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Cannot undo Faraday rotation and update reference frequency simultaneously.", original->filename);
    return 0;
  }
  if(original->isDeFarad == 1 && undo == 0) {
    if(update == 0) {
      if(verbose.verbose) {
 for(i = 0; i < verbose.indent; i++)
   printf(" ");
 printf("  Data in %s is already de-Faraday rotated, will not correct data.\n", original->filename);
      }
      return 1;
    }else {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Updating reference frequency from ");
      if(inffreq_old)
 printf("infinity");
      else
 printf("%f MHz", freq_ref);
      if(inffreq)
 printf(" to infinity\n");
      else
 printf(" to %f MHz\n", original->freq_ref);
      if((inffreq == 1 && inffreq_old == 1) || fabs(original->freq_ref - freq_ref) < 0.001) {
 if(verbose.verbose) {
   for(i = 0; i < verbose.indent; i++)
     printf(" ");
   printf("  Reference frequency is the same, nothing done\n");
 }
 return 1;
      }
    }
  }
  if(original->isDeFarad == 0 && update) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Data in %s is not de-Faraday rotated, cannot update reference frequency.\n", original->filename);
      return 1;
    }
  }
  if(original->isDeFarad == 0 && undo) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
    }
    printf("  Data in %s is not de-Faraday rotated, will not undo correction.\n", original->filename);
    return 1;
  }
  if(original->isDeFarad == -1) {
    if(undo == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_deFaraday (%s): unknown de-Faraday rotation state. Data will not be de-Faraday rotated.", original->filename);
    }else {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_deFaraday (%s): unknown de-Faraday rotation state. De-Faraday rotation will not be undone.", original->filename);
    }
    return 1;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Remove Faraday rotation from data which still contains the four Stokes parameters.", original->filename);
    return 0;
  }
  if(original->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Expected four polarization channels.", original->filename);
    return 0;
  }
  if(original->poltype != POLTYPE_STOKES) {
    fflush(stdout);
    if(preprocess_stokes(original, verbose2) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Conversion into Stokes parameters failed.", original->filename);
      return 0;
    }
    if(original->poltype != POLTYPE_STOKES) {
      printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Conversion into Stokes parameters failed.", original->filename);
      return 0;
    }
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): only works if data is loaded into memory.", original->filename);
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Cannot handle PA data.", original->filename);
    return 0;
  }
  if(fabs(original->rm) < 1e-3) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_deFaraday (%s): RM does not appear to be set", original->filename);
    return 1;
  }
  if(verbose.debug) {
    for(n = 0; n < verbose.indent; n++)
      printf(" ");
    if(original->isDebase != 1) {
      printf("WARNING preprocess_deFaraday: Baseline variations in time/freq will be introduced if baseline is not removed\n");
    }
  }
  dphi = 0;
  if(update) {
    double freq;
    freq = get_channel_freq(*original, 0, verbose);
    dphi = -calcRMAngle(freq, freq_ref, inffreq_old, original->rm);
    dphi += calcRMAngle(freq, original->freq_ref, inffreq, original->rm);
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    if(update == 0)
      dphi = calcRMAngle(get_channel_freq(*original, original->NrFreqChan/2, verbose), original->freq_ref, inffreq, original->rm);
    printf("  Rotating Q&U of central frequency channel (%f MHz) by %f rad = %f deg\n", get_channel_freq(*original, original->NrFreqChan/2, verbose), 2.0*dphi, 2.0*dphi*180.0/M_PI);
  }
  pulseQ = (float *)malloc(original->NrBins*sizeof(float));
  pulseU = (float *)malloc(original->NrBins*sizeof(float));
  if(pulseQ == NULL || pulseU == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Cannot allocate memory.", original->filename);
    return 0;
  }
  for(f = 0; f < original->NrFreqChan; f++) {
    if(update == 0) {
      dphi = calcRMAngle(get_channel_freq(*original, f, verbose), original->freq_ref, inffreq, original->rm);
      if(undo)
 dphi *= -1.0;
    }
    for(n = 0; n < original->NrSubints; n++) {
      if(readPulsePSRData(original, n, 1, f, 0, original->NrBins, pulseQ, verbose) != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Read error.", original->filename);
 return 0;
      }
      if(readPulsePSRData(original, n, 2, f, 0, original->NrBins, pulseU, verbose) != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Read error.", original->filename);
 return 0;
      }
      for(b = 0; b < original->NrBins; b++) {
 L = sqrt(pulseQ[b]*pulseQ[b]+pulseU[b]*pulseU[b]);
 phi = atan2(pulseU[b], pulseQ[b]);
 phi -= 2.0*dphi;
 pulseQ[b] = L*cos(phi);
 pulseU[b] = L*sin(phi);
      }
      if(writePulsePSRData(original, n, 1, f, 0, original->NrBins, pulseQ, verbose) != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Write error.", original->filename);
 return 0;
      }
      if(writePulsePSRData(original, n, 2, f, 0, original->NrBins, pulseU, verbose) != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Write error.", original->filename);
 return 0;
      }
    }
  }
  if(undo == 0)
    original->isDeFarad = 1;
  else
    original->isDeFarad = 0;
  free(pulseQ);
  free(pulseU);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done              \n");
  }
  return 2;
}
int preprocess_corrParAng(datafile_definition *original, datafile_definition *clone, int undo, verbose_definition verbose)
{
  int i;
  long n;
  double parang, parang2;
  verbose_definition verbose2;
  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("De-rotate data to ");
    if(undo)
      printf("add the effect of parallactic angle\n");
    else
      printf("remove the effect of parallactic angle\n");
  }
  if(original->isDePar == 1 && undo == 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
    }
    printf("  Data in %s is already parallactic angle corrected, so cannot do correction.\n", original->filename);
    return 1;
  }
  if(original->isDePar == 0 && undo) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
    }
    printf("  Data in %s is not parallactic angle corrected, so cannot undo correction.\n", original->filename);
    return 1;
  }
  if(original->isDePar == -1 && undo == 0) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): unknown parallactic angle state. Parallactic angle correction is ignored.", original->filename);
    return 2;
  }
  if(original->isDePar == -1 && undo) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): unknown parallactic angle state. Will not undo Parallactic angle correction.", original->filename);
    return 2;
  }
  if(verbose.verbose) {
    if(data_parang(*original, -1, &parang, verbose) == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): Cannot calculate parallactic angle, ignoring parallactic angle correction.", original->filename);
      return 2;
    }
    parang *= 180.0/M_PI;
    if(undo)
      parang *= -1.0;
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  Correcting %s for a parallactic angle (which is %f deg for midpoint of observation)\n", original->filename, parang);
    if(original->NrSubints > 1) {
      if(data_parang(*original, 0, &parang, verbose) == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): Cannot calculate parallactic angle, ignoring parallactic angle correction.", original->filename);
 return 2;
      }
      if(data_parang(*original, original->NrSubints-1, &parang2, verbose) == 0) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): Cannot calculate parallactic angle, ignoring parallactic angle correction.", original->filename);
 return 2;
      }
      parang = parang2 - parang;
      parang *= 180.0/M_PI;
      if(undo)
 parang *= -1.0;
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Total change in parallactic angle throughout observation is %f deg.\n", parang);
    }
  }
  if(clone != NULL) {
    if(make_clone(*original, clone, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_corrParAng: Cannot make clone of data, cannot remove effect of parallactic angle.");
      return 0;
    }
  }
  if(original->poltype != POLTYPE_STOKES) {
    if(clone != NULL) {
      if(preprocess_stokes(clone, verbose2) == 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR preprocess_corrParAng (%s): Cannot form the Stokes parameters.", original->filename);
 return 0;
      }
    }else {
      if(preprocess_stokes(original, verbose2) == 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR preprocess_corrParAng (%s): Cannot form the Stokes parameters.", original->filename);
 return 0;
      }
    }
  }
  for(n = 0; n < original->NrSubints; n++) {
    if(data_parang(*original, n, &parang, verbose) == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): Cannot calculate parallactic angle, ignoring parallactic angle correction.", original->filename);
      return 2;
    }
    parang *= 180.0/M_PI;
    if(undo)
      parang *= -1.0;
    if(verbose.verbose && original->NrSubints > 1 && n == 0) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("For first subint:\n");
    }
    if(n != 0) {
      verbose2.verbose = 0;
      verbose2.nocounters = 0;
    }
    if(clone != NULL) {
      if(preprocess_rotateStokes(clone, NULL, 1, n, 2.0*parang, 1, 2, verbose2) != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR preprocess_corrParAng (%s): Cannot rotate Stokes Q and U.", original->filename);
 return 0;
      }
    }else {
      if(preprocess_rotateStokes(original, NULL, 1, n, 2.0*parang, 1, 2, verbose2) != 1) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR preprocess_corrParAng (%s): Cannot rotate Stokes Q and U.", original->filename);
 return 0;
      }
    }
  }
  if(undo) {
    if(clone != NULL)
      clone->isDePar = 0;
    else
      original->isDePar = 0;
  }else {
    if(clone != NULL)
      clone->isDePar = 1;
    else
      original->isDePar = 1;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  done          \n");
  }
  return 3;
}
int preprocess_make_profile(datafile_definition original, datafile_definition *profile, int stokesI, verbose_definition verbose)
{
  int i;
  datafile_definition clone, clone2;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Constructing ");
    if(stokesI)
      printf("Stokes I ");
    printf("profile\n");
    verbose.indent += 2;
  }
  if(make_clone(original, &clone, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_make_profile: Cannot make clone of data, cannot construct profile");
    return 0;
  }
  if(original.NrPols == 4 && stokesI) {
    if(original.poltype != POLTYPE_STOKES) {
      if(preprocess_stokes(&clone, verbose) == 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR preprocess_make_profile: forming Stokes parameters failed, cannot construct profile");
 return 0;
      }
    }
    if(original.NrPols > 1) {
      if(preprocess_polselect(clone, &clone2, 0, verbose) == 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR preprocess_make_profile: Selecting Stokes I failed, cannot construct profile");
 return 0;
      }
      swap_orig_clone(&clone, &clone2, verbose);
    }
  }
  if(clone.NrFreqChan > 1) {
    if(preprocess_dedisperse(&clone, 0, 0, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_make_profile: de-dispersing failed, cannot construct profile");
      return 0;
    }
    if(preprocess_addsuccessiveFreqChans(clone, &clone2, clone.NrFreqChan, NULL, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_make_profile: summing frequency channels failed, cannot construct profile");
      return 0;
    }
    swap_orig_clone(&clone, &clone2, verbose);
  }
  if(preprocess_addsuccessivepulses(clone, profile, clone.NrSubints, 0, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_make_profile: summing subints failed, cannot construct profile");
    return 0;
  }
  closePSRData(&clone, verbose);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Constructing profile done\n");
  }
  return 1;
}
int preprocess_changeRefFreq(datafile_definition *original, double freq_ref_new, verbose_definition verbose)
{
  double freq_ref_cur;
  int i, inffreq_cur, inffreq_new;
  verbose_definition verbose2;
  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 4;
  freq_ref_cur = original->freq_ref;
  inffreq_cur = 0;
  if(original->freq_ref > -1.1 && original->freq_ref < -0.9) {
    inffreq_cur = 1;
  }else if(original->freq_ref < 0) {
    if(original->isDeDisp || original->isDeFarad) {
      printwarning(verbose.debug, "WARNING preprocess_changeRefFreq: Current reference frequency is unknown, it is assumed it was infinite frequency");
    }
    inffreq_cur = 1;
    freq_ref_cur = -1;
  }
  if(freq_ref_new > -1.1 && freq_ref_new < -0.9) {
    inffreq_new = 1;
  }else {
    inffreq_new = 0;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Updating reference frequency of %s from ", original->filename);
    if(inffreq_cur)
      printf("infinity");
    else
      printf("%f MHz", freq_ref_cur);
    if(inffreq_new)
      printf(" to infinity\n");
    else
      printf(" to %f MHz\n", freq_ref_new);
  }
  original->freq_ref = freq_ref_new;
  if(original->isDeDisp == 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Data not yet dispersed, so no dedispersion done.\n");
    }
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Re-dedispersing data.\n");
    }
    if(preprocess_dedisperse(original, 1, freq_ref_cur, verbose2) == 0) {
      printerror(verbose.debug, "WARNING preprocess_changeRefFreq: Re-dedispersion failed");
      return 0;
    }
  }
  if(original->isDeFarad == 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Data not yet de-Faraday rotated, so no de-Faraday rotation done.\n");
    }
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      printf("  Re-de-Farady rotating data.\n");
    }
    if(preprocess_deFaraday(original, 0, 1, freq_ref_cur, verbose2) == 0) {
      printerror(verbose.debug, "WARNING preprocess_changeRefFreq: Re-de-Faraday rotating data failed");
      return 0;
    }
  }
  return 1;
}

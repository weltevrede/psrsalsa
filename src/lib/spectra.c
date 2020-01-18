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
#include "psrsalsa.h"
#define USEFFTW3 1
#ifdef USEFFTW3
  #include <complex.h>
  #include <fftw3.h>
#else
  #include "nr.h"
  #include "nrutil.h"
#endif
int calc2DFS(float *data, long nry, long nrx, unsigned long fft_size, float *twodfs, pulselongitude_regions_definition *onpulse, int region, verbose_definition verbose)
{
  unsigned long nr_fftblocks;
  float junk_float;
  long junk_int, i, nf, nb, nb2, np, nrx2, bin_offpulse_left;
  int ok;
  #ifdef USEFFTW3
    float *inputdata, pwr;
    fftwf_complex *fftdata;
    fftwf_plan plan;
  #else
    float ***inputdata, **speq;
  #endif
  ok = 0;
  pwr = 0;
  junk_float = log(fft_size)/log(2);
  junk_int = junk_float;
  junk_float = pow(2, junk_int);
  if(fabs(junk_float-fft_size) > 0.1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: fft length is not a power of two!");
    return 0;
  }
  if(onpulse->nrRegions < region) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: selected region is not defined");
    return 0;
  }
  if(onpulse->bins_defined[region] == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: Region is not defined in bins");
    return 0;
  }
  nrx2 = onpulse->right_bin[region]-onpulse->left_bin[region]+1;
  #ifdef USEFFTW3
    junk_int = nrx2 / 2;
    junk_int *= 2;
    junk_int -= nrx2;
    if(junk_int != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calc2DFS: Length onpulse region should be an even number.");
      return 0;
    }
  #else
    junk_float = log(nrx2)/log(2);
    junk_int = junk_float;
    junk_float = pow(2, junk_int);
    if(fabs(junk_float-nrx2) > 0.1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calc2DFS: Length onpulse region is not a power of two. Please recompile using fftw3 support.");
      return 0;
    }
  #endif
  nr_fftblocks = nry/fft_size;
  if(nr_fftblocks == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: Cannot calculate 2dfs for %ld pulses (fft size = %ld)", nry, fft_size);
    return 0;
  }
  if(verbose.verbose) printf("Calculating 2DFS (%ld blocks)\n", nr_fftblocks);
  if(onpulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: Onpulse region is undefined.");
    return 0;
  }
  if(onpulse->nrRegions == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: Onpulse region is undefined.");
    return 0;
  }
  #ifdef USEFFTW3
    inputdata = (float *)malloc(nrx2*fft_size*sizeof(float));
    fftdata = (fftwf_complex *)fftwf_malloc(nrx2*(1+fft_size/2)*sizeof(fftwf_complex));
    if(fftdata == NULL || inputdata == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calc2DFS: fftwf_malloc failed.");
      return 0;
    }
    plan = fftwf_plan_dft_r2c_2d(nrx2, fft_size, inputdata, fftdata, FFTW_ESTIMATE);
  #else
    inputdata = f3tensor(1,1,1,nrx2,1,fft_size);
    speq = matrix(1,1,1,2*nrx2);
    if(inputdata == NULL || speq == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calc2DFS: Cannot allocate memory");
      return 0;
    }
  #endif
  bin_offpulse_left = -1;
  if(onpulse != NULL) {
    for(bin_offpulse_left = 0; bin_offpulse_left < nrx; bin_offpulse_left++) {
      ok = 1;
      for(i = 0; i < nrx2; i++) {
 if(i+bin_offpulse_left >= nrx) {
   ok = 0;
   break;
 }
 if(checkRegions(i+bin_offpulse_left, onpulse, 0, verbose) != 0) {
   ok = 0;
   break;
 }
      }
      if(ok) {
 break;
      }
    }
    if(ok == 0)
      bin_offpulse_left = -1;
  }
  if(bin_offpulse_left >= 0) {
    if(verbose.verbose) printf("  Found suitable offpulse region (%ld %ld).\n", bin_offpulse_left, bin_offpulse_left+nrx2-1);
  }else {
    if(verbose.verbose) printwarning(verbose.debug, "  WARNING calc2DFS: Didn't found suitable offpulse region.");
  }
  for(i = 0; i < nrx2*(1+fft_size/2); i++)
    twodfs[i] = 0;
  for(nf = 0; nf < nr_fftblocks; nf++) {
    for(nb = 0; nb < nrx2; nb++) {
      for(np = 0; np < fft_size; np++) {
#ifdef USEFFTW3
 inputdata[nb*fft_size + np] = data[(nf*fft_size+np)*nrx + nb + onpulse->left_bin[region]];
#else
 inputdata[1][nb+1][np+1] = data[(nf*fft_size+np)*nrx + nb + onpulse->left_bin[region]];
#endif
      }
    }
#ifdef USEFFTW3
    fftwf_execute(plan);
    for(nb = 0; nb < nrx2; nb++) {
      nb2 = nb+nrx2/2;
      if(nb2 >= nrx2)
 nb2 -= nrx2;
      for(np = 0; np < fft_size/2+1; np++) {
 if(np != 0) {
   pwr = cabs(fftdata[nb*(fft_size/2+1)+np]);
   twodfs[np*nrx2+nb2] += pwr*pwr;
 }
      }
    }
#else
    rlft3(inputdata, speq, 1, nrx2, fft_size, 1);
    for(nb = 0; nb < nrx2; nb++) {
      nb2 = nb+nrx2/2;
      if(nb2 >= nrx2)
 nb2 -= nrx2;
      for(np = 0; np < fft_size/2; np++) {
 if(np != 0)
   twodfs[np*nrx2+nb2] += inputdata[1][nb+1][2*np+1]*inputdata[1][nb+1][2*np+1]+inputdata[1][nb+1][2*np+2]*inputdata[1][nb+1][2*np+2];
      }
      twodfs[(fft_size/2)*nrx2+nb2] += speq[1][2*nb+1]*speq[1][2*nb+1]+speq[1][2*nb+2]*speq[1][2*nb+2];
    }
#endif
    if(bin_offpulse_left >= 0) {
      for(nb = 0; nb < nrx2; nb++) {
 for(np = 0; np < fft_size; np++) {
#ifdef USEFFTW3
   inputdata[nb*fft_size + np] = data[(nf*fft_size+np)*nrx + nb + bin_offpulse_left];
#else
   inputdata[1][nb+1][np+1] = data[(nf*fft_size+np)*nrx + nb + bin_offpulse_left];
#endif
 }
      }
#ifdef USEFFTW3
    fftwf_execute(plan);
    for(nb = 0; nb < nrx2; nb++) {
      nb2 = nb+nrx2/2;
      if(nb2 >= nrx2)
 nb2 -= nrx2;
      for(np = 0; np < fft_size/2+1; np++) {
 if(np != 0)
   pwr = cabs(fftdata[nb*(fft_size/2+1)+np]);
   twodfs[np*nrx2+nb2] -= pwr*pwr;
      }
    }
#else
      rlft3(inputdata, speq, 1, nrx2, fft_size, 1);
      for(nb = 0; nb < nrx2; nb++) {
       nb2 = nb+nrx2/2;
 if(nb2 >= nrx2)
   nb2 -= nrx2;
 for(np = 0; np < fft_size/2; np++) {
   if(np != 0)
     twodfs[np*nrx2+nb2] -= inputdata[1][nb+1][2*np+1]*inputdata[1][nb+1][2*np+1]+inputdata[1][nb+1][2*np+2]*inputdata[1][nb+1][2*np+2];
 }
 twodfs[(fft_size/2)*nrx2+nb2] -= speq[1][2*nb+1]*speq[1][2*nb+1]+speq[1][2*nb+2]*speq[1][2*nb+2];
      }
#endif
    }
    if(nr_fftblocks > 1 && verbose.nocounters == 0) {
      printf("Block %ld of the %ld     \r", nf+1, nr_fftblocks);
      fflush(stdout);
    }
  }
  if(nr_fftblocks > 1 && verbose.nocounters == 0)
    printf("Done                       \n");
#ifdef USEFFTW3
  fftwf_destroy_plan(plan);
  fftwf_free(fftdata);
  free(inputdata);
#else
  free_matrix(speq, 1,1,1,2*nrx2);
  free_f3tensor(inputdata,1,1,1,nrx2,1,fft_size);
#endif
  return 1;
}
int calcLRFS(float *data, long nry, long nrx, unsigned long fft_size, float *lrfs, int subtractDC, float *avrg_offpulse_lrfs_power, float *phase_track, float *phase_track_phases, int calcPhaseTrack, float freq_min, float freq_max, int track_only_first_region, float *subpulseAmplitude, int calcsubpulseAmplitude, int mask_freqs, int inverseFFT, pulselongitude_regions_definition *regions, float *var_rms, int argc, char **argv, verbose_definition verbose)
{
  long fftblock, binnr, pulsenr;
  unsigned long nr_fftblocks;
  float *data1, *lrfs_tmp, pwr, pwrtot, freq, var_mean;
  double avr_power_subtracted_from_lrfs;
  long i, j, n, tot_var_rms_samples;
  float zapmin, zapmax, p3;
  int k;
  float *phase_track_complex_template, *phase_track_complex;
  long ok, itteration, nrphasetracks, nspecbins;
#ifdef USEFFTW3
  fftwf_plan plan1;
  fftwf_plan plan2;
#endif
  if(regions != NULL) {
    *var_rms = 0;
    var_mean = 0;
    tot_var_rms_samples = 0;
    avr_power_subtracted_from_lrfs = 0;
  }
  phase_track_complex_template = NULL;
  phase_track_complex = NULL;
  nspecbins = 0;
  nr_fftblocks = nry/fft_size;
  if(verbose.verbose) printf("Calculating LRFS (%ld blocks)\n", nr_fftblocks);
  if(nr_fftblocks == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calcLRFS: Cannot calculate lrfs for %ld pulses (smaller than fft size = %ld)", nry, fft_size);
    return 0;
  }
  data1 = (float *)malloc((fft_size+2)*sizeof(float));
  lrfs_tmp = (float *)malloc(nrx*(fft_size/2+1)*sizeof(float));
  if(data1 == NULL || lrfs_tmp == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calcLRFS: Cannot allocate memory");
    return 0;
  }
#ifdef USEFFTW3
  plan1 = fftwf_plan_dft_r2c_1d(fft_size, data1, (fftwf_complex *)data1, FFTW_ESTIMATE);
  plan2 = fftwf_plan_dft_c2r_1d(fft_size, (fftwf_complex *)data1, data1, FFTW_ESTIMATE);
  if(fft_size > 2147483640) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calcLRFS: requested fft too long.");
    return 0;
  }
#endif
  if(calcPhaseTrack || inverseFFT || calcsubpulseAmplitude) {
    nspecbins = 0;
    for(i = 0; i <= fft_size/2; i++) {
      freq = i/(float)fft_size;
      if(freq >= freq_min && freq <= freq_max) {
 nspecbins++;
      }
    }
    if(verbose.verbose) {
      printf("  Using %ld spectral bins for the calculation of the subpulse track.\n", nspecbins);
    }
    if(nspecbins == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calcLRFS: Cannot calculate phase track for %ld frequency bins", nspecbins);
      return 0;
    }else if(verbose.debug) {
      for(i = 0; i <= fft_size/2; i++) {
 freq = i/(float)fft_size;
 if(freq >= freq_min && freq <= freq_max) {
   printf("    This is frequency bin with a frequency of %ld/%ld=%f cpp (P3 = %f P1).\n", i, fft_size, freq, 1.0/freq);
 }
      }
    }
  }
  if(calcPhaseTrack || calcsubpulseAmplitude) {
    phase_track_complex = (float *)malloc(2*nrx*nspecbins*nr_fftblocks*sizeof(float));
    phase_track_complex_template = (float *)malloc(2*nrx*sizeof(float));
    if(phase_track_complex == NULL || phase_track_complex_template == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calcLRFS: Cannot allocate memory");
      return 0;
    }
  }
  for(i = 0; i < nrx*(1+fft_size/2); i++)
    lrfs[i] = 0;
  nrphasetracks = 0;
  for(fftblock = 0; fftblock < nr_fftblocks; fftblock++) {
    for(i = 0; i < nrx*(fft_size/2+1); i++)
      lrfs_tmp[i] = 0;
    for(binnr = 0; binnr < nrx; binnr++) {
      pwrtot = 0;
      for(pulsenr = 0; pulsenr < fft_size; pulsenr++) {
 data1[pulsenr] = data[(fftblock*fft_size+pulsenr)*nrx + binnr];
 pwrtot += data1[pulsenr];
      }
#ifdef USEFFTW3
      fftwf_execute(plan1);
#else
      realft(data1-1, fft_size, 1);
      data1[2*(fft_size/2)] = data1[1];
      data1[2*(fft_size/2)+1] = 0;
      data1[1] = 0;
#endif
      if(calcPhaseTrack) {
 phase_track_complex[2*(nrx*nrphasetracks + binnr)] = 0;
 phase_track_complex[2*(nrx*nrphasetracks + binnr)+1] = 0;
      }
      for(i = 0; i <= fft_size/2; i++) {
 freq = i/(float)fft_size;
 pwr = data1[2*i]*data1[2*i] + data1[2*i+1]*data1[2*i+1];
 if(i == 0 && subtractDC) {
   pwr -= pwrtot*pwrtot;
 }
 if(mask_freqs == 0 || (freq >= freq_min && freq <= freq_max)) {
   lrfs_tmp[i*nrx+binnr] = pwr;
   if(calcPhaseTrack) {
     if(freq >= freq_min && freq <= freq_max) {
       phase_track_complex[2*(nrx*nrphasetracks + binnr)] += data1[2*i];
       phase_track_complex[2*(nrx*nrphasetracks + binnr)+1] += data1[2*i+1];
       nrphasetracks++;
     }
   }
 }
      }
      nrphasetracks -= nspecbins;
    }
    if(regions != NULL) {
      if(regions->nrRegions > 0) {
 for(i = 0; i <= fft_size/2; i++) {
   pwrtot = 0;
   n = 0;
   for(binnr = 0; binnr < nrx; binnr++) {
     if(checkRegions(binnr, regions, 0, verbose) == 0) {
       pwrtot += lrfs_tmp[i*nrx+binnr];
       n++;
       *var_rms += lrfs_tmp[i*nrx+binnr]*lrfs_tmp[i*nrx+binnr];
       var_mean += lrfs_tmp[i*nrx+binnr];
       tot_var_rms_samples++;
     }
   }
   if(n > 0) {
     pwrtot /= (float)n;
   }else if(i == 0 && fftblock == 0) {
     fflush(stdout);
     printwarning(verbose.debug, "WARNING: onpulse region is defined, but no offpulse region is available!");
   }
   for(binnr = 0; binnr < nrx; binnr++) {
     lrfs_tmp[i*nrx+binnr] -= pwrtot;
   }
   avr_power_subtracted_from_lrfs += pwrtot;
 }
      }
    }
    for(binnr = 0; binnr < nrx; binnr++) {
      for(i = 0; i <= fft_size/2; i++) {
 lrfs[i*nrx+binnr] += lrfs_tmp[i*nrx+binnr];
      }
    }
    nrphasetracks += nspecbins;
    if(verbose.verbose && verbose.nocounters == 0) {
      printf("  Block %ld of the %ld     \r", fftblock+1, nr_fftblocks);
      fflush(stdout);
    }
  }
  if(regions != NULL) {
    *var_rms /= (float)tot_var_rms_samples;
    var_mean /= (float)tot_var_rms_samples;
    *var_rms = sqrt((*var_rms)-var_mean*var_mean);
  }
  if(calcPhaseTrack || calcsubpulseAmplitude) {
    int itteration_max = 100;
    double complex corr, phase, *phases;
    float total_phase_offset;
    if(verbose.verbose) printf("  Coherently add phase tracks\n");
    phases = (double complex *)malloc(nrphasetracks*sizeof(double complex));
    if(phases == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calcLRFS: Cannot allocate memory");
      return 0;
    }
    for(i = 0; i < nrphasetracks; i++) {
      phases[i] = 1;
    }
    for(itteration = 0; itteration < itteration_max+1; itteration++) {
      total_phase_offset = 0;
      for(binnr = 0; binnr < 2*nrx; binnr++) {
 phase_track_complex_template[binnr] = 0;
 for(i = 0; i < nrphasetracks; i++) {
   phase_track_complex_template[binnr] += phase_track_complex[binnr+i*2*nrx];
 }
      }
      if(itteration < itteration_max) {
 for(i = 0; i < nrphasetracks; i++) {
   corr = 0;
   for(binnr = 0; binnr < nrx; binnr++) {
     ok = 1;
     if(regions != NULL) {
       if(regions->nrRegions > 0) {
  if(track_only_first_region) {
    if(checkRegions(binnr, regions, 1, verbose) == 0) {
      ok = 0;
    }
  }else {
    if(checkRegions(binnr, regions, 0, verbose) == 0) {
      ok = 0;
    }
  }
       }
     }
     if(ok) {
       corr += (phase_track_complex_template[2*binnr]+I*phase_track_complex_template[2*binnr+1])*(phase_track_complex[2*binnr+i*2*nrx]-I*phase_track_complex[2*binnr+1+i*2*nrx]);
     }
   }
   if(cabs(corr) > 0) {
     corr /= cabs(corr);
     total_phase_offset += fabs(carg(corr));
     for(binnr = 0; binnr < nrx; binnr++) {
       phase = phase_track_complex[2*binnr+i*2*nrx] + I*phase_track_complex[2*binnr+1+i*2*nrx];
       phase *= corr;
       phase_track_complex[2*binnr+i*2*nrx] = creal(phase);
       phase_track_complex[2*binnr+1+i*2*nrx] = cimag(phase);
     }
     phases[i] *= corr;
   }
 }
 total_phase_offset *= 180.0/(float)(M_PI*nrphasetracks);
 if(verbose.debug) printf("  Phase offset: %f\n", total_phase_offset);
 if(total_phase_offset <= 0.00001) {
   itteration = itteration_max-1;
 }
      }
    }
    for(binnr = 0; binnr < nrx; binnr++) {
      phase_track[binnr] = polar_angle_rad(phase_track_complex_template[2*binnr], phase_track_complex_template[2*binnr+1])*180.0/M_PI;
      phase_track[binnr] = derotate_deg(phase_track[binnr]);
    }
    if(phase_track_phases != NULL) {
      if(nspecbins > 1) {
 printwarning(verbose.debug, "  WARNING: More than one P3 modulation frequency channel is selected. The subpulse phase tracks for the individual frequency bins will be averaged, but not with an optimal weighting. You probably want to use shorter fft transforms to ensure that only one frequency bin is selected for the computation of the phase track.");
      }
      for(i = 0; i < nry/fft_size; i++) {
 phase_track_phases[i] = 0;
 for(j = 0; j < nspecbins; j++) {
   phase_track_phases[i] += carg(phases[i*nspecbins+j])*180.0/(M_PI*(float)nspecbins);
 }
      }
    }
    free(phases);
  }
  if(argc > 0 && argv != NULL) {
    for(i = 0; i < argc-1; i++) {
      if(strcmp(argv[i], "-p3zap") == 0) {
 j = sscanf(argv[i+1], "%f %f", &zapmin, &zapmax);
 if(j != 2) {
   fflush(stdout);
   printerror(verbose.debug, "Cannot parse -p3zap option.");
   exit(0);
 }
 for(j = 0; j <= (fft_size/2); j++) {
   if(zapmin > 0.9 || zapmax > 0.9) {
     p3 = j;
   }else {
     p3 = j/(float)(fft_size);
   }
   if(p3 >= zapmin && p3 <= zapmax) {
     for(k = 0; k < nrx; k++) {
       lrfs[j*nrx+k] = 0;
     }
   }
 }
      }
    }
  }
  if(calcsubpulseAmplitude) {
    for(binnr = 0; binnr < nrx; binnr++)
      subpulseAmplitude[binnr] = 0;
    float peakpwr;
    peakpwr = 0;
    for(binnr = 0; binnr < nrx; binnr++) {
      subpulseAmplitude[binnr] = sqrt(phase_track_complex_template[2*binnr]*phase_track_complex_template[2*binnr]+phase_track_complex_template[2*binnr+1]*phase_track_complex_template[2*binnr+1])/sqrt(nspecbins);
      pwrtot = 0;
      for(i = 0; i < nr_fftblocks*fft_size; i++) {
 pwrtot += data[i*nrx+binnr];
      }
      if(binnr == 0 || pwrtot > peakpwr)
 peakpwr = pwrtot;
    }
    for(binnr = 0; binnr < nrx; binnr++) {
      subpulseAmplitude[binnr] /= peakpwr;
      subpulseAmplitude[binnr] *= 2;
    }
  }
  if(verbose.verbose && verbose.nocounters == 0)
    printf("Done                       \n");
  avr_power_subtracted_from_lrfs /= (double)(fft_size/2+1);
  if(avrg_offpulse_lrfs_power != NULL) {
    *avrg_offpulse_lrfs_power = avr_power_subtracted_from_lrfs;
  }
  if(verbose.debug) {
    if(regions != NULL) {
      if(regions->nrRegions > 0) {
 printf("  Average (off-pulse) power subtracted from LRFS per frequency/pulse longitude bin: %e\n", avr_power_subtracted_from_lrfs);
      }
    }
  }
  free(data1);
  free(lrfs_tmp);
  if(calcPhaseTrack || calcsubpulseAmplitude) {
    free(phase_track_complex);
    free(phase_track_complex_template);
  }
#ifdef USEFFTW3
  fftwf_destroy_plan(plan1);
  fftwf_destroy_plan(plan2);
#endif
  return 1;
}
void calcModindex(float *lrfs, float *profile, long nrx, unsigned long fft_size, unsigned long nrpulses, float *sigma, float *rms_sigma, float *modind, float *rms_modind, pulselongitude_regions_definition *regions, float var_rms, float *avrg_offpulse_lrfs_power, verbose_definition verbose)
{
  long b, f, n, nrblocks;
  float var, max, rms, av_sigma, rms_var;
  for(b = 0; b < nrx; b++) {
    var = 0;
    for(f = 0; f <= fft_size/2; f++) {
      var += lrfs[f*nrx+b];
    }
    sigma[b] = var;
  }
  nrblocks = nrpulses / (fft_size);
  max = profile[0];
  for(b = 0; b < nrx; b++) {
    if(max < profile[b])
      max = profile[b];
  }
  if(verbose.debug) {
    printf("Modulation index calculation: profile normalisation = %e\n", max);
  }
  if(verbose.verbose) {
    if(avrg_offpulse_lrfs_power != NULL) {
      float offpulse_stddev;
      int nrspecchannels = 1+fft_size/2;
      offpulse_stddev = (*avrg_offpulse_lrfs_power)*nrspecchannels;
      offpulse_stddev *= nrpulses*2.0/(float)fft_size;
      offpulse_stddev = sqrt(offpulse_stddev);
      offpulse_stddev /= max;
      printf("Spectrally determined normalised off-pulse standard deviation = %e\n", offpulse_stddev);
    }
  }
  for(b = 0; b < nrx; b++) {
    profile[b] /= max;
    sigma[b] /= max*max;
    sigma[b] = sqrt(fabs(sigma[b])*nrpulses*2.0/(float)fft_size);
  }
  rms = 0;
  rms_var = 0;
  av_sigma = 0;
  n = 0;
  if(regions != NULL) {
    if(regions->nrRegions > 0) {
      for(b = 0; b < nrx; b++) {
 if(checkRegions(b, regions, 0, verbose) == 0) {
   rms_var += (sigma[b])*(sigma[b])*(sigma[b])*(sigma[b]);
   rms += profile[b]*profile[b];
   av_sigma += sigma[b]*sigma[b];
   n++;
 }
      }
      if(n > 0) {
 av_sigma /= (float)n;
 rms_var /= (float)n;
 rms /= (float)n;
 for(b = 0; b < nrx; b++) {
   rms_sigma[b] = 0.5*sqrt(rms_var-av_sigma*av_sigma)/sigma[b];
 }
 rms = sqrt(rms);
      }
    }
    if(var_rms > 0) {
      var_rms /= max*max;
      var_rms *= sqrt(n*(fft_size/2)*nrblocks);
      *rms_sigma = sqrt(fabs(var_rms*var_rms)*nrpulses*2.0/(float)fft_size);
    }
  }
  for(b = 0; b < nrx; b++) {
    if(profile[b] != 0)
      modind[b] = sigma[b]/profile[b];
    else
      modind[b] = -1e10;
    if(regions != NULL) {
      if(regions->nrRegions > 0) {
 rms_modind[b] = sqrt(fabs(modind[b]*modind[b]*(pow(*rms_sigma/sigma[b], 2.0) + pow(rms/profile[b], 2.0))));
      }
    }
  }
}

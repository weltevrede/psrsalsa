/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "psrsalsa.h"

#define MaxNrPolarizations 5

int pgplotProfile(char *plotDevice, int windowwidth, int windowheight, float *profile, float *stddev, float *rms_stddev, float *modindex, float *rms_modindex, int nrx, float xmin, float xmax, char *xlabel, char *ylabel, char *title, int stddev_flag, int mod_flag, int zoom_flag, float xmin_zoom, float xmax_zoom, verbose_definition verbose);
int main(int argc, char **argv)
{
  int fft_size, index, originalNrPols, selectMoreOnpulseRegions, powertwo, track_only_first_region;
  int profile_flag, lrfs_flag, stddev_flag, mod_flag, twodfs_flag, bootstrap, subtractDC, track_flag, amplitude_flag, ftrack_mask, inverseFFT, write_flag, modSimple_flag, zoom_flag, zoom_flag1, p2range_set, regionnr, p3_fold_flag, p3_fold_refine, p3_fold_cpb, p3_fold_nbin, p3_fold_onpulse_flag, originalNrPolsP3, p3fold_nosmooth, s2dfs_p3_flag, s2dfs_p2_flag;
  long fft_blocks, junk_int;
  long i, j, k, l, p, nrpointsrms;
  float xmin, xmax, xmin_zoom, xmax_zoom, mod_sigma, stddev_sigma, sampleI, freq_min, freq_max, var_rms, p3_fold, p3_fold_smoothWidth, p3fold_dphase;
  float *profileI, *lrfs, *lrfs2, *stddev, *modindex, *rms_sigma, *rms_modindex, *twodfs, *clone_profileI, *phase_track, *phase_track_phases, *amplitude_profile, slope, track_dphase;
  float *p3foldmap, *p3foldmap2;
  float zapmin, zapmax, p2min, p2max, p2, p3, junk_float;
  double *stddev_av, *modindex_av, *stddev_square, *modindex_square, rms, avrg;
  char lrfsdevice[1000], onpulseselectdevice[1000], profiledevice[1000], trackdevice[1000], amplitudedevice[1000], twodfsdevice[1000], outputname[1000], txt[1000], p3fold_device[1000], s2dfs_p3_device[1000], s2dfs_p2_device[1000];
  FILE *fout_ascii;
  psrsalsaApplication application;
  datafile_definition fin[MaxNrPolarizations], clone, fout;
  pgplot_options_definition pgplot_options;
  verbose_definition noverbose;
  initApplication(&application, "pspec", "[options] inputfile");
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_formatlist = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.switch_nocounters = 1;
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_onpulse = 1;
  application.switch_onpulsef = 1;
  application.switch_itf = 1;
  application.switch_rebin = 1;
  application.switch_nskip = 1;
  application.switch_nread = 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_conshift= 1;
  application.switch_circshift= 1;
  application.switch_shuffle = 1;
  application.switch_libversions = 1;
  fft_size = 512;
  powertwo = 0;
  lrfs_flag = 0;
  profile_flag = 0;
  stddev_flag = 0;
  mod_flag = 0;
  mod_sigma = -1;
  stddev_sigma = -1;
  bootstrap = 0;
  subtractDC = 1;
  track_flag = 0;
  amplitude_flag = 0;
  freq_min = 0;
  freq_max = 0.5;
  ftrack_mask = 0;
  inverseFFT = 0;
  write_flag = 0;
  modSimple_flag = 0;
  zoom_flag = 0;
  zoom_flag1 = 0;
  slope = 0;
  track_dphase = 0;
  twodfs_flag = 0;
  p2range_set = 0;
  selectMoreOnpulseRegions = 0;
  pgplot_clear_options(&pgplot_options);
  sprintf(lrfsdevice, "?");
  sprintf(onpulseselectdevice, "?");
  sprintf(profiledevice, "?");
  sprintf(trackdevice, "?");
  sprintf(amplitudedevice, "?");
  sprintf(twodfsdevice, "?");
  sprintf(p3fold_device, "?");
  sprintf(s2dfs_p3_device, "?");
  sprintf(s2dfs_p2_device, "?");
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;
  track_only_first_region = 0;
  p3_fold_flag = 0;
  p3_fold_cpb = 1;
  p3_fold_refine = 1;
  p3_fold_smoothWidth = -1;
  p3fold_dphase = 0;
  p3_fold_onpulse_flag = 1;
  p3fold_nosmooth = 0;
  s2dfs_p3_flag = 0;
  s2dfs_p2_flag = 0;

  application.oformat = FITS_format;
  if(argv[argc-1][0] == '-' && strcmp(argv[argc-1], "-formatlist") != 0 && strcmp(argv[argc-1], "-headerlist") != 0) {
    printerror(application.verbose_state.debug, "pspec: Last command line option is expected to be a file name.\nRun pspec without command line arguments to show help");
    terminateApplication(&application);
    return 0;
  }

  if(argc < 2) {
    printf("Program to analyse (folded) single pulse data using mostly Fourier techniques.\nIt is assumed the data contains a single pulse in each subint and that the\nbaseline is subtracted (use pmod -debase).\n\n");
    printApplicationHelp(&application);
    printf("General options:\n");
    printf("  -nfft               Set size of fft's [default=%d].\n", fft_size);
    printf("  -powertwo           When manually selecting onpulse regions, they are forced\n");
    printf("                      to be a power of two bins wide.\n");
    printf("  -w                  Write out the results to files.\n");
    printf("  -bootstrap          Find error bars on the standard deviation, modulation\n");
    printf("                      index and subpulse phase by random adding noise to the\n");
    printf("                      data. This will be done for the specified number of times\n");
    printf("                      (larger value will be more precise, but takes longer). The\n");
    printf("                      error bars (although somewhat overestimated) are more\n");
    printf("                      accurate than the analytic approximation used by default.\n");
    printf("\nOutput options:\n");
    printf("  -prof               Compute pulse profile.\n");
    printf("  -lrfs               Compute LRFS.\n");
    printf("  -DC                 Leave the DC channel in the LRFS.\n");
    printf("  -stddev             Compute standard deviation profile.\n");
    printf("  -stddev_sigma       Specify sigma threshold for the stddev output to file.\n");
    printf("                      The plot (shown with -prof) only shows 3 sigma detections.\n");
    printf("  -mod                Compute modulation index profile.\n");
    printf("  -mod_sigma          Specify sigma threshold for the mod. index output to file.\n");
    printf("                      The plot (shown with -prof) only shows 3 sigma detections.\n");
    printf("  -track              Compute subpulse phase track (use with -freq).\n");
    printf("  -track_dphase       Add specified offset (in deg) to the subpulse phase track.\n");
    printf("  -track_firstregion  Only use the first selected onpulse region to find the\n");
    printf("                      alignments of the phases of the different fft blocks.\n");
    printf("                      The other onpulse regions are still used to subtract from\n");
    printf("                      the LRFS from which the phases are derived.\n");
    printf("  -slope              Subtract slope from subpulse phases (in degrees subpulse\n");
    printf("                      phase per degree pulse longitude).\n");
    printf("  -amplitude          Compute modulation amplitude (use with -freq).\n");
    printf("  -2dfs               Compute 2DFS.\n");
    printf("  -s2dfs_p3           Compute S2DFS (sliding 2DFS P3 map).\n");
    printf("  -s2dfs_p2           Compute S2DFS (sliding 2DFS P2 map)\n");
    printf("                      (for first selected region only).\n");
    printf("  -freq               Define which fluctuation frequencies (in cpp) are used for\n");
    printf("                      the subpulse phase track/amplitude calculation\n");
    printf("  -p3fold             \"P3 n\": Fold the data using this P3 value in pulse\n");
    printf("                      periods and the P3 cycle is divided in n bins in the final\n");
    printf("                      result. If n > P3, the different bins are not independent\n");
    printf("                      and you might want to use -p3fold_smooth option to make\n");
    printf("                      the effective resolution equal to P3.\n");
    printf("  -p3fold_dphase      Add this subpulse phase in degrees (can also use -slope).\n");
    printf("  -p3fold_norefine    Do not attemt to align subsequent blocks, i.e. fixed\n");
    printf("                      period folding\n");
    printf("  -p3fold_nritt       Set the number of itterations, which produces a template\n");
    printf("                      first thereby producing better results. Default is %d.\n", p3_fold_refine);
    printf("  -p3fold_cpb         Set the number of cycles per block used in the cross\n");
    printf("                      correlation used to compensate for P3 variations. More\n");
    printf("                      means more signal to correlate (more precise alignment of\n");
    printf("                      the blocks, less means less smearing in each block because\n");
    printf("                      of P3 variation within the block. Default is %d.\n", p3_fold_cpb);
    printf("  -p3fold_smooth      Replace the tophat weight used to assign power to the P3\n");
    printf("                      bins with a Gausian weight with this width in pulse\n");
    printf("                      periods. This could make oversampling look nicer and\n");
    printf("                      reduce the effective resolution. Example: if P3=10P, you\n");
    printf("                      could set n in the -p3fold option to 20, resulting in\n");
    printf("                      oversampling with a factor 2. By setting -p3fold_smooth\n");
    printf("                      to 2, the effective resolution is reduced by a factor 2\n");
    printf("                      because each input pulse is smeared out with this width\n");
    printf("                      in pulse periods.\n");
    printf("  -p3fold_noonpulse   Ignore selected pulse longitude range, but use the full\n");
    printf("                      range when doing the cross correlations\n");
    printf("  -p2zap              \"P2min P2max\" Zap fluctuations in this P2 range in cpp.\n");
    printf("  -p3zap              \"P3min P3max\" Zap fluctuations in this P3 range.\n");
    printf("                      P3min and P3max can be specified as bins or in cpp.\n");
    printf("\nGraphics options:\n");
    printf("  -onpulsed           Set pgplot device for the selection of the onpulse region.\n");
    printf("  -profd              Set pgplot device for the pulse profile.\n");
    printf("  -lrfsd              Set pgplot device for the LRFS.\n");
    printf("  -trackd             Set pgplot device for the subpulse phase.\n");
    printf("  -amplituded         Set pgplot device for the modulation amplitude.\n");
    printf("  -2dfsd              Set pgplot device for the 2DFS.\n");
    printf("  -s2dfs_p3d          Set pgplot device for the S2DFS (P3 map).\n");
    printf("  -s2dfs_p2d          Set pgplot device for the S2DFS (P2 map).\n");
    printf("  -p3foldd            Set pgplot device for the P3 fold map.\n");
    printf("  -onpulsegr          Enables graphical selection of additional on-pulse regions\n");
    printf("                      to those defined with the -onpulse option.\n");

    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about the lrfs/2dfs/modulation index can be found in:\n");
    printf(" - Weltevrede et al. 2006, A&A, 445, 243\n");
    printf(" - Weltevrede et al. 2007, A&A, 469, 607\n");
    printf("More information about bootstrap/subpulse phase track & amplitude can be found in:\n");
    printf(" - Weltevrede et al. 2012, MNRAS, 424, 843\n");
    printf("More information about the sliding 2dfs can be found in:\n");
    printf(" - Serylak et al. 2009, A&A, 506, 865\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else if(argc > 2 || strcmp(argv[argc-1], "-formatlist") == 0 || strcmp(argv[argc-1], "-headerlist") == 0) {
    int lastindex;
    lastindex = argc-1;
    if(strcmp(argv[argc-1], "-formatlist") == 0 || strcmp(argv[argc-1], "-headerlist") == 0)
      lastindex++;
    for(i = 1; i < lastindex; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
 i = index;
      }else if(strcmp(argv[i], "-nfft") == 0 || strcmp(argv[i], "-fft_size") == 0 || strcmp(argv[i], "-fftsize") == 0 || strcmp(argv[i], "-fft_length") == 0 || strcmp(argv[i], "-fftlength") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &fft_size, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-bootstrap") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &bootstrap, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-onpulsed") == 0) {
 strcpy(onpulseselectdevice, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-profd") == 0) {
 strcpy(profiledevice, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-p3foldd") == 0) {
 strcpy(p3fold_device, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-lrfsd") == 0) {
 strcpy(lrfsdevice, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-trackd") == 0) {
 strcpy(trackdevice, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-amplituded") == 0) {
 strcpy(amplitudedevice, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-2dfsd") == 0) {
 strcpy(twodfsdevice, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-s2dfs_p3d") == 0) {
 strcpy(s2dfs_p3_device, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-s2dfs_p2d") == 0) {
 strcpy(s2dfs_p2_device, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-w") == 0) {
 write_flag = 1;
      }else if(strcmp(argv[i], "-prof") == 0) {
 profile_flag = 1;
      }else if(strcmp(argv[i], "-stddev") == 0) {
 stddev_flag = 1;
      }else if(strcmp(argv[i], "-mod") == 0) {
 mod_flag = 1;
      }else if(strcmp(argv[i], "-lrfs") == 0) {
 lrfs_flag = 1;
      }else if(strcmp(argv[i], "-track") == 0) {
 track_flag = 1;
      }else if(strcmp(argv[i], "-track_firstregion") == 0) {
 track_only_first_region = 1;
      }else if(strcmp(argv[i], "-amplitude") == 0) {
 amplitude_flag = 1;
      }else if(strcmp(argv[i], "-2dfs") == 0) {
 twodfs_flag = 1;
      }else if(strcmp(argv[i], "-s2dfs_p3") == 0) {
 s2dfs_p3_flag = 1;
      }else if(strcmp(argv[i], "-s2dfs_p2") == 0) {
 s2dfs_p2_flag = 1;
      }else if(strcmp(argv[i], "-p2zap") == 0) {
 i++;
      }else if(strcmp(argv[i], "-p3zap") == 0) {
 i++;
      }else if(strcmp(argv[i], "-powertwo") == 0) {
 powertwo = 1;
      }else if(strcmp(argv[i], "-DC") == 0) {
 subtractDC = 0;
      }else if(strcmp(argv[i], "-onpulsegr") == 0) {
 selectMoreOnpulseRegions = 1;
      }else if(strcmp(argv[i], "-p3fold") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %d", &p3_fold, &p3_fold_nbin, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 p3_fold_flag = 1;
 i++;
      }else if(strcmp(argv[i], "-p3fold_nritt") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &p3_fold_refine, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-p3fold_cpb") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &p3_fold_cpb, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-p3fold_smooth") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &p3_fold_smoothWidth, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-p3fold_norefine") == 0) {
 p3_fold_refine = 0;
      }else if(strcmp(argv[i], "-p3fold_noonpulse") == 0) {
 p3_fold_onpulse_flag = 0;
      }else if(strcmp(argv[i], "-track_dphase") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &track_dphase, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-p3fold_dphase") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &p3fold_dphase, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-slope") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &slope, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-freq") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &freq_min, &freq_max, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-mod_sigma") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &mod_sigma, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-stddev_sigma") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &stddev_sigma, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else {
        printerror(application.verbose_state.debug, "Unknown option: %s", argv[i]);
 terminateApplication(&application);
 return 0;
      }
    }
  }
  junk_float = log(fft_size)/log(2);
  junk_int = junk_float;
  junk_float = pow(2, junk_int);
  if(fabs(junk_float-fft_size) > 0.1) {
    printerror(application.verbose_state.debug, "ERROR pspec: fft length is not a power of two.");
    return 0;
  }




  for(i = 0; i < MaxNrPolarizations; i++)
    cleanPSRData(&fin[i], application.verbose_state);






  if(application.iformat <= 0)
    application.iformat = guessPSRData_format(argv[argc-1], 0, application.verbose_state);
  if(application.iformat < 1) {
    printerror(application.verbose_state.debug, "ERROR pspec: Unknown input file format.\n");
    return 0;
  }


  closePSRData(&fin[0], 0, application.verbose_state);


  if(!openPSRData(&fin[0], argv[argc-1], application.iformat, 0, 1, 0, application.verbose_state))
    return 0;



  if(PSRDataHeader_parse_commandline(&fin[0], argc, argv, application.verbose_state) == 0)
    return 0;
  for(i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-header") == 0) {
      printwarning(application.verbose_state.debug, "WARNING pspec: If using the -header option, be aware it applied BEFORE the preprocessing.");
    }
  }
  if(preprocessApplication(&application, &fin[0]) == 0) {
    return 0;
  }

  double period;
  int ret;
  ret = get_period(fin[0], 0, &period, application.verbose_state);
  if(ret == 2) {
    printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
    return 0;
  }
  if(period < 0.001) {
    printerror(application.verbose_state.debug, "ERROR pspec: The period does not appear to be set in the header. Consider using the -header option.");
    return 0;
  }
  if(get_tsamp(fin[0], 0, application.verbose_state) < 0.0000001 || get_tsamp(fin[0], 0, application.verbose_state) >= 10) {
    printerror(application.verbose_state.debug, "ERROR pspec: The sampling time does not appear to be set correctly in the header. Consider using the -header option.");
    return 0;
  }
  if(fin[0].isDebase == 0) {
    printerror(application.verbose_state.debug, "ERROR pspec: Baseline is not subtracted. Use pmod -debase first.\n");
    return 0;
  }else if(fin[0].isDebase != 1) {
    printwarning(application.verbose_state.debug, "WARNING pspec:  It is not known if baseline is already subtracted. Use pmod -debase first.\n");
  }
  if(check_baseline_subtracted(fin[0], application.verbose_state) == 0) {
    printwarning(application.verbose_state.debug, "WARNING pspec: Baseline does not appear to be subtracted. Use pmod -debase first.\n");
  }



  region_frac_to_int(&(application.onpulse), fin[0].NrBins, 0);
  if(fin[0].NrPols > MaxNrPolarizations) {
    printerror(application.verbose_state.debug, "ERROR pspec: Maximum supported input parameters is exceeded.\n");
    return 0;
  }


  originalNrPols = fin[0].NrPols;

  if(fin[0].NrPols > 0) {
    for(i = fin[0].NrPols-1; i >= 0; i--) {
      if(preprocess_polselect(fin[0], &clone, i, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspec: Error selecting Stokes parameter %ld.", i);
 return 0;
      }
      swap_orig_clone(&fin[i], &clone, application.verbose_state);
    }
  }



  cleanPSRData(&fout, application.verbose_state);
  copy_params_PSRData(fin[0], &fout, application.verbose_state);
  profileI = (float *)malloc(fin[0].NrBins*sizeof(float));
  if(profileI == NULL) {
    printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
    return 0;
  }
  fft_blocks = fin[0].NrSubints/fft_size;
  junk_int = fft_blocks*fft_size;
  if(application.verbose_state.verbose && (lrfs_flag || stddev_flag || mod_flag || track_flag || amplitude_flag || modSimple_flag || profile_flag || twodfs_flag || s2dfs_p3_flag || s2dfs_p2_flag
))
    printf("Only using %ld of the %ld pulses for the spectra being generated (%ld blocks with fft size %d).\n", junk_int, fin[0].NrSubints, fft_blocks, fft_size);
  if(junk_int == 0) {
    printerror(application.verbose_state.debug, "ERROR pspec: Not enough pulses, try a shorter fft length.");
    return 0;
  }
  read_partprofilePSRData(fin[0], profileI, NULL, 0, 0, fft_blocks*fft_size, noverbose);
  xmin = 0;
  ret = get_period(fin[0], 0, &period, application.verbose_state);
  if(ret == 2) {
    printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
    return 0;
  }
  xmax = 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period;
  if(application.onpulse.nrRegions == 0 || selectMoreOnpulseRegions) {
    strcpy(pgplot_options.viewport.plotDevice, onpulseselectdevice);
    strcpy(pgplot_options.box.xlabel, "Bin");
    strcpy(pgplot_options.box.ylabel, "Intensity");
    strcpy(pgplot_options.box.title, "Select on-pulse region");
    selectRegions(profileI, fin[0].NrBins, &pgplot_options, 0, powertwo, 1, &application.onpulse, application.verbose_state);
  }else {
    if(strcmp(profiledevice, "?") == 0)
      printf("Specify plotting device to show the profile showing the selected regions: \n  ");
    strcpy(pgplot_options.viewport.plotDevice, profiledevice);
    strcpy(pgplot_options.box.xlabel, "Phase[deg]");
    strcpy(pgplot_options.box.ylabel, "Intensity");
    strcpy(pgplot_options.box.title, fin[0].psrname);
    if(pgplotGraph1(&pgplot_options, profileI, NULL, NULL, fin[0].NrBins, xmin, xmax, 0, xmin, xmax, 0, 0, 0, 1, 0, 0, 1, 1, &application.onpulse, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR pspec: Unable to open plotdevice.\n");
      return 0;
    }
  }
  region_int_to_frac(&(application.onpulse), 1.0/(float)fin[0].NrBins, 0);
  regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);
  xmin_zoom = xmin;
  xmax_zoom = xmax;
  if(zoom_flag) {
    if(application.onpulse.nrRegions > 0) {
      if(zoom_flag1) {
 if(application.onpulse.bins_defined[0] == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: region not defined in bins");
   return 0;
 }
 xmin_zoom = application.onpulse.left_bin[0];
 xmax_zoom = application.onpulse.right_bin[0];
      }else {
 xmin_zoom = 0;
 for(i = 0; i < fin[0].NrBins; i++) {
   if(checkRegions(i, &application.onpulse, 0, application.verbose_state) != 0) {
     xmin_zoom = i;
     break;
   }
 }
 xmax_zoom = fin[0].NrBins-1;
 for(i = fin[0].NrBins-1; i >= 0; i--) {
   if(checkRegions(i, &application.onpulse, 0, application.verbose_state) != 0) {
     xmax_zoom = i;
     break;
   }
 }
      }
      if(zoom_flag == 2) {
 i = xmax_zoom - xmin_zoom;
 xmin_zoom -= i;
 xmax_zoom += i;
 if(xmin_zoom < 0)
   xmin_zoom = 0;
 if(xmax_zoom >= fin[0].NrBins)
   xmax_zoom = fin[0].NrBins-1;
      }
      xmin_zoom = (xmax-xmin)*xmin_zoom/(float)fin[0].NrBins+xmin;
      xmax_zoom = (xmax-xmin)*xmax_zoom/(float)fin[0].NrBins+xmin;
    }
  }
  stddev_av = modindex_av = stddev_square = modindex_square = NULL;
  if(lrfs_flag || stddev_flag || mod_flag || track_flag || amplitude_flag || modSimple_flag || profile_flag) {
    lrfs = (float *)malloc(originalNrPols*(fft_size/2+1)*fin[0].NrBins*sizeof(float));
    stddev = (float *)calloc(fin[0].NrBins, sizeof(float));
    rms_sigma = (float *)calloc(fin[0].NrBins, sizeof(float));
    modindex = (float *)calloc(fin[0].NrBins, sizeof(float));
    rms_modindex = (float *)malloc(fin[0].NrBins*sizeof(float));
    phase_track = (float *)malloc(fin[0].NrBins*(bootstrap+2)*sizeof(float));
    phase_track_phases = (float *)malloc((fin[0].NrSubints/fft_size)*sizeof(float));
    amplitude_profile = (float *)malloc(fin[0].NrBins*sizeof(float));
    if(lrfs == NULL || stddev == NULL || rms_sigma == NULL || modindex == NULL || rms_modindex == NULL || phase_track == NULL || phase_track_phases == NULL || amplitude_profile == NULL) {
      printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
      return 0;
    }
    if(bootstrap > 0) {
      stddev_av = (double *)calloc(fin[0].NrBins, sizeof(double));
      modindex_av = (double *)calloc(fin[0].NrBins, sizeof(double));
      stddev_square = (double *)calloc(fin[0].NrBins, sizeof(double));
      modindex_square = (double *)calloc(fin[0].NrBins, sizeof(double));
      clone_profileI = (float *)malloc(fin[0].NrBins*sizeof(float));
      if(stddev_av == NULL || modindex_av == NULL || stddev_square == NULL || modindex_square == NULL || clone_profileI == NULL) {
 printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
 return 0;
      }
      printf("Starting bootstrap\n");
      if(application.onpulse.nrRegions == 0) {
 printerror(application.verbose_state.debug, "ERROR pspec: Cannot bootstrap without a selected region");
 return 0;
      }
      nrpointsrms = 0;
      rms = 0;
      avrg = 0;
      for(i = 0; i < fin[0].NrSubints; i++) {
 for(j = 0; j < fin[0].NrBins; j++) {
   if(checkRegions(j, &application.onpulse, 0, application.verbose_state) == 0) {
     nrpointsrms++;
     if(readPulsePSRData(&fin[0], i, 0, 0, j, 1, &sampleI, application.verbose_state) != 1) {
       printerror(application.verbose_state.debug, "ERROR pspec: Error reading data");
       return 0;
     }
     rms += sampleI*sampleI;
     avrg += sampleI;
   }
 }
      }
      rms /= (double)nrpointsrms;
      avrg /= (double)nrpointsrms;
      rms -= avrg*avrg;
      rms = sqrt(rms);
      if(application.verbose_state.verbose)
 printf("  Average off-pulse intensity = %e, rms = %e based on %ld points\n", avrg, rms, nrpointsrms);
      for(j = 0; j < fin[0].NrBins; j++) {
 modindex_av[j] = 0;
 modindex_square[j] = 0;
 stddev_av[j] = 0;
 stddev_square[j] = 0;
      }
      for(i = 0; i < bootstrap; i++) {
 if(application.verbose_state.verbose && application.verbose_state.nocounters == 0) {
   printf("\r  bootstrap step %ld/%d         ", i+1, bootstrap);
   fflush(stdout);
 }
 if(preprocess_addNoise(fin[0], &clone, rms, noverbose) != 1) {
   printerror(application.verbose_state.debug, "ERROR pspec: Adding noise failed");
   return 0;
 }
 verbose_definition verbose2;
 copyVerboseState(application.verbose_state, &verbose2);
 if(application.verbose_state.verbose && (!application.verbose_state.nocounters))
   verbose2.verbose = 1;
 else
   verbose2.verbose = 0;
 if(calcLRFS(clone.data, fin[0].NrSubints, fin[0].NrBins, fft_size, lrfs, subtractDC, &phase_track[(i+2)*fin[0].NrBins], NULL, track_flag, freq_min, freq_max, track_only_first_region, NULL, 0, 0, 0, &application.onpulse, &var_rms, argc, argv, verbose2) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate LRFS");
   return 0;
 }
 read_partprofilePSRData(clone, clone_profileI, NULL, 0, 0, fft_blocks*fft_size, noverbose);
 calcModindex(lrfs, clone_profileI, fin[0].NrBins, fft_size, fft_blocks*fft_size, stddev, rms_sigma, modindex, rms_modindex, &application.onpulse, var_rms, application.verbose_state);
 for(j = 0; j < fin[0].NrBins; j++) {
   modindex_av[j] += modindex[j];
   modindex_square[j] += modindex[j]*modindex[j];
   stddev_av[j] += stddev[j];
   stddev_square[j] += stddev[j]*stddev[j];
 }
 clone.opened_flag = 1;
 closePSRData(&clone, 0, application.verbose_state);
      }
      printf("Bootstrap finished\n");
    }
    for(i = 0; i < originalNrPols; i++) {
      int track_flag_pol, amplitude_flag_pol;
      track_flag_pol = track_flag;
      amplitude_flag_pol = amplitude_flag;
      if(i != 0) {
 track_flag_pol = 0;
 amplitude_flag_pol = 0;
      }
      if(calcLRFS(fin[i].data, fin[i].NrSubints, fin[i].NrBins, fft_size, &lrfs[i*(fft_size/2+1)*fin[0].NrBins], subtractDC, phase_track, phase_track_phases, track_flag_pol, freq_min, freq_max, track_only_first_region, amplitude_profile, amplitude_flag_pol, ftrack_mask, inverseFFT, &application.onpulse, &var_rms, argc, argv, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate LRFS");
 return 0;
      }
    }
    if(profile_flag || stddev_flag || mod_flag || modSimple_flag) {
      calcModindex(lrfs, profileI, fin[0].NrBins, fft_size, fft_blocks*fft_size, stddev, rms_sigma, modindex, rms_modindex, &application.onpulse, var_rms, application.verbose_state);
      if(bootstrap > 0) {
 for(j = 0; j < fin[0].NrBins; j++) {
   stddev_av[j] /= (double)bootstrap;
   stddev_square[j] /= (double)bootstrap;
   modindex_av[j] /= (double)bootstrap;
   modindex_square[j] /= (double)bootstrap;
   rms_sigma[j] = sqrt(stddev_square[j] - stddev_av[j]*stddev_av[j]);
   rms_modindex[j] = sqrt(modindex_square[j] - modindex_av[j]*modindex_av[j]);
   modindex[j] = modindex_av[j];
   stddev[j] = stddev_av[j];
 }
      }
    }
    if(profile_flag) {
      if(strcmp(profiledevice, "?") == 0)
 printf("Specify plotting device to show the profile: \n  ");
      ret = get_period(fin[0], 0, &period, application.verbose_state);
      if(ret == 2) {
 printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
 return 0;
      }
      if(pgplotProfile(profiledevice, -1, -1, profileI, stddev, rms_sigma, modindex, rms_modindex, fin[0].NrBins, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period, "Phase[deg]", "Intensity", fin[0].psrname, stddev_flag, mod_flag|modSimple_flag, zoom_flag, xmin_zoom, xmax_zoom, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspec: Unable to open plotdevice.\n");
 return 0;
      }
      if(write_flag) {
  if(change_filename_extension(argv[argc-1], outputname, "profile", 1000, application.verbose_state) == 0)
    return 0;
  fout_ascii = fopen(outputname, "w");
  if(fout_ascii == NULL) {
    printerror(application.verbose_state.debug, "ERROR pspec: Unable to open %s.", outputname);
    return 0;
  }
  if(application.verbose_state.verbose) {
    printf("Writing bin nr, intensity, stddev, rms stddev, modindex, rms modindex to %s\n", outputname);
  }
  for(i = 0; i < fin[0].NrBins; i++) {
    fprintf(fout_ascii, "%ld %e", i, profileI[i]);
    if(stddev_sigma <= 0 || stddev[i]/rms_sigma[i] >= stddev_sigma)
      fprintf(fout_ascii, " %e %e", stddev[i], rms_sigma[i]);
    else
      fprintf(fout_ascii, " %e %e", -1.0, 0.0);
    if(mod_sigma <= 0 || modindex[i]/rms_modindex[i] >= mod_sigma)
      fprintf(fout_ascii, " %e %e", modindex[i], rms_modindex[i]);
    else
      fprintf(fout_ascii, " %e %e", -1.0, 0.0);
    fprintf(fout_ascii, "\n");
  }
  fclose(fout_ascii);
       }
    }
    if(lrfs_flag) {
      if(strcmp(lrfsdevice, "?") == 0)
 printf("Specify plotting device to show the LRFS: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, lrfsdevice);
      strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.ylabel, "P3 [cpp]");
      strcpy(pgplot_options.box.title, "LRFS");
      ret = get_period(fin[0], 0, &period, application.verbose_state);
      if(ret == 2) {
 printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
 return 0;
      }
      pgplotMap(&pgplot_options, lrfs, fin[0].NrBins, 1+fft_size/2, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period, xmin_zoom, xmax_zoom, 0, 0.5, 0, 0.5, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, application.verbose_state);
      if(write_flag) {
 fout.NrSubints = 1+fft_size/2;
 fout.NrBins = fin[0].NrBins;
 fout.NrPols = originalNrPols;
 fout.gentype = GENTYPE_LRFS;
 fout.yrangeset = 1;
 fout.yrange[0] = 0;
 fout.yrange[1] = 0.5;
 fout.tsubMode = TSUBMODE_FIXEDTSUB;
 if(fout.tsub_list != NULL) {
   free(fout.tsub_list);
 }
 fout.tsub_list = (double *)malloc(sizeof(double));
 if(fout.tsub_list == NULL) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
   return 0;
 }
 fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
 if(change_filename_extension(argv[argc-1], outputname, "lrfs", 1000, application.verbose_state) == 0) {
   return 0;
 }
 if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
   return 0;
 if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
   return 0;
 }
 lrfs2 = (float *)malloc(originalNrPols*(fft_size/2+1)*fin[0].NrBins*sizeof(float));
 if(lrfs2 == NULL) {
   printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
   return 0;
 }
 for(i = 0; i < fft_size/2+1; i++) {
   for(p = 0; p < originalNrPols; p++) {
     memcpy(&lrfs2[(originalNrPols*i+p)*fin[0].NrBins], &lrfs[(p*(fft_size/2+1)+i)*fin[0].NrBins], fin[0].NrBins*sizeof(float));
   }
 }
 if(writePSRData(&fout, lrfs2, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
   return 0;
 }
 free(lrfs2);
 closePSRData(&fout, 1, application.verbose_state);
 fout.gentype = GENTYPE_UNDEFINED;
 fout.yrangeset = 0;
 fout.xrangeset = 0;
 fout.NrPols = 1;
      }
    }
    if(track_flag) {
      for(j = 0; j < fin[0].NrBins; j++) {
 phase_track[fin[0].NrBins+j] = -1;
      }
      if(bootstrap > 0) {
 for(k = 0; k < 100; k++) {
   for(j = 0; j < fin[0].NrBins; j++) {
     phase_track[fin[0].NrBins+j] = 0;
   }
   for(i = 0; i < bootstrap; i++) {
     for(j = 0; j < fin[0].NrBins; j++) {
       float xdeg;
       xdeg = phase_track[j] - phase_track[(i+2)*fin[0].NrBins+j];
       if(fabs(xdeg+360) < fabs(xdeg))
  xdeg += 360;
       if(fabs(xdeg-360) < fabs(xdeg))
  xdeg -= 360;
       phase_track[fin[0].NrBins+j] += xdeg*xdeg;
     }
   }
   for(j = 0; j < fin[0].NrBins; j++) {
     phase_track[fin[0].NrBins+j] = sqrt(phase_track[fin[0].NrBins+j]/(float)bootstrap);
   }
   float offset, weight;
   for(i = 0; i < bootstrap; i++) {
     offset = 0;
     weight = 0;
     for(j = 0; j < fin[0].NrBins; j++) {
       if(checkRegions(j, &application.onpulse, 0, application.verbose_state) != 0) {
  float xdeg;
  xdeg = phase_track[(i+2)*fin[0].NrBins+j] - phase_track[j];
  if(fabs(xdeg+360) < fabs(xdeg))
    xdeg += 360;
  if(fabs(xdeg-360) < fabs(xdeg))
    xdeg -= 360;
  offset += (xdeg)/phase_track[fin[0].NrBins+j];
  weight += 1.0/phase_track[fin[0].NrBins+j];
       }
     }
     offset /= weight;
     if(application.verbose_state.verbose && application.verbose_state.nocounters == 0) {
       printf("Itteration %9ld: Offset track %9ld: %20e\n", k+1, i, offset);
     }
     for(j = 0; j < fin[0].NrBins; j++) {
       phase_track[(i+2)*fin[0].NrBins+j] -= offset;
       if(phase_track[(i+2)*fin[0].NrBins+j] < 0)
  phase_track[(i+2)*fin[0].NrBins+j] += 360;
       if(phase_track[(i+2)*fin[0].NrBins+j] >= 360)
  phase_track[(i+2)*fin[0].NrBins+j] -= 360;
     }
   }
 }
      }
      for(i = 0; i < fin[0].NrBins; i++) {
 float xdeg;
 xdeg = (xmax-xmin)*i/(float)(fin[0].NrBins-1);
 phase_track[i] -= slope*xdeg - track_dphase;
 phase_track[i] = derotate_deg(phase_track[i]);
      }
      if(bootstrap > 0) {
 for(j = 0; j < bootstrap; j++) {
   for(i = 0; i < fin[0].NrBins; i++) {
     float xdeg;
     xdeg = (xmax-xmin)*i/(float)(fin[0].NrBins-1);
     phase_track[(j+2)*fin[0].NrBins+i] -= slope*xdeg - track_dphase;
     phase_track[(j+2)*fin[0].NrBins+i] = derotate_deg(phase_track[(j+2)*fin[0].NrBins+i]);
   }
 }
      }
      if(strcmp(trackdevice, "?") == 0)
 printf("Specify plotting device to show the subpulse phase track: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, trackdevice);
      strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.ylabel, "Subpulse phase");
      strcpy(pgplot_options.box.title, "Subpulse phase track");
      if(bootstrap > 0) {
 pgplotGraph1(&pgplot_options, phase_track, NULL, &phase_track[fin[0].NrBins], fin[0].NrBins, xmin, xmax, 0, xmin_zoom, xmax_zoom, 0, 0, 0, 0, 0, 0, 1, 1, NULL, application.verbose_state);
      }else {
 pgplotGraph1(&pgplot_options, phase_track, NULL, NULL, fin[0].NrBins, xmin, xmax, 0, xmin_zoom, xmax_zoom, 0, 0, 0, 0, 0, 0, 1, 1, NULL, application.verbose_state);
      }
      if(write_flag) {
  if(change_filename_extension(argv[argc-1], outputname, "track", 1000, application.verbose_state) == 0)
    return 0;
  fout_ascii = fopen(outputname, "w");
  if(fout_ascii == NULL) {
    printerror(application.verbose_state.debug, "ERROR pspec: Unable to open %s.", outputname);
    return 0;
  }
  for(i = 0; i < fin[0].NrBins; i++) {
    fprintf(fout_ascii, "%ld %f %e\n", i, phase_track[i], phase_track[fin[0].NrBins+i]);
  }
  fclose(fout_ascii);
      }
    }
    if(amplitude_flag) {
      if(strcmp(amplitudedevice, "?") == 0)
 printf("Specify plotting device to show the subpulse amplitude profile: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, amplitudedevice);
      pgplot_options.viewport.dontclose = 1;
      strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.ylabel, "Subpulse amplitude");
      strcpy(pgplot_options.box.title, "Subpulse amplitude profile");
      if(!(profile_flag || stddev_flag || mod_flag || modSimple_flag)) {
 float imax;
 imax = 1;
 for(i = 0; i < fin[0].NrBins; i++) {
   if(i == 0 || profileI[i] > imax)
     imax = profileI[i];
 }
 for(i = 0; i < fin[0].NrBins; i++) {
   profileI[i] /= imax;
 }
      }
      pgplotGraph1(&pgplot_options, profileI, NULL, NULL, fin[0].NrBins, xmin, xmax, 0, xmin_zoom, xmax_zoom, 0, 0, 0, 0, 0, 0, 1, 1, NULL, application.verbose_state);
      ppgsci(2);
      for(i = 0; i < fin[0].NrBins; i++) {
 float x;
 x = (xmax-xmin)*i/(float)(fin[0].NrBins-1);
 if(i == 0) {
   ppgmove(x, amplitude_profile[i]);
 }else {
   ppgdraw(x, amplitude_profile[i]);
 }
      }
      ppgsci(1);
      ppgclos();
      if(write_flag) {
 if(change_filename_extension(argv[argc-1], outputname, "amplitude", 1000, application.verbose_state) == 0)
   return 0;
 fout_ascii = fopen(outputname, "w");
 if(fout_ascii == NULL) {
   printerror(application.verbose_state.debug, "ERROR pspec: Unable to open %s.", outputname);
   return 0;
 }
 for(i = 0; i < fin[0].NrBins; i++) {
   fprintf(fout_ascii, "%ld %f\n", i, amplitude_profile[i]);
 }
 fclose(fout_ascii);
      }
    }
    free(lrfs);
    free(stddev);
    free(modindex);
    free(rms_sigma);
    free(rms_modindex);
    free(phase_track);
    free(phase_track_phases);
    free(amplitude_profile);
  }
  if(twodfs_flag) {
    for(regionnr = 0; regionnr < application.onpulse.nrRegions; regionnr++) {
      if(application.onpulse.bins_defined[regionnr] == 0) {
 printerror(application.verbose_state.debug, "ERROR pspec: region not defined in bins");
 return 0;
      }
      twodfs = (float *)malloc((1+fft_size/2)*(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1)*sizeof(float));
      if(twodfs == NULL) {
 printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
 return 0;
      }
      if(calc2DFS(fin[0].data, fin[0].NrSubints, fin[0].NrBins, fft_size, twodfs, &application.onpulse, regionnr, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate 2DFS");
 return 0;
      }
      for(i = 1; i < argc-1; i++) {
 if(strcmp(argv[i], "-p3zap") == 0) {
   if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &zapmin, &zapmax, NULL) == 0) {
     printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
     return 0;
   }
   for(j = 0; j < (1+fft_size/2); j++) {
     if(zapmin > 0.9 || zapmax > 0.9) {
       p3 = j;
     }else {
       p3 = j/(float)fft_size;
     }
     if(p3 >= zapmin && p3 <= zapmax) {
       for(k = application.onpulse.left_bin[regionnr]; k <= application.onpulse.right_bin[regionnr]; k++) {
  twodfs[j*(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1)+(k-application.onpulse.left_bin[regionnr])] = 0;
       }
     }
   }
 }
 if(strcmp(argv[i], "-p2zap") == 0) {
   if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &zapmin, &zapmax, NULL) == 0) {
     printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
     return 0;
   }
   for(j = 0; j < (application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1); j++) {
     p2 = j*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1)-fin[0].NrBins/2-0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1);
     if(p2 >= zapmin && p2 <= zapmax) {
       for(k = 0; k < (1+fft_size/2); k++) {
  twodfs[k*(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1)+(j)] = 0;
       }
     }
   }
 }
      }
      if(p2range_set == 0) {
 p2min = -fin[0].NrBins/2-0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1);
 p2max = +fin[0].NrBins/2-0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1);
      }
      if(strcmp(twodfsdevice, "?") == 0)
 printf("Specify plotting device to show the 2DFS: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, twodfsdevice);
      strcpy(pgplot_options.box.xlabel, "P2 [cpp]");
      strcpy(pgplot_options.box.ylabel, "P3 [cpp]");
      strcpy(pgplot_options.box.title, "2DFS");
      pgplotMap(&pgplot_options, twodfs, application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1, fft_size/2+1, -fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1), fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1), p2min, p2max, 0, 0.5, 0, 0.5, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state);
      if(write_flag) {
 fout.NrSubints = fft_size/2+1;
 fout.NrBins = application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1;
 fout.gentype = GENTYPE_2DFS;
 fout.tsubMode = TSUBMODE_FIXEDTSUB;
 if(fout.tsub_list != NULL)
   free(fout.tsub_list);
 fout.tsub_list = (double *)malloc(sizeof(double));
 if(fout.tsub_list == NULL) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
   return 0;
 }
 fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
 fout.yrangeset = 1;
 fout.yrange[0] = 0;
 fout.yrange[1] = 0.5;
 fout.xrangeset = 1;
 fout.xrange[0] = -fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1);
 fout.xrange[1] = fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1);
 sprintf(txt, "%d.2dfs", regionnr+1);
 if(change_filename_extension(argv[argc-1], outputname, txt, 1000, application.verbose_state) == 0) {
   return 0;
 }
 if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
   return 0;
 if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
   return 0;
 }
 if(writePSRData(&fout, twodfs, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
   return 0;
 }
 closePSRData(&fout, 1, application.verbose_state);
 fout.gentype = GENTYPE_UNDEFINED;
 fout.yrangeset = 0;
 fout.xrangeset = 0;
      }
      free(twodfs);
    }
  }
  if(s2dfs_p3_flag || s2dfs_p2_flag) {
    float *s2dfs_p3, *s2dfs_p2;
    printf("Calculating S2DFS\n");
    if(application.onpulse.nrRegions < 1) {
      printerror(application.verbose_state.debug, "ERROR pspec: region is not defined");
      return 0;
    }
    if(application.onpulse.bins_defined[0] == 0) {
      printerror(application.verbose_state.debug, "ERROR pspec: region not defined in bins");
      return 0;
    }
    twodfs = (float *)malloc((1+fft_size/2)*(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1)*sizeof(float));
    s2dfs_p3 = (float *)malloc((1+fft_size/2)*(fin[0].NrSubints-fft_size+1)*sizeof(float));
    s2dfs_p2 = (float *)malloc((application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1)*(fin[0].NrSubints-fft_size+1)*sizeof(float));
    if(twodfs == NULL || s2dfs_p3 == NULL || s2dfs_p2 == NULL) {
      printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
      return 0;
    }
    verbose_definition noverbose;
    copyVerboseState(application.verbose_state, &noverbose);
    noverbose.verbose = 0;
    noverbose.nocounters = 1;
    for(i = 0; i < fin[0].NrSubints-fft_size+1; i++) {
      if(calc2DFS(&fin[0].data[i*fin[0].NrBins], fft_size, fin[0].NrBins, fft_size, twodfs, &application.onpulse, 0, noverbose) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate 2DFS");
 return 0;
      }
      for(l = 1; l < argc-1; l++) {
 if(strcmp(argv[l], "-p3zap") == 0) {
   if(parse_command_string(application.verbose_state, argc, argv, l+1, 0, -1, "%f %f", &zapmin, &zapmax, NULL) == 0) {
     printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
     return 0;
   }
   for(j = 0; j < (1+fft_size/2); j++) {
     p3 = j/(float)fft_size;
     if(p3 >= zapmin && p3 <= zapmax) {
       for(k = 0; k < (application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1); k++) {
  twodfs[j*(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1)+k] = 0;
       }
     }
   }
 }
      }
      if(s2dfs_p3_flag) {
 for(j = 0; j < (1+fft_size/2); j++) {
   s2dfs_p3[j*(fin[0].NrSubints-fft_size+1)+i] = 0;
   for(k = 0; k < (application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1); k++) {
     s2dfs_p3[j*(fin[0].NrSubints-fft_size+1)+i] += twodfs[j*(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1)+k];
   }
 }
      }
      if(s2dfs_p2_flag) {
 for(k = 0; k < (application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1); k++) {
   s2dfs_p2[k*(fin[0].NrSubints-fft_size+1)+i] = 0;
   for(j = 0; j < (1+fft_size/2); j++) {
     s2dfs_p2[k*(fin[0].NrSubints-fft_size+1)+i] += twodfs[j*(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1)+k];
   }
 }
      }
      if(!application.verbose_state.nocounters)
 printf("  Block %ld of the %ld  (%.2f%%)      \r", i, fin[0].NrSubints-fft_size+1, 100*(i+1)/(float)(fin[0].NrSubints-fft_size+1));
    }
    printf("  done                                             \n");
    if(s2dfs_p3_flag) {
      if(strcmp(s2dfs_p3_device, "?") == 0)
 printf("Specify plotting device to show the S2DFS P3 map: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, s2dfs_p3_device);
      strcpy(pgplot_options.box.xlabel, "Block number");
      strcpy(pgplot_options.box.ylabel, "P3 [cpp]");
      strcpy(pgplot_options.box.title, "S2DFS");
      pgplotMap(&pgplot_options, s2dfs_p3, fin[0].NrSubints-fft_size+1, fft_size/2+1, 0, fin[0].NrSubints-fft_size+1, 0, fin[0].NrSubints-fft_size+1, 0, 0.5, 0, 0.5, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state);
      if(write_flag) {
 fout.NrSubints = fft_size/2+1;
 fout.NrBins = fin[0].NrSubints-fft_size+1;
 fout.gentype = GENTYPE_S2DFSP3;
 fout.tsubMode = TSUBMODE_FIXEDTSUB;
 if(fout.tsub_list != NULL)
   free(fout.tsub_list);
 fout.tsub_list = (double *)malloc(sizeof(double));
 if(fout.tsub_list == NULL) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
   return 0;
 }
 fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
 fout.yrangeset = 1;
 fout.yrange[0] = 0;
 fout.yrange[1] = 0.5;
 if(change_filename_extension(argv[argc-1], outputname, "s2dfs_p3", 1000, application.verbose_state) == 0) {
   return 0;
 }
 if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
   return 0;
 if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
   return 0;
 }
 if(writePSRData(&fout, s2dfs_p3, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
   return 0;
 }
 closePSRData(&fout, 1, application.verbose_state);
      }
    }
    if(s2dfs_p2_flag) {
      if(strcmp(s2dfs_p2_device, "?") == 0)
 printf("Specify plotting device to show the S2DFS P2 map: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, s2dfs_p2_device);
      strcpy(pgplot_options.box.xlabel, "Block number");
      strcpy(pgplot_options.box.ylabel, "P2 [cpp]");
      strcpy(pgplot_options.box.title, "S2DFS");
      if(p2range_set == 0) {
 p2min = -fin[0].NrBins/2-0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1);
 p2max = +fin[0].NrBins/2-0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1);
      }
      pgplotMap(&pgplot_options, s2dfs_p2, fin[0].NrSubints-fft_size+1, (application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1), 0, fin[0].NrSubints-fft_size+1, 0, fin[0].NrSubints-fft_size+1, -fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1), fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1), p2min, p2max, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state);
      if(write_flag) {
 fout.NrSubints = (application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1);
 fout.NrBins = fin[0].NrSubints-fft_size+1;
 fout.gentype = GENTYPE_S2DFSP2;
 fout.tsubMode = TSUBMODE_FIXEDTSUB;
 if(fout.tsub_list != NULL)
   free(fout.tsub_list);
 fout.tsub_list = (double *)malloc(sizeof(double));
 if(fout.tsub_list == NULL) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
   return 0;
 }
 fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
 fout.yrangeset = 1;
 fout.yrange[0] = -fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1);
 fout.yrange[1] = fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1);
 if(change_filename_extension(argv[argc-1], outputname, "s2dfs_p2", 1000, application.verbose_state) == 0) {
   return 0;
 }
 if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
   return 0;
 if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
   return 0;
 }
 if(writePSRData(&fout, s2dfs_p2, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
   return 0;
 }
 closePSRData(&fout, 1, application.verbose_state);
      }
    }
    free(twodfs);
    free(s2dfs_p3);
    free(s2dfs_p2);
  }
  if(p3_fold_flag) {
    p3foldmap = malloc(originalNrPols*fin[0].NrBins * p3_fold_nbin * sizeof(float));
    if(p3foldmap == NULL) {
      printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
      return 0;
    }
    originalNrPolsP3 = originalNrPols;
      if(originalNrPols > 1) {
 printwarning(application.verbose_state.debug, "WARNING pspec: Only first polarization chanel is folded.");
 originalNrPolsP3 = 1;
      }
    for(i = 0; i < originalNrPolsP3; i++) {
      if(p3_fold_onpulse_flag) {
 if(foldP3(fin[i].data, fin[i].NrSubints, fin[i].NrBins, &p3foldmap[i*fin[0].NrBins * p3_fold_nbin], p3_fold_nbin, p3_fold, p3_fold_refine, p3_fold_cpb, p3fold_nosmooth, p3_fold_smoothWidth, slope*360.0/(float)fin[i].NrBins, p3fold_dphase, &application.onpulse
, application.verbose_state) == 0) {
   return 0;
 }
      }else {
 if(foldP3(fin[i].data, fin[i].NrSubints, fin[i].NrBins, &p3foldmap[i*fin[0].NrBins * p3_fold_nbin], p3_fold_nbin, p3_fold, p3_fold_refine, p3_fold_cpb, p3fold_nosmooth, p3_fold_smoothWidth, slope*360.0/(float)fin[i].NrBins, p3fold_dphase, NULL
, application.verbose_state) == 0) {
   return 0;
 }
      }
    }
    if(strcmp(p3fold_device, "?") == 0)
      printf("Specify plotting device to show the P3 fold map: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, p3fold_device);
      strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.ylabel, "P3 [pulse periods]");
      strcpy(pgplot_options.box.title, "P3 fold");
      ret = get_period(fin[0], 0, &period, application.verbose_state);
      if(ret == 2) {
 printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
 return 0;
      }
      pgplotMap(&pgplot_options, p3foldmap, fin[0].NrBins, p3_fold_nbin, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period, xmin_zoom, xmax_zoom, 0.5*p3_fold/(float)p3_fold_nbin, 0.5*p3_fold/(float)p3_fold_nbin + p3_fold*(p3_fold_nbin-1)/(float)p3_fold_nbin, 0, p3_fold, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, application.verbose_state);
    if(write_flag) {
      fout.NrSubints = p3_fold_nbin;
      fout.NrBins = fin[0].NrBins;
      fout.NrPols = originalNrPolsP3;
      fout.gentype = GENTYPE_P3FOLD;
      fout.tsubMode = TSUBMODE_FIXEDTSUB;
      if(fout.tsub_list != NULL)
 free(fout.tsub_list);
      fout.tsub_list = (double *)malloc(sizeof(double));
      if(fout.tsub_list == NULL) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
 return 0;
      }
      fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
      fout.yrangeset = 1;
      fout.yrange[0] = 0.5*p3_fold/(float)p3_fold_nbin;
      fout.yrange[1] = 0.5*p3_fold/(float)p3_fold_nbin + p3_fold*(p3_fold_nbin-1)/(float)p3_fold_nbin;
      if(change_filename_extension(argv[argc-1], outputname, "p3fold", 1000, application.verbose_state) == 0) {
 return 0;
      }
      if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state)) {
 printerror(application.verbose_state.debug, "ERROR pspec: Unable to open file for writing.\n");
 return 0;
      }
      if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
 return 0;
      }
      p3foldmap2 = malloc(originalNrPols*fin[0].NrBins * p3_fold_nbin * sizeof(float));
      if(p3foldmap2 == NULL) {
 printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
 return 0;
      }
      for(i = 0; i < p3_fold_nbin; i++) {
 for(p = 0; p < originalNrPolsP3; p++) {
   memcpy(&p3foldmap2[(originalNrPolsP3*i+p)*fin[0].NrBins], &p3foldmap[(p*p3_fold_nbin+i)*fin[0].NrBins], fin[0].NrBins*sizeof(float));
 }
      }
      if(writePSRData(&fout, p3foldmap2, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
 return 0;
      }
      free(p3foldmap2);
      closePSRData(&fout, 1, application.verbose_state);
      fout.gentype = GENTYPE_UNDEFINED;
      fout.yrangeset = 0;
      fout.xrangeset = 0;
      fout.NrPols = 1;
    }
    free(p3foldmap);
  }
  closePSRData(&fout, 0, application.verbose_state);
  for(i = 0; i < MaxNrPolarizations; i++)
    closePSRData(&fin[i], 0, application.verbose_state);
  free(profileI);
  ppgend();
  if(bootstrap > 0) {
    free(stddev_av);
    free(modindex_av);
    free(stddev_square);
    free(modindex_square);
    free(clone_profileI);
  }
  terminateApplication(&application);
  return 0;
}
int pgplotProfile(char *plotDevice, int windowwidth, int windowheight, float *profile, float *stddev, float *rms_stddev, float *modindex, float *rms_modindex, int nrx, float xmin, float xmax, char *xlabel, char *ylabel, char *title, int stddev_flag, int mod_flag, int zoom_flag, float xmin_zoom, float xmax_zoom, verbose_definition verbose)
{
  float min, max;
  float x, y, dy;
  int i, deviceID;
  deviceID = ppgopen(plotDevice);
  if(deviceID <= 0) {
    printerror(verbose.debug, "ERROR pspec: Cannot open plot device");
    return 0;
  }
  if(windowwidth > 0 && windowheight > 0) {
    x = windowwidth*0.01175548589341692789994673739445429916373;
    y = (windowheight-1)/(float)windowwidth;
    ppgpap(x,y);
  }
  ppgask(0);
  ppgslw(1);
  ppgpage();
  ppgsvp(0.1, 0.9, 0.1, 0.9);
  min = max = profile[0];
  for(i = 1; i < nrx; i++) {
    if(profile[i] > max)
      max = profile[i];
    if(profile[i] < min)
      min = profile[i];
    if(stddev_flag) {
      if(stddev[i] > max)
 max = stddev[i];
    }
    if(mod_flag) {
      if(modindex[i] > 3*rms_modindex[i]) {
 if(modindex[i]+rms_modindex[i] > max)
   max = modindex[i]+rms_modindex[i];
      }
    }
  }
  ppgsci(1);
  ppgswin(xmin_zoom, xmax_zoom, min, max*1.03);
  ppgbox("bcnsti",0.0,0,"bcntsi",0.0,0);
  ppglab(xlabel, ylabel, title);
  ppgsci(1);
  ppgslw(5);
  for(i = 0; i < nrx; i++) {
    x = i*(xmax-xmin)/(float)nrx + xmin;
    y = profile[i];
    if(i == 0)
      ppgmove(x, y);
    else
      ppgdraw(x, y);
  }
  ppgslw(1);
  if(stddev_flag) {
    for(i = 0; i < nrx; i++) {
      x = i*(xmax-xmin)/(float)nrx + xmin;
      y = stddev[i];
      if(y > 3*rms_stddev[i]) {
 ppgpt1(x, y, 4);
      }
    }
  }
  if(mod_flag) {
    for(i = 0; i < nrx; i++) {
      x = i*(xmax-xmin)/(float)nrx + xmin;
      y = modindex[i];
      dy = rms_modindex[i];
      if(y > 3*dy) {
 ppgslw(5);
 ppgpt1(x, y, -1);
 ppgslw(1);
 ppgerr1(6, x, y, dy, 3);
      }
    }
  }
  ppgclos();
  return 1;
}

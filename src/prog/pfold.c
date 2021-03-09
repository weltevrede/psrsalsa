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
int main(int argc, char **argv)
{
  int index, originalNrPols, originalNrPolsP3;
  long i, p;
  int p3_fold_flag, p3_fold_refine, p3_fold_cpb, p3_fold_nbin, p3_fold_onpulse_flag;
  int write_flag, zoom_flag, zoom_flag1, selectMoreOnpulseRegions;
  float p3_fold, p3_fold_smoothWidth, p3fold_dphase, p3fold_nosmooth, slope;
  float xmin, xmax, xmin_zoom, xmax_zoom, *profileI, *p3foldmap, *p3foldmap2;
  char onpulseselectdevice[1000], p3fold_device[1000], outputname[1000];
  psrsalsaApplication application;
  pgplot_options_definition pgplot_options;
  verbose_definition noverbose;
  datafile_definition fin[MaxNrPolarizations], clone, fout;
  initApplication(&application, "pfold", "[options] inputfile");
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
  application.switch_rebin = 1;
  application.switch_nskip = 1;
  application.switch_nread = 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_conshift= 1;
  application.switch_circshift= 1;
  application.switch_libversions = 1;
  application.switch_history_cmd_only = 1;
  write_flag = 0;
  zoom_flag = 0;
  zoom_flag1 = 0;
  selectMoreOnpulseRegions = 0;
  slope = 0;
  p3fold_dphase = 0;
  p3fold_nosmooth = 0;
  p3_fold_flag = 0;
  p3_fold_cpb = 1;
  p3_fold_refine = 1;
  p3_fold_smoothWidth = -1;
  p3_fold_onpulse_flag = 1;
  sprintf(onpulseselectdevice, "?");
  sprintf(p3fold_device, "?");
  pgplot_clear_options(&pgplot_options);
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;
  application.oformat = FITS_format;
  if(argv[argc-1][0] == '-' && strcmp(argv[argc-1], "-formatlist") != 0 && strcmp(argv[argc-1], "-headerlist") != 0) {
    printerror(application.verbose_state.debug, "pfold: Last command line option is expected to be a file name.\nRun pfold without command line arguments to show help");
    terminateApplication(&application);
    return 0;
  }
  if(argc < 2) {
    printf("Program to fold (folded) single pulse data thereby visualising the subpulse modulation cycle.\n\n");
    printApplicationHelp(&application);
    printf("General options:\n");
    printf("  -w                  Write out the results to files.\n");
    printf("\nOutput options related to p3 folding:\n");
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
    printf("  -slope              Subtract slope from subpulse phases (in degrees subpulse\n");
    printf("                      phase per degree pulse longitude).\n");
    printf("\nGraphics options:\n");
    printf("  -onpulsed           Set pgplot device for the selection of the onpulse region.\n");
    printf("  -p3foldd            Set pgplot device for the P3 fold map.\n");
    printf("  -onpulsegr          Enables graphical selection of additional on-pulse regions\n");
    printf("                      to those defined with the -onpulse option.\n");
    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else if(argc > 2 || strcmp(argv[argc-1], "-formatlist") == 0 || strcmp(argv[argc-1], "-headerlist") == 0) {
    int lastindex;
    lastindex = argc-1;
    if(strcmp(argv[argc-1], "-formatlist") == 0 || strcmp(argv[argc-1], "-headerlist") == 0)
      lastindex++;
    for(i = 1; i <= lastindex; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
 i = index;
      }else if(strcmp(argv[i], "-onpulsed") == 0) {
 strcpy(onpulseselectdevice, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-p3foldd") == 0) {
 strcpy(p3fold_device, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-w") == 0) {
 write_flag = 1;
      }else if(strcmp(argv[i], "-onpulsegr") == 0) {
 selectMoreOnpulseRegions = 1;
      }else if(strcmp(argv[i], "-p3fold") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %d", &p3_fold, &p3_fold_nbin, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 p3_fold_flag = 1;
 i++;
      }else if(strcmp(argv[i], "-p3fold_nritt") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &p3_fold_refine, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-p3fold_cpb") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &p3_fold_cpb, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-p3fold_smooth") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &p3_fold_smoothWidth, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-p3fold_norefine") == 0) {
 p3_fold_refine = 0;
      }else if(strcmp(argv[i], "-p3fold_noonpulse") == 0) {
 p3_fold_onpulse_flag = 0;
      }else if(strcmp(argv[i], "-p3fold_dphase") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &p3fold_dphase, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-slope") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &slope, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else {
 if(argv[i][0] == '-') {
   printerror(application.verbose_state.debug, "ERROR pfold: Unknown option: %s\nRun 'pfold' without command line arguments to show help", argv[i]);
   return 0;
 }
 else {
   if(applicationAddFilename(i, application.verbose_state) == 0) {
     printerror(application.verbose_state.debug, "ERROR pfold: applicationAddFilename() failed", argv[i]);
     terminateApplication(&application);
     return 0;
   }
 }
      }
    }
  }
  char *filename_ptr;
  filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state);
  if(filename_ptr == NULL) {
    printerror(application.verbose_state.debug, "ERROR pfold: Cannot obtain next file name to process.");
    return 0;
  }
  for(i = 0; i < MaxNrPolarizations; i++)
    cleanPSRData(&fin[i], application.verbose_state);
  if(application.iformat <= 0) {
    application.iformat = guessPSRData_format(argv[argc-1], 0, application.verbose_state);
    if(application.iformat == -2 || application.iformat == -3)
      return 0;
  }
  if(isValidPSRDATA_format(application.iformat) == 0) {
    printerror(application.verbose_state.debug, "ERROR pfold: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", filename_ptr);
    return 0;
  }
  closePSRData(&fin[0], 0, application.verbose_state);
  if(!openPSRData(&fin[0], argv[argc-1], application.iformat, 0, 1, 0, application.verbose_state))
    return 0;
  if(PSRDataHeader_parse_commandline(&fin[0], argc, argv, application.verbose_state) == 0)
    return 0;
  for(i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-header") == 0) {
      printwarning(application.verbose_state.debug, "WARNING pfold: If using the -header option, be aware it applied BEFORE the preprocessing.");
    }
  }
  if(preprocessApplication(&application, &fin[0]) == 0) {
    return 0;
  }
  double period;
  int ret;
  ret = get_period(fin[0], 0, &period, application.verbose_state);
  if(ret == 2) {
    printerror(application.verbose_state.debug, "ERROR pfold (%s): Cannot obtain period", fin[0].filename);
    return 0;
  }
  if(period < 0.001) {
    printerror(application.verbose_state.debug, "ERROR pfold: The period does not appear to be set in the header. Consider using the -header option.");
    return 0;
  }
  if(get_tsamp(fin[0], 0, application.verbose_state) < 0.0000001 || get_tsamp(fin[0], 0, application.verbose_state) >= 10) {
    printerror(application.verbose_state.debug, "ERROR pfold: The sampling time does not appear to be set correctly in the header. Consider using the -header option.");
    return 0;
  }
  if(fin[0].isDebase == 0) {
    printerror(application.verbose_state.debug, "ERROR pfold: Baseline is not subtracted. Use pmod -debase first.");
    return 0;
  }else if(fin[0].isDebase != 1) {
    printwarning(application.verbose_state.debug, "WARNING pfold:  It is not known if baseline is already subtracted. Use pmod -debase first.");
  }
  if(check_baseline_subtracted(fin[0], application.verbose_state) == 0) {
    printwarning(application.verbose_state.debug, "WARNING pfold: Baseline does not appear to be subtracted. Use pmod -debase first.");
  }
  region_frac_to_int(&(application.onpulse), fin[0].NrBins, 0);
  if(p3_fold_flag) {
    if(fin[0].NrPols > MaxNrPolarizations) {
      printerror(application.verbose_state.debug, "ERROR pfold: Maximum supported input parameters is exceeded.\n");
      return 0;
    }
    originalNrPols = fin[0].NrPols;
    if(fin[0].NrPols > 0) {
      for(i = fin[0].NrPols-1; i >= 0; i--) {
 if(preprocess_polselect(fin[0], &clone, i, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR pfold: Error selecting Stokes parameter %ld.", i);
   return 0;
 }
 swap_orig_clone(&fin[i], &clone, application.verbose_state);
      }
    }
  }
  cleanPSRData(&fout, application.verbose_state);
  copy_params_PSRData(fin[0], &fout, application.verbose_state);
  if(p3_fold_flag) {
    profileI = (float *)malloc(fin[0].NrBins*sizeof(float));
    if(profileI == NULL) {
      printerror(application.verbose_state.debug, "ERROR pfold: Cannot allocate memory");
      return 0;
    }
    read_profilePSRData(fin[0], profileI, NULL, 0, noverbose);
    xmin = 0;
    ret = get_period(fin[0], 0, &period, application.verbose_state);
    if(ret == 2) {
      printerror(application.verbose_state.debug, "ERROR pfold (%s): Cannot obtain period", fin[0].filename);
      return 0;
    }
    xmax = 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period;
    if(application.onpulse.nrRegions == 0 || selectMoreOnpulseRegions) {
      strcpy(pgplot_options.viewport.plotDevice, onpulseselectdevice);
      strcpy(pgplot_options.box.xlabel, "Bin");
      strcpy(pgplot_options.box.ylabel, "Intensity");
      strcpy(pgplot_options.box.title, "Select on-pulse region");
      selectRegions(profileI, fin[0].NrBins, &pgplot_options, 0, 0, 1, &application.onpulse, application.verbose_state);
    }else {
      if(strcmp(onpulseselectdevice, "?") == 0)
 printf("Specify plotting device to show the profile showing the selected regions: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, onpulseselectdevice);
      strcpy(pgplot_options.box.xlabel, "Phase[deg]");
      strcpy(pgplot_options.box.ylabel, "Intensity");
      strcpy(pgplot_options.box.title, fin[0].psrname);
      if(pgplotGraph1(&pgplot_options, profileI, NULL, NULL, fin[0].NrBins, xmin, xmax, 0, xmin, xmax, 0, 0, 0, 1, 0, 1, 0, 1, 1, &application.onpulse, -1, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pfold: Unable to open plotdevice.\n");
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
     printerror(application.verbose_state.debug, "ERROR pfold: region not defined in bins");
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
  }
  if(p3_fold_flag) {
    p3foldmap = malloc(originalNrPols*fin[0].NrBins * p3_fold_nbin * sizeof(float));
    if(p3foldmap == NULL) {
      printerror(application.verbose_state.debug, "ERROR pfold: Cannot allocate memory");
      return 0;
    }
    originalNrPolsP3 = originalNrPols;
      if(originalNrPols > 1) {
 printwarning(application.verbose_state.debug, "WARNING pfold: Only first polarization chanel is folded.");
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
    if(strcmp(p3fold_device, "?") == 0) {
      printf("Specify plotting device to show the P3 fold map: \n  ");
    }
    strcpy(pgplot_options.viewport.plotDevice, p3fold_device);
    strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
    strcpy(pgplot_options.box.ylabel, "P3 [pulse periods]");
    strcpy(pgplot_options.box.title, "P3 fold");
    ret = get_period(fin[0], 0, &period, application.verbose_state);
    if(ret == 2) {
      printerror(application.verbose_state.debug, "ERROR pfold (%s): Cannot obtain period", fin[0].filename);
      return 0;
    }
    pgplotMap(&pgplot_options, p3foldmap, fin[0].NrBins, p3_fold_nbin, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period, xmin_zoom, xmax_zoom, 0.5*p3_fold/(float)p3_fold_nbin, 0.5*p3_fold/(float)p3_fold_nbin + p3_fold*(p3_fold_nbin-1)/(float)p3_fold_nbin, 0, p3_fold, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, application.verbose_state);
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
 printerror(application.verbose_state.debug, "ERROR pfold: Memory allocation error");
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
 printerror(application.verbose_state.debug, "ERROR pfold: Unable to open file for writing.\n");
 return 0;
      }
      if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pfold: Unable to write header.\n");
 return 0;
      }
      p3foldmap2 = malloc(originalNrPols*fin[0].NrBins * p3_fold_nbin * sizeof(float));
      if(p3foldmap2 == NULL) {
 printerror(application.verbose_state.debug, "ERROR pfold: Cannot allocate memory");
 return 0;
      }
      for(i = 0; i < p3_fold_nbin; i++) {
 for(p = 0; p < originalNrPolsP3; p++) {
   memcpy(&p3foldmap2[(originalNrPolsP3*i+p)*fin[0].NrBins], &p3foldmap[(p*p3_fold_nbin+i)*fin[0].NrBins], fin[0].NrBins*sizeof(float));
 }
      }
      if(writePSRData(&fout, p3foldmap2, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pfold: Unable to write data.\n");
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
  if(p3_fold_flag) {
    free(profileI);
  }
  ppgend();
  terminateApplication(&application);
  return 0;
}

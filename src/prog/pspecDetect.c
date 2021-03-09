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
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "psrsalsa.h"
double f2_min, f2_max, f3_min, f3_max, fl_min, fl_max, f2n_min, f2n_max, f3n_min, f3n_max;
double twodsf_Imin;
double I_noise_max, I_noise_min;
float centroid_x, centroid_y;
float centroid_err_x, centroid_err_y;
int centroid_calculated;
int SelectedPlot;
int PlotAvrgProfile;
int SelectedComponent;
int Centered;
int KeyCode;
datafile_definition twodfs, lrfs, AverageProfile, noise;
#define Max_nr_noise_patches 100
#define MaxNrNoiseBins 100
int f2npatch_min[Max_nr_noise_patches],f2npatch_max[Max_nr_noise_patches],f3npatch_min[Max_nr_noise_patches],f3npatch_max[Max_nr_noise_patches];
int nr_noise_patches;
void PlotWindow(verbose_definition verbose);
int MakeSelection(double sigma_noise);
void SelectRegion();
double calculate_noise_sigma(void);
void calculate_2dfs_Centroid(double sigma_noise, FILE *writesel_fptr, int SelectedComponent, verbose_definition verbose);
int main(int argc, char **argv)
{
  int i, xi, yi, index, extprefix;
  double I;
  char filename[MaxFilenameLength+1], txt[MaxFilenameLength+1], *input_filename_ptr;
  double sigma_noise;
  psrsalsaApplication application;
  initApplication(&application, "pspecDetect", "[options] pulse_stack,\nwhere pulse_stack is the file name of the pulse stack that has been processed by\npspec to produce the 2DFS and LRFS.");
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_device = 1;
  cleanPSRData(&twodfs, application.verbose_state);
  cleanPSRData(&lrfs, application.verbose_state);
  cleanPSRData(&AverageProfile, application.verbose_state);
  cleanPSRData(&noise, application.verbose_state);
  closePSRData(&twodfs, 0, application.verbose_state);
  closePSRData(&lrfs, 0, application.verbose_state);
  closePSRData(&AverageProfile, 0, application.verbose_state);
  closePSRData(&noise, 0, application.verbose_state);
  SelectedPlot = 0;
  sigma_noise = 0;
  centroid_calculated = 0;
  nr_noise_patches = 0;
  SelectedComponent = 1;
  Centered = 1;
  PlotAvrgProfile = 1;
  extprefix = 0;
  if(argc < 2) {
    printf("Interactive program designed to analyse features in the 2DFS to obtain\ncentroid P2 and P3 values and corresponding error-bars.\n\n");
    printApplicationHelp(&application);
    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about the how to use the centroid information can be found in:\n");
    printf(" - Weltevrede et al. 2006, A&A, 445, 243\n");
    printf(" - Weltevrede et al. 2007, A&A, 469, 607.\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
 i = index;
      }else {
 if(argv[i][0] == '-') {
   printerror(application.verbose_state.debug, "pspecDetect: Unknown option: %s\n\nRun pspecDetect without command line arguments to show help", argv[i]);
   terminateApplication(&application);
   return 0;
 }else {
   if(applicationAddFilename(i, application.verbose_state) == 0)
     return 0;
 }
      }
    }
  }
  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pspecDetect: No files specified");
    return 0;
  }
  FILE *writesel_fptr;
  writesel_fptr = NULL;
  input_filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state);
  if(extprefix == 0) {
    sprintf(txt, "%d.2dfs", 1);
  }
  if(change_filename_extension(input_filename_ptr, filename, txt, 1000, application.verbose_state) == 0) {
    return 0;
  }
  if(application.verbose_state.verbose)
    printf("Reading %s\n", filename);
  if(!openPSRData(&twodfs, filename, 0, 0, 1, 0, application.verbose_state))
    return 0;
  if(twodfs.NrPols > 1) {
    datafile_definition clone;
    if(preprocess_polselect(twodfs, &clone, 0, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "Cannot select first polarization channel\n");
      return 0;
    }
    swap_orig_clone(&twodfs, &clone, application.verbose_state);
  }
  twodsf_Imin = 0;
  for(xi = 0; xi < twodfs.NrBins; xi++) {
    for(yi = 0; yi < twodfs.NrSubints; yi++) {
      I = twodfs.data[yi*twodfs.NrBins+xi];
      if(I < twodsf_Imin)
 twodsf_Imin = I;
    }
  }
  if(application.verbose_state.verbose)
    printf("%ldx%ld points read from 2dfs\n", twodfs.NrBins, twodfs.NrSubints);
  if(preprocess_polselect(twodfs, &noise, 0, application.verbose_state) != 1)
    exit(0);
  if(extprefix == 0) {
    sprintf(txt, "lrfs");
  }
  if(change_filename_extension(input_filename_ptr, filename, txt, 1000, application.verbose_state) == 0) {
    return 0;
  }
  if(application.verbose_state.verbose)
    printf("Reading %s\n", filename);
  if(!openPSRData(&lrfs, filename, 0, 0, 1, 0, application.verbose_state))
    return 0;
  if(application.verbose_state.verbose)
    printf("%ldx%ld points read from lrfs\n", lrfs.NrBins, lrfs.NrSubints);
  if(lrfs.NrPols > 1) {
    datafile_definition clone;
    if(preprocess_polselect(lrfs, &clone, 0, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "Cannot select first polarization channel\n");
      return 0;
    }
    swap_orig_clone(&lrfs, &clone, application.verbose_state);
  }
  fl_min = 0;
  fl_max = lrfs.NrBins-1;
  if(application.verbose_state.verbose)
    printf("Reading %s\n", filename);
  if(!openPSRData(&AverageProfile, input_filename_ptr, 0, 0, 0, 0, application.verbose_state))
    return 0;
  if(!readHeaderPSRData(&AverageProfile, 0, 0, application.verbose_state))
    return 0;
  AverageProfile.data = malloc(AverageProfile.NrBins*sizeof(float));
  if(AverageProfile.data == NULL) {
    printerror(application.verbose_state.debug, "Memory allocation error\n");
    return 0;
  }
  if(!read_profilePSRData(AverageProfile, AverageProfile.data, NULL, 0, application.verbose_state))
    return 0;
  if(AverageProfile.NrBins != lrfs.NrBins) {
    printwarning(application.verbose_state.debug, "WARNING: It looks like data is rebinned? Check the units.");
    datafile_definition clone;
    AverageProfile.format = MEMORY_format;
    AverageProfile.NrPols = 1;
    AverageProfile.NrFreqChan = 1;
    AverageProfile.NrSubints = 1;
    if(preprocess_rebin(AverageProfile, &clone, lrfs.NrBins, application.verbose_state) == 0) {
      printwarning(application.verbose_state.debug, "WARNING: Rebinning of profile failed.");
      return 0;
    }
    swap_orig_clone(&AverageProfile, &clone, application.verbose_state);
    printwarning(application.verbose_state.debug, "WARNING: Assuming the number of bins = %ld and the sampling time = %lf s.", AverageProfile.NrBins, AverageProfile.fixedtsamp);
  }
  if(application.verbose_state.verbose)
    printf("%ld points read from avgprof\n\n", AverageProfile.NrBins);
  if(twodfs.xrangeset) {
    f2_min = f2n_min = twodfs.xrange[0];
    f2_max = f2n_max = twodfs.xrange[1];
  }else {
    f2_min = f2n_min = -AverageProfile.NrBins/2.0;
    f2_max = f2n_max = -AverageProfile.NrBins/2.0 + AverageProfile.NrBins*(twodfs.NrBins-1.0)/(float)twodfs.NrBins;
  }
  f3_min = 0;
  f3n_min = 0;
  f3_max = f3n_max = 0.5;
  ppgopen(application.pgplotdevice);
  ppgask(0);
  ppgslw(1);
    printf("Press h for help\n");
  SelectRegion();
  PlotWindow(application.verbose_state);
  do {
      MakeSelection(sigma_noise);
    switch(KeyCode) {
    case 0:
      PlotWindow(application.verbose_state);
      break;
    case 113:
    case 81:
    case 27: KeyCode = 27; break;
    case 13:
 sigma_noise = calculate_noise_sigma();
      calculate_2dfs_Centroid(sigma_noise, writesel_fptr, SelectedComponent, application.verbose_state);
      PlotWindow(application.verbose_state);
      break;
    case 67:
    case 99:
      if(Centered == 0) {
 printf("P2 Centring on\n");
 Centered = 1;
      }else {
 Centered = 0;
 printf("P2 Centring off\n");
      }
      printf("\nBe aware: P2 centring makes it less likely to detect a significant P2 offset just because the region included in the centroid calculation is offset. However, for clearly significant features you might want to turn P2 centring off to ensure the region included in the centroid calculation is centred on the feature, avoiding the corresponding bias.\n");
      break;
    case 82:
    case 114:
      sigma_noise = 0;
      centroid_calculated = 0;
      if(twodfs.xrangeset) {
 f2_min = f2n_min = twodfs.xrange[0];
 f2_max = f2n_max = twodfs.xrange[1];
      }else {
 f2_min = f2n_min = -AverageProfile.NrBins/2.0;
 f2_max = f2n_max = -AverageProfile.NrBins/2.0 + AverageProfile.NrBins*(twodfs.NrBins-1.0)/(float)twodfs.NrBins;
      }
      f3_min = 0;
      f3n_min = 0;
      f3_max = f3n_max = 0.5;
      fl_min = 0;
      fl_max = lrfs.NrBins-1;
      if(noise.opened_flag)
        closePSRData(&noise, 0, application.verbose_state);
      cleanPSRData(&noise, application.verbose_state);
      if(preprocess_polselect(twodfs, &noise, 0, application.verbose_state) != 1)
        exit(0);
      nr_noise_patches = 0;
      SelectRegion();
      PlotWindow(application.verbose_state);
      break;
    case 70:
    case 102:
      if(twodfs.xrangeset) {
 f2_min = twodfs.xrange[0];
 f2_max = twodfs.xrange[1];
      }else {
 f2_min = -AverageProfile.NrBins/2.0;
 f2_max = -AverageProfile.NrBins/2.0 + AverageProfile.NrBins*(twodfs.NrBins-1.0)/(float)twodfs.NrBins;
      }
      f3_min = 0;
      f3_max = 0.5;
      fl_min = 0;
      fl_max = lrfs.NrBins-1;
      SelectRegion();
      PlotWindow(application.verbose_state);
      break;
    case 49:
    case 50:
    case 51:
    case 52:
    case 53:
    case 54:
    case 55:
    case 56:
    case 57:
      SelectedComponent = KeyCode - 48;
      strcpy(filename, argv[argc - 1]);
      if(extprefix == 0) {
 sprintf(txt, "%d.2dfs", SelectedComponent);
      }
      if(change_filename_extension(input_filename_ptr, filename, txt, 1000, application.verbose_state) == 0) {
 return 0;
      }
      if(application.verbose_state.verbose)
 printf("Reading %s\n", filename);
      if(closePSRData(&twodfs, 0, application.verbose_state)) {
 printerror(0, "Closing file failed\n");
 return 0;
      }
      if(!openPSRData(&twodfs, filename, 0, 0, 1, 0, application.verbose_state))
 return 0;
      if(twodfs.NrPols > 1) {
 datafile_definition clone;
 if(preprocess_polselect(twodfs, &clone, 0, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "Cannot select first polarization channel\n");
   return 0;
 }
 swap_orig_clone(&twodfs, &clone, application.verbose_state);
      }
      if(application.verbose_state.verbose)
        printf("%ldx%ld points read from 2dfs\n", twodfs.NrBins, twodfs.NrSubints);
      twodsf_Imin = 0;
      for(xi = 0; xi < twodfs.NrBins; xi++) {
 for(yi = 0; yi < twodfs.NrSubints; yi++) {
   I = twodfs.data[yi*twodfs.NrBins+xi];
   if(I < twodsf_Imin)
     twodsf_Imin = I;
 }
      }
      sigma_noise = 0;
      centroid_calculated = 0;
      nr_noise_patches = 0;
      if(twodfs.xrangeset) {
 f2_min = f2n_min = twodfs.xrange[0];
 f2_max = f2n_max = twodfs.xrange[1];
      }else {
 f2_min = f2n_min = -AverageProfile.NrBins/2.0;
 f2_max = f2n_max = -AverageProfile.NrBins/2.0 + AverageProfile.NrBins*(twodfs.NrBins-1.0)/(float)twodfs.NrBins;
      }
      f3_min = 0;
      f3n_min = 0;
      f3_max = f3n_max = 0.5;
      if(noise.opened_flag)
        closePSRData(&noise, 0, application.verbose_state);
      cleanPSRData(&noise, application.verbose_state);
      if(preprocess_polselect(twodfs, &noise, 0, application.verbose_state) != 1)
        exit(0);
      SelectRegion();
      PlotWindow(application.verbose_state);
      break;
    case 32:
      SelectedPlot++;
      if(SelectedPlot > 2
  )
 SelectedPlot = 0;
      SelectRegion();
      PlotWindow(application.verbose_state);
      break;
    case 97:
      PlotAvrgProfile++;
      if(PlotAvrgProfile > 1)
 PlotAvrgProfile = 0;
      SelectRegion();
      PlotWindow(application.verbose_state);
      break;
    case 72:
    case 104:
    case 63:
      printf("\nGeneral strategy: Use space to show the 2dfs with the title \"For noise calculation ...\". Flag all signal only noise is visible. This sets the rms used in the error calculation. Use space to select the 2dfs and zoom in on the feature of interest. Depending on the situation, toggling 'c' would be beneficial (see help when using this option). Press return to calculate centroid of the zoomed in area (indicated by cross) and 1-sigma error box. The error is purely statistical, the systematic error resulting from the decision of what area to include in the centroid calculation can be assessed by selecting slightly different regions and repeat the calculation.\n\n");
      printf("h/?   = Help\n");
      printf("1..9  = Switch to component 1..9 (load new 2dfs)\n");
      printf("R     = Reset (zoom + noise-flagging settings)\n");
      printf("F     = Select new feature (but keep noise-flagging settings)\n");
      printf("C     = Toggle P2 centering (allows selection to be non-symmetric)\n");
      printf("A     = Toggle showing average profile superimposed over the LRFS\n");
      printf("SPACE = Switch between the 2dfs, lrfs, and the 2dfs samples included in the noise calculation\n");
      printf("ENTER = Calculate centroid\n");
      printf("Q/ESC = Quit\n");
      break;
    default: printf("Unknown key: %d\n", KeyCode); break;
    }
  }while(KeyCode != 27
  );
  ppgend();
  closePSRData(&twodfs, 0, application.verbose_state);
  closePSRData(&lrfs, 0, application.verbose_state);
  closePSRData(&AverageProfile, 0, application.verbose_state);
  if(noise.opened_flag)
    closePSRData(&noise, 0, application.verbose_state);
  terminateApplication(&application);
  return 0;
}
void SelectRegion()
{
  int i, xi, yi, offset;
  if(SelectedPlot == 2 || SelectedPlot == 4) {
    if(nr_noise_patches > 0) {
      for(i = 0; i < nr_noise_patches; i++) {
 for(xi = f2npatch_min[i]; xi <= f2npatch_max[i]; xi++) {
   for(yi = f3npatch_min[i]; yi <= f3npatch_max[i]; yi++) {
     offset = yi*noise.NrBins+xi;
     if(offset < 0 || offset >= noise.NrBins*noise.NrSubints) {
       printerror(0, "Bug!!\n");
     }else {
       noise.data[offset] = -10*fabs(twodsf_Imin);
     }
   }
 }
      }
    }
  }
}
void PlotWindow(verbose_definition verbose)
{
  char txt[100];
  double Imax;
  double I;
  int i, xi, yi;
    ppgpage();
    if(SelectedPlot == 0) {
      ppgsvp(0.1, 0.9, 0.1, 0.9);
      ppgswin(f2_min,f2_max,f3_min,f3_max);
      sprintf(txt, "2dfs feature component %d", SelectedComponent);
      pgplot_options_definition pgplot_options;
      pgplot_clear_options(&pgplot_options);
      pgplot_options.viewport.dontopen = 1;
      pgplot_options.viewport.dontclose = 1;
      strcpy(pgplot_options.box.xlabel, "Fluctuation frequency (cycles/period)");
      strcpy(pgplot_options.box.ylabel, "Fluctuation frequency (cycles/period)");
      strcpy(pgplot_options.box.title, txt);
      float xrange0, xrange1;
      if(twodfs.xrangeset) {
 xrange0 = twodfs.xrange[0];
 xrange1 = twodfs.xrange[1];
      }else {
 xrange0 = -AverageProfile.NrBins/2.0;
 xrange1 = -AverageProfile.NrBins/2.0 + AverageProfile.NrBins*(twodfs.NrBins-1.0)/(float)twodfs.NrBins;
      }
      pgplotMap(&pgplot_options, twodfs.data, twodfs.NrBins, twodfs.NrSubints, xrange0, xrange1, f2_min, f2_max, 0, 0.5, f3_min, f3_max, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, verbose);
      if(centroid_calculated) {
 ppgsci(2);
 ppgslw(1);
 ppgmove(centroid_x,0.5);
 ppgdraw(centroid_x, -0.5);
        ppgmove(f2_min,centroid_y);
        ppgdraw(f2_max,centroid_y);
        ppgmove(centroid_x-centroid_err_x,centroid_y-centroid_err_y);
        ppgdraw(centroid_x+centroid_err_x,centroid_y-centroid_err_y);
        ppgdraw(centroid_x+centroid_err_x,centroid_y+centroid_err_y);
        ppgdraw(centroid_x-centroid_err_x,centroid_y+centroid_err_y);
        ppgdraw(centroid_x-centroid_err_x,centroid_y-centroid_err_y);
        ppgsci(1);
      }
    }else if(SelectedPlot == 1) {
      ppgsvp(0.1, 0.9, 0.1, 0.9);
      ppgswin(fl_min,fl_max,f3_min,f3_max);
      pgplot_options_definition pgplot_options;
      pgplot_clear_options(&pgplot_options);
      pgplot_options.viewport.dontopen = 1;
      pgplot_options.viewport.dontclose = 1;
      strcpy(pgplot_options.box.xlabel, "bins");
      strcpy(pgplot_options.box.ylabel, "Fluctuation frequency (cycles/period)");
      strcpy(pgplot_options.box.title, "lrfs feature");
      pgplotMap(&pgplot_options, lrfs.data, lrfs.NrBins, lrfs.NrSubints, 0, lrfs.NrBins-1, fl_min, fl_max, 0, 0.5, f3_min, f3_max, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, verbose);
      if(PlotAvrgProfile != 0) {
        i = 0;
 Imax = AverageProfile.data[0];
 for(xi=0; xi < AverageProfile.NrBins; xi++) {
   if(AverageProfile.data[xi] > Imax)
     Imax = AverageProfile.data[xi];
 }
 for(xi=0; xi < AverageProfile.NrBins; xi++) {
   ppgsci(2
);
   if(xi >= fl_min && xi <= fl_max) {
     I = AverageProfile.data[xi]*(f3_max-f3_min)/Imax + f3_min;
     if(I < f3_min)
       I = f3_min;
     if(i == 0) {
       ppgmove(xi,I);
       i = 1;
     }else {
       ppgdraw(xi,I);
     }
   }
 }
 ppgsci(1);
      }
      if(centroid_calculated) {
        ppgsci(2);
        ppgmove(0,centroid_y);
        ppgdraw(AverageProfile.NrBins,centroid_y);
        ppgsci(1);
      }
    }else if(SelectedPlot == 2) {
      for(xi = 0; xi < noise.NrBins; xi++) {
        for(yi = 0; yi < noise.NrSubints; yi++) {
          I = noise.data[yi*noise.NrBins+xi];
          if(I < twodsf_Imin)
            noise.data[yi*noise.NrBins+xi] = 0;
        }
      }
      pgplot_options_definition pgplot_options;
      pgplot_clear_options(&pgplot_options);
      pgplot_options.viewport.dontopen = 1;
      pgplot_options.viewport.dontclose = 1;
      strcpy(pgplot_options.box.xlabel, "Fluctuation frequency (cycles/period)");
      strcpy(pgplot_options.box.ylabel, "Fluctuation frequency (cycles/period)");
      strcpy(pgplot_options.box.title, "For noise calulation: All signal should be flagged in this 2dfs plot");
      pgplotMap(&pgplot_options, noise.data, noise.NrBins, noise.NrSubints,
  -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)noise.NrBins, +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)noise.NrBins, f2_min, f2_max, 0, 0.5, f3_min, f3_max, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, verbose);
      SelectRegion();
    }
}
int MakeSelection(double sigma_noise)
{
  float x0, x1, y0, y1, dummy;
  char c;
  x0 = y0 = x1 = y1 = c = 0;
  KeyCode = 0;
  ppgsci(2);
  ppgband(0, 0, 0.0, 0.0, &x0, &y0, &c);
  if(c != 65) {
    ppgsci(1);
    KeyCode = c;
    return 0;
  }
  if(SelectedPlot == 3) {
    y0 = x0;
    x0 = (y0 - f3_min)/(f3_max-f3_min)*(f2_max-f2_min)+f2_min;
    printf("P3[P0]  = %lf\n", 1/y0);
    printf("P2[deg] = %lf\n", 360/x0);
    KeyCode = 0;
    ppgsci(1);
    return 0;
  }
  ppgband(2, 0, x0, y0, &x1, &y1, &c);
  ppgsci(1);
  if(c != 65)
    return 0;
  if(y0 > y1) {
    dummy = y0;
    y0 = y1;
    y1 = dummy;
  }
  if(x0 > x1) {
    dummy = x0;
    x0 = x1;
    x1 = dummy;
  }
  if(SelectedPlot == 0) {
    if(Centered) {
      if(fabs(x0) > fabs(x1)) {
 x1 = fabs(x0);
 x0 = -fabs(x0);
      }else {
 x0 = -fabs(x1);
 x1 = fabs(x0);
      }
    }
    f2_min = x0;
    f2_max = x1;
    f3_min = y0;
    f3_max = y1;
    int nx, ny;
    float dx, dy;
    pgplotMapCoordinate_dbl(f2_min, f3_min, &nx, &ny);
    pgplotMapCoordinateInverse_dbl(&f2_min, &f3_min, nx, ny);
    pgplotMapCoordinate_dbl(f2_max, f3_max, &nx, &ny);
    pgplotMapCoordinateInverse_dbl(&f2_max, &f3_max, nx, ny);
    pgplotMapCoordinateBinSize(&dx, &dy);
    f2_min -= 0.49*dx;
    f2_max += 0.49*dx;
    f3_min -= 0.49*dy;
    f3_max += 0.49*dy;
  }else if(SelectedPlot == 1) {
    fl_min = x0;
    fl_max = x1;
    f3_min = y0;
    f3_max = y1;
  }else if(SelectedPlot == 2) {
    if(nr_noise_patches == Max_nr_noise_patches) {
      printf("Too many patches\n");
      nr_noise_patches--;
    }
    pgplotMapCoordinate(x0, y0, &(f2npatch_min[nr_noise_patches]), &(f3npatch_min[nr_noise_patches]));
    pgplotMapCoordinate(x1, y1, &(f2npatch_max[nr_noise_patches]), &(f3npatch_max[nr_noise_patches]));
    nr_noise_patches++;
    SelectRegion();
      sigma_noise = calculate_noise_sigma();
  }
  return 1;
}
double calculate_noise_sigma(void)
{
  double I, sigma_noise;
  int xi, yi;
  int nrpoints_flagged;
  sigma_noise = 0;
  nrpoints_flagged = 0;
  for(yi = 0; yi < noise.NrSubints; yi++) {
    for(xi = 0; xi < noise.NrBins; xi++) {
      I = noise.data[yi*noise.NrBins+xi];
      if(I < twodsf_Imin)
        nrpoints_flagged++;
      else
        sigma_noise += I*I;
    }
  }
  sigma_noise = sqrt(sigma_noise/(double)(noise.NrBins*noise.NrSubints-nrpoints_flagged));
  printf("Sigma noise = %e (%d points flagged)\n", sigma_noise, nrpoints_flagged);
  if(nrpoints_flagged == 0) {
    printf("\nNo signal is flagged, so the current error-bar is based on the rms of the noise + pulsar signal, overestimating the actual error. Press 'h' for a general help.\n\n");
  }
  return sigma_noise;
}
void calculate_2dfs_Centroid(double sigma_noise, FILE *writesel_fptr, int SelectedComponent, verbose_definition verbose)
{
  double I, Itot;
  int xi, yi, xstart, xend, ystart, yend, nrbins;
  double x, y, xcent, ycent, xerr, yerr;
  float xf, yf, binsizex, binsizey;
  if(SelectedPlot == 0) {
    pgplotMapCoordinate(f2_min, f3_min, &xstart, &ystart);
    pgplotMapCoordinate(f2_max, f3_max, &xend, &yend);
    if(verbose.verbose)
      printf("Centroid:                            ");
    xcent = ycent = Itot = 0;
    nrbins = 0;
    for(yi = ystart; yi <= yend; yi++) {
      for(xi = xstart; xi <= xend; xi++) {
 I = fabs(twodfs.data[yi*twodfs.NrBins+xi]);
 pgplotMapCoordinateInverse(&xf, &yf, xi, yi);
 x = xf;
 y = yf;
 Itot += I;
 xcent += I*x;
 ycent += I*y;
 nrbins++;
      }
    }
    xcent /= Itot;
    ycent /= Itot;
    if(verbose.verbose) {
      printf("(%lf, %lf) cpp", xcent, ycent);
      printf(" based on %d selected bins\n", nrbins);
      printf("Statistical error caused by noise:   ");
    }
    xerr = yerr = 0;
    for(yi = ystart; yi <= yend; yi++) {
      for(xi = xstart; xi <= xend; xi++) {
 I = twodfs.data[yi*twodfs.NrBins+xi];
 pgplotMapCoordinateInverse(&xf, &yf, xi, yi);
 x = xf;
 y = yf;
 xerr += (x - xcent)*(x - xcent);
 yerr += (y - ycent)*(y - ycent);
      }
    }
    xerr *= sigma_noise*sigma_noise/(Itot*Itot);
    yerr *= sigma_noise*sigma_noise/(Itot*Itot);
    pgplotMapCoordinateBinSize(&binsizex, &binsizey);
    xerr = sqrt(xerr);
    yerr = sqrt(yerr);
    if(verbose.verbose) {
      printf("(%e, %e) cpp\n", xerr, yerr);
      printf("Half bin size:                       (%e, %e) cpp\n", 0.5*binsizex, 0.5*binsizey);
    }
    centroid_x = xcent;
    centroid_y = ycent;
    centroid_calculated = 1;
    centroid_err_x = xerr;
    centroid_err_y = yerr;
      printf("P3[P0]  = %lf +- %lf\n", 1/ycent, (centroid_err_y)/(ycent*ycent));
      printf("P2[deg] = %lf +- %lf\n", 360/xcent, 360*(centroid_err_x)/(xcent*xcent));
  }
}

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
#include "psrsalsa.h"

int SelectP2Region, SelectP3Integrate, SelectLRegion, doFlip, ExtraVerticalMaximaSkip, noside, overlaypp, inside, noylabels, noxlabels, notop;
double f2_min, f2_max, f3_min, f3_max, fl_min, fl_max, oversaturize, oversaturize2, oversaturizel, dl;
double f2_min2, f2_max2, f3_min2, f3_max2, Imin, Imax;
double P3IntegrateLow, P3IntegrateHigh, P2RegionLow, P2RegionHigh, LRegionLow, LRegionHigh;
double P3IntegrateLow2, P3IntegrateHigh2, P2RegionLow2, P2RegionHigh2;
float labelscale;
char title2[100];
datafile_definition twodfs, twodfs2, lrfs, AverageProfile, VarianceProfile, ModProfile, VarianceProfileErr, ModProfileErr;

void Plot2dfs(int Number, int plot_xlabel, int plot_ylabel, int noside, int plot_ylabeltop, int nomain, double scaleFig_x, double scaleFig_y, int showwedge, int normaliseSide, int nointegrateNumbers, verbose_definition verbose);
void Plot2dfsOnly(int Number, int plot_xlabel, int plot_ylabel, int noside, int nomain, int showwedge, char *title, double scaleFig_x, double scaleFig_y, int normaliseSide, int nointegrateNumbers, verbose_definition verbose);
void PlotLRFS(int plot_xlabel, int plot_ylabel, int noside, int plot_ylabeltop, int lineColour, int showwedge, char *title, double scaleFig_x, double scaleFig_y, int normaliseSide, int nointegrateNumbers, int usephase, verbose_definition verbose);

int main(int argc, char **argv)
{
  int i, j, xi, SelectP3Region, LoadTwo, NrSelectedOverSaturize, plot_xlabel, plot_ylabel, nomain, normaliseSide, nointegrateNumbers;
  int plot_ylabeltop, file_number, ImaxSet, IminSet, type_of_plots, maxSubpulsePhaseSet, minSubpulsePhaseSet, do_phase_slope, showwedge, normalise_spectra;
  char PlotDevice[100], filename[1000], txt[1000];
  int title_index, Load2dfs, LoadLRFS, plotvariance, plotmodindex, ok_flag, altProf, lineStyle, lineColour, ret, usephase;
  double P3RegionLow, P3RegionHigh, maxSubpulsePhase, minSubpulsePhase, scaleFig_x, scaleFig_y;
  float I, x, profilescale;
  float maxvalue_mod, maxvalue_stddev, maxsigma_stddev, maxsigma_mod, ImaxValue, IminValue, phase_slope_g, phase_slope_o;
  FILE *fout_ascii;
  long k;
  datafile_definition clone, subpulseTrackProfile, subpulseTrackProfileErr, subpulseAmpProfile;
  verbose_definition noverbose;
  psrsalsaApplication application;

  initApplication(&application, "pspecFig", "");





  type_of_plots = 0;
  file_number = 1;

  scaleFig_x = 1;
  scaleFig_y = 1;
  do_phase_slope = 0;

  if(argc < 2) {
    printf("Program to plot some of the pspec output\n\nUsage: pspecFig [options] stack_file (i.e. the name of the pulse stack that has been processed by pspec). By default the profile/modulation index/standard deviation profile/lrfs/2dfs are combined in a single plot (mode A). When -phaseplot is specified, a plot of the profile/subpulse amplitude and subpulse phase is produced (mode B).\n\n");
    printf("Where optional options are:\n\n");
    printf("General options:\n");
    printf("-headerlist       Show options of the -header option.\n");
    printf("-header           Change this header parameter.\n");
    printf("                  Example -header 'name J0123+4567'\n");
    printf("-v                Verbose mode (to get a better idea what is happening)\n");
    printf("-debug            Enable more output (where implemented)\n");
    printf("-nocounters       Don't show counters etc (useful when generating log files)\n");
    printf("-phaseplot        Enable mode B operation, see above.\n");

    printf("\nGeneral file/panel options:\n");
    printf("-altProf          Load this alternative .profile file (generated by pspec).\n");
    printf("                  Useful if profile/modindex/stddev curves are calculated from\n");
    printf("                  a different file than the spectra.\n");
    printf("-notop            Don't show top plot\n");
    printf("\nMode A specific file/panel options:\n");
    printf("-2                Load two 2dfs's corresponding to two pulse longitude ranges\n");
    printf("-2dfsnr           Override 2dfs file number of first 2dfs. Default is %d.\n", file_number);
    printf("-nolrfs           Do not load lrfs\n");
    printf("-no2dfs           Do not load 2dfs's\n");
    printf("-noside           Don't show side panels\n");
    printf("-normspectra      Normalise the spectra (peak value=1)\n");
    printf("-normside         Normalise the side panels (peak value=1)\n");

    printf("\nGeneral range options:\n");
    printf("-Imax             Set maximum value of the y-range of the top plot\n");
    printf("-Imin             Set minimum value of the y-range of the top plot\n");
    printf("-l  \"low high\"    Set horizontal range shown lrfs/profile.\n");
    printf("-phase            Use pulse longitude in phase rather than degrees.\n");
    printf("-dl               Shift lrfs/profile by this amount of degrees or phase.\n");
    printf("\nMode A specific range options:\n");
    printf("-p3 \"low high\"    Set vertical range shown in lrfs/2dfs.\n");
    printf("-p2 \"low high\"    Set horizontal range shown in 2dfs.\n");
    printf("-int  \"low high\"  Select vertical integration range in 2dfs\n");
    printf("                  (affects the side panels).\n");
    printf("-modsigma         Set minimum significance for the modulation index (def. is 3).\n");
    printf("-stddevsigma      Same for the stddev values (default is 3).\n");
    printf("-modmax           Set the max allowed value for the modulation index\n");
    printf("                  (other values are ignored).\n");
    printf("-stddevmax        Same for the stddev values\n");
    printf("\nMode B specific range options:\n");
    printf("-spmax            Set maximum value of the y-range of the subpulse phase plot\n");
    printf("-spmin            Set minimum value of the y-range of the subpulse phase plot\n");
    printf("\nOther graphics options:\n");
    printf("-device  \"...\"    Plot device\n");
    printf("-inside           Place tick marks inside\n");
    printf("-labelscale       Set size of labels (default is 1)\n");
    printf("-noylabels        Don't show ylabels\n");
    printf("-linestyle        Set the PGPLOT line style of the pulse profile\n");
    printf("-linecolor        Set the PGPLOT line color of the pulse profile\n");
    printf("-scalefig \"x y\"   Scale size of panel with factors x and y\n");
    printf("-title            Set the title\n");
    printf("-ytop             Show y-label of top plot\n");
    printf("\nMode A specific graphics options:\n");
    printf("-intnrs           Show nrs along side panels axis, instead of just a tick\n");
    printf("-scalel           Specify scale, default is 1. This option multiplies the lrfs\n");
    printf("                  powers with this factor. This results in clipping, making\n");
    printf("                  weak features clearer\n");
    printf("-scale2           Similar to -scalel, but for 2dfs instead of lrfs\n\n");
    printf("                  When this option is used twice, the 2nd time the\n");
    printf("                  option is used applies to the second 2dfs shown.\n");
    printf("-scalep           Scale profile by this factor\n");
    printf("-f                Do not flip 2DFS horizontally. If specified positive drift\n");
    printf("                  corresponds to power in the left-hand side of the diagram.\n");
    printf("-overlay          Overlay pulse profile over LRFS\n");


    printf("-noxlabels        Don't show xlabels on 2dfs bottom integration panel.\n");
    printf("-nomod            Do not plot a modulation index profile\n");
    printf("-nostddev         Do not plot a standard deviation profile\n");
    printf("-showwedge        Plot an annotated wedge to show color scale\n");
    printf("-textside         Set the text printed in top left corner of the graph\n");
    printf("-xlabel           Show x-label of LRFS and 2DFS in cpp\n");
    printf("-xlabel2          Show x-label of LRFS and 2DFS in P0/P2\n");
    printf("-xlabel3          Show x-label of LRFS and 2DFS in P/P2\n");
    printf("-ylabel           Show y-label of LRFS in cpp\n");
    printf("-ylabel2          Show y-label of LRFS in P0/P3\n");
    printf("-ylabel3          Show y-label of LRFS in P/P3\n");
    printf("-s                Skip the specified number of P3 bins when determining\n");
    printf("                  the range in left-hand side integrations panels. By default\n");
    printf("                  the first bin is skipped.\n");
    printf("\nMode B specific graphics options:\n");
    printf("-f                Do not change sign of subpulse phase. If specified positive\n");
    printf("                  drift corresponds to a decrease subpulse phase as function of\n");
    printf("                  pulse longitude.\n");
    printf("-xlabel           Show x-label in pulse longitude\n");
    printf("-ylabel           Show y-label in subpulse phase\n");
    printf("-phaseslope \"g o\" In the subpulse phase plot, add a line with gradient g (in deg\n");
    printf("                  per deg) and offset o in deg.\n");
    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about the lrfs/2dfs/modulation index can be found in:\n");
    printf(" - Weltevrede et al. 2006, A&A, 445, 243\n");
    printf(" - Weltevrede et al. 2007, A&A, 469, 607\n");
    printf("More information about bootstrap/subpulse phase track & amplitude can be found in:\n");
    printf(" - Weltevrede et al. 2012, MNRAS, 424, 843\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }

  plot_ylabeltop = 0;
  nomain = 0;
  SelectP3Region = 0;
  SelectP2Region = 0;
  SelectLRegion = 0;
  SelectP3Integrate = 0;
  LoadTwo = 0;
  strcpy(PlotDevice, "?");
  oversaturize = 1;
  oversaturize2 = 1;
  oversaturizel = 1;
  NrSelectedOverSaturize = 0;
  maxvalue_mod = -1;
  maxsigma_mod = 3;
  maxvalue_stddev = -1;
  maxsigma_stddev = 3;
  title_index = -1;
  title2[0] = 0;
  doFlip = 1;
  ExtraVerticalMaximaSkip = 1;
  Load2dfs = 1;
  LoadLRFS = 1;
  plotvariance = 1;
  plotmodindex = 1;
  profilescale = 1;
  plot_xlabel = 0;
  notop = 0;
  noside = 0;
  overlaypp = 0;
  inside = 0;
  labelscale = 1.0;
  noylabels = 0;
  noxlabels = 0;

  plot_ylabel = 0;
  altProf = 0;
  dl = 0;
  ImaxSet = 0;
  IminSet = 0;
  lineStyle = 1;
  lineColour = 1;
  maxSubpulsePhaseSet = 0;
  minSubpulsePhaseSet = 0;
  maxSubpulsePhase = 0;
  minSubpulsePhase = 0;
  showwedge = 0;
  normalise_spectra = 0;
  normaliseSide = 0;
  nointegrateNumbers = 1;
  usephase = 0;


  for(i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-headerlist") == 0) {
      printHeaderCommandlineOptions(stdout);
      terminateApplication(&application);
      return 0;
    }
  }
  int argclast;
  argclast = argc-1;
  if(argv[argc-1][0] == '-')
    argclast += 1;
  for(i = 1; i < argclast; i++) {
    if(strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "-D") == 0 || strcmp(argv[i], "-device") == 0 || strcmp(argv[i], "-dev") == 0) {
      strcpy(PlotDevice, argv[i+1]);
      i++;
    }else if(strcmp(argv[i], "-v") == 0) {
      application.verbose_state.verbose = 1;
    }else if(strcmp(argv[i], "-phaseplot") == 0) {
      type_of_plots = 1;
    }else if(strcmp(argv[i], "-debug") == 0) {
      application.verbose_state.debug = 1;
    }else if(strcmp(argv[i], "-nocounters") == 0) {
      application.verbose_state.nocounters = 1;
    }else if(strcmp(argv[i], "-xlabel") == 0
      ) {
      plot_xlabel = 1;
    }else if(strcmp(argv[i], "-xlabel2") == 0
      ) {
      plot_xlabel = 2;
    }else if(strcmp(argv[i], "-xlabel3") == 0
      ) {
      plot_xlabel = 3;
    }else if(strcmp(argv[i], "-ylabel") == 0
      ) {
      plot_ylabel = 1;
    }else if(strcmp(argv[i], "-ylabel2") == 0
      ) {
      plot_ylabel = 2;
    }else if(strcmp(argv[i], "-ylabel3") == 0
      ) {
      plot_ylabel = 3;
    }else if(strcmp(argv[i], "-nolrfs") == 0
      ) {
      LoadLRFS = 0;
    }else if(strcmp(argv[i], "-2") == 0) {
      LoadTwo = 1;
    }else if(strcmp(argv[i], "-no2dfs") == 0
      ) {
      Load2dfs = 0;
    }else if(strcmp(argv[i], "-inside") == 0) {
      inside = 1;
    }else if(strcmp(argv[i], "-nostddev") == 0
      ) {
      plotvariance = 0;
    }else if(strcmp(argv[i], "-nomod") == 0
      ) {
      plotmodindex = 0;
    }else if(strcmp(argv[i], "-f") == 0) {
      doFlip = 0;
    }else if(strcmp(argv[i], "-notop") == 0) {
      notop = 1;
    }else if(strcmp(argv[i], "-noside") == 0) {
      noside = 1;
    }else if(strcmp(argv[i], "-nomain") == 0) {
      nomain = 1;
    }else if(strcmp(argv[i], "-overlay") == 0) {
      overlaypp = 1;
    }else if(strcmp(argv[i], "-noylabels") == 0) {
      noylabels = 1;
    }else if(strcmp(argv[i], "-noxlabels") == 0) {
      noxlabels = 1;
    }else if(strcmp(argv[i], "-showwedge") == 0) {
      showwedge = 1;
    }else if(strcmp(argv[i], "-normspectra") == 0) {
      normalise_spectra = 1;
    }else if(strcmp(argv[i], "-normside") == 0) {
      normaliseSide = 1;
    }else if(strcmp(argv[i], "-intnrs") == 0) {
      nointegrateNumbers = 0;
    }else if(strcmp(argv[i], "-phase") == 0) {
      usephase = 1;
    }else if(strcmp(argv[i], "-linestyle") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &lineStyle, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-linecolor") == 0 || strcmp(argv[i], "-linecolour") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &lineColour, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-scalep") == 0
      ) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &profilescale, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-dl") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &dl, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-labelscale") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &labelscale, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-Imax") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &ImaxValue, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      ImaxSet = 1;
      i++;
    }else if(strcmp(argv[i], "-Imin") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &IminValue, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      IminSet = 1;
      i++;
    }else if(strcmp(argv[i], "-spmax") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &maxSubpulsePhase, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      maxSubpulsePhaseSet = 1;
      i++;
    }else if(strcmp(argv[i], "-spmin") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &minSubpulsePhase, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      minSubpulsePhaseSet = 1;
      i++;
    }else if(strcmp(argv[i], "-header") == 0) {
      i++;
    }else if(strcmp(argv[i], "-altProf") == 0) {
      altProf = i+1;
      i++;
    }else if(strcmp(argv[i], "-p3") == 0) {
      SelectP3Region = 1;
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &P3RegionLow, &P3RegionHigh, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-p2") == 0) {
      if(SelectP2Region == 0) {
 SelectP2Region = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &P2RegionLow, &P2RegionHigh, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 P2RegionLow2 = P2RegionLow;
 P2RegionHigh2 = P2RegionHigh;
      }else {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &P2RegionLow2, &P2RegionHigh2, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
      }
      i++;
    }else if(strcmp(argv[i], "-l") == 0) {
      SelectLRegion = 1;
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &LRegionLow, &LRegionHigh, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-phaseslope") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &phase_slope_g, &phase_slope_o, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      do_phase_slope = 1;
      i++;
    }else if(strcmp(argv[i], "-scalel") == 0
      ) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &oversaturizel, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-scale2") == 0
      ) {
      if(NrSelectedOverSaturize == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &oversaturize, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 oversaturize2 = oversaturize;
 NrSelectedOverSaturize = 1;
      }else {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &oversaturize2, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
      }
      i++;
    }else if(strcmp(argv[i], "-modsigma") == 0
      ) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0,-1, "%f", &maxsigma_mod, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-stddevsigma") == 0
      ) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &maxsigma_stddev, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-modmax") == 0
      ) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &maxvalue_mod, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-stddevmax") == 0
      ) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &maxvalue_stddev, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-s") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &ExtraVerticalMaximaSkip, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-title") == 0
      ) {
      title_index = i+1;
      i++;
    }else if(strcmp(argv[i], "-textside") == 0
      ) {
      sprintf(title2, "%s", argv[i+1]);
      i++;
    }else if(strcmp(argv[i], "-2dfsnr") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &file_number, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-int") == 0
      ) {
      if(SelectP3Integrate == 0) {
 SelectP3Integrate = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &P3IntegrateLow, &P3IntegrateHigh, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 P3IntegrateLow2 = P3IntegrateLow;
 P3IntegrateHigh2 = P3IntegrateHigh;
      }else {
 SelectP3Integrate = 2;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &P3IntegrateLow2, &P3IntegrateHigh2, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
      }
      i++;
    }else if(strcmp(argv[i], "-ytop") == 0
      ) {
      plot_ylabeltop = 1;
    }else if(strcmp(argv[i], "-scalefig") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &scaleFig_x, &scaleFig_y, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else {
      printerror(application.verbose_state.debug, "Unknown option: %s\nRun pspecFig without command line options for a help", argv[i]);
      terminateApplication(&application);
      return 0;
    }
  }
  cleanPSRData(&twodfs, application.verbose_state);
  cleanPSRData(&lrfs, application.verbose_state);
  cleanPSRData(&AverageProfile, application.verbose_state);
  cleanPSRData(&VarianceProfile, application.verbose_state);
  cleanPSRData(&VarianceProfileErr, application.verbose_state);
  cleanPSRData(&ModProfile, application.verbose_state);
  cleanPSRData(&subpulseTrackProfile, application.verbose_state);
  cleanPSRData(&subpulseTrackProfileErr, application.verbose_state);
  cleanPSRData(&subpulseAmpProfile, application.verbose_state);
  if(type_of_plots == 0) {
    sprintf(txt, "lrfs");
    if(change_filename_extension(argv[argc-1], filename, txt, 1000, application.verbose_state) == 0) {
      return 0;
    }
    if(application.verbose_state.verbose)
      printf("Reading %s\n", filename);
    closePSRData(&lrfs, 0, application.verbose_state);
    if(!openPSRData(&lrfs, filename, 0, 0, 1, 0, application.verbose_state))
      return 0;
    long i;
    double max;
    if(normalise_spectra) {
      for(i = 0; i < lrfs.NrSubints * lrfs.NrBins; i++) {
 if(i == 0 || fabs(lrfs.data[i]) > max) {
   max = fabs(lrfs.data[i]);
 }
      }
      for(i = 0; i < lrfs.NrSubints * lrfs.NrBins; i++) {
 lrfs.data[i] /= 0.999*max;
      }
    }
    if(PSRDataHeader_parse_commandline(&lrfs, argc, argv, application.verbose_state) == 0)
      return 0;
    double period;
    int ret;
    ret = get_period(lrfs, 0, &period, application.verbose_state);
    if(ret == 2) {
      printerror(application.verbose_state.debug, "ERROR pspecFig (%s): Cannot obtain period", lrfs.filename);
      return 0;
    }
    if(period < 0.001) {
      printerror(application.verbose_state.debug, "The period does not appear to be set in the header. Consider using the -header option.");
      return 0;
    }
    if(get_tsamp(lrfs, 0, application.verbose_state) < 0.0000001) {
      printerror(application.verbose_state.debug, "The sampling time does not appear to be set in the header. Consider using the -header option.");
      return 0;
    }
    if(application.verbose_state.verbose)
      printf("%ldx%ld points read from lrfs\n", lrfs.NrBins, lrfs.NrSubints);
    if(lrfs.NrPols > 1) {
      if(preprocess_polselect(lrfs, &clone, 0, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Error selecting polarization channel 0.");
 return 0;
      }
      swap_orig_clone(&lrfs, &clone, application.verbose_state);
    }
  }
  copyVerboseState(application.verbose_state, &noverbose);
  noverbose.verbose = 0;
  if(application.verbose_state.verbose)
    printf("Reading %s\n", argv[argc-1]);
  closePSRData(&AverageProfile, 0, application.verbose_state);
  if(!openPSRData(&AverageProfile, argv[argc-1], 0, 0, 0, 0, application.verbose_state))
    return 0;
  if(!readHeaderPSRData(&AverageProfile, 0, 0, application.verbose_state))
    return 0;
  if(PSRDataHeader_parse_commandline(&AverageProfile, argc, argv, application.verbose_state) == 0)
    return 0;
  if(type_of_plots == 0) {
    if(AverageProfile.NrBins != lrfs.NrBins) {
      printwarning(application.verbose_state.debug, "WARNING: It looks like data is rebinned? Check the units.");
      AverageProfile.fixedtsamp *= AverageProfile.NrBins/(double)lrfs.NrBins;
      AverageProfile.tsampMode = TSAMPMODE_FIXEDTSAMP;
      AverageProfile.NrBins = lrfs.NrBins;
      printwarning(application.verbose_state.debug, "WARNING: Assuming the number of bins = %ld and the sampling time = %lf s.", AverageProfile.NrBins, AverageProfile.fixedtsamp);
    }
  }
  double period;
  int ret_prd;
  ret_prd = get_period(AverageProfile, 0, &period, application.verbose_state);
  if(ret_prd == 2) {
    printerror(application.verbose_state.debug, "ERROR pspecFig (%s): Cannot obtain period", AverageProfile.filename);
    return 0;
  }
  if(period < 0.001) {
    printerror(application.verbose_state.debug, "The period does not appear to be set in the header. Consider using the -header option.");
    return 0;
  }
  if(get_tsamp(AverageProfile, 0, application.verbose_state) < 0.0000001) {
    printerror(application.verbose_state.debug, "The sampling time does not appear to be set in the header. Consider using the -header option.");
    return 0;
  }
  AverageProfile.format = MEMORY_format;
  AverageProfile.NrSubints = 1;
  AverageProfile.NrFreqChan = 1;
  AverageProfile.NrPols = 1;
  AverageProfile.data = malloc(AverageProfile.NrBins*sizeof(float));
  if(AverageProfile.data == NULL) {
    printerror(application.verbose_state.debug, "Memory allocation error");
    return 0;
  }
  if(change_filename_extension(argv[argc-1], filename, "profile", 1000, application.verbose_state) == 0)
    return 0;
  if(altProf > 0) {
    strncpy(filename, argv[altProf], 1000);
  }
  if(set_filename_PSRData(&AverageProfile, filename, application.verbose_state) == 0) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR psecFig: Setting file name failed.");
    return 0;
  }
  if(type_of_plots == 0) {
    copy_params_PSRData(AverageProfile, &VarianceProfile, application.verbose_state);
    copy_params_PSRData(AverageProfile, &ModProfile, application.verbose_state);
    copy_params_PSRData(AverageProfile, &VarianceProfileErr, application.verbose_state);
    copy_params_PSRData(AverageProfile, &ModProfileErr, application.verbose_state);
    VarianceProfile.data = malloc(AverageProfile.NrBins*sizeof(float));
    ModProfile.data = malloc(AverageProfile.NrBins*sizeof(float));
    VarianceProfileErr.data = malloc(AverageProfile.NrBins*sizeof(float));
    ModProfileErr.data = malloc(AverageProfile.NrBins*sizeof(float));
    if(VarianceProfile.data == NULL || ModProfile.data == NULL || VarianceProfileErr.data == NULL || ModProfileErr.data == NULL) {
      printerror(application.verbose_state.debug, "Memory allocation error");
      return 0;
    }
  }else {
    copy_params_PSRData(AverageProfile, &subpulseTrackProfile, application.verbose_state);
    copy_params_PSRData(AverageProfile, &subpulseTrackProfileErr, application.verbose_state);
    copy_params_PSRData(AverageProfile, &subpulseAmpProfile, application.verbose_state);
    subpulseTrackProfile.data = malloc(AverageProfile.NrBins*sizeof(float));
    subpulseTrackProfileErr.data = malloc(AverageProfile.NrBins*sizeof(float));
    subpulseAmpProfile.data = malloc(AverageProfile.NrBins*sizeof(float));
    if(subpulseTrackProfile.data == NULL || subpulseTrackProfileErr.data == NULL || subpulseAmpProfile.data == NULL) {
      printerror(application.verbose_state.debug, "Memory allocation error");
      return 0;
    }
  }
  if(application.verbose_state.verbose)
    printf("Reading %s\n", filename);
  fout_ascii = fopen(filename, "r");
  if(fout_ascii == NULL) {
    printerror(application.verbose_state.debug, "ERROR pspecFig: Unable to open %s.", filename);
    return 0;
  }
  for(i = 0; i < AverageProfile.NrBins; i++) {
    if(type_of_plots == 0) {
      j = fscanf(fout_ascii, "%ld %f %f %f %f %f", &k, &AverageProfile.data[i], &VarianceProfile.data[i], &VarianceProfileErr.data[i], &ModProfile.data[i], &ModProfileErr.data[i]);
    }else {
      float junk;
      j = fscanf(fout_ascii, "%ld %f %f %f %f %f", &k, &AverageProfile.data[i], &junk, &junk, &junk, &junk);
    }
    if(j != 6) {
      if(type_of_plots == 0) {
 printerror(application.verbose_state.debug, "Unexpected end of file (bin %d). Resolution in profile doesn't match resolution in lrfs?", i+1);
 return 0;
      }else {
 printwarning(application.verbose_state.debug, "WARNING: It looks like data is rebinned? Check the units.");
 AverageProfile.fixedtsamp *= AverageProfile.NrBins/(double)(i);
 AverageProfile.tsampMode = TSAMPMODE_FIXEDTSAMP;
 AverageProfile.NrBins = i;
 printwarning(application.verbose_state.debug, "WARNING: Assuming the number of bins = %ld and the sampling time = %lf s.", AverageProfile.NrBins, AverageProfile.fixedtsamp);
 copy_params_PSRData(AverageProfile, &subpulseTrackProfile, application.verbose_state);
 copy_params_PSRData(AverageProfile, &subpulseTrackProfileErr, application.verbose_state);
 copy_params_PSRData(AverageProfile, &subpulseAmpProfile, application.verbose_state);
 break;
      }
    }
    if(k != i) {
      printerror(application.verbose_state.debug, "Unexpected bin number");
      return 0;
    }
  }
  fclose(fout_ascii);
  if(application.verbose_state.verbose)
    printf("%ld points read from %s\n", AverageProfile.NrBins, filename);
  ret = get_period(AverageProfile, 0, &period, application.verbose_state);
  if(ret == 2) {
    printerror(application.verbose_state.debug, "ERROR pspecFig (%s): Cannot obtain period", AverageProfile.filename);
    return 0;
  }
  fl_min = 0 + dl;
  fl_max = 360*(AverageProfile.NrBins-1)*get_tsamp(AverageProfile, 0, application.verbose_state)/period + dl;
  if(usephase) {
    fl_min /= 360.0;
    fl_max /= 360.0;
  }
  if(type_of_plots == 1) {
    if(change_filename_extension(argv[argc-1], filename, "amplitude", 1000, application.verbose_state) == 0)
      return 0;
    if(application.verbose_state.verbose)
      printf("Reading %s\n", filename);
    fout_ascii = fopen(filename, "r");
    if(fout_ascii == NULL) {
      printerror(application.verbose_state.debug, "ERROR pspecFig: Unable to open %s.", filename);
      return 0;
    }
    for(i = 0; i < AverageProfile.NrBins; i++) {
      j = fscanf(fout_ascii, "%ld %f", &k, &subpulseAmpProfile.data[i]);
      if(j != 2) {
 printerror(application.verbose_state.debug, "Unexpected end of file (bin %d). Resolution in subpulse phase track doesn't match resolution in profile file?", i+1);
 return 0;
      }
      if(k != i) {
 printerror(application.verbose_state.debug, "Unexpected bin number");
 return 0;
      }
    }
    fclose(fout_ascii);
    if(application.verbose_state.verbose)
      printf("%ld points read from %s\n", AverageProfile.NrBins, filename);
  }
  if(type_of_plots == 1) {
    if(change_filename_extension(argv[argc-1], filename, "track", 1000, application.verbose_state) == 0)
      return 0;
    if(application.verbose_state.verbose)
      printf("Reading %s\n", filename);
    fout_ascii = fopen(filename, "r");
    if(fout_ascii == NULL) {
      printerror(application.verbose_state.debug, "ERROR pspecFig: Unable to open %s.", filename);
      return 0;
    }
    for(i = 0; i < AverageProfile.NrBins; i++) {
      j = fscanf(fout_ascii, "%ld %f %f", &k, &subpulseTrackProfile.data[i], &subpulseTrackProfileErr.data[i]);
      if(j != 3) {
 printerror(application.verbose_state.debug, "Unexpected end of file (bin %d). Resolution in subpulse phase track doesn't match resolution in profile file?", i+1);
 return 0;
      }
      if(k != i) {
 printerror(application.verbose_state.debug, "Unexpected bin number");
 return 0;
      }
    }
    fclose(fout_ascii);
    if(application.verbose_state.verbose)
      printf("%ld points read from %s\n", AverageProfile.NrBins, filename);
  }
  if(type_of_plots == 0) {
    sprintf(txt, "%d.2dfs", file_number);
    if(change_filename_extension(argv[argc-1], filename, txt, 1000, application.verbose_state) == 0) {
      return 0;
    }
    if(application.verbose_state.verbose)
      printf("Reading %s\n", filename);
    closePSRData(&twodfs, 0, application.verbose_state);
    if(!openPSRData(&twodfs, filename, 0, 0, 1, 0, application.verbose_state))
      return 0;
    long i;
    double max;
    if(normalise_spectra) {
      for(i = 0; i < twodfs.NrSubints * twodfs.NrBins; i++) {
 if(i == 0 || fabs(twodfs.data[i]) > max) {
   max = fabs(twodfs.data[i]);
 }
      }
      for(i = 0; i < twodfs.NrSubints * twodfs.NrBins; i++) {
 twodfs.data[i] /= 0.999*max;
      }
    }
    if(PSRDataHeader_parse_commandline(&twodfs, argc, argv, application.verbose_state) == 0)
      return 0;
    if(application.verbose_state.verbose)
      printf("%ldx%ld points read from 2dfs\n", twodfs.NrBins, twodfs.NrSubints);
    if(twodfs.NrPols > 1) {
      if(preprocess_polselect(twodfs, &clone, 0, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Error selecting polarization channel 0.");
 return 0;
      }
      swap_orig_clone(&twodfs, &clone, application.verbose_state);
    }
    f2_min = -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins;
    f2_max = +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins;
    f3_min = 0;
    f3_max = 0.5;
    if(P3IntegrateHigh > f3_max)
      P3IntegrateHigh = f3_max;
    if(LoadTwo != 0) {
      sprintf(txt, "%d.2dfs", file_number+1);
      if(change_filename_extension(argv[argc-1], filename, txt, 1000, application.verbose_state) == 0) {
 return 0;
      }
      if(application.verbose_state.verbose)
 printf("Reading %s\n", filename);
      if(!openPSRData(&twodfs2, filename, 0, 0, 1, 0, application.verbose_state))
 return 0;
      long i;
      double max;
      if(normalise_spectra) {
 for(i = 0; i < twodfs2.NrSubints * twodfs2.NrBins; i++) {
   if(i == 0 || fabs(twodfs2.data[i]) > max) {
     max = fabs(twodfs2.data[i]);
   }
 }
 for(i = 0; i < twodfs2.NrSubints * twodfs2.NrBins; i++) {
   twodfs2.data[i] /= 0.999*max;
 }
      }
      if(PSRDataHeader_parse_commandline(&twodfs2, argc, argv, application.verbose_state) == 0)
 return 0;
      if(application.verbose_state.verbose)
 printf("%ldx%ld points read from 2dfs\n", twodfs2.NrBins, twodfs2.NrSubints);
      if(twodfs2.NrPols > 1) {
 if(preprocess_polselect(twodfs2, &clone, 0, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Error selecting polarization channel 0.");
   return 0;
 }
 swap_orig_clone(&twodfs2, &clone, application.verbose_state);
      }
      f2_min2 = -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs2.NrBins;
      f2_max2 = +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs2.NrBins;
      f3_min2 = 0;
      f3_max2 = 0.5;
      if(P3IntegrateHigh2 > f3_max2)
 P3IntegrateHigh2 = f3_max2;
    }
    if(f3_min < 0)
      f3_min = 0;
    if(f3_min2 < 0)
      f3_min2 = 0;
    if(SelectP3Region != 0) {
      f3_min = P3RegionLow;
      f3_max = P3RegionHigh;
      f3_min2 = P3RegionLow;
      f3_max2 = P3RegionHigh;
    }
    if(SelectP3Integrate == 0) {
      P3IntegrateLow = f3_min;
      P3IntegrateHigh = f3_max;
      P3IntegrateLow2 = f3_min;
      P3IntegrateHigh2 = f3_max;
    }
  }
  ppgopen(PlotDevice);
  ppgask(0);
  ppgslw(1);
  ppgpage();
  ppgslw(1);
  ppgscf(1);
  ppgsch(0.38*labelscale);
  if(notop == 0) {
    ppgsch(0.55*labelscale);
    ppgscf(2);
    ppgslw(2);
    ppgsvp(0.2, 0.2+0.11*scaleFig_x, 0.95-0.15*scaleFig_y, 0.95);
    if(title_index > 0) {
      ppgmtxt("t",1,0.5,0.5,argv[title_index]);
    }
    ppgslw(1);
    ppgscf(1);
    ppgsch(0.38*labelscale);
    Imin = 0;
    Imax = 0;
    ret = get_period(AverageProfile, 0, &period, application.verbose_state);
    if(ret == 2) {
      printerror(application.verbose_state.debug, "ERROR pspecFig (%s): Cannot obtain period", AverageProfile.filename);
      return 0;
    }
    for(xi=0; xi < AverageProfile.NrBins; xi++) {
      double xpos;
      xpos = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
      if(usephase)
 xpos /= 360.0;
      if(xpos >= fl_min && xpos <= fl_max) {
 I = AverageProfile.data[xi];
 if(I > Imax)
   Imax = I;
 if(I < Imin)
   Imin = I;
      }
    }
    if(SelectLRegion != 0) {
      fl_min = LRegionLow;
      fl_max = LRegionHigh;
    }
    if(type_of_plots == 0) {
      for(xi=0; xi < ModProfile.NrBins; xi++) {
 double xpos;
 xpos = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
 if(usephase)
   xpos /= 360.0;
 if(xpos+dl >= fl_min && xpos+dl <= fl_max) {
   ok_flag = 1;
   if(maxsigma_mod > 0 && ModProfile.data[xi]/ModProfileErr.data[xi] < maxsigma_mod)
     ok_flag = 0;
   if(maxvalue_mod > 0 && ModProfile.data[xi] > maxvalue_mod)
     ok_flag = 0;
   if(ok_flag) {
     I = ModProfile.data[xi];
     if(I+ModProfileErr.data[xi] > Imax)
       Imax = I+ModProfileErr.data[xi];
     if(I-ModProfileErr.data[xi] < Imin)
       Imin = I-ModProfileErr.data[xi];
   }
 }
      }
      for(xi=0; xi < VarianceProfile.NrBins; xi++) {
 double xpos;
 xpos = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
 if(usephase)
   xpos /= 360.0;
 if(xpos + dl >= fl_min && xpos + dl <= fl_max) {
   I = VarianceProfile.data[xi];
   if(I > Imax)
     Imax = I;
   if(I < Imin)
     Imin = I;
 }
      }
    }else {
      for(xi=0; xi < subpulseAmpProfile.NrBins; xi++) {
 double xpos;
 xpos = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
 if(usephase)
   xpos /= 360.0;
 if(xpos + dl >= fl_min && xpos + dl <= fl_max) {
   I = subpulseAmpProfile.data[xi];
   if(I > Imax)
     Imax = I;
   if(I < Imin)
     Imin = I;
 }
      }
    }
    if(-0.05*Imax < Imin)
      Imin = -0.05*Imax;
    if(ImaxSet)
      Imax = ImaxValue/1.05;
    if(IminSet)
      Imin = IminValue/1.05;
    ppgswin(fl_min, fl_max, Imin, 1.05*Imax);
    if(LoadLRFS == 0) {
      if(inside) {
 if(noylabels) {
   ppgbox("bcnst",0.0,0,"bcst",0.0,0);
 }else {
   ppgbox("bcnst",0.0,0,"bcnst",0.0,0);
 }
      }else {
 if(noylabels) {
   ppgbox("bcnst",0.0,0,"bcsti",0.0,0);
 }else {
   ppgbox("bcnst",0.0,0,"bcnsti",0.0,0);
 }
      }
    }else {
      if(inside) {
 if(noylabels) {
   ppgbox("cst",0.0,0,"bcst",0.0,0);
 }else {
   ppgbox("cst",0.0,0,"bcnst",0.0,0);
 }
      }else {
 if(noylabels) {
   ppgbox("cst",0.0,0,"bcsti",0.0,0);
 }else {
   ppgbox("cst",0.0,0,"bcnsti",0.0,0);
 }
      }
    }
    if(plot_xlabel != 0 && LoadLRFS == 0) {
      ppgsch(0.3*labelscale);
      if(usephase)
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (phase)");
      else
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (deg)");
      ppgsch(0.38*labelscale);
    }
    ppgsls(lineStyle);
    ppgslw(2);
    ppgsci(lineColour);
    i = 0;
    for(xi=0; xi < AverageProfile.NrBins; xi++) {
      double xpos;
      xpos = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
      if(usephase)
 xpos /= 360.0;
      if(xpos+dl >= fl_min && xpos+dl <= fl_max) {
 I = AverageProfile.data[xi]*profilescale;
 if(i == 0) {
   ppgmove(xpos+dl,I);
   i = 1;
 }else {
   ppgdraw(xpos+dl,I);
 }
      }
    }
    ppgsci(1);
    if(type_of_plots == 0) {
      if(plotvariance) {
 i = 0;
 ppgsls(1);
 ppgslw(1);
 ppgscr(20, 0.5, 0.5, 0.5);
 ppgsci(20);
 for(xi=0; xi < AverageProfile.NrBins; xi++) {
   x = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
   if(usephase)
     x /= 360.0;
   x += dl;
   if(x >= fl_min && x <= fl_max) {
     I = VarianceProfile.data[xi];
     ok_flag = 1;
     if(maxsigma_stddev > 0 && VarianceProfile.data[xi]/VarianceProfileErr.data[xi] < maxsigma_stddev)
       ok_flag = 0;
     if(maxvalue_stddev > 0 && VarianceProfile.data[xi] > maxvalue_stddev)
       ok_flag = 0;
     if(ok_flag == 0) {
       I = 0;
       i = 0;
     }else {
       if(i == 0) {
  ppgmove(x,I);
  i = 1;
       }else {
  ppgsci(20);
  ppgdraw(x,I);
       }
       ppgsci(1);
       ppgpt1(x, I, 4);
     }
   }
 }
 ppgsci(1);
      }
      if(plotmodindex) {
 i = 0;
 ppgsls(1);
 ppgslw(1);
 for(xi=0; xi < ModProfile.NrBins; xi++) {
   x = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
   if(usephase)
     x /= 360.0;
   x += dl;
   if(x >= fl_min && x <= fl_max) {
     I = ModProfile.data[xi];
     ok_flag = 1;
     if(maxsigma_mod > 0 && ModProfile.data[xi]/ModProfileErr.data[xi] < maxsigma_mod)
       ok_flag = 0;
     if(maxvalue_mod > 0 && ModProfile.data[xi] > maxvalue_mod)
       ok_flag = 0;
     if(ok_flag == 0) {
       I = 0;
       i = 0;
     }else {
       if(i == 0) {
  ppgmove(x,I);
  i = 1;
       }else {
  ppgdraw(x,I);
       }
     }
   }
 }
 ppgsls(1);
 ppgslw(1);
 for(xi=0; xi < ModProfile.NrBins; xi++) {
   x = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
   if(usephase)
     x /= 360.0;
   x += dl;
   if(x >= fl_min && x <= fl_max) {
     I = ModProfile.data[xi];
     ok_flag = 1;
     if(maxsigma_mod > 0 && ModProfile.data[xi]/ModProfileErr.data[xi] < maxsigma_mod)
       ok_flag = 0;
     if(maxvalue_mod > 0 && ModProfile.data[xi] > maxvalue_mod)
       ok_flag = 0;
     if(ok_flag != 0) {
       ppgpt1(x, I, -1);
       ppgerr1(6,x,I,ModProfileErr.data[xi], 1);
     }
   }
 }
      }
    }else {
      i = 0;
      ppgsls(4);
      ppgslw(1);
      for(xi=0; xi < AverageProfile.NrBins; xi++) {
 x = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
 if(usephase)
   x /= 360.0;
 x += dl;
 if(x >= fl_min && x <= fl_max) {
   I = subpulseAmpProfile.data[xi];
   if(i == 0) {
     ppgmove(x, I);
     i = 1;
   }else {
     ppgdraw(x, I);
   }
 }
      }
      ppgsls(1);
    }
    if(plot_ylabel != 0) {
      ppgsch(0.3*labelscale);
      if(plot_ylabeltop) {
 if(type_of_plots == 0) {
   ppgmtxt("l",2.8,0.5,0.5,"Intensity/Modulation index");
 }else {
   ppgmtxt("l",2.8,0.5,0.5,"Intensity");
 }
      }
      ppgsch(0.38*labelscale);
    }
  }
  if(type_of_plots == 0) {
    if(Load2dfs != 0 && LoadLRFS == 0) {
      if(title_index > 0 && notop && LoadLRFS == 0) {
 Plot2dfsOnly(0, plot_xlabel, plot_ylabel, noside, nomain, showwedge, argv[title_index], scaleFig_x, scaleFig_y, normaliseSide, nointegrateNumbers, application.verbose_state);
      }else {
 Plot2dfsOnly(0, plot_xlabel, plot_ylabel, noside, nomain, showwedge, NULL, scaleFig_x, scaleFig_y, normaliseSide, nointegrateNumbers, application.verbose_state);
      }
      if (LoadTwo != 0)
 Plot2dfsOnly(1, plot_xlabel, plot_ylabel, noside, nomain, showwedge, NULL, scaleFig_x, scaleFig_y, normaliseSide, nointegrateNumbers, application.verbose_state);
    }
    else {
      if(LoadLRFS != 0) {
 if(title_index > 0 && notop) {
   PlotLRFS(plot_xlabel, plot_ylabel, noside, plot_ylabeltop, lineColour, showwedge, argv[title_index], scaleFig_x, scaleFig_y, normaliseSide, nointegrateNumbers, usephase, application.verbose_state);
 }else {
   PlotLRFS(plot_xlabel, plot_ylabel, noside, plot_ylabeltop, lineColour, showwedge, NULL, scaleFig_x, scaleFig_y, normaliseSide, nointegrateNumbers, usephase, application.verbose_state);
 }
      }
      if(Load2dfs != 0)
 Plot2dfs(0, plot_xlabel, plot_ylabel, noside, plot_ylabeltop, nomain, scaleFig_x, scaleFig_y, showwedge, normaliseSide, nointegrateNumbers, application.verbose_state);
      if(LoadTwo != 0)
 Plot2dfs(1, plot_xlabel, plot_ylabel, noside, plot_ylabeltop, nomain, scaleFig_x, scaleFig_y, showwedge, normaliseSide, nointegrateNumbers, application.verbose_state);
    }
  }else {
    ppgsvp(0.2, 0.2+0.11*scaleFig_x, 0.95-0.3*scaleFig_y, 0.95-0.15*scaleFig_y);
    ppgslw(1);
    ppgscf(1);
    ppgsch(0.38*labelscale);
    Imin = 0;
    Imax = 0;
    for(xi=0; xi < subpulseTrackProfile.NrBins; xi++) {
      double xpos;
      xpos = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
      if(usephase)
 xpos /= 360.0;
      if(xpos+dl >= fl_min && xpos+dl <= fl_max) {
 I = subpulseTrackProfile.data[xi];
 if(doFlip)
   I *= -1;
 if(I+subpulseTrackProfileErr.data[xi] > Imax)
   Imax = I+subpulseTrackProfileErr.data[xi];
 if(I-subpulseTrackProfileErr.data[xi] < Imin)
   Imin = I-subpulseTrackProfileErr.data[xi];
      }
    }
    if(SelectLRegion != 0) {
      fl_min = LRegionLow;
      fl_max = LRegionHigh;
    }
    if(maxSubpulsePhaseSet)
      Imax = maxSubpulsePhase/1.05;
    if(minSubpulsePhaseSet)
      Imin = minSubpulsePhase/1.05;
    ppgswin(fl_min, fl_max, Imin, 1.05*Imax);
    if(inside) {
      if(noylabels) {
 ppgbox("bcnst",0.0,0,"bcst",0.0,0);
      }else {
 ppgbox("bcnst",0.0,0,"bcnst",0.0,0);
      }
    }else {
      if(noylabels) {
 ppgbox("bcnst",0.0,0,"bcsti",0.0,0);
      }else {
 ppgbox("bcnst",0.0,0,"bcnsti",0.0,0);
      }
    }
    if(plot_xlabel != 0 && LoadLRFS == 0) {
      ppgsch(0.3*labelscale);
      if(usephase)
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (phase)");
      else
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (deg)");
      ppgsch(0.38*labelscale);
    }
    if(do_phase_slope) {
      ppgsls(4);
      if(usephase == 0) {
 x = dl;
 I = x*phase_slope_g + phase_slope_o;
 I = derotate_deg(I);
 ppgmove(x, I);
 ppgdraw(x+360, I+phase_slope_g*360);
 ppgmove(x, I-360);
 ppgdraw(x+360, I-360+phase_slope_g*360);
 ppgmove(x, I-2*360);
 ppgdraw(x+360, I-2*360+phase_slope_g*360);
 ppgmove(x, I+360);
 ppgdraw(x+360, I+360+phase_slope_g*360);
 ppgmove(x, I+2*360);
 ppgdraw(x+360, I+2*360+phase_slope_g*360);
      }else {
 x = dl;
 I = x*360.0*phase_slope_g + phase_slope_o;
 I = derotate_deg(I);
 ppgmove(x, I);
 ppgdraw(x+1.0, I+phase_slope_g*360);
 ppgmove(x, I-360);
 ppgdraw(x+1.0, I-360+phase_slope_g*360);
 ppgmove(x, I-2*360);
 ppgdraw(x+1.0, I-2*360+phase_slope_g*360);
 ppgmove(x, I+360);
 ppgdraw(x+1.0, I+360+phase_slope_g*360);
 ppgmove(x, I+2*360);
 ppgdraw(x+1.0, I+2*360+phase_slope_g*360);
      }
      ppgsls(1);
    }
    ppgscr(20, 0.5, 0.5, 0.5);
    ppgsci(20);
    ppgsls(1);
    ppgslw(1);
    for(xi=0; xi < subpulseTrackProfile.NrBins; xi++) {
      x = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
      if(usephase)
 x /= 360.0;
      x += dl;
      if(x >= fl_min && x <= fl_max) {
 I = subpulseTrackProfile.data[xi];
 if(doFlip)
   I *= -1;
 ppgerr1(6,x,I,subpulseTrackProfileErr.data[xi], 1);
 ppgerr1(6,x,I+360,subpulseTrackProfileErr.data[xi], 1);
 ppgerr1(6,x,I-360,subpulseTrackProfileErr.data[xi], 1);
      }
    }
    ppgsci(1);
    ppgslw(3);
    for(xi=0; xi < subpulseTrackProfile.NrBins; xi++) {
      x = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
      if(usephase)
 x /= 360.0;
      x += dl;
      if(x >= fl_min && x <= fl_max) {
 I = subpulseTrackProfile.data[xi];
 if(doFlip)
   I *= -1;
 ppgpt1(x, I, -1);
 ppgpt1(x, I+360, -1);
 ppgpt1(x, I-360, -1);
      }
    }
    ppgslw(1);
    if(plot_xlabel != 0) {
      ppgsch(0.3*labelscale);
      if(usephase)
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (phase)");
      else
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (deg)");
      ppgsch(0.38*labelscale);
    }
    if(plot_ylabel != 0) {
      ppgsch(0.3*labelscale);
      if(usephase)
 ppgmtxt("l",2.8,0.5,0.5,"Subpulse phase (phase)");
      else
 ppgmtxt("l",2.8,0.5,0.5,"Subpulse phase (deg)");
      ppgsch(0.38*labelscale);
    }
  }
  if(preprocess_checknan(AverageProfile, 1, noverbose)) {
    printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have NAN's in them. Artifacts can be expected in the plot.");
  }
  if(preprocess_checkinf(AverageProfile, 1, noverbose)) {
    printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have INF's in them. Artifacts can be expected in the plot.");
  }
  if(type_of_plots == 0) {
    if(preprocess_checknan(VarianceProfile, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have NAN's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checkinf(VarianceProfile, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have INF's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checknan(ModProfile, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have NAN's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checkinf(ModProfile, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have INF's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checknan(VarianceProfileErr, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have NAN's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checkinf(VarianceProfileErr, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have INF's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checknan(ModProfileErr, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have NAN's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checkinf(ModProfileErr, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have INF's in them. Artifacts can be expected in the plot.");
    }
  }else {
    if(preprocess_checknan(subpulseTrackProfile, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have NAN's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checkinf(subpulseTrackProfile, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have INF's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checknan(subpulseTrackProfileErr, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have NAN's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checkinf(subpulseTrackProfileErr, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have INF's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checknan(subpulseAmpProfile, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have NAN's in them. Artifacts can be expected in the plot.");
    }
    if(preprocess_checkinf(subpulseAmpProfile, 1, noverbose)) {
      printwarning(application.verbose_state.debug, "WARNING: The profile data appears to have INF's in them. Artifacts can be expected in the plot.");
    }
  }
  ppgend();
  closePSRData(&AverageProfile, 0, application.verbose_state);
    closePSRData(&twodfs, 0, application.verbose_state);
    closePSRData(&lrfs, 0, application.verbose_state);
    closePSRData(&VarianceProfile, 0, application.verbose_state);
    closePSRData(&ModProfile, 0, application.verbose_state);
    closePSRData(&VarianceProfileErr, 0, application.verbose_state);
    closePSRData(&ModProfileErr, 0, application.verbose_state);
    closePSRData(&subpulseTrackProfile, 0, application.verbose_state);
    closePSRData(&subpulseTrackProfileErr, 0, application.verbose_state);
    closePSRData(&subpulseAmpProfile, 0, application.verbose_state);
  terminateApplication(&application);
  return 0;
}
void GetExtremesSubset(float *Imin, float *Imax)
{
  fprintf(stderr, "GetExtremesSubset function disabled\n");
}
void GetExtremesSubsetVertical(float *Imin, float *Imax, int Number)
{
  float I, x, y;
  int xi, yi;
  *Imin = 0;
  *Imax = 0;
  if(Number == 0) {
    for(xi = 0; xi < twodfs.NrBins; xi++) {
      I = 0;
      for(yi = 0; yi < twodfs.NrSubints; yi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(y >= P3IntegrateLow && y <= P3IntegrateHigh)
   I += twodfs.data[yi*twodfs.NrBins+xi];
      }
      if(I < *Imin)
 *Imin = I;
      if(I > *Imax)
 *Imax = I;
    }
  }else if(Number == 1) {
    for(xi = 0; xi < twodfs2.NrBins; xi++) {
      I = 0;
      for(yi = 0; yi < twodfs2.NrSubints; yi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(y >= P3IntegrateLow2 && y <= P3IntegrateHigh2)
   I += twodfs2.data[yi*twodfs2.NrBins+xi];
      }
      if(I < *Imin)
 *Imin = I;
      if(I > *Imax)
 *Imax = I;
    }
  }
  *Imin *= 2.0;
  *Imax *= 2.0;
}
void GetExtremesSubsetHorizontal(float *Imin, float *Imax, int Number)
{
  float I, x, y;
  int xi, yi;
  *Imin = 0;
  *Imax = 0;
  if(Number == -1) {
    for(yi = ExtraVerticalMaximaSkip; yi < lrfs.NrSubints; yi++) {
      I = 0;
      for(xi = 0; xi < lrfs.NrBins; xi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(x >= fl_min && x <= fl_max && y >= f3_min && y <= f3_max)
   I += lrfs.data[yi*lrfs.NrBins+xi];
      }
      if(I < *Imin)
 *Imin = I;
      if(I > *Imax)
 *Imax = I;
    }
  }else if(Number == 0) {
    for(yi = ExtraVerticalMaximaSkip; yi < twodfs.NrSubints; yi++) {
      I = 0;
      for(xi = 0; xi < twodfs.NrBins; xi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(x >= f2_min && x <= f2_max && y >= f3_min && y <= f3_max)
   I += twodfs.data[yi*twodfs.NrBins+xi];
      }
      if(I < *Imin)
 *Imin = I;
      if(I > *Imax)
 *Imax = I;
    }
  }else if(Number == 1) {
    for(yi = ExtraVerticalMaximaSkip; yi < twodfs2.NrSubints; yi++) {
      I = 0;
      for(xi = 0; xi < twodfs2.NrBins; xi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(x >= f2_min && x <= f2_max && y >= f3_min && y <= f3_max)
   I += twodfs2.data[yi*twodfs2.NrBins+xi];
      }
      if(I < *Imin)
 *Imin = I;
      if(I > *Imax)
 *Imax = I;
    }
  }
  *Imin *= 2.0;
  *Imax *= 2.0;
}
void IntegrateSubsetHorizontal(int Number, double scaleFig_x, double scaleFig_y, int normalise, int nointegrateNumbers)
{
  float I, Imin, Imax, x, y;
  int xi, yi;
  float offset;
  offset = 0+0.03*(labelscale-1.0)*((float)Number+1.0);
  if(Number == -1) {
    ppgsvp(0.2-0.04*scaleFig_x, 0.2, 0.95-0.3*scaleFig_y, 0.95-0.15*scaleFig_y);
  }else {
    ppgsvp(0.2-0.04*scaleFig_x, 0.2, 0.95-0.48*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y, 0.95-0.33*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y);
  }
  GetExtremesSubsetHorizontal(&Imin, &Imax, Number);
  double scale;
  scale = fabs(Imax);
  if(fabs(Imin) > scale)
    scale = fabs(Imin);
  if(normalise == 0)
    scale = 1.0;
  ppgswin(Imin/scale,1.05*Imax/scale, f3_min,f3_max);
  ppgbox("bcv",0.0,0,"bnsti",0.0,0);
  ppgbox("",0.0,0,"c",0.0,0);
  if(Number == -1) {
    if(title2[0] != 0) {
      ppgsch(0.55*labelscale);
      ppgscf(2);
      ppgslw(2);
      if (notop == 1) {
 ppgtext(Imin, (f3_max-f3_min)*1.15, title2);
      } else {
        ppgtext(Imin, (f3_max-f3_min)*2.095, title2);
      }
      ppgslw(1);
      ppgscf(1);
      ppgsch(0.38*labelscale);
    }
  }
  y = floor(log10(Imax*0.33/scale));
  x = floor(Imax*0.33/(pow(10,y)*scale));
  x = x*pow(10,y);
  char labelnumbers[3];
  if(nointegrateNumbers == 0) {
    strcpy(labelnumbers, "n");
  }else {
    strcpy(labelnumbers, "");
  }
  ppgsch(0.38*labelscale*0.66);
  ppgaxis(labelnumbers,0,f3_max,Imax*0.33/scale,f3_max,0,Imax*0.33/scale,x,1,0.3,0,0,-0.5,90);
  ppgsch(0.38*labelscale);
  if(Number == -1) {
    for(yi = 0; yi < lrfs.NrSubints; yi++) {
      I = 0;
      for(xi = 0; xi < lrfs.NrBins; xi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(x >= fl_min && x <= fl_max && y >= f3_min && y <= f3_max) {
   I += lrfs.data[yi*lrfs.NrBins+xi];
 }
      }
      if(yi == 0) {
 ppgmove(2.0*I/scale, y);
      }else {
 ppgdraw(2.0*I/scale, y);
      }
    }
  }else if(Number == 0) {
    for(yi = 0; yi < twodfs.NrSubints; yi++) {
      I = 0;
      for(xi = 0; xi < twodfs.NrBins; xi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(x >= f2_min && x <= f2_max && y >= f3_min && y <= f3_max)
   I += twodfs.data[yi*twodfs.NrBins+xi];
      }
      if(yi == 0) {
 ppgmove(2.0*I/scale, y);
      }else {
 ppgdraw(2.0*I/scale, y);
      }
    }
  }else if(Number == 1) {
    for(yi = 0; yi < twodfs2.NrSubints; yi++) {
      I = 0;
      for(xi = 0; xi < twodfs2.NrBins; xi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(x >= f2_min2 && x <= f2_max2 && y >= f3_min2 && y <= f3_max2)
   I += twodfs2.data[yi*twodfs2.NrBins+xi];
      }
      if(yi == 0) {
 ppgmove(2.0*I/scale, y);
      }else {
 ppgdraw(2.0*I/scale, y);
      }
    }
  }
}
void IntegrateSubsetHorizontal_2dfsOnly(int Number, double scaleFig_x, double scaleFig_y, int normalise, int nointegrateNumbers)
{
  float I, Imin, Imax, x, y, offset;
  int xi, yi;
  offset = 0+0.03*(labelscale-1.0)*((float)Number+1.0);
  if(Number == -1) {
    ppgsvp(0.2-0.04*scaleFig_x, 0.2, 0.95-0.3*scaleFig_y, 0.95-0.15*scaleFig_y);
  }else {
    ppgsvp(0.2-0.04*scaleFig_x, 0.2, 0.95-0.33*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y, 0.95-0.18*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y);
  }
  GetExtremesSubsetHorizontal(&Imin, &Imax, Number);
  double scale;
  scale = fabs(Imax);
  if(fabs(Imin) > scale)
    scale = fabs(Imin);
  if(normalise == 0)
    scale = 1.0;
  ppgswin(Imin/scale,1.05*Imax/scale, f3_min,f3_max);
  ppgbox("bcv",0.0,0,"bnsti",0.0,0);
  ppgbox("",0.0,0,"c",0.0,0);
  if(Number == 0) {
    if(strlen(title2) != 0 && notop == 0) {
      ppgsch(0.55*labelscale);
      ppgscf(2);
      ppgslw(2);
      ppgtext(Imin, (f3_max-f3_min)+0.65, title2);
      ppgslw(1);
      ppgscf(1);
      ppgsch(0.38*labelscale);
    }else if (strlen(title2) != 0 && notop == 1) {
      ppgsch(0.55*labelscale);
      ppgscf(2);
      ppgslw(2);
      ppgtext(Imin, (f3_max-f3_min)+0.1, title2);
      ppgslw(1);
      ppgscf(1);
      ppgsch(0.38*labelscale);
    }
  }
  y = floor(log10(Imax*0.33/scale));
  x = floor(Imax*0.33/(pow(10,y)*scale));
  x = x*pow(10,y);
  char labelnumbers[3];
  if(nointegrateNumbers == 0) {
    strcpy(labelnumbers, "n");
  }else {
    strcpy(labelnumbers, "");
  }
  ppgsch(0.38*labelscale*0.66);
  ppgaxis(labelnumbers,0,f3_max,Imax*0.33/scale,f3_max,0,Imax*0.33/scale,x,1,0.3,0,0,-0.5,90);
  ppgsch(0.38*labelscale);
  if(Number == 0) {
    for(yi = 0; yi < twodfs.NrSubints; yi++) {
      I = 0;
      for(xi = 0; xi < twodfs.NrBins; xi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 I += twodfs.data[yi*twodfs.NrBins+xi];
    }
      if(yi == 0) {
 ppgmove(2.0*I/scale, y);
      }else {
 ppgdraw(2.0*I/scale, y);
      }
    }
  }else if(Number == 1) {
    for(yi = 0; yi < twodfs2.NrSubints; yi++) {
      I = 0;
      for(xi = 0; xi < twodfs2.NrBins; xi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 I += twodfs2.data[yi*twodfs2.NrBins+xi];
    }
      if(yi == 0) {
 ppgmove(2.0*I/scale, y);
      }else {
 ppgdraw(2.0*I/scale, y);
      }
    }
  }
}
void PlotLRFS(int plot_xlabel, int plot_ylabel, int noside, int plot_ylabeltop, int lineColour, int showwedge, char *title, double scaleFig_x, double scaleFig_y, int normaliseSide, int nointegrateNumbers, int usephase, verbose_definition verbose)
{
  float I;
  int xi, i;
  ppgsvp(0.2, 0.2+0.11*scaleFig_x, 0.95-0.3*scaleFig_y, 0.95-0.15*scaleFig_y);
  if(title != NULL) {
    ppgsch(0.55*labelscale);
    ppgscf(2);
    ppgslw(2);
    ppgmtxt("t",1,0.5,0.5,title);
    ppgslw(1);
    ppgscf(1);
    ppgsch(0.38*labelscale);
  }
  if(SelectLRegion != 0) {
    fl_min = LRegionLow;
    fl_max = LRegionHigh;
  }
  ppgswin(fl_min,fl_max,f3_min,f3_max);
  if(plot_xlabel != 0) {
    ppgsch(0.3*labelscale);
    if(usephase)
      ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (phase)");
    else
      ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (deg)");
    ppgsch(0.38*labelscale);
  }
  if(plot_ylabel != 0) {
    ppgsch(0.3*labelscale);
    if(noside) {
      if(plot_ylabel == 1)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plot_ylabel == 2)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plot_ylabel == 3)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }else {
      float offset;
      offset = 10 - (labelscale-1.0)*(7-2.8);
      offset += 7*(scaleFig_x-1.0);
      if(plot_ylabel == 1)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plot_ylabel == 2)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plot_ylabel == 3)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }
    ppgsch(0.38*labelscale);
  }
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);
  pgplot_options.box.box_labelsize = 0.3*labelscale;
  pgplot_options.viewport.dontopen = 1;
  pgplot_options.viewport.dontclose = 1;
  pgplot_options.viewport.noclear = 1;
  double period;
  int ret;
  ret = get_period(lrfs, 0, &period, verbose);
  if(ret == 2) {
    printerror(verbose.debug, "ERROR pspecFig (%s): Cannot obtain period", lrfs.filename);
    exit(0);
  }
  double xright;
  xright = 360*(lrfs.NrBins-1)*get_tsamp(lrfs, 0, verbose)/period;
  if(usephase)
    xright /= 360.0;
  pgplotMap(&pgplot_options, lrfs.data, lrfs.NrBins, lrfs.NrSubints, 0+dl, xright+dl, fl_min, fl_max, 0, 0.5, f3_min, f3_max, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, oversaturizel, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, showwedge, 0, 0, verbose);
  ret = get_period(AverageProfile, 0, &period, verbose);
  if(ret == 2) {
    printerror(verbose.debug, "ERROR pspecFig (%s): Cannot obtain period", AverageProfile.filename);
    exit(0);
  }
  if(overlaypp) {
    ppgsci(lineColour);
    ppgsls(1);
    ppgslw(1);
    i = 0;
    for(xi=0; xi < AverageProfile.NrBins; xi++) {
      double xpos;
      xpos = xi*get_tsamp(AverageProfile, 0, verbose)*360.0/period;
      if(usephase)
 xpos /= 360.0;
      if(xpos+dl >= fl_min && xpos+dl <= fl_max) {
 I = AverageProfile.data[xi]*0.5;
 if(i == 0) {
   ppgmove(xpos+dl,I);
   i = 1;
 }else {
   ppgdraw(xpos+dl,I);
 }
      }
    }
    ppgsci(1);
  }
  ppgswin(fl_min,fl_max,f3_min,f3_max);
  ppgsch(0.38*labelscale);
  if(noside) {
    if(inside) {
      if(noylabels)
 ppgbox("bcnst",0.0,0,"bcst",0.0,0);
      else
 ppgbox("bcnst",0.0,0,"bcnst",0.0,0);
    }else {
      if(noylabels)
 ppgbox("bcnsti",0.0,0,"bcsti",0.0,0);
      else
 ppgbox("bcnsti",0.0,0,"bcnsti",0.0,0);
    }
  }else {
    if(inside) {
      ppgbox("bcnst",0.0,0,"cst",0.0,0);
    }else {
      ppgbox("bcnsti",0.0,0,"csti",0.0,0);
    }
  }
  if(noside == 0)
    IntegrateSubsetHorizontal(-1, scaleFig_x, scaleFig_y, normaliseSide, nointegrateNumbers);
  ppgsch(0.38*labelscale);
}
void IntegrateSubsetVertical(int Number, double scaleFig_x, double scaleFig_y, int normalise, int nointegrateNumbers)
{
  float I, Imin, Imax, x, y, offset;
  int xi, yi;
  offset = 0+0.03*(labelscale-1.0)*((float)Number+1.0);
  ppgsvp(0.2, 0.2+0.11*scaleFig_x, 0.95-0.52*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y, 0.95-0.48*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y);
  GetExtremesSubsetVertical(&Imin, &Imax, Number);
  double scale;
  scale = fabs(Imax);
  if(fabs(Imin) > scale)
    scale = fabs(Imin);
  if(normalise == 0)
    scale = 1.0;
  if(Number == 0)
    ppgswin(f2_min,f2_max,Imin/scale,1.05*Imax/scale);
  else
    ppgswin(f2_min2,f2_max2,Imin/scale,1.05*Imax/scale);
  if(noxlabels == 0)
    ppgbox("bnst",0.0,0,"bcvi",0.0,0);
  else
    ppgbox("bst",0.0,0,"bcvi",0.0,0);
  ppgbox("c",0.0,0,"",0.0,0);
  y = floor(log10(Imax*0.7/scale));
  x = floor(Imax*0.7/(pow(10,y)*scale));
  x = x*pow(10,y);
  char labelnumbers[3];
  if(nointegrateNumbers == 0) {
    strcpy(labelnumbers, "n");
  }else {
    strcpy(labelnumbers, "");
  }
  ppgsch(0.38*labelscale*0.66);
  if(Number == 0)
    ppgaxis(labelnumbers,f2_min, 0,f2_min, Imax*.7/scale,0,Imax*.7/scale,x,1,0.3,0,0,-0.8,90);
  else
    ppgaxis(labelnumbers,f2_min2,0,f2_min2,Imax*.7/scale,0,Imax*.7/scale,x,1,0.3,0,0,-0.8,90);
  ppgsch(0.38*labelscale);
  if(Number == 0) {
    for(xi = 0; xi < twodfs.NrBins; xi++) {
      I = 0;
      for(yi = 0; yi < twodfs.NrSubints; yi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(y >= P3IntegrateLow && y <= P3IntegrateHigh)
   I += twodfs.data[yi*twodfs.NrBins+xi];
      }
      if(xi == 0) {
 ppgmove(x,2.0*I/scale);
      }else {
 ppgdraw(x,2.0*I/scale);
      }
    }
  }else if(Number == 1) {
    for(xi = 0; xi < twodfs2.NrBins; xi++) {
      I = 0;
      for(yi = 0; yi < twodfs2.NrSubints; yi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(y >= P3IntegrateLow2 && y <= P3IntegrateHigh2)
   I += twodfs2.data[yi*twodfs2.NrBins+xi];
      }
      if(xi == 0) {
 ppgmove(x,2.0*I/scale);
      }else {
 ppgdraw(x,2.0*I/scale);
      }
    }
  }
}
void IntegrateSubsetVertical_2dfsOnly(int Number, double scaleFig_x, double scaleFig_y, int normalise, int nointegrateNumbers, verbose_definition verbose)
{
  float I, Imin, Imax, x, y, offset;
  int xi, yi;
  offset = 0+0.03*(labelscale-1.0)*((float)Number+1.0);
  ppgsvp(0.2, 0.2+0.11*scaleFig_x, 0.95-0.37*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y, 0.95-0.33*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y);
  GetExtremesSubsetVertical(&Imin, &Imax, Number);
  double scale;
  scale = fabs(Imax);
  if(fabs(Imin) > scale)
    scale = fabs(Imin);
  if(normalise == 0)
    scale = 1.0;
  if(Number == 0)
    ppgswin(f2_min,f2_max,Imin/scale,1.05*Imax/scale);
  else
    ppgswin(f2_min2,f2_max2,Imin/scale,1.05*Imax/scale);
  if(noxlabels == 0)
    ppgbox("bnst",0.0,0,"bcvi",0.0,0);
  else
    ppgbox("bst",0.0,0,"bcvi",0.0,0);
  ppgbox("c",0.0,0,"",0.0,0);
  y = floor(log10(Imax*0.7/scale));
  x = floor(Imax*0.7/(pow(10,y)*scale));
  x = x*pow(10,y);
  char labelnumbers[3];
  if(nointegrateNumbers == 0) {
    strcpy(labelnumbers, "n");
  }else {
    strcpy(labelnumbers, "");
  }
  ppgsch(0.38*labelscale*0.66);
  if(Number == 0)
    ppgaxis(labelnumbers,f2_min, 0,f2_min, Imax*.7/scale,0,Imax*.7/scale,x,1,0.3,0,0,-0.8,90);
  else
    ppgaxis(labelnumbers,f2_min2,0,f2_min2,Imax*.7/scale,0,Imax*.7/scale,x,1,0.3,0,0,-0.8,90);
  ppgsch(0.38*labelscale);
  if(Number == 0) {
    for(xi = 0; xi < twodfs.NrBins; xi++) {
      I = 0;
      for(yi = 0; yi < twodfs.NrSubints; yi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(y >= P3IntegrateLow && y <= P3IntegrateHigh)
   I += twodfs.data[yi*twodfs.NrBins+xi];
      }
      if(xi == 0) {
 ppgmove(x,2.0*I);
      }else {
 ppgdraw(x,2.0*I);
      }
    }
  }else if(Number == 1) {
    for(xi = 0; xi < twodfs2.NrBins; xi++) {
      I = 0;
      for(yi = 0; yi < twodfs2.NrSubints; yi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(y >= P3IntegrateLow2 && y <= P3IntegrateHigh2)
   I += twodfs2.data[yi*twodfs2.NrBins+xi];
      }
      if(xi == 0) {
 ppgmove(x,2.0*I);
      }else {
 ppgdraw(x,2.0*I);
      }
    }
  }
}
void Plot2dfs(int Number, int plot_xlabel, int plot_ylabel, int noside, int plot_ylabeltop, int nomain, double scaleFig_x, double scaleFig_y, int showwedge, int normaliseSide, int nointegrateNumbers, verbose_definition verbose)
{
  float offset;
  offset = 0+0.03*(labelscale-1.0)*((float)Number+1.0);
  ppgsvp(0.2, 0.2+0.11*scaleFig_x, 0.95-0.48*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y, 0.95-0.33*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y);
  if(SelectP2Region != 0) {
    if(Number == 0) {
      f2_min = P2RegionLow;
      f2_max = P2RegionHigh;
    }else {
      f2_min2 = P2RegionLow2;
      f2_max2 = P2RegionHigh2;
    }
  }
  ppgswin(f2_min,f2_max,f3_min,f3_max);
  if(plot_ylabel != 0) {
    ppgsch(0.3*labelscale);
    if(noside) {
      if(plot_ylabel == 1)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plot_ylabel == 2)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plot_ylabel == 3)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }else {
      float offset;
      offset = 10 - (labelscale-1.0)*(7-2.8);
      offset += 7*(scaleFig_x-1.0);
      if(plot_ylabel == 1)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plot_ylabel == 2)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plot_ylabel == 3)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }
    ppgsch(0.38*labelscale);
  }
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);
  pgplot_options.box.box_labelsize = 0.3*labelscale;
  pgplot_options.viewport.noclear = 1;
  pgplot_options.viewport.dontopen = 1;
  pgplot_options.viewport.dontclose = 1;
  if(!nomain) {
    if(Number == 0) {
      if(doFlip) {
 pgplotMap(&pgplot_options, twodfs.data, twodfs.NrBins, twodfs.NrSubints, +AverageProfile.NrBins/2.0+0.5*AverageProfile.NrBins/(float)twodfs.NrBins, -AverageProfile.NrBins/2.0+0.5*AverageProfile.NrBins/(float)twodfs.NrBins, f2_min, f2_max, 0, 0.5, f3_min, f3_max, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, oversaturize, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, showwedge, 0, 0, verbose);
      }else {
 pgplotMap(&pgplot_options, twodfs.data, twodfs.NrBins, twodfs.NrSubints, -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins, +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins, f2_min, f2_max, 0, 0.5, f3_min, f3_max, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, oversaturize, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, showwedge, 0, 0, verbose);
      }
    }else {
      if(doFlip)
 pgplotMap(&pgplot_options, twodfs2.data, twodfs2.NrBins, twodfs2.NrSubints, +AverageProfile.NrBins/2.0+0.5*AverageProfile.NrBins/(float)twodfs2.NrBins, -AverageProfile.NrBins/2.0+0.5*AverageProfile.NrBins/(float)twodfs2.NrBins, f2_min2, f2_max2, 0, 0.5, f3_min2, f3_max2, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, oversaturize2, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, showwedge, 0, 0, verbose);
      else
 pgplotMap(&pgplot_options, twodfs2.data, twodfs2.NrBins, twodfs2.NrSubints, -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs2.NrBins, +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs2.NrBins, f2_min2, f2_max2, 0, 0.5, f3_min2, f3_max2, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, oversaturize2, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, showwedge, 0, 0, verbose);
    }
    ppgsch(0.38*labelscale);
    if(noside) {
      if(inside) {
 if(noylabels)
   ppgbox("cst",0.0,0,"bcst",0.0,0);
 else
   ppgbox("cst",0.0,0,"bcnst",0.0,0);
      }else {
 if(noylabels)
   ppgbox("csti",0.0,0,"bcsti",0.0,0);
 else
   ppgbox("csti",0.0,0,"bcnsti",0.0,0);
      }
    }else {
      if(inside) {
 ppgbox("cst",0.0,0,"cst",0.0,0);
      }else {
 ppgbox("csti",0.0,0,"csti",0.0,0);
      }
    }
    ppgbox("c",0.0,0,"c",0.0,0);
    if(SelectP3Integrate != 0) {
      ppgsls(2);
      ppgslw(1);
      if(Number == 0) {
 ppgmove(f2_min, P3IntegrateLow);
 ppgdraw(f2_max, P3IntegrateLow);
      }else {
 ppgmove(f2_min, P3IntegrateLow2);
 ppgdraw(f2_max, P3IntegrateLow2);
      }
      ppgsls(2);
      if(Number == 0) {
 ppgmove(f2_min, P3IntegrateHigh);
 ppgdraw(f2_max, P3IntegrateHigh);
      }else {
 ppgmove(f2_min, P3IntegrateHigh2);
 ppgdraw(f2_max, P3IntegrateHigh2);
      }
      ppgsls(1);
      ppgslw(1);
    }
  }
  IntegrateSubsetVertical(Number, scaleFig_x, scaleFig_y, normaliseSide, nointegrateNumbers);
  if(plot_xlabel == 1) {
    ppgsch(0.3*labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (cpp)");
    ppgsch(0.38*labelscale);
  }else if(plot_xlabel == 2) {
    ppgsch(0.3*labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d2\\u)");
    ppgsch(0.38*labelscale);
  }else if(plot_xlabel == 3) {
    ppgsch(0.3*labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (P/P\\d2\\u)");
    ppgsch(0.38*labelscale);
  }
  if(!noside)
    IntegrateSubsetHorizontal(Number, scaleFig_x, scaleFig_y, normaliseSide, nointegrateNumbers);
}
void Plot2dfsOnly(int Number, int plot_xlabel, int plot_ylabel, int noside, int nomain, int showwedge, char *title, double scaleFig_x, double scaleFig_y, int normaliseSide, int nointegrateNumbers, verbose_definition verbose)
{
  float offset;
  offset = 0+0.03*(labelscale-1.0)*((float)Number+1.0);
  ppgsvp(0.2, 0.2+0.11*scaleFig_x, 0.95-0.33*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y, 0.95-0.18*scaleFig_y-0.22*Number*scaleFig_y-offset*scaleFig_y);
  if(title != NULL) {
    ppgsch(0.55*labelscale);
    ppgscf(2);
    ppgslw(2);
    ppgmtxt("t",1,0.5,0.5,title);
    ppgslw(1);
    ppgscf(1);
    ppgsch(0.38*labelscale);
  }
  if(SelectP2Region != 0) {
    if(Number == 0) {
      f2_min = P2RegionLow;
      f2_max = P2RegionHigh;
    }else {
      f2_min2 = P2RegionLow2;
      f2_max2 = P2RegionHigh2;
    }
  }
  ppgswin(f2_min,f2_max,f3_min,f3_max);
  if(plot_ylabel != 0) {
    ppgsch(0.3*labelscale);
    if(noside) {
      if(plot_ylabel == 1)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plot_ylabel == 2)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plot_ylabel == 3)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }else {
      float offset;
      offset = 10 - (labelscale-1.0)*(7-2.8);
      offset += 7*(scaleFig_x-1.0);
      if(plot_ylabel == 1)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plot_ylabel == 2)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plot_ylabel == 3)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }
    ppgsch(0.38*labelscale);
  }
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);
  pgplot_options.box.box_labelsize = 0.3*labelscale;
  pgplot_options.viewport.noclear = 1;
  pgplot_options.viewport.dontopen = 1;
  pgplot_options.viewport.dontclose = 1;
  if(!nomain) {
    if(Number == 0) {
      if(doFlip)
 pgplotMap(&pgplot_options, twodfs.data, twodfs.NrBins, twodfs.NrSubints, +AverageProfile.NrBins/2.0+0.5*AverageProfile.NrBins/(float)twodfs.NrBins, -AverageProfile.NrBins/2.0+0.5*AverageProfile.NrBins/(float)twodfs.NrBins, f2_min, f2_max, 0, 0.5, f3_min, f3_max, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, oversaturize, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, showwedge, 0, 0, verbose);
      else
 pgplotMap(&pgplot_options, twodfs.data, twodfs.NrBins, twodfs.NrSubints, -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins, +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins, f2_min, f2_max, 0, 0.5, f3_min, f3_max, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, oversaturize, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, showwedge, 0, 0, verbose);
    }else {
      if(doFlip)
 pgplotMap(&pgplot_options, twodfs2.data, twodfs2.NrBins, twodfs2.NrSubints, +AverageProfile.NrBins/2.0+0.5*AverageProfile.NrBins/(float)twodfs2.NrBins, -AverageProfile.NrBins/2.0+0.5*AverageProfile.NrBins/(float)twodfs2.NrBins, f2_min2, f2_max2, 0, 0.5, f3_min2, f3_max2, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, oversaturize2, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, showwedge, 0, 0, verbose);
      else
 pgplotMap(&pgplot_options, twodfs2.data, twodfs2.NrBins, twodfs2.NrSubints, -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs2.NrBins, +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs2.NrBins, f2_min2, f2_max2, 0, 0.5, f3_min2, f3_max2, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, oversaturize2, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, showwedge, 0, 0, verbose);
    }
    ppgsch(0.38*labelscale);
    if(noside) {
      if(inside) {
 if(noylabels)
   ppgbox("cst",0.0,0,"bcst",0.0,0);
 else
   ppgbox("cst",0.0,0,"bcnst",0.0,0);
      }else {
 if(noylabels)
   ppgbox("csti",0.0,0,"bcsti",0.0,0);
 else
   ppgbox("csti",0.0,0,"bcnsti",0.0,0);
      }
    }else {
      if(inside) {
 ppgbox("cst",0.0,0,"cst",0.0,0);
      }else {
 ppgbox("csti",0.0,0,"csti",0.0,0);
      }
    }
    if(SelectP3Integrate != 0) {
      ppgsls(2);
      ppgslw(1);
      if(Number == 0) {
 ppgmove(f2_min, P3IntegrateLow);
 ppgdraw(f2_max, P3IntegrateLow);
      }else {
 ppgmove(f2_min, P3IntegrateLow2);
 ppgdraw(f2_max, P3IntegrateLow2);
      }
      ppgsls(2);
      if(Number == 0) {
 ppgmove(f2_min, P3IntegrateHigh);
 ppgdraw(f2_max, P3IntegrateHigh);
      }else {
 ppgmove(f2_min, P3IntegrateHigh2);
 ppgdraw(f2_max, P3IntegrateHigh2);
      }
      ppgsls(1);
      ppgslw(1);
    }
  }
  IntegrateSubsetVertical_2dfsOnly(Number, scaleFig_x, scaleFig_y, normaliseSide, nointegrateNumbers, verbose);
  if(plot_xlabel == 1) {
    ppgsch(0.3*labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (cpp)");
    ppgsch(0.38*labelscale);
  }else if(plot_xlabel == 2) {
    ppgsch(0.3*labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d2\\u)");
    ppgsch(0.38*labelscale);
  }else if(plot_xlabel == 3) {
    ppgsch(0.3*labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (P/P\\d2\\u)");
    ppgsch(0.38*labelscale);
  }
  if(!noside)
    IntegrateSubsetHorizontal_2dfsOnly(Number, scaleFig_x, scaleFig_y, normaliseSide, nointegrateNumbers);
}

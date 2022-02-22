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
#include <gsl/gsl_sort.h>
#include "psrsalsa.h"
typedef struct {
  datafile_definition datafile;
  double f2_min, f2_max;
  double f3_min, f3_max;
  double P3IntegrateLow, P3IntegrateHigh;
  int SelectP3Integrate;
  double percentage_power_drift, percentage_power_drift_error;
  double oversaturize;
}twodfs_def;
typedef struct {
  char *textside;
  double fl_min, fl_max;
  double dl;
  double LRegionLow, LRegionHigh;
  int SelectLRegion;
  int doFlip;
  int ExtraVerticalMaximaSkip;
  int noside;
  int overlaypp;
  int inside;
  int noylabels;
  int noxlabels;
  int notop;
  int autooversaturizel, autooversaturize2, autozoomP2;
  float labelscale, titlescale, titledy, textsidescale;
  double oversaturizel;
  double Imin, Imax;
  double scaleFig_x, scaleFig_y;
  double markmodindex, markmodindex_error, markmodindex_shift;
  int plot_xlabel, plot_ylabel, plot_ylabeltop;
  int nomain;
  int showwedge;
  int showwedge_max;
  int normaliseSide;
  int nointegrateNumbers;
  int intflip;
  int usephase;
  int lineColour;
  int normalise_spectra;
  int SelectP3Region;
  int domarkdriftpower;
}plotoptions_def;
int loadLRFS(datafile_definition *lrfs, int extprefix, int longsnap, int argc, char **argv, plotoptions_def *plotoptions, verbose_definition verbose);
int load2dfs(twodfs_def *twodfs_allinfo, datafile_definition AverageProfile, int file_number, int extprefix, int altname, int argc, char **argv, plotoptions_def plotoptions, verbose_definition verbose);
void Plot2dfs(twodfs_def twodfs_allinfo, twodfs_def twodfs2_allinfo, datafile_definition AverageProfile, int twodfsonly, int Number, char *title, plotoptions_def plotoptions, int argc, char **argv, verbose_definition verbose);
void PlotLRFS(datafile_definition lrfs, datafile_definition AverageProfile, twodfs_def twodfs_allinfo, char *title, plotoptions_def plotoptions, verbose_definition verbose);
int loadHeaderPulseStack(datafile_definition *AverageProfile, datafile_definition lrfs, double *period, int type_of_plots, int argc, char **argv, verbose_definition verbose);
void autozoomP2(twodfs_def *twodfs_allinfo, datafile_definition AverageProfile, plotoptions_def plotoptions, verbose_definition verbose);
int main(int argc, char **argv)
{
  int i, j, xi, SelectP2Region, LoadTwo, NrSelectedOverSaturize;
  int file_number, ImaxSet, IminSet, type_of_plots, maxSubpulsePhaseSet, minSubpulsePhaseSet, do_phase_slope, extprefix, domarkavmod;
  char PlotDevice[100], filename[1000], txt[1000];
  int title_index, Load2dfs, LoadLRFS, plotvariance, plotmodindex, ok_flag, altProf, lineStyle, ret, longsnap;
  double P3RegionLow, P3RegionHigh, maxSubpulsePhase, minSubpulsePhase;
  float I, x, profilescale;
  float maxvalue_mod, maxvalue_stddev, maxsigma_stddev, maxsigma_mod, ImaxValue, IminValue, phase_slope_g, phase_slope_o;
  FILE *fout_ascii;
  long k;
  datafile_definition subpulseTrackProfile, subpulseTrackProfileErr, subpulseAmpProfile;
  datafile_definition AverageProfile, VarianceProfile, ModProfile, VarianceProfileErr, ModProfileErr, lrfs;
  verbose_definition noverbose;
  psrsalsaApplication application;
  plotoptions_def plotoptions;
  twodfs_def twodfs_allinfo, twodfs2_allinfo;
  initApplication(&application, "pspecFig", "");
  type_of_plots = 0;
  file_number = 1;
  plotoptions.scaleFig_x = 1;
  plotoptions.scaleFig_y = 1;
  do_phase_slope = 0;
  plotoptions.autooversaturizel = 0;
  plotoptions.autooversaturize2 = 0;
  plotoptions.autozoomP2 = 0;
  plotoptions.markmodindex = -1;
  plotoptions.markmodindex_error = -1;
  plotoptions.markmodindex_shift = 0.0;
  extprefix = 0;
  longsnap = 0;
  domarkavmod = 0;
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
    printf("-long  \"low high\" Set horizontal range shown lrfs/profile.\n");
    printf("-phase            Use pulse longitude in phase rather than degrees.\n");
    printf("-dlong            Shift lrfs/profile by this amount of degrees or phase.\n");
    printf("\nMode A specific range options:\n");
    printf("-p3 \"low high\"    Set vertical range shown in lrfs/2dfs in cpp.\n");
    printf("-p2 \"low high\"    Set horizontal range as shown in 2dfs (can use this option\n");
    printf("                  twice if -2 is used) in cpp.\n");
    printf("-int  \"low high\"  Select vertical integration range in 2dfs in cpp\n");
    printf("                  (affects the side panels, can be use twice if -2 is used).\n");
    printf("-noylabels        Don't show ylabels\n");
    printf("-linestyle        Set the PGPLOT line style of the pulse profile\n");
    printf("-linecolor        Set the PGPLOT line color of the pulse profile\n");
    printf("-scalefig \"x y\"   Scale size of panel with factors x and y\n");
    printf("-title            Set the title\n");
    printf("-ytop             Show y-label of top plot\n");
    printf("\nMode A specific graphics options:\n");
    printf("-intflip          Add a flipped version of the vertical integration of the 2dfs\n");
    printf("-intnrs           Show nrs along side panels axis, instead of just a tick\n");
    printf("-scalel scale     The values in the lrfs are boosted by factor scale, resulting\n");
    printf("                  in clipping, which can highlight weaker features.\n");
    printf("-scale2 scale     As -scalel, but for 2dfs instead of lrfs.\n");
    printf("                  When this option is used twice, the 2nd time the\n");
    printf("                  option is used applies to the second 2dfs shown.\n");
    printf("-scalep           Scale profile by this factor\n");
    printf("-noflip           Do not flip 2DFS horizontally. If specified positive drift\n");
    printf("                  corresponds to power in the negative side of the diagram.\n");
    printf("-overlay          Overlay pulse profile over LRFS\n");
    printf("-noxlabels        Don't show xlabels on 2dfs bottom integration panel.\n");
    printf("-nomod            Do not plot a modulation index profile\n");
    printf("-nostddev         Do not plot a standard deviation profile\n");
    printf("-showwedge        Plot an annotated wedge to show color scale\n");
    printf("-showwedge_showmax  Plot a label indicating the maximum value in the data.\n");
    printf("                  With the -scalel or -scale2 option this is not necessarily\n");
    printf("                  the maximum of the colour scale.\n");
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
    printf("-noflip           Do not change sign of subpulse phase. If specified positive\n");
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
  plotoptions.plot_ylabeltop = 0;
  plotoptions.nomain = 0;
  plotoptions.SelectP3Region = 0;
  plotoptions.domarkdriftpower = 0;
  SelectP2Region = 0;
  plotoptions.SelectLRegion = 0;
  twodfs_allinfo.SelectP3Integrate = 0;
  twodfs2_allinfo.SelectP3Integrate = 0;
  LoadTwo = 0;
  strcpy(PlotDevice, "?");
  twodfs_allinfo.oversaturize = 1;
  twodfs2_allinfo.oversaturize = 1;
  twodfs_allinfo.percentage_power_drift = 0;
  twodfs_allinfo.percentage_power_drift_error = -1;
  twodfs2_allinfo.percentage_power_drift = 0;
  twodfs2_allinfo.percentage_power_drift_error = -1;
  plotoptions.oversaturizel = 1;
  NrSelectedOverSaturize = 0;
  maxvalue_mod = -1;
  maxsigma_mod = 3;
  maxvalue_stddev = -1;
  maxsigma_stddev = 3;
  title_index = -1;
  plotoptions.textside = NULL;
  plotoptions.doFlip = 1;
  plotoptions.ExtraVerticalMaximaSkip = 1;
  Load2dfs = 1;
  LoadLRFS = 1;
  plotvariance = 1;
  plotmodindex = 1;
  profilescale = 1;
  plotoptions.plot_xlabel = 0;
  plotoptions.plot_ylabel = 0;
  plotoptions.notop = 0;
  plotoptions.noside = 0;
  plotoptions.overlaypp = 0;
  plotoptions.inside = 0;
  plotoptions.labelscale = 1.0;
  plotoptions.titlescale = 1.0;
  plotoptions.textsidescale = 1.0;
  plotoptions.noylabels = 0;
  plotoptions.noxlabels = 0;
  plotoptions.titledy = 0;
  altProf = 0;
  plotoptions.dl = 0;
  ImaxSet = 0;
  IminSet = 0;
  lineStyle = 1;
  plotoptions.lineColour = 1;
  maxSubpulsePhaseSet = 0;
  minSubpulsePhaseSet = 0;
  maxSubpulsePhase = 0;
  minSubpulsePhase = 0;
  plotoptions.showwedge = 0;
  plotoptions.showwedge_max = 0;
  plotoptions.normalise_spectra = 0;
  plotoptions.normaliseSide = 0;
  plotoptions.nointegrateNumbers = 1;
  plotoptions.intflip = 0;
  plotoptions.usephase = 0;
  twodfs_allinfo.f2_min = 0.0;
  twodfs_allinfo.f2_max = 0.0;
  twodfs2_allinfo.f2_min = 0.0;
  twodfs2_allinfo.f2_max = 0.0;
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
      plotoptions.plot_xlabel = 1;
    }else if(strcmp(argv[i], "-xlabel2") == 0
      ) {
      plotoptions.plot_xlabel = 2;
    }else if(strcmp(argv[i], "-xlabel3") == 0
      ) {
      plotoptions.plot_xlabel = 3;
    }else if(strcmp(argv[i], "-ylabel") == 0
      ) {
      plotoptions.plot_ylabel = 1;
    }else if(strcmp(argv[i], "-ylabel2") == 0
      ) {
      plotoptions.plot_ylabel = 2;
    }else if(strcmp(argv[i], "-ylabel3") == 0
      ) {
      plotoptions.plot_ylabel = 3;
    }else if(strcmp(argv[i], "-nolrfs") == 0
      ) {
      LoadLRFS = 0;
    }else if(strcmp(argv[i], "-2") == 0) {
      LoadTwo = 1;
    }else if(strcmp(argv[i], "-no2dfs") == 0
      ) {
      Load2dfs = 0;
    }else if(strcmp(argv[i], "-inside") == 0) {
      plotoptions.inside = 1;
    }else if(strcmp(argv[i], "-nostddev") == 0
      ) {
      plotvariance = 0;
    }else if(strcmp(argv[i], "-nomod") == 0
      ) {
      plotmodindex = 0;
    }else if(strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "-noflip") == 0) {
      plotoptions.doFlip = 0;
    }else if(strcmp(argv[i], "-notop") == 0) {
      plotoptions.notop = 1;
    }else if(strcmp(argv[i], "-noside") == 0) {
      plotoptions.noside = 1;
    }else if(strcmp(argv[i], "-nomain") == 0) {
      plotoptions.nomain = 1;
    }else if(strcmp(argv[i], "-overlay") == 0) {
      plotoptions.overlaypp = 1;
    }else if(strcmp(argv[i], "-noylabels") == 0) {
      plotoptions.noylabels = 1;
    }else if(strcmp(argv[i], "-noxlabels") == 0) {
      plotoptions.noxlabels = 1;
    }else if(strcmp(argv[i], "-showwedge") == 0) {
      plotoptions.showwedge = 1;
    }else if(strcmp(argv[i], "-showwedge_showmax") == 0) {
      plotoptions.showwedge_max = 1;
    }else if(strcmp(argv[i], "-normspectra") == 0) {
      plotoptions.normalise_spectra = 1;
    }else if(strcmp(argv[i], "-normside") == 0) {
      plotoptions.normaliseSide = 1;
    }else if(strcmp(argv[i], "-intnrs") == 0) {
      plotoptions.nointegrateNumbers = 0;
    }else if(strcmp(argv[i], "-intflip") == 0) {
      plotoptions.intflip = 1;
    }else if(strcmp(argv[i], "-phase") == 0) {
      plotoptions.usephase = 1;
    }else if(strcmp(argv[i], "-linestyle") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &lineStyle, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-linecolor") == 0 || strcmp(argv[i], "-linecolour") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &(plotoptions.lineColour), NULL) == 0) {
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
    }else if(strcmp(argv[i], "-dlong") == 0 || strcmp(argv[i], "-dl") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &(plotoptions.dl), NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-labelscale") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &(plotoptions.labelscale), NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-titlescale") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &(plotoptions.titlescale), NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-titledy") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &(plotoptions.titledy), NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-textsidescale") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &(plotoptions.textsidescale), NULL) == 0) {
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
      plotoptions.SelectP3Region = 1;
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &P3RegionLow, &P3RegionHigh, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-p2") == 0) {
      if(SelectP2Region == 0) {
 SelectP2Region = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &(twodfs_allinfo.f2_min), &(twodfs_allinfo.f2_max), NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 twodfs2_allinfo.f2_min = twodfs_allinfo.f2_min;
 twodfs2_allinfo.f2_max = twodfs_allinfo.f2_max;
      }else {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &(twodfs2_allinfo.f2_min), &(twodfs2_allinfo.f2_max), NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
      }
      i++;
    }else if(strcmp(argv[i], "-long") == 0 || strcmp(argv[i], "-l") == 0) {
      plotoptions.SelectLRegion = 1;
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &(plotoptions.LRegionLow), &(plotoptions.LRegionHigh), NULL) == 0) {
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
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &(plotoptions.oversaturizel), NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-scale2") == 0
      ) {
      if(NrSelectedOverSaturize == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &(twodfs_allinfo.oversaturize), NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 twodfs2_allinfo.oversaturize = twodfs_allinfo.oversaturize;
 NrSelectedOverSaturize = 1;
      }else {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &(twodfs2_allinfo.oversaturize), NULL) == 0) {
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
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &(plotoptions.ExtraVerticalMaximaSkip), NULL) == 0) {
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
      plotoptions.textside = argv[i+1];
      i++;
    }else if(strcmp(argv[i], "-2dfsnr") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &file_number, NULL) == 0) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-int") == 0
      ) {
      if(twodfs_allinfo.SelectP3Integrate == 0) {
 twodfs_allinfo.SelectP3Integrate = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &(twodfs_allinfo.P3IntegrateLow), &(twodfs_allinfo.P3IntegrateHigh), NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 twodfs2_allinfo.SelectP3Integrate = 1;
 twodfs2_allinfo.P3IntegrateLow = twodfs_allinfo.P3IntegrateLow;
 twodfs2_allinfo.P3IntegrateHigh = twodfs_allinfo.P3IntegrateHigh;
      }else {
 twodfs2_allinfo.SelectP3Integrate = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &(twodfs2_allinfo.P3IntegrateLow), &(twodfs2_allinfo.P3IntegrateHigh), NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
      }
      i++;
    }else if(strcmp(argv[i], "-ytop") == 0
      ) {
      plotoptions.plot_ylabeltop = 1;
    }else if(strcmp(argv[i], "-scalefig") == 0) {
      if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &(plotoptions.scaleFig_x), &(plotoptions.scaleFig_y), NULL) == 0) {
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
  copyVerboseState(application.verbose_state, &noverbose);
  noverbose.verbose = 0;
  cleanPSRData(&(twodfs_allinfo.datafile), application.verbose_state);
  cleanPSRData(&(twodfs2_allinfo.datafile), application.verbose_state);
  cleanPSRData(&lrfs, application.verbose_state);
  cleanPSRData(&AverageProfile, application.verbose_state);
  cleanPSRData(&VarianceProfile, application.verbose_state);
  cleanPSRData(&VarianceProfileErr, application.verbose_state);
  cleanPSRData(&ModProfile, application.verbose_state);
  cleanPSRData(&ModProfileErr, application.verbose_state);
  cleanPSRData(&subpulseTrackProfile, application.verbose_state);
  cleanPSRData(&subpulseTrackProfileErr, application.verbose_state);
  cleanPSRData(&subpulseAmpProfile, application.verbose_state);
  if(type_of_plots == 0) {
    if(loadLRFS(&lrfs, extprefix, longsnap, argc, argv, &plotoptions, application.verbose_state) == 0) {
      closePSRData(&AverageProfile, 0, 0, application.verbose_state);
      closePSRData(&(twodfs_allinfo.datafile), 0, 0, application.verbose_state);
      closePSRData(&(twodfs2_allinfo.datafile), 0, 0, application.verbose_state);
      closePSRData(&lrfs, 0, 0, application.verbose_state);
      closePSRData(&VarianceProfile, 0, 0, application.verbose_state);
      closePSRData(&ModProfile, 0, 0, application.verbose_state);
      closePSRData(&VarianceProfileErr, 0, 0, application.verbose_state);
      closePSRData(&ModProfileErr, 0, 0, application.verbose_state);
      closePSRData(&subpulseTrackProfile, 0, 0, application.verbose_state);
      closePSRData(&subpulseTrackProfileErr, 0, 0, application.verbose_state);
      closePSRData(&subpulseAmpProfile, 0, 0, application.verbose_state);
      terminateApplication(&application);
      return 0;
    }
  }
  double period;
  if(loadHeaderPulseStack(&AverageProfile, lrfs, &period, type_of_plots, argc, argv, application.verbose_state) == 0) {
    return 0;
  }
  if(extprefix == 0) {
    sprintf(txt, "profile");
  }
  if(change_filename_extension(argv[argc-1], filename, txt, 1000, application.verbose_state) == 0) {
   terminateApplication(&application);
   return 0;
  }
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
  char *line;
  line = malloc(10000);
  if(line == NULL) {
    printerror(application.verbose_state.debug, "ERROR pspecFig: Memory allocation error.");
    return 0;
  }
  if(application.verbose_state.verbose)
    printf("Reading %s\n", filename);
  fout_ascii = fopen(filename, "r");
  if(fout_ascii == NULL) {
    printerror(application.verbose_state.debug, "ERROR pspecFig: Unable to open %s.", filename);
    return 0;
  }
  for(i = 0; i < AverageProfile.NrBins; i++) {
    if(ascii_file_get_next_line(fout_ascii, line, 9999, '#', application.verbose_state) == 0) {
      printwarning(application.verbose_state.debug, "WARNING pspecFig: Read error while reading bin %d from %s.", i+1, filename);
      j = 0;
    }else {
      if(type_of_plots == 0) {
 j = sscanf(line, "%ld %f %f %f %f %f", &k, &AverageProfile.data[i], &VarianceProfile.data[i], &VarianceProfileErr.data[i], &ModProfile.data[i], &ModProfileErr.data[i]);
      }else {
 float junk;
 j = sscanf(line, "%ld %f %f %f %f %f", &k, &AverageProfile.data[i], &junk, &junk, &junk, &junk);
      }
    }
    if(j != 6) {
 printwarning(application.verbose_state.debug, "WARNING: It looks like profile data is rebinned? Check the units.");
 AverageProfile.fixedtsamp *= AverageProfile.NrBins/(double)(i);
 AverageProfile.tsampMode = TSAMPMODE_FIXEDTSAMP;
 AverageProfile.NrBins = i;
 printwarning(application.verbose_state.debug, "WARNING: Assuming the number of bins = %ld and the sampling time = %lf s.", AverageProfile.NrBins, AverageProfile.fixedtsamp);
 if(type_of_plots == 0) {
   copy_params_PSRData(AverageProfile, &VarianceProfile, application.verbose_state);
   copy_params_PSRData(AverageProfile, &ModProfile, application.verbose_state);
   copy_params_PSRData(AverageProfile, &VarianceProfileErr, application.verbose_state);
   copy_params_PSRData(AverageProfile, &ModProfileErr, application.verbose_state);
 }else {
   copy_params_PSRData(AverageProfile, &subpulseTrackProfile, application.verbose_state);
   copy_params_PSRData(AverageProfile, &subpulseTrackProfileErr, application.verbose_state);
   copy_params_PSRData(AverageProfile, &subpulseAmpProfile, application.verbose_state);
 }
 break;
    }
    if(k != i) {
      printerror(application.verbose_state.debug, "Unexpected bin number");
      return 0;
    }
  }
  rewind(fout_ascii);
  int foundavmod;
  foundavmod = 0;
  while(fgets(line, 9999, fout_ascii) != NULL) {
    if(domarkavmod) {
      double avmod, avmoderr;
      if(sscanf(line, "#Average modulation index = %lf +- %lf", &avmod, &avmoderr) == 2) {
 plotoptions.markmodindex = avmod;
 plotoptions.markmodindex_error = avmoderr;
 foundavmod = 1;
      }else {
 if(sscanf(line, "#Average modulation index = %lf", &avmod) == 1) {
   plotoptions.markmodindex = avmod;
   plotoptions.markmodindex_error = -1;
   foundavmod = 1;
 }
      }
    }
  }
  if(domarkavmod && foundavmod == 0) {
    printwarning(application.verbose_state.debug, "WARNING pspecFig (%s): -markavmod is specified, but the average modulation index was not found in data file with modulation index points", AverageProfile.filename);
  }
  fclose(fout_ascii);
  free(line);
  if(application.verbose_state.verbose)
    printf("%ld points read from %s\n", AverageProfile.NrBins, filename);
  ret = get_period(AverageProfile, 0, &period, application.verbose_state);
  if(ret == 2) {
    printerror(application.verbose_state.debug, "ERROR pspecFig (%s): Cannot obtain period", AverageProfile.filename);
    return 0;
  }
  plotoptions.fl_min = 0 + plotoptions.dl;
  plotoptions.fl_max = 360*(AverageProfile.NrBins-1)*get_tsamp(AverageProfile, 0, application.verbose_state)/period + plotoptions.dl;
  if(plotoptions.usephase) {
    plotoptions.fl_min /= 360.0;
    plotoptions.fl_max /= 360.0;
  }
  if(type_of_plots == 1) {
    if(extprefix == 0) {
      sprintf(txt, "amplitude");
    }
    if(change_filename_extension(argv[argc-1], filename, txt, 1000, application.verbose_state) == 0)
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
    if(extprefix == 0) {
      sprintf(txt, "track");
    }
    if(change_filename_extension(argv[argc-1], filename, txt, 1000, application.verbose_state) == 0)
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
  twodfs_allinfo.f3_min = 0;
  twodfs_allinfo.f3_max = 0.5;
  if(type_of_plots == 0 && Load2dfs != 0) {
    if(load2dfs(&twodfs_allinfo, AverageProfile, file_number, extprefix, 0, argc, argv, plotoptions, application.verbose_state) == 0) {
      return 0;
    }
    if(LoadTwo != 0) {
      file_number++;
      if(load2dfs(&twodfs2_allinfo, AverageProfile, file_number, extprefix, 0, argc, argv, plotoptions, application.verbose_state) == 0) {
 return 0;
      }
    }
    if(plotoptions.SelectP3Region != 0) {
      twodfs_allinfo.f3_min = P3RegionLow;
      twodfs_allinfo.f3_max = P3RegionHigh;
      twodfs2_allinfo.f3_min = P3RegionLow;
      twodfs2_allinfo.f3_max = P3RegionHigh;
    }
    if(twodfs_allinfo.SelectP3Integrate == 0) {
      twodfs_allinfo.P3IntegrateLow = twodfs_allinfo.f3_min;
      twodfs_allinfo.P3IntegrateHigh = twodfs_allinfo.f3_max;
      twodfs2_allinfo.P3IntegrateLow = twodfs_allinfo.f3_min;
      twodfs2_allinfo.P3IntegrateHigh = twodfs_allinfo.f3_max;
    }
  }
  ppgopen(PlotDevice);
  ppgask(0);
  ppgslw(1);
  ppgpage();
  ppgslw(1);
  ppgscf(1);
  ppgsch(0.38*plotoptions.labelscale);
  if(plotoptions.notop == 0) {
    ppgsvp(0.2, 0.2+0.11*plotoptions.scaleFig_x, 0.95-0.15*plotoptions.scaleFig_y, 0.95);
    ppgslw(1);
    ppgscf(1);
    ppgsch(0.38*plotoptions.labelscale);
    plotoptions.Imin = 0;
    plotoptions.Imax = 0;
    ret = get_period(AverageProfile, 0, &period, application.verbose_state);
    if(ret == 2) {
      printerror(application.verbose_state.debug, "ERROR pspecFig (%s): Cannot obtain period", AverageProfile.filename);
      return 0;
    }
    for(xi=0; xi < AverageProfile.NrBins; xi++) {
      double xpos;
      xpos = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
      if(plotoptions.usephase)
 xpos /= 360.0;
      if(xpos >= plotoptions.fl_min && xpos <= plotoptions.fl_max) {
 I = AverageProfile.data[xi];
 if(I > plotoptions.Imax)
   plotoptions.Imax = I;
 if(I < plotoptions.Imin)
   plotoptions.Imin = I;
      }
    }
    if(plotoptions.SelectLRegion != 0) {
      plotoptions.fl_min = plotoptions.LRegionLow;
      plotoptions.fl_max = plotoptions.LRegionHigh;
    }
    if(type_of_plots == 0) {
      for(xi=0; xi < ModProfile.NrBins; xi++) {
 double xpos;
 xpos = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
 if(plotoptions.usephase)
   xpos /= 360.0;
 if(xpos+plotoptions.dl >= plotoptions.fl_min && xpos+plotoptions.dl <= plotoptions.fl_max) {
   ok_flag = 1;
   if(maxsigma_mod > 0 && ModProfile.data[xi]/ModProfileErr.data[xi] < maxsigma_mod)
     ok_flag = 0;
   if(maxvalue_mod > 0 && ModProfile.data[xi] > maxvalue_mod)
     ok_flag = 0;
   if(ok_flag) {
     I = ModProfile.data[xi];
     if(I+ModProfileErr.data[xi] > plotoptions.Imax)
       plotoptions.Imax = I+ModProfileErr.data[xi];
     if(I-ModProfileErr.data[xi] < plotoptions.Imin)
       plotoptions.Imin = I-ModProfileErr.data[xi];
   }
 }
      }
      for(xi=0; xi < VarianceProfile.NrBins; xi++) {
 double xpos;
 xpos = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
 if(plotoptions.usephase)
   xpos /= 360.0;
 if(xpos + plotoptions.dl >= plotoptions.fl_min && xpos + plotoptions.dl <= plotoptions.fl_max) {
   I = VarianceProfile.data[xi];
   if(I > plotoptions.Imax)
     plotoptions.Imax = I;
   if(I < plotoptions.Imin)
     plotoptions.Imin = I;
 }
      }
    }else {
      for(xi=0; xi < subpulseAmpProfile.NrBins; xi++) {
 double xpos;
 xpos = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
 if(plotoptions.usephase)
   xpos /= 360.0;
 if(xpos + plotoptions.dl >= plotoptions.fl_min && xpos + plotoptions.dl <= plotoptions.fl_max) {
   I = subpulseAmpProfile.data[xi];
   if(I > plotoptions.Imax)
     plotoptions.Imax = I;
   if(I < plotoptions.Imin)
     plotoptions.Imin = I;
 }
      }
    }
    if(-0.05*plotoptions.Imax < plotoptions.Imin)
      plotoptions.Imin = -0.05*plotoptions.Imax;
    if(ImaxSet)
      plotoptions.Imax = ImaxValue/1.05;
    if(IminSet)
      plotoptions.Imin = IminValue/1.05;
    ppgswin(plotoptions.fl_min, plotoptions.fl_max, plotoptions.Imin, 1.05*plotoptions.Imax);
    ppgsch(0.55*plotoptions.labelscale*plotoptions.titlescale);
    ppgscf(2);
    ppgslw(2);
    if(title_index > 0) {
      ppgptxt(plotoptions.fl_min + 0.5*(plotoptions.fl_max-plotoptions.fl_min), 1.05*plotoptions.Imax + (0.10+plotoptions.titledy)*(1.05*plotoptions.Imax-plotoptions.Imin), 0.0, 0.5, argv[title_index]);
    }
    ppgslw(1);
    ppgscf(1);
    ppgsch(0.38*plotoptions.labelscale);
    if(LoadLRFS == 0) {
      if(plotoptions.inside) {
 if(plotoptions.noylabels) {
   ppgbox("bcnst",0.0,0,"bcst",0.0,0);
 }else {
   ppgbox("bcnst",0.0,0,"bcnst",0.0,0);
 }
      }else {
 if(plotoptions.noylabels) {
   ppgbox("bcnst",0.0,0,"bcsti",0.0,0);
 }else {
   ppgbox("bcnst",0.0,0,"bcnsti",0.0,0);
 }
      }
    }else {
      if(plotoptions.inside) {
 if(plotoptions.noylabels) {
   ppgbox("cst",0.0,0,"bcst",0.0,0);
 }else {
   ppgbox("cst",0.0,0,"bcnst",0.0,0);
 }
      }else {
 if(plotoptions.noylabels) {
   ppgbox("cst",0.0,0,"bcsti",0.0,0);
 }else {
   ppgbox("cst",0.0,0,"bcnsti",0.0,0);
 }
      }
    }
    if(plotoptions.plot_xlabel != 0 && LoadLRFS == 0) {
      ppgsch(0.3*plotoptions.labelscale);
      if(plotoptions.usephase)
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (phase)");
      else
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (deg)");
      ppgsch(0.38*plotoptions.labelscale);
    }
    if(application.onpulse.nrRegions > 0) {
      region_frac_to_int(&(application.onpulse), AverageProfile.NrBins, 0);
      int region_nr;
      float *xarray;
      float *yarray;
      xarray = malloc((AverageProfile.NrBins+2)*sizeof(float));
      yarray = malloc((AverageProfile.NrBins+2)*sizeof(float));
      if(xarray == NULL || yarray == NULL) {
 printerror(application.verbose_state.debug, "ERROR pspecFig: Memory allocation error");
 return 0;
      }
      for(region_nr = 0; region_nr < application.onpulse.nrRegions; region_nr++) {
 ppgsls(1);
 ppgslw(1);
 double deg_per_sample;
 deg_per_sample = get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
 int cur_index;
 cur_index = 0;
 for(xi=0; xi < AverageProfile.NrBins; xi++) {
   double xpos;
   xpos = xi*deg_per_sample;
   if(plotoptions.usephase)
     xpos /= 360.0;
   if(xpos+plotoptions.dl >= plotoptions.fl_min-deg_per_sample && xpos+plotoptions.dl <= plotoptions.fl_max+deg_per_sample) {
     if(checkRegions(xi, &(application.onpulse), 1+region_nr, application.verbose_state) != 0) {
       if(cur_index == 0) {
  xarray[cur_index] = xpos;
  yarray[cur_index] = plotoptions.Imin-plotoptions.Imax;
  cur_index++;
       }
       xarray[cur_index] = xpos;
       yarray[cur_index] = AverageProfile.data[xi]*profilescale;
       cur_index++;
     }
   }
 }
 xarray[cur_index] = xarray[cur_index-1];
 yarray[cur_index] = yarray[0];
 cur_index++;
 ppgsfs(1);
 ppgscr(16, 0.85, 0.85, 0.85);
 ppgsci(16);
 ppgpoly(cur_index, xarray, yarray);
 ppgsci(15);
 ppgshs(45+90*region_nr, 1.0, 0.0);
 ppgsfs(3);
 ppgpoly(cur_index, xarray, yarray);
 ppgsfs(2);
 ppgpoly(cur_index, xarray, yarray);
      }
      free(xarray);
      free(yarray);
      ppgsci(1);
    }
    ppgsls(lineStyle);
    ppgslw(2);
    ppgsci(plotoptions.lineColour);
    i = 0;
    double deg_per_sample;
    deg_per_sample = get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
    for(xi=0; xi < AverageProfile.NrBins; xi++) {
      double xpos;
      xpos = xi*deg_per_sample;
      if(plotoptions.usephase)
 xpos /= 360.0;
      if(xpos+plotoptions.dl >= plotoptions.fl_min-deg_per_sample && xpos+plotoptions.dl <= plotoptions.fl_max+deg_per_sample) {
 I = AverageProfile.data[xi]*profilescale;
 if(i == 0) {
   ppgmove(xpos+plotoptions.dl,I);
   i = 1;
 }else {
   ppgdraw(xpos+plotoptions.dl,I);
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
   if(plotoptions.usephase)
     x /= 360.0;
   x += plotoptions.dl;
   if(x >= plotoptions.fl_min-deg_per_sample && x <= plotoptions.fl_max+deg_per_sample) {
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
   if(plotoptions.usephase)
     x /= 360.0;
   x += plotoptions.dl;
   if(x >= plotoptions.fl_min-deg_per_sample && x <= plotoptions.fl_max+deg_per_sample) {
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
   if(plotoptions.usephase)
     x /= 360.0;
   x += plotoptions.dl;
   if(x >= plotoptions.fl_min-deg_per_sample && x <= plotoptions.fl_max+deg_per_sample) {
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
 if(plotoptions.usephase)
   x /= 360.0;
 x += plotoptions.dl;
 if(x >= plotoptions.fl_min-deg_per_sample && x <= plotoptions.fl_max+deg_per_sample) {
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
    if(plotoptions.plot_ylabel != 0) {
      ppgsch(0.3*plotoptions.labelscale);
      if(plotoptions.plot_ylabeltop) {
 if(type_of_plots == 0) {
   ppgmtxt("l",2.8,0.5,0.5,"Intensity/Modulation index");
 }else {
   ppgmtxt("l",2.8,0.5,0.5,"Intensity");
 }
      }
      ppgsch(0.38*plotoptions.labelscale);
    }
  }
  if(type_of_plots == 0) {
    if(Load2dfs != 0 && LoadLRFS == 0) {
      char *title;
      title = NULL;
      if(title_index > 0 && plotoptions.notop && LoadLRFS == 0) {
 title = argv[title_index];
      }
      Plot2dfs(twodfs_allinfo, twodfs2_allinfo, AverageProfile, 1, 0, title, plotoptions, argc, argv, application.verbose_state);
      if (LoadTwo != 0) {
 Plot2dfs(twodfs_allinfo, twodfs2_allinfo, AverageProfile, 1, 1, NULL, plotoptions, argc, argv, application.verbose_state);
      }
    }else {
      if(LoadLRFS != 0) {
 char *title;
 title = NULL;
 if(title_index > 0 && plotoptions.notop) {
   title = argv[title_index];
 }
 PlotLRFS(lrfs, AverageProfile, twodfs_allinfo, title, plotoptions, application.verbose_state);
      }
      if(Load2dfs != 0) {
 Plot2dfs(twodfs_allinfo, twodfs2_allinfo, AverageProfile, 0, 0, NULL, plotoptions, argc, argv, application.verbose_state);
      }
      if(LoadTwo != 0) {
 Plot2dfs(twodfs_allinfo, twodfs2_allinfo, AverageProfile, 0, 1, NULL, plotoptions, argc, argv, application.verbose_state);
      }
    }
  }else {
    ppgsvp(0.2, 0.2+0.11*plotoptions.scaleFig_x, 0.95-0.3*plotoptions.scaleFig_y, 0.95-0.15*plotoptions.scaleFig_y);
    ppgslw(1);
    ppgscf(1);
    ppgsch(0.38*plotoptions.labelscale);
    plotoptions.Imin = 0;
    plotoptions.Imax = 0;
    for(xi=0; xi < subpulseTrackProfile.NrBins; xi++) {
      double xpos;
      xpos = xi*get_tsamp(AverageProfile, 0, application.verbose_state)*360.0/period;
      if(plotoptions.usephase)
 xpos /= 360.0;
      if(xpos+plotoptions.dl >= plotoptions.fl_min && xpos+plotoptions.dl <= plotoptions.fl_max) {
 I = subpulseTrackProfile.data[xi];
 if(plotoptions.doFlip)
   I *= -1;
 if(I+subpulseTrackProfileErr.data[xi] > plotoptions.Imax)
   plotoptions.Imax = I+subpulseTrackProfileErr.data[xi];
 if(I-subpulseTrackProfileErr.data[xi] < plotoptions.Imin)
   plotoptions.Imin = I-subpulseTrackProfileErr.data[xi];
      }
    }
    if(plotoptions.SelectLRegion != 0) {
      plotoptions.fl_min = plotoptions.LRegionLow;
      plotoptions.fl_max = plotoptions.LRegionHigh;
    }
    if(maxSubpulsePhaseSet)
      plotoptions.Imax = maxSubpulsePhase/1.05;
    if(minSubpulsePhaseSet)
      plotoptions.Imin = minSubpulsePhase/1.05;
    ppgswin(plotoptions.fl_min, plotoptions.fl_max, plotoptions.Imin, 1.05*plotoptions.Imax);
    if(plotoptions.inside) {
      if(plotoptions.noylabels) {
 ppgbox("bcnst",0.0,0,"bcst",0.0,0);
      }else {
 ppgbox("bcnst",0.0,0,"bcnst",0.0,0);
      }
    }else {
      if(plotoptions.noylabels) {
 ppgbox("bcnst",0.0,0,"bcsti",0.0,0);
      }else {
 ppgbox("bcnst",0.0,0,"bcnsti",0.0,0);
      }
    }
    if(plotoptions.plot_xlabel != 0 && LoadLRFS == 0) {
      ppgsch(0.3*plotoptions.labelscale);
      if(plotoptions.usephase)
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (phase)");
      else
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (deg)");
      ppgsch(0.38*plotoptions.labelscale);
    }
    if(do_phase_slope) {
      ppgsls(4);
      if(plotoptions.usephase == 0) {
 x = plotoptions.dl;
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
 x = plotoptions.dl;
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
      if(plotoptions.usephase)
 x /= 360.0;
      x += plotoptions.dl;
      if(x >= plotoptions.fl_min && x <= plotoptions.fl_max) {
 I = subpulseTrackProfile.data[xi];
 if(plotoptions.doFlip)
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
      if(plotoptions.usephase)
 x /= 360.0;
      x += plotoptions.dl;
      if(x >= plotoptions.fl_min && x <= plotoptions.fl_max) {
 I = subpulseTrackProfile.data[xi];
 if(plotoptions.doFlip)
   I *= -1;
 ppgpt1(x, I, -1);
 ppgpt1(x, I+360, -1);
 ppgpt1(x, I-360, -1);
      }
    }
    ppgslw(1);
    if(plotoptions.plot_xlabel != 0) {
      ppgsch(0.3*plotoptions.labelscale);
      if(plotoptions.usephase)
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (phase)");
      else
 ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (deg)");
      ppgsch(0.38*plotoptions.labelscale);
    }
    if(plotoptions.plot_ylabel != 0) {
      ppgsch(0.3*plotoptions.labelscale);
      if(plotoptions.usephase)
 ppgmtxt("l",2.8,0.5,0.5,"Subpulse phase (phase)");
      else
 ppgmtxt("l",2.8,0.5,0.5,"Subpulse phase (deg)");
      ppgsch(0.38*plotoptions.labelscale);
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
  closePSRData(&AverageProfile, 0, 0, application.verbose_state);
  closePSRData(&(twodfs_allinfo.datafile), 0, 0, application.verbose_state);
  closePSRData(&(twodfs2_allinfo.datafile), 0, 0, application.verbose_state);
  closePSRData(&lrfs, 0, 0, application.verbose_state);
  closePSRData(&VarianceProfile, 0, 0, application.verbose_state);
  closePSRData(&ModProfile, 0, 0, application.verbose_state);
  closePSRData(&VarianceProfileErr, 0, 0, application.verbose_state);
  closePSRData(&ModProfileErr, 0, 0, application.verbose_state);
  closePSRData(&subpulseTrackProfile, 0, 0, application.verbose_state);
  closePSRData(&subpulseTrackProfileErr, 0, 0, application.verbose_state);
  closePSRData(&subpulseAmpProfile, 0, 0, application.verbose_state);
  terminateApplication(&application);
  return 0;
}
void GetExtremesSubsetVertical(twodfs_def twodfs_allinfo, float *Imin, float *Imax)
{
  float I, x, y;
  int xi, yi;
  *Imin = 0;
  *Imax = 0;
  for(xi = 0; xi < twodfs_allinfo.datafile.NrBins; xi++) {
    I = 0;
    for(yi = 0; yi < twodfs_allinfo.datafile.NrSubints; yi++) {
      pgplotMapCoordinateInverse(&x, &y, xi, yi);
      if(y >= twodfs_allinfo.P3IntegrateLow && y <= twodfs_allinfo.P3IntegrateHigh) {
 if(x >= twodfs_allinfo.f2_min && x <= twodfs_allinfo.f2_max) {
   I += twodfs_allinfo.datafile.data[yi*twodfs_allinfo.datafile.NrBins+xi];
 }
      }
    }
    if(I < *Imin)
      *Imin = I;
    if(I > *Imax)
      *Imax = I;
  }
  *Imin *= 2.0;
  *Imax *= 2.0;
}
void GetExtremesSubsetHorizontal(datafile_definition spectrum, int istwodfs, double f2_min, double f2_max, double f3_min, double f3_max, float *Imin, float *Imax, plotoptions_def plotoptions)
{
  float I, x, y;
  int xi, yi, ok;
  *Imin = 0;
  *Imax = 0;
  for(yi = plotoptions.ExtraVerticalMaximaSkip; yi < spectrum.NrSubints; yi++) {
    I = 0;
    for(xi = 0; xi < spectrum.NrBins; xi++) {
      pgplotMapCoordinateInverse(&x, &y, xi, yi);
      ok = 0;
      if(istwodfs == 0) {
 if(x >= plotoptions.fl_min && x <= plotoptions.fl_max && y >= f3_min && y <= f3_max) {
   ok = 1;
 }
      }else {
 if(x >= f2_min && x <= f2_max && y >= f3_min && y <= f3_max) {
   ok = 1;
 }
      }
      if(ok) {
 I += spectrum.data[yi*spectrum.NrBins+xi];
      }
    }
    if(I < *Imin)
      *Imin = I;
    if(I > *Imax)
      *Imax = I;
  }
  *Imin *= 2.0;
  *Imax *= 2.0;
}
void set_color_featurenumber(int featurenr, int passnr)
{
 if(featurenr == 0) {
   if(passnr == 0) {
     ppgscr(16+featurenr, 1.0, 0.9, 0.8);
   }else {
     ppgscr(16+featurenr, 1.0, 0.8, 0.6);
   }
   ppgsci(16+featurenr);
 }else if(featurenr == 1) {
   if(passnr == 0) {
     ppgscr(16+featurenr, 0.9, 0.9, 1.0);
   }else {
     ppgscr(16+featurenr, 0.8, 0.8, 1.0);
   }
   ppgsci(16+featurenr);
 }else if(featurenr == 2) {
   if(passnr == 0) {
     ppgscr(16+featurenr, 0.9, 1.0, 0.9);
   }else {
     ppgscr(16+featurenr, 0.8, 1.0, 0.7);
   }
   ppgsci(16+featurenr);
 }else if(featurenr == 3) {
   if(passnr == 0) {
     ppgscr(16+featurenr, 1.0, 1.0, 0.85);
   }else {
     ppgscr(16+featurenr, 1.0, 1.0, 0.5);
   }
   ppgsci(16+featurenr);
 }else if(featurenr == 4) {
   if(passnr == 0) {
     ppgscr(16+featurenr, 1.0, 0.85, 1.0);
   }else {
     ppgscr(16+featurenr, 1.0, 0.7, 1.0);
   }
   ppgsci(16+featurenr);
 }else {
   ppgsci(featurenr+2);
 }
}
void IntegrateSubsetHorizontal(datafile_definition lrfs, twodfs_def twodfs_allinfo, twodfs_def twodfs2_allinfo, int twodfsonly, int Number, int normalise, plotoptions_def plotoptions)
{
  float I, Imin, Imax, x, y;
  int xi, yi;
  float offset;
  offset = 0+0.03*(plotoptions.labelscale-1.0)*((float)Number+1.0);
  if(Number == -1) {
    ppgsvp(0.2-0.04*plotoptions.scaleFig_x, 0.2, 0.95-0.3*plotoptions.scaleFig_y, 0.95-0.15*plotoptions.scaleFig_y);
  }else {
    if(twodfsonly == 0) {
      ppgsvp(0.2-0.04*plotoptions.scaleFig_x, 0.2, 0.95-0.48*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y, 0.95-0.33*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y);
    }else {
      ppgsvp(0.2-0.04*plotoptions.scaleFig_x, 0.2, 0.95-0.33*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y, 0.95-0.18*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y);
    }
  }
  if(Number == -1) {
    GetExtremesSubsetHorizontal(lrfs, 0, 0.0, 0.0, twodfs_allinfo.f3_min, twodfs_allinfo.f3_max, &Imin, &Imax, plotoptions);
  }else if(Number == 0) {
    GetExtremesSubsetHorizontal(twodfs_allinfo.datafile, 1, twodfs_allinfo.f2_min, twodfs_allinfo.f2_max, twodfs_allinfo.f3_min, twodfs_allinfo.f3_max, &Imin, &Imax, plotoptions);
  }else {
    GetExtremesSubsetHorizontal(twodfs2_allinfo.datafile, 1, twodfs2_allinfo.f2_min, twodfs2_allinfo.f2_max, twodfs2_allinfo.f3_min, twodfs2_allinfo.f3_max, &Imin, &Imax, plotoptions);
  }
  double scale;
  scale = fabs(Imax);
  if(fabs(Imin) > scale)
    scale = fabs(Imin);
  if(normalise == 0)
    scale = 1.0;
  ppgswin(Imin/scale,1.05*Imax/scale, twodfs_allinfo.f3_min,twodfs_allinfo.f3_max);
  if(Number == 0 && twodfs_allinfo.SelectP3Integrate) {
    ppgsci(4);
    ppgsls(4);
    if(twodfs_allinfo.SelectP3Integrate == 2) {
      ppgslw(3);
    }
    ppgmove(Imin/scale, twodfs_allinfo.P3IntegrateHigh);
    ppgdraw(1.05*Imax/scale, twodfs_allinfo.P3IntegrateHigh);
    ppgmove(Imin/scale, twodfs_allinfo.P3IntegrateLow);
    ppgdraw(1.05*Imax/scale, twodfs_allinfo.P3IntegrateLow);
    ppgsci(1);
    ppgsls(1);
    ppgslw(1);
  }else if(Number == 1 && twodfs2_allinfo.SelectP3Integrate) {
    ppgsci(4);
    ppgsls(4);
    if(twodfs2_allinfo.SelectP3Integrate == 2) {
      ppgslw(3);
    }
    ppgmove(Imin/scale, twodfs2_allinfo.P3IntegrateHigh);
    ppgdraw(1.05*Imax/scale, twodfs2_allinfo.P3IntegrateHigh);
    ppgmove(Imin/scale, twodfs2_allinfo.P3IntegrateLow);
    ppgdraw(1.05*Imax/scale, twodfs2_allinfo.P3IntegrateLow);
    ppgsci(1);
    ppgsls(1);
    ppgslw(1);
  }
  ppgbox("bcv",0.0,0,"bnsti",0.0,0);
  ppgbox("",0.0,0,"c",0.0,0);
  if(plotoptions.textside != NULL) {
    if(twodfsonly == 0) {
      if(Number == -1) {
 ppgsch(0.55*plotoptions.labelscale*plotoptions.textsidescale);
 ppgscf(2);
 ppgslw(2);
 if (plotoptions.notop == 1) {
   ppgtext(Imin/scale, twodfs_allinfo.f3_min+(twodfs_allinfo.f3_max-twodfs_allinfo.f3_min)*1.15, plotoptions.textside);
 }else {
   ppgtext(Imin/scale, twodfs_allinfo.f3_min+(twodfs_allinfo.f3_max-twodfs_allinfo.f3_min)*2.095, plotoptions.textside);
 }
 ppgslw(1);
 ppgscf(1);
 ppgsch(0.38*plotoptions.labelscale);
      }
    }else {
      if(strlen(plotoptions.textside) != 0 && plotoptions.notop == 0) {
 ppgsch(0.55*plotoptions.labelscale*plotoptions.textsidescale);
 ppgscf(2);
 ppgslw(2);
 ppgtext(Imin/scale, twodfs_allinfo.f3_min+(twodfs_allinfo.f3_max-twodfs_allinfo.f3_min)+0.65, plotoptions.textside);
 ppgslw(1);
 ppgscf(1);
 ppgsch(0.38*plotoptions.labelscale);
      }else if (strlen(plotoptions.textside) != 0 && plotoptions.notop == 1) {
 ppgsch(0.55*plotoptions.labelscale*plotoptions.textsidescale);
 ppgscf(2);
 ppgslw(2);
 ppgtext(Imin/scale, twodfs_allinfo.f3_min+(twodfs_allinfo.f3_max-twodfs_allinfo.f3_min)+0.1, plotoptions.textside);
 ppgslw(1);
 ppgscf(1);
 ppgsch(0.38*plotoptions.labelscale);
      }
    }
  }
  y = floor(log10(Imax*0.33/scale));
  x = floor(Imax*0.33/(pow(10,y)*scale));
  x = x*pow(10,y);
  char labelnumbers[3];
  if(plotoptions.nointegrateNumbers == 0) {
    strcpy(labelnumbers, "n");
  }else {
    strcpy(labelnumbers, "");
  }
  ppgsch(0.38*plotoptions.labelscale*0.66);
  ppgaxis(labelnumbers,0,twodfs_allinfo.f3_max,Imax*0.33/scale,twodfs_allinfo.f3_max,0,Imax*0.33/scale,x,1,0.3,0,0,-0.5,90);
  ppgsch(0.38*plotoptions.labelscale);
  datafile_definition *spectrum;
  double xmin, xmax, ymin, ymax;
  if(Number == -1) {
    spectrum = &lrfs;
    xmin = plotoptions.fl_min;
    xmax = plotoptions.fl_max;
    ymin = twodfs_allinfo.f3_min;
    ymax = twodfs_allinfo.f3_max;
  }else if(Number == 0) {
    spectrum = &(twodfs_allinfo.datafile);
    xmin = twodfs_allinfo.f2_min;
    xmax = twodfs_allinfo.f2_max;
    ymin = twodfs_allinfo.f3_min;
    ymax = twodfs_allinfo.f3_max;
  }else if(Number == 1) {
    spectrum = &(twodfs2_allinfo.datafile);
    xmin = twodfs2_allinfo.f2_min;
    xmax = twodfs2_allinfo.f2_max;
    ymin = twodfs2_allinfo.f3_min;
    ymax = twodfs2_allinfo.f3_max;
  }
  for(yi = 0; yi < spectrum->NrSubints; yi++) {
    I = 0;
    for(xi = 0; xi < spectrum->NrBins; xi++) {
      pgplotMapCoordinateInverse(&x, &y, xi, yi);
      if(x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
 I += spectrum->data[yi*spectrum->NrBins+xi];
      }
    }
    if(yi == 0) {
      ppgmove(2.0*I/scale, y);
    }else {
      ppgdraw(2.0*I/scale, y);
    }
  }
  ppgslw(1);
}
void PlotLRFS(datafile_definition lrfs, datafile_definition AverageProfile, twodfs_def twodfs_allinfo, char *title, plotoptions_def plotoptions, verbose_definition verbose)
{
  float I;
  int xi, i;
  ppgsvp(0.2, 0.2+0.11*plotoptions.scaleFig_x, 0.95-0.3*plotoptions.scaleFig_y, 0.95-0.15*plotoptions.scaleFig_y);
  if(title != NULL) {
    ppgsch(0.55*plotoptions.labelscale*plotoptions.titlescale);
    ppgscf(2);
    ppgslw(2);
    ppgmtxt("t",1,0.5,0.5,title);
    ppgslw(1);
    ppgscf(1);
    ppgsch(0.38*plotoptions.labelscale);
  }
  if(plotoptions.SelectLRegion != 0) {
    plotoptions.fl_min = plotoptions.LRegionLow;
    plotoptions.fl_max = plotoptions.LRegionHigh;
  }
  ppgswin(plotoptions.fl_min,plotoptions.fl_max,twodfs_allinfo.f3_min,twodfs_allinfo.f3_max);
  if(plotoptions.plot_xlabel != 0) {
    ppgsch(0.3*plotoptions.labelscale);
    if(plotoptions.usephase)
      ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (phase)");
    else
      ppgmtxt("b",3.0,0.5,0.5,"Pulse longitude (deg)");
    ppgsch(0.38*plotoptions.labelscale);
  }
  if(plotoptions.plot_ylabel != 0) {
    ppgsch(0.3*plotoptions.labelscale);
    if(plotoptions.noside) {
      if(plotoptions.plot_ylabel == 1)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plotoptions.plot_ylabel == 2)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plotoptions.plot_ylabel == 3)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }else {
      float offset;
      offset = 10 - (plotoptions.labelscale-1.0)*(7-2.8);
      offset += 7*(plotoptions.scaleFig_x-1.0);
      if(plotoptions.plot_ylabel == 1)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plotoptions.plot_ylabel == 2)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plotoptions.plot_ylabel == 3)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }
    ppgsch(0.38*plotoptions.labelscale);
  }
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);
  pgplot_options.box.box_labelsize = 0.3*plotoptions.labelscale;
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
  if(plotoptions.usephase)
    xright /= 360.0;
  pgplotMap(&pgplot_options, lrfs.data, lrfs.NrBins, lrfs.NrSubints, 0+plotoptions.dl, xright+plotoptions.dl, plotoptions.fl_min, plotoptions.fl_max, 0, 0.5, twodfs_allinfo.f3_min, twodfs_allinfo.f3_max, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, plotoptions.oversaturizel, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, plotoptions.showwedge, plotoptions.showwedge_max, 0, 0, verbose);
  ret = get_period(AverageProfile, 0, &period, verbose);
  if(ret == 2) {
    printerror(verbose.debug, "ERROR pspecFig (%s): Cannot obtain period", AverageProfile.filename);
    exit(0);
  }
  if(plotoptions.overlaypp) {
    ppgsci(plotoptions.lineColour);
    ppgsls(1);
    ppgslw(1);
    i = 0;
    for(xi=0; xi < AverageProfile.NrBins; xi++) {
      double xpos;
      xpos = xi*get_tsamp(AverageProfile, 0, verbose)*360.0/period;
      if(plotoptions.usephase)
 xpos /= 360.0;
      if(xpos+plotoptions.dl >= plotoptions.fl_min && xpos+plotoptions.dl <= plotoptions.fl_max) {
 I = AverageProfile.data[xi]*0.5;
 if(i == 0) {
   ppgmove(xpos+plotoptions.dl,I);
   i = 1;
 }else {
   ppgdraw(xpos+plotoptions.dl,I);
 }
      }
    }
    ppgsci(1);
  }
  ppgswin(plotoptions.fl_min,plotoptions.fl_max,twodfs_allinfo.f3_min,twodfs_allinfo.f3_max);
  ppgsch(0.38*plotoptions.labelscale);
  if(plotoptions.noside) {
    if(plotoptions.inside) {
      if(plotoptions.noylabels)
 ppgbox("bcnst",0.0,0,"bcst",0.0,0);
      else
 ppgbox("bcnst",0.0,0,"bcnst",0.0,0);
    }else {
      if(plotoptions.noylabels)
 ppgbox("bcnsti",0.0,0,"bcsti",0.0,0);
      else
 ppgbox("bcnsti",0.0,0,"bcnsti",0.0,0);
    }
  }else {
    if(plotoptions.inside) {
      ppgbox("bcnst",0.0,0,"cst",0.0,0);
    }else {
      ppgbox("bcnsti",0.0,0,"csti",0.0,0);
    }
  }
  if(plotoptions.noside == 0) {
    IntegrateSubsetHorizontal(lrfs, twodfs_allinfo, twodfs_allinfo, 0, -1, plotoptions.normaliseSide, plotoptions);
  }
  ppgsch(0.38*plotoptions.labelscale);
}
void IntegrateSubsetVertical(twodfs_def twodfs_allinfo, int twodfsonly, int Number, int normalise, int nointegrateNumbers, plotoptions_def plotoptions)
{
  float I, Imin, Imax, x, y, offset;
  int xi, yi;
  offset = 0+0.03*(plotoptions.labelscale-1.0)*((float)Number+1.0);
  if(twodfsonly == 0) {
    ppgsvp(0.2, 0.2+0.11*plotoptions.scaleFig_x, 0.95-0.52*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y, 0.95-0.48*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y);
  }else {
    ppgsvp(0.2, 0.2+0.11*plotoptions.scaleFig_x, 0.95-0.37*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y, 0.95-0.33*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y);
  }
  GetExtremesSubsetVertical(twodfs_allinfo, &Imin, &Imax);
  double scale;
  scale = fabs(Imax);
  if(fabs(Imin) > scale)
    scale = fabs(Imin);
  if(normalise == 0)
    scale = 1.0;
  ppgswin(twodfs_allinfo.f2_min,twodfs_allinfo.f2_max,Imin/scale,1.05*Imax/scale);
  if(plotoptions.noxlabels == 0)
    ppgbox("bnst",0.0,0,"bcvi",0.0,0);
  else
    ppgbox("bst",0.0,0,"bcvi",0.0,0);
  ppgbox("c",0.0,0,"",0.0,0);
  y = floor(log10(Imax*0.7/scale));
  x = floor(Imax*0.7/(pow(10,y)*scale));
  x = x*pow(10,y);
  char labelnumbers[3];
  if(plotoptions.nointegrateNumbers == 0) {
    strcpy(labelnumbers, "n");
  }else {
    strcpy(labelnumbers, "");
  }
  ppgsch(0.38*plotoptions.labelscale*0.66);
  ppgaxis(labelnumbers,twodfs_allinfo.f2_min, 0,twodfs_allinfo.f2_min, Imax*.7/scale,0,Imax*.7/scale,x,1,0.3,0,0,-0.8,90);
  ppgsch(0.38*plotoptions.labelscale);
  int direction;
  for(direction = 1-plotoptions.intflip; direction < 2; direction++) {
    if(direction == 0) {
      ppgsci(8);
      ppgsls(2);
    }else {
      ppgsci(1);
      ppgsls(1);
    }
    for(xi = 0; xi < twodfs_allinfo.datafile.NrBins; xi++) {
      I = 0;
      for(yi = 0; yi < twodfs_allinfo.datafile.NrSubints; yi++) {
 pgplotMapCoordinateInverse(&x, &y, xi, yi);
 if(y >= twodfs_allinfo.P3IntegrateLow && y <= twodfs_allinfo.P3IntegrateHigh)
   I += twodfs_allinfo.datafile.data[yi*twodfs_allinfo.datafile.NrBins+xi];
      }
      if(direction == 0)
 x *= -1.0;
      if(xi == 0) {
 ppgmove(x,2.0*I/scale);
      }else {
 ppgdraw(x,2.0*I/scale);
      }
    }
  }
  ppgslw(1);
}
void ppgplot_set_internal_mapping_coordinated(double xmin, double xmax, double ymin, double ymax, int nrx, int nry);
void Plot2dfs(twodfs_def twodfs_allinfo, twodfs_def twodfs2_allinfo, datafile_definition AverageProfile, int twodfsonly, int Number, char *title, plotoptions_def plotoptions, int argc, char **argv, verbose_definition verbose)
{
  datafile_definition *curtwodfs;
  if(Number == 0) {
    curtwodfs = &(twodfs_allinfo.datafile);
  }else {
    curtwodfs = &(twodfs2_allinfo.datafile);
  }
  float cpp_min;
  float cpp_max;
  if(curtwodfs->xrangeset) {
    cpp_min = curtwodfs->xrange[0];
    cpp_max = curtwodfs->xrange[1];
  }else {
    cpp_min = -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)curtwodfs->NrBins;
    cpp_max = +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)curtwodfs->NrBins;
  }
  if(plotoptions.doFlip) {
    cpp_min *= -1.0;
    cpp_max *= -1.0;
  }
  ppgplot_set_internal_mapping_coordinated(cpp_min, cpp_max, 0.0, 0.5, curtwodfs->NrBins, curtwodfs->NrSubints);
  if(Number == 0) {
    IntegrateSubsetVertical(twodfs_allinfo, twodfsonly, Number, plotoptions.normaliseSide, plotoptions.nointegrateNumbers, plotoptions);
  }else {
    IntegrateSubsetVertical(twodfs2_allinfo, twodfsonly, Number, plotoptions.normaliseSide, plotoptions.nointegrateNumbers, plotoptions);
  }
  if(plotoptions.plot_xlabel == 1) {
    ppgsch(0.3*plotoptions.labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (cpp)");
    ppgsch(0.38*plotoptions.labelscale);
  }else if(plotoptions.plot_xlabel == 2) {
    ppgsch(0.3*plotoptions.labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d2\\u)");
    ppgsch(0.38*plotoptions.labelscale);
  }else if(plotoptions.plot_xlabel == 3) {
    ppgsch(0.3*plotoptions.labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (P/P\\d2\\u)");
    ppgsch(0.38*plotoptions.labelscale);
  }
  if(!(plotoptions.noside)) {
    IntegrateSubsetHorizontal(twodfs_allinfo.datafile, twodfs_allinfo, twodfs2_allinfo, twodfsonly, Number, plotoptions.normaliseSide, plotoptions);
  }
  float offset;
  offset = 0+0.03*(plotoptions.labelscale-1.0)*((float)Number+1.0);
  if(twodfsonly == 0) {
    ppgsvp(0.2, 0.2+0.11*plotoptions.scaleFig_x, 0.95-0.48*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y, 0.95-0.33*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y);
  }else {
    ppgsvp(0.2, 0.2+0.11*plotoptions.scaleFig_x, 0.95-0.33*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y, 0.95-0.18*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y);
  }
  if(twodfsonly && title != NULL) {
    ppgsch(0.55*plotoptions.labelscale*plotoptions.titlescale);
    ppgscf(2);
    ppgslw(2);
    ppgmtxt("t",1,0.5,0.5,title);
    ppgslw(1);
    ppgscf(1);
    ppgsch(0.38*plotoptions.labelscale);
  }
  ppgswin(twodfs_allinfo.f2_min,twodfs_allinfo.f2_max,twodfs_allinfo.f3_min,twodfs_allinfo.f3_max);
  if(plotoptions.plot_ylabel != 0) {
    ppgsch(0.3*plotoptions.labelscale);
    if(plotoptions.noside) {
      if(plotoptions.plot_ylabel == 1)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plotoptions.plot_ylabel == 2)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plotoptions.plot_ylabel == 3)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }else {
      float offset;
      offset = 10 - (plotoptions.labelscale-1.0)*(7-2.8);
      offset += 7*(plotoptions.scaleFig_x-1.0);
      if(plotoptions.plot_ylabel == 1)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plotoptions.plot_ylabel == 2)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plotoptions.plot_ylabel == 3)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }
    ppgsch(0.38*plotoptions.labelscale);
  }
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);
  pgplot_options.box.box_labelsize = 0.3*plotoptions.labelscale;
  pgplot_options.viewport.noclear = 1;
  pgplot_options.viewport.dontopen = 1;
  pgplot_options.viewport.dontclose = 1;
  ppgsci(3);
  pgplot_options.box.drawbox = 1;
  pgplot_options.box.drawlabels = 1;
  sprintf(pgplot_options.box.xlabel, "Hello");
  if(!plotoptions.nomain) {
    float cur_f2_min, cur_f2_max;
    float cur_f3_min, cur_f3_max;
    float cur_oversaturize;
    if(Number == 0) {
      cur_f2_min = twodfs_allinfo.f2_min;
      cur_f2_max = twodfs_allinfo.f2_max;
      cur_f3_min = twodfs_allinfo.f3_min;
      cur_f3_max = twodfs_allinfo.f3_max;
      cur_oversaturize = twodfs_allinfo.oversaturize;
    }else {
      cur_f2_min = twodfs2_allinfo.f2_min;
      cur_f2_max = twodfs2_allinfo.f2_max;
      cur_f3_min = twodfs2_allinfo.f3_min;
      cur_f3_max = twodfs2_allinfo.f3_max;
      cur_oversaturize = twodfs2_allinfo.oversaturize;
    }
    pgplotMap(&pgplot_options, curtwodfs->data, curtwodfs->NrBins, curtwodfs->NrSubints, cpp_min, cpp_max, cur_f2_min, cur_f2_max, 0, 0.5, cur_f3_min, cur_f3_max, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, cur_oversaturize, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, plotoptions.showwedge, plotoptions.showwedge_max, 0, 0, verbose);
    ppgsch(0.38*plotoptions.labelscale);
    if(plotoptions.noside) {
      if(plotoptions.inside) {
 if(plotoptions.noylabels)
   ppgbox("cst",0.0,0,"bcst",0.0,0);
 else
   ppgbox("cst",0.0,0,"bcnst",0.0,0);
      }else {
 if(plotoptions.noylabels)
   ppgbox("csti",0.0,0,"bcsti",0.0,0);
 else
   ppgbox("csti",0.0,0,"bcnsti",0.0,0);
      }
    }else {
      if(plotoptions.inside) {
 ppgbox("cst",0.0,0,"cst",0.0,0);
      }else {
 ppgbox("csti",0.0,0,"csti",0.0,0);
      }
    }
    ppgbox("c",0.0,0,"c",0.0,0);
    ppgbox("b",-1.0,0,"b",-1.0,0);
    if(twodfs_allinfo.SelectP3Integrate != 0) {
      int loworhigh;
      ppgsls(1);
      ppgslw(1);
      for(loworhigh = 0; loworhigh < 2; loworhigh++) {
 int pointnr;
 int nrpoints;
 nrpoints = 10;
 for(pointnr = -1; pointnr < nrpoints; pointnr++) {
   float x, y, dx;
   if(Number == 0) {
     dx = (twodfs_allinfo.f2_max - twodfs_allinfo.f2_min)/(float)(nrpoints);
     x = twodfs_allinfo.f2_min + (pointnr + 0.5)*dx;
     if(loworhigh == 0) {
       y = twodfs_allinfo.P3IntegrateLow;
     }else {
       y = twodfs_allinfo.P3IntegrateHigh;
     }
   }else {
     dx = (twodfs2_allinfo.f2_max - twodfs2_allinfo.f2_min)/(float)(nrpoints);
     x = twodfs2_allinfo.f2_min + (pointnr + 0.5)*dx;
     if(loworhigh == 0) {
       y = twodfs2_allinfo.P3IntegrateLow;
     }else {
       y = twodfs2_allinfo.P3IntegrateHigh;
     }
   }
   ppgsci(7);
   ppgmove(x-0.05*dx, y);
   ppgdraw(x+0.05*dx, y);
   ppgsci(4);
   ppgmove(x-0.05*dx+0.5*dx, y);
   ppgdraw(x+0.05*dx+0.5*dx, y);
 }
      }
      ppgsci(1);
      ppgsls(1);
      ppgslw(1);
    }
  }
  int i, featurenr;
  featurenr = 0;
  for(i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-markP2P3") == 0) {
      double markp2, markp2errpos, markp2errneg, markp3, markp3err;
      int marktwodfsnr;
      if(parse_command_string(verbose, argc, argv, i+1, 0, -1, "%d %lf %lf %lf %lf %lf", &marktwodfsnr, &markp2, &markp2errpos, &markp2errneg, &markp3, &markp3err, NULL) == 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR pspecFig: Cannot parse '%s' option.", argv[i]);
 exit(0);
      }
      if(marktwodfsnr == Number + 1) {
 set_color_featurenumber(featurenr, 1);
 ppgslw(13);
 ppgerr1(1, markp2, markp3, markp2errpos, 1.0);
 ppgerr1(3, markp2, markp3, markp2errneg, 1.0);
 ppgerr1(6, markp2, markp3, markp3err, 1.0);
 ppgslw(1);
 ppgsci(1);
 ppgerr1(1, markp2, markp3, markp2errpos, 1.0);
 ppgerr1(3, markp2, markp3, markp2errneg, 1.0);
 ppgerr1(6, markp2, markp3, markp3err, 1.0);
 featurenr++;
      }
      i++;
    }
  }
}
void Plot2dfs_old(twodfs_def twodfs_allinfo, twodfs_def twodfs2_allinfo, datafile_definition AverageProfile, int twodfsonly, int Number, char *title, plotoptions_def plotoptions, verbose_definition verbose)
{
  float offset;
  offset = 0+0.03*(plotoptions.labelscale-1.0)*((float)Number+1.0);
  if(twodfsonly == 0) {
    ppgsvp(0.2, 0.2+0.11*plotoptions.scaleFig_x, 0.95-0.48*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y, 0.95-0.33*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y);
  }else {
    ppgsvp(0.2, 0.2+0.11*plotoptions.scaleFig_x, 0.95-0.33*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y, 0.95-0.18*plotoptions.scaleFig_y-0.22*Number*plotoptions.scaleFig_y-offset*plotoptions.scaleFig_y);
  }
  if(twodfsonly && title != NULL) {
    ppgsch(0.55*plotoptions.labelscale*plotoptions.titlescale);
    ppgscf(2);
    ppgslw(2);
    ppgmtxt("t",1,0.5,0.5,title);
    ppgslw(1);
    ppgscf(1);
    ppgsch(0.38*plotoptions.labelscale);
  }
  ppgswin(twodfs_allinfo.f2_min,twodfs_allinfo.f2_max,twodfs_allinfo.f3_min,twodfs_allinfo.f3_max);
  if(plotoptions.plot_ylabel != 0) {
    ppgsch(0.3*plotoptions.labelscale);
    if(plotoptions.noside) {
      if(plotoptions.plot_ylabel == 1)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plotoptions.plot_ylabel == 2)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plotoptions.plot_ylabel == 3)
 ppgmtxt("l",2.8,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }else {
      float offset;
      offset = 10 - (plotoptions.labelscale-1.0)*(7-2.8);
      offset += 7*(plotoptions.scaleFig_x-1.0);
      if(plotoptions.plot_ylabel == 1)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (cpp)");
      else if(plotoptions.plot_ylabel == 2)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d3\\u)");
      else if(plotoptions.plot_ylabel == 3)
 ppgmtxt("l",offset,0.5,0.5,"Fluctuation frequency (P/P\\d3\\u)");
    }
    ppgsch(0.38*plotoptions.labelscale);
  }
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);
  pgplot_options.box.box_labelsize = 0.3*plotoptions.labelscale;
  pgplot_options.viewport.noclear = 1;
  pgplot_options.viewport.dontopen = 1;
  pgplot_options.viewport.dontclose = 1;
  if(!plotoptions.nomain) {
    datafile_definition *curtwodfs;
    float cur_f2_min, cur_f2_max;
    float cur_f3_min, cur_f3_max;
    float cur_oversaturize;
    if(Number == 0) {
      curtwodfs = &(twodfs_allinfo.datafile);
      cur_f2_min = twodfs_allinfo.f2_min;
      cur_f2_max = twodfs_allinfo.f2_max;
      cur_f3_min = twodfs_allinfo.f3_min;
      cur_f3_max = twodfs_allinfo.f3_max;
      cur_oversaturize = twodfs_allinfo.oversaturize;
    }else {
      curtwodfs = &(twodfs2_allinfo.datafile);
      cur_f2_min = twodfs2_allinfo.f2_min;
      cur_f2_max = twodfs2_allinfo.f2_max;
      cur_f3_min = twodfs2_allinfo.f3_min;
      cur_f3_max = twodfs2_allinfo.f3_max;
      cur_oversaturize = twodfs2_allinfo.oversaturize;
    }
    float cpp_min;
    float cpp_max;
    if(curtwodfs->xrangeset) {
      cpp_min = curtwodfs->xrange[0];
      cpp_max = curtwodfs->xrange[1];
    }else {
      cpp_min = -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)curtwodfs->NrBins;
      cpp_max = +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)curtwodfs->NrBins;
    }
    if(plotoptions.doFlip) {
      cpp_min *= -1.0;
      cpp_max *= -1.0;
    }
    pgplotMap(&pgplot_options, curtwodfs->data, curtwodfs->NrBins, curtwodfs->NrSubints, cpp_min, cpp_max, cur_f2_min, cur_f2_max, 0, 0.5, cur_f3_min, cur_f3_max, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, cur_oversaturize, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, plotoptions.showwedge, plotoptions.showwedge_max, 0, 0, verbose);
    ppgsch(0.38*plotoptions.labelscale);
    if(plotoptions.noside) {
      if(plotoptions.inside) {
 if(plotoptions.noylabels)
   ppgbox("cst",0.0,0,"bcst",0.0,0);
 else
   ppgbox("cst",0.0,0,"bcnst",0.0,0);
      }else {
 if(plotoptions.noylabels)
   ppgbox("csti",0.0,0,"bcsti",0.0,0);
 else
   ppgbox("csti",0.0,0,"bcnsti",0.0,0);
      }
    }else {
      if(plotoptions.inside) {
 ppgbox("cst",0.0,0,"cst",0.0,0);
      }else {
 ppgbox("csti",0.0,0,"csti",0.0,0);
      }
    }
    ppgbox("c",0.0,0,"c",0.0,0);
    if(twodfs_allinfo.SelectP3Integrate != 0) {
      int loworhigh;
      ppgsls(1);
      ppgslw(1);
      for(loworhigh = 0; loworhigh < 2; loworhigh++) {
 int pointnr;
 int nrpoints;
 nrpoints = 10;
 for(pointnr = -1; pointnr < nrpoints; pointnr++) {
   float x, y, dx;
   if(Number == 0) {
     dx = (twodfs_allinfo.f2_max - twodfs_allinfo.f2_min)/(float)(nrpoints);
     x = twodfs_allinfo.f2_min + (pointnr + 0.5)*dx;
     if(loworhigh == 0) {
       y = twodfs_allinfo.P3IntegrateLow;
     }else {
       y = twodfs_allinfo.P3IntegrateHigh;
     }
   }else {
     dx = (twodfs2_allinfo.f2_max - twodfs2_allinfo.f2_min)/(float)(nrpoints);
     x = twodfs2_allinfo.f2_min + (pointnr + 0.5)*dx;
     if(loworhigh == 0) {
       y = twodfs2_allinfo.P3IntegrateLow;
     }else {
       y = twodfs2_allinfo.P3IntegrateHigh;
     }
   }
   ppgsci(7);
   ppgmove(x-0.05*dx, y);
   ppgdraw(x+0.05*dx, y);
   ppgsci(4);
   ppgmove(x-0.05*dx+0.5*dx, y);
   ppgdraw(x+0.05*dx+0.5*dx, y);
 }
      }
      ppgsci(1);
      ppgsls(1);
      ppgslw(1);
    }
  }
  if(Number == 0) {
    IntegrateSubsetVertical(twodfs_allinfo, twodfsonly, Number, plotoptions.normaliseSide, plotoptions.nointegrateNumbers, plotoptions);
  }else {
    IntegrateSubsetVertical(twodfs2_allinfo, twodfsonly, Number, plotoptions.normaliseSide, plotoptions.nointegrateNumbers, plotoptions);
  }
  if(plotoptions.plot_xlabel == 1) {
    ppgsch(0.3*plotoptions.labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (cpp)");
    ppgsch(0.38*plotoptions.labelscale);
  }else if(plotoptions.plot_xlabel == 2) {
    ppgsch(0.3*plotoptions.labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (P\\d0\\u/P\\d2\\u)");
    ppgsch(0.38*plotoptions.labelscale);
  }else if(plotoptions.plot_xlabel == 3) {
    ppgsch(0.3*plotoptions.labelscale);
    ppgmtxt("b",3.0,0.5,0.5,"Fluctuation frequency (P/P\\d2\\u)");
    ppgsch(0.38*plotoptions.labelscale);
  }
  if(!(plotoptions.noside)) {
    IntegrateSubsetHorizontal(twodfs_allinfo.datafile, twodfs_allinfo, twodfs2_allinfo, twodfsonly, Number, plotoptions.normaliseSide, plotoptions);
  }
}
int loadLRFS(datafile_definition *lrfs, int extprefix, int longsnap, int argc, char **argv, plotoptions_def *plotoptions, verbose_definition verbose)
{
  char filename[MaxFilenameLength], txt[MaxFilenameLength];
  datafile_definition clone;
  if(extprefix == 0) {
    sprintf(txt, "lrfs");
  }
  if(change_filename_extension(argv[argc-1], filename, txt, MaxFilenameLength, verbose) == 0) {
    return 0;
  }
  if(verbose.verbose) {
    printf("Reading %s\n", filename);
  }
  closePSRData(lrfs, 0, 0, verbose);
  if(!openPSRData(lrfs, filename, 0, 0, 1, 0, verbose)) {
    return 0;
  }
  if(PSRDataHeader_parse_commandline(lrfs, argc, argv, verbose) == 0) {
    return 0;
  }
  double period;
  int ret;
  ret = get_period(*lrfs, 0, &period, verbose);
  if(ret == 2) {
    printerror(verbose.debug, "ERROR pspecFig (%s): Cannot obtain period", lrfs->filename);
    return 0;
  }
  if(period < 0.001) {
    printerror(verbose.debug, "ERROR pspecFig (%s): The period does not appear to be set in the header. Consider using the -header option.", lrfs->filename);
    return 0;
  }
  if(get_tsamp(*lrfs, 0, verbose) < 0.0000001) {
    printerror(verbose.debug, "ERROR pspecFig (%s): The sampling time does not appear to be set in the header. Consider using the -header option.", lrfs->filename);
    return 0;
  }
  if(verbose.verbose) {
    printf("%ldx%ld points read from lrfs\n", lrfs->NrBins, lrfs->NrSubints);
  }
  if(lrfs->NrPols > 1) {
    if(preprocess_polselect(*lrfs, &clone, 0, verbose) == 0) {
      printerror(verbose.debug, "ERROR pspecFig (%s): Error selecting polarization channel 0.", lrfs->filename);
      return 0;
    }
    swap_orig_clone(lrfs, &clone, verbose);
  }
  long i;
  double max;
  if(plotoptions->normalise_spectra) {
    for(i = 0; i < lrfs->NrSubints * lrfs->NrBins; i++) {
      if(i == 0 || fabs(lrfs->data[i]) > max) {
 max = fabs(lrfs->data[i]);
      }
    }
    for(i = 0; i < lrfs->NrSubints * lrfs->NrBins; i++) {
      lrfs->data[i] /= 0.999*max;
    }
  }
  return 1;
}
int loadHeaderPulseStack(datafile_definition *AverageProfile, datafile_definition lrfs, double *period, int type_of_plots, int argc, char **argv, verbose_definition verbose)
{
  if(verbose.verbose)
    printf("Reading %s\n", argv[argc-1]);
  closePSRData(AverageProfile, 0, 0, verbose);
  if(!openPSRData(AverageProfile, argv[argc-1], 0, 0, 0, 0, verbose))
    return 0;
  if(!readHeaderPSRData(AverageProfile, 0, 0, verbose))
    return 0;
  if(PSRDataHeader_parse_commandline(AverageProfile, argc, argv, verbose) == 0)
    return 0;
  if(type_of_plots == 0) {
    if(AverageProfile->NrBins != lrfs.NrBins) {
      printwarning(verbose.debug, "WARNING: Nr of bins in pulse stack and LRFS do not match. It looks like data is rebinned? Check the units.");
    }
  }
  int ret_prd;
  ret_prd = get_period(*AverageProfile, 0, period, verbose);
  if(ret_prd == 2) {
    printerror(verbose.debug, "ERROR pspecFig (%s): Cannot obtain period", AverageProfile->filename);
    return 0;
  }
  if(*period < 0.001) {
    printerror(verbose.debug, "ERROR pspecFig (%s): The period does not appear to be set in the header. Consider using the -header option.", AverageProfile->filename);
    return 0;
  }
  if(get_tsamp(*AverageProfile, 0, verbose) < 0.0000001) {
    printerror(verbose.debug, "ERROR pspecFig (%s): The sampling time does not appear to be set in the header. Consider using the -header option.", AverageProfile->filename);
    return 0;
  }
  closePSRData(AverageProfile, 1, 0, verbose);
  AverageProfile->format = MEMORY_format;
  AverageProfile->NrSubints = 1;
  AverageProfile->NrFreqChan = 1;
  AverageProfile->NrPols = 1;
  AverageProfile->data = malloc(AverageProfile->NrBins*sizeof(float));
  if(AverageProfile->data == NULL) {
    printerror(verbose.debug, "ERROR pspecFig (%s): Memory allocation error", AverageProfile->filename);
    return 0;
  }
  return 1;
}
int load2dfs(twodfs_def *twodfs_allinfo, datafile_definition AverageProfile, int file_number, int extprefix, int altname, int argc, char **argv, plotoptions_def plotoptions, verbose_definition verbose)
{
  char filename[MaxFilenameLength], txt[MaxFilenameLength];
  datafile_definition clone;
  if(altname) {
    strcpy(filename, argv[altname]);
  }else {
    if(extprefix == 0) {
      sprintf(txt, "%d.2dfs", file_number);
    }
    if(change_filename_extension(argv[argc-1], filename, txt, MaxFilenameLength, verbose) == 0) {
      return 0;
    }
  }
  if(verbose.verbose) {
    printf("Reading %s\n", filename);
  }
  closePSRData(&(twodfs_allinfo->datafile), 0, 0, verbose);
  if(!openPSRData(&(twodfs_allinfo->datafile), filename, 0, 0, 1, 0, verbose)) {
    return 0;
  }
  if(PSRDataHeader_parse_commandline(&(twodfs_allinfo->datafile), argc, argv, verbose) == 0) {
    return 0;
  }
  if(verbose.verbose) {
    printf("%ldx%ld points read from 2dfs\n", twodfs_allinfo->datafile.NrBins, twodfs_allinfo->datafile.NrSubints);
  }
  if(twodfs_allinfo->datafile.NrPols > 1) {
    if(preprocess_polselect(twodfs_allinfo->datafile, &clone, 0, verbose) == 0) {
      printerror(verbose.debug, "ERROR pspecFig: Error selecting polarization channel 0.");
      return 0;
    }
    swap_orig_clone(&(twodfs_allinfo->datafile), &clone, verbose);
  }
  long i;
  double max;
  if(plotoptions.normalise_spectra) {
    for(i = 0; i < twodfs_allinfo->datafile.NrSubints * twodfs_allinfo->datafile.NrBins; i++) {
      if(i == 0 || fabs(twodfs_allinfo->datafile.data[i]) > max) {
 max = fabs(twodfs_allinfo->datafile.data[i]);
      }
    }
    for(i = 0; i < twodfs_allinfo->datafile.NrSubints * twodfs_allinfo->datafile.NrBins; i++) {
      twodfs_allinfo->datafile.data[i] /= 0.999*max;
    }
  }
  if(twodfs_allinfo->f2_min == 0.0 && twodfs_allinfo->f2_min == 0.0) {
    if(twodfs_allinfo->datafile.xrangeset) {
      twodfs_allinfo->f2_min = twodfs_allinfo->datafile.xrange[0];
      twodfs_allinfo->f2_max = twodfs_allinfo->datafile.xrange[1];
    }else {
      twodfs_allinfo->f2_min = -AverageProfile.NrBins/2.0;
      twodfs_allinfo->f2_max = -AverageProfile.NrBins/2.0 + AverageProfile.NrBins*(twodfs_allinfo->datafile.NrBins-1.0)/(float)twodfs_allinfo->datafile.NrBins;
    }
    float f2_resolution;
    f2_resolution = (twodfs_allinfo->f2_max - twodfs_allinfo->f2_min)/(float)(twodfs_allinfo->datafile.NrBins-1);
    twodfs_allinfo->f2_min -= 0.5*f2_resolution;
    twodfs_allinfo->f2_max += 0.5*f2_resolution;
    if(plotoptions.doFlip) {
      float swap;
      swap = twodfs_allinfo->f2_min;
      twodfs_allinfo->f2_min = -twodfs_allinfo->f2_max;
      twodfs_allinfo->f2_max = -swap;
    }
  }
  twodfs_allinfo->f3_min = 0;
  twodfs_allinfo->f3_max = 0.5;
  if(twodfs_allinfo->SelectP3Integrate) {
    if(twodfs_allinfo->P3IntegrateHigh > twodfs_allinfo->f3_max) {
      twodfs_allinfo->P3IntegrateHigh = twodfs_allinfo->f3_max;
    }
  }
  if(twodfs_allinfo->f3_min < 0) {
    twodfs_allinfo->f3_min = 0;
  }
  char *notes;
  notes = get_history_notes_last(&(twodfs_allinfo->datafile));
  if(notes != NULL) {
    int ret, nrbins_from_2dfs_header, onpulse_start_from_2dfs_header, onpulse_end_from_2dfs_header, nfft_from_2dfs_header;
    double power_drift, power_tot, percentage_power_drift, percentage_power_drift_error;
    ret = sscanf(notes, "2DFS generated from pulse stack with %d bins, using bins %d to %d and nfft=%d. Asymmetric power is %lf / %lf = (%lf +- %lf)%%'", &nrbins_from_2dfs_header, &onpulse_start_from_2dfs_header, &onpulse_end_from_2dfs_header, &nfft_from_2dfs_header, &power_drift, &power_tot, &percentage_power_drift, &percentage_power_drift_error);
    if(ret == 8) {
      twodfs_allinfo->percentage_power_drift = percentage_power_drift;
      twodfs_allinfo->percentage_power_drift_error = percentage_power_drift_error;
    }
  }
  return 1;
}

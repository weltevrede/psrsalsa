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
#define max_nr_zap 10250
#define max_nr_stack 1000
#define MaxNrPointsInPolygon 20
enum xUnitsSwitch_values {
  XUNIT_BINS,
  XUNIT_DEG,
  XUNIT_PHASE,
  XUNIT_TIME,
  XUNIT_ENDOFLIST
};
float ypos(float *stackI, int bin, long PulseNr, float scale, int NrBins, long pulse_bottom, long pulse_top, long subint_start, float yUnitCmdLine, float dyshift);
int setBaselineParams(datafile_definition fin, float *baseline, float *dxshift, int *xUnitsSwitch, verbose_definition verbose);
int zapVectors(int nrZappedVectors, int *zappedVectors, long subint_start, datafile_definition fin, float *dataSubset, int onlysinglechannel, verbose_definition verbose);
int zapSubints(int nrZappedSubints, int *zappedSubints, long subint_start, int didtranspose_orig_nrbin, datafile_definition fin, float *dataSubset, int onlysinglesubint, verbose_definition verbose);
void yvec2unit(datafile_definition fin, int yUnitsSwitch, int vectorMinMaxSpecified, float vectorMin, float vectorMax, int didtranspose_orig_nrbin, int type, float valuein, float *valueout, int mapmode, int inverse, verbose_definition verbose);
typedef struct {
  long subint_start, subint_end;
  long x1, x2;
  int grayscalemode;
  int nrZappedVectors, nrZappedSubints;
}plot_state_def;
void copystackstate(plot_state_def state1, plot_state_def *state2)
{
  state2->subint_start = state1.subint_start;
  state2->subint_end = state1.subint_end;
  state2->x1 = state1.x1;
  state2->x2 = state1.x2;
  state2->grayscalemode = state1.grayscalemode;
  state2->nrZappedVectors = state1.nrZappedVectors;
  state2->nrZappedSubints = state1.nrZappedSubints;
}
typedef struct {
  float characterheight;
  int linewidth;
  int font;
}pgplot_text_state_def;
int main(int argc, char **argv)
{
  long subint_start, subint_end, subint_range_defined;
  int bin1_start, bin2_start;
  float dxshift_start;
  float dyshift;
  int change_pulse_numbering_flag, change_pulse_numbering_value;
  int longitudeRangeSet_flag;
  float longitude_start, longitude_end;
  int xUnitsSwitch;
  int yUnitsSwitch;
  float yUnitCmdLine;
  int vectorMinMaxSpecified;
  float vectorMin, vectorMax;
  int yUnitsTobs;
  int yUnitsMHz;
  int fixverticalscale_flag;
  int current_polnr;
  int interactive_flag;
  int *zappedVectors, *zappedSubints;
  char *inputfilename;
  int viewportOptionsSet;
  float viewport_startx, viewport_starty, viewport_endx, viewport_endy;
  float ymargin_both, ymargin_top;
  int nrpanelsx, nrpanelsy;
  int appendframes_flag;
  int nokeypresses_flag;
  float scale_start, scale2_start;
  plot_state_def *stack_state;
  int current_stack_pos;
  int grayscalemode_start;
  int showTwice_flag;
  int nonumside_flag;
  int disable_x_numbers, disable_y_numbers;
  int showtop, showright;
  int showwedge;
  int histswitch;
  int polymode_flag;
  int noboxx, noboxy;
  int heading_string_index;
  pgplot_text_state_def heading_font;
  pgplot_frame_def_internal pgplot_frame;
  pgplot_text_state_def title_font;
  pgplot_text_state_def label_font;
  pgplot_text_state_def box_font;
  pgplot_text_state_def copy_font;
  int plotlw;
  int notitleset;
  char title[1000];
  int xtitle_set, ytitleset;
  char xtitle[1000], ytitle[1000];
  int wedgelabel_set;
  char wedgelabel[1000];
  long i, j, k;
  int ignorebins, ignorebins2;
  int ignorelastpulses, scalerange;
  float xmax_oscil, xmin_oscil, scalerange_min, scalerange_max;
  psrsalsaApplication application;
  datafile_definition fin;
  initApplication(&application, "pplot", "[options] inputfile(s)");
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_rebin = 1;
  application.switch_device = 1;
  application.switch_tscr = 1;
  application.switch_tscr_complete = 1;
  application.switch_TSCR = 1;
  application.switch_fscr = 1;
  application.switch_FSCR = 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_debase = 1;
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_itf = 1;
  application.switch_formatlist = 1;
  application.switch_iformat = 1;
  application.switch_noweights = 1;
  application.switch_useweights = 1;
  application.switch_uniformweights = 1;
  application.switch_polselect = 1;
  application.switch_nocounters = 1;
  application.switch_conshift= 1;
  application.switch_circshift= 1;
  application.switch_dedisperse= 1;
  application.switch_deFaraday= 1;
  application.switch_stokes = 1;
  application.switch_coherence = 1;
  application.switch_filelist = 1;
  application.switch_size = 1;
  application.switch_macro = 1;
  application.switch_cmaplist = 1;
  application.switch_cmap = 1;
  application.switch_onpulse = 1;
  application.switch_onpulsef = 1;
  application.switch_onpulsegr = 1;
  application.switch_deparang = 1;
  application.switch_changeRefFreq = 1;
  application.switch_norm = 1;
  application.switch_normglobal = 1;
  application.switch_clip = 1;
  application.switch_align = 1;
  application.switch_templatedata = 1;
  application.switch_template = 1;
  application.switch_libversions = 1;
  application.cmap = PPGPLOT_INVERTED_HEAT;
  strcpy(application.pgplotdevice, "/xs");
  interactive_flag = 0;
  dxshift_start = 0;
  dyshift = 0;
  xUnitsSwitch = XUNIT_BINS;
  yUnitsSwitch = 0;
  yUnitCmdLine = 1;
  longitudeRangeSet_flag = 0;
  current_polnr = 0;
  subint_range_defined = 0;
  heading_font.characterheight = 1;
  heading_font.linewidth = 1;
  heading_font.font = 1;
  heading_string_index = 0;
  title_font.characterheight = 1;
  title_font.linewidth = 1;
  title_font.font = 1;
  box_font.characterheight = 1;
  box_font.linewidth = 1;
  box_font.font = 1;
  label_font.characterheight = 1;
  label_font.linewidth = 1;
  label_font.font = 1;
  copy_font.characterheight = 1.3;
  copy_font.linewidth = 1;
  copy_font.font = 3;
  notitleset = 1;
  title[0] = 0;
  xtitle_set = 0;
  ytitleset = 0;
  wedgelabel_set = 0;
  viewportOptionsSet = 0;
  viewport_startx = 0.15;
  viewport_starty = 0.15;
  viewport_endx = 0.9;
  viewport_endy = 0.9;
  ymargin_both = 0;
  ymargin_top = 0;
  appendframes_flag = 0;
  showTwice_flag = 0;
  fixverticalscale_flag = 0;
  nokeypresses_flag = 0;
  grayscalemode_start = -1;
  scale_start = 1;
  scale2_start = 0;
  nrpanelsx = 1;
  nrpanelsy = 1;
  vectorMinMaxSpecified = 0;
  bin1_start = -1;
  bin2_start = -1;
  polymode_flag = 1;
  ignorebins = 0;
  ignorebins2 = 0;
  showtop = 0;
  showright = 0;
  showwedge = 0;
  change_pulse_numbering_flag = 0;
  change_pulse_numbering_value = 0;
  plotlw = 1;
  ignorelastpulses = 0;
  histswitch = 0;
  xmax_oscil = 1e10;
  xmin_oscil = -1e10;
  noboxx = noboxy = 0;
  disable_x_numbers = 0;
  disable_y_numbers = 0;
  nonumside_flag = 0;
  yUnitsTobs = 0;
  yUnitsMHz = 0;
  clear_pgplot_frame(&pgplot_frame);
  scalerange = 0;
  if(argc < 2) {
    printf("Program to plot pulsar data in various ways.\n");
    printApplicationHelp(&application);
    printf("Other options:\n");
    printf("-interactive    Turn interactive mode on.\n");
    printf("-ia             short for -interactive.\n");
    printf("\nData selection options:\n");
    printf("-prange     Specify subint range to be plotted, default is all.\n");
    printf("-b          Specify bin range to be plotted, default is all.\n");
    printf("-long       Specify longitude range in degrees.\n");
    printf("-time       Specify longitude range in seconds.\n");
    printf("-phase      Specify longitude range in phase.\n");
    printf("\nPanel options:\n");
    printf("-appendframes Remove space between the frames of multiple adjacent panels\n");
    printf("-N            \"nrx nry\" Create nrx by nry panels, rather than a single plot\n");
    printf("              per page\n");
    printf("\nPlot options:\n");
    printf("-hist         Turn on histogram mode. The bin centre is plotted at an integer\n");
    printf("              bin nr. If plotting versus pulse longitude/time/phase, you can\n");
    printf("              use the -dl option.\n");
    printf("-lw           set line width of line plot (not used in map mode).\n");
    printf("-lp           Show graphics as a \"joy division\" line plot rather than a colour\n");
    printf("              map (-map). For a single row of data -lp is the default.\n");
    printf("-map          Show graphics as a colour map (default if multiple rows in data).\n");
    printf("              The bin centre is plotted at an integer bin nr. If plotting versus\n");
    printf("              pulse longitude/time/phase, you can use the -dl option.\n");
    printf("-nokeypress   Do not wait for key presses to go to next page in output plot\n");
    printf("-nopoly       Do not use filled polygons but only line drawing, which sometimes\n");
    printf("              could give better results (or worse).\n");
    printf("-showwedge    Plot an annotated wedge to show color scale.\n");
    printf("-showtop      Show a panel at top of the map (use with -map).\n");
    printf("-showright    Show a panel at right of the map (use with -map).\n");
    printf("-showtwice    Plot the map twice above each other (use with -map).\n");
    printf("-textkeywords Lists of keywords you can use with");
    printf(" -title.\n");
    printf("-dx           Set viewport x-range (start) to this value (default is %.2f).\n", viewport_startx);
    printf("-x            Set viewport x-range (end) to this value (default is %.2f).\n", viewport_endx);
    printf("-ys           Set viewport y-range (start) to this value (default is %.2f).\n", viewport_starty);
    printf("-y            Set viewport y-range (end) to this value (default is %.2f).\n", viewport_endy);
    printf("\nScaling and offset of data:\n");
    printf("-scale      Specify scale, default is 1. This option multiplies the intensity by\n");
    printf("            this factor. In -map mode, a value > 1 results in clipping, making\n");
    printf("            weak features clearer.\n");
    printf("-scale2     Specify scale2 (a value between 0 and 1), default is 0, only used in\n");
    printf("            map mode. This clips the low intensity values to emphasize the\n");
    printf("            bright features.\n");
    printf("-scalerange \"min max\" In map mode: The used color range is between min and max.\n");
    printf("-dl         Shift plot in pulse longitude.\n");
    printf("-half       Shift pulses by half a pulse.\n");
    printf("-0          Force first pulse to be plotted as pulse zero (ignored with -map).\n");
    printf("-1          Force first pulse to be plotted as pulse one (ignored with -map).\n");
    printf("-yunit      Multiply the y-axis with this number\n");
    printf("-vecRange   Specify value corresponding to first and last vector (vertical axis).\n");
    printf("-tobs       For subint/pulse phase plots, use observing time as y-axis.\n");
    printf("-MHz        For frequency/pulse phase plots, use MHz as y-axis.\n");
    printf("\nLabel options:\n");
    printf("-title      \"...\" Specify title (supports keywords listed with -textkeywords)\n");
    printf("-heading    \"...\" Specify heading.\n");
    printf("-xtitle     \"...\" Set title of xaxis.\n");
    printf("-ytitle     \"...\" Set title of yaxis.\n");
    printf("-wedgetitle \"...\" Set title of wedge generated with -showwedge.\n");
    printf("-nonumx           Disable numbering x-as.\n");
    printf("-nonumy           Disable numbering y-as.\n");
    printf("-nonumside        Disable numbering on the side panels.\n");
    printf("-noboxx           Disable drawing the x-axis.\n");
    printf("-noboxy           Disable drawing the y-axis.\n");
    printf("-labels           \"heading_ch heading_lw heading_f title_ch title_lw title_f label_ch label_lw label_f box_ch box_lw box_f");
    printf("\".\n");
    printf("                  default is \"%.1f %d %d %.1f %d %d %.1f %d %d %.1f %d %d", heading_font.characterheight, heading_font.linewidth, heading_font.font, title_font.characterheight, title_font.linewidth, title_font.font, label_font.characterheight, label_font.linewidth, label_font.font, box_font.characterheight, box_font.linewidth, box_font.font);
    printf("\".\n");
    printf("\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      int index;
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
 i = index;
      }else if(strcmp(argv[i], "-nokeypress") == 0) {
 nokeypresses_flag = 1;
      }else if(strcmp(argv[i], "-prange") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld %ld", &subint_start, &subint_end, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 subint_range_defined = 1;
        i++;
      }else if(strcmp(argv[i], "-appendframes") == 0) {
 appendframes_flag = 1;
      }else if(strcmp(argv[i], "-lw") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &plotlw, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-N") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &nrpanelsx, &nrpanelsy, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-dl") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &dxshift_start, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-vecRange") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &vectorMin, &vectorMax, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 vectorMinMaxSpecified = 1;
 yUnitsSwitch = 1;
        i++;
      }else if(strcmp(argv[i], "-tobs") == 0) {
 yUnitsTobs = 1;
      }else if(strcmp(argv[i], "-MHz") == 0) {
 yUnitsMHz = 1;
      }else if(strcmp(argv[i], "-b") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &bin1_start, &bin2_start, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-textkeywords") == 0) {
 str_list_replace_keys(0);
 terminateApplication(&application);
 return 0;
      }else if(strcmp(argv[i], "-yunit") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &yUnitCmdLine, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-ytitle") == 0) {
 strcpy(ytitle, argv[i+1]);
 ytitleset = 1;
 i++;
      }else if(strcmp(argv[i], "-xtitle") == 0) {
 strcpy(xtitle, argv[i+1]);
 xtitle_set = 1;
 i++;
      }else if(strcmp(argv[i], "-wedgetitle") == 0) {
 strcpy(wedgelabel, argv[i+1]);
 wedgelabel_set = 1;
 i++;
      }else if(strcmp(argv[i], "-nonumside") == 0) {
 nonumside_flag = 1;
      }else if(strcmp(argv[i], "-nonumx") == 0) {
 disable_x_numbers = 1;
      }else if(strcmp(argv[i], "-nonumy") == 0) {
 disable_y_numbers = 1;
      }else if(strcmp(argv[i], "-showtwice") == 0) {
 showTwice_flag = 1;
      }else if(strcmp(argv[i], "-labels") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %d %d %f %d %d %f %d %d %f %d %d", &heading_font.characterheight, &heading_font.linewidth, &heading_font.font, &title_font.characterheight, &title_font.linewidth, &title_font.font, &label_font.characterheight, &label_font.linewidth, &label_font.font, &box_font.characterheight, &box_font.linewidth, &box_font.font, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-scale") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &scale_start, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-scale2") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &scale2_start, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-scalerange") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &scalerange_min, &scalerange_max, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 scalerange = 1;
 i++;
      }else if(strcmp(argv[i], "-title") == 0) {
 strcpy(title, argv[i+1]);
 notitleset = 0;
 i++;
      }else if(strcmp(argv[i], "-heading") == 0) {
 heading_string_index = i+1;
 i++;
      }else if(strcmp(argv[i], "-x") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &viewport_endx, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 viewportOptionsSet = 1;
 i++;
      }else if(strcmp(argv[i], "-dx") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &viewport_startx, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 viewportOptionsSet = 1;
 i++;
      }else if(strcmp(argv[i], "-ys") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &viewport_starty, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 viewportOptionsSet = 1;
 i++;
      }else if(strcmp(argv[i], "-y") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &viewport_endy, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 viewportOptionsSet = 1;
 i++;
      }else if(strcmp(argv[i], "-long") == 0 || strcmp(argv[i], "-time") == 0 || strcmp(argv[i], "-phase") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &longitude_start, &longitude_end, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 xUnitsSwitch = XUNIT_DEG;
 if(strcmp(argv[i], "-phase") == 0) {
   xUnitsSwitch = XUNIT_PHASE;
 }else if(strcmp(argv[i], "-time") == 0) {
   xUnitsSwitch = XUNIT_TIME;
 }
 longitudeRangeSet_flag = 1;
 i++;
      }else if(strcmp(argv[i], "-0") == 0) {
 change_pulse_numbering_flag = 1;
 change_pulse_numbering_value = 0;
      }else if(strcmp(argv[i], "-1") == 0) {
 change_pulse_numbering_flag = 1;
 change_pulse_numbering_value = 1;
      }else if(strcmp(argv[i], "-half") == 0) {
 dyshift = 0.5;
      }else if(strcmp(argv[i], "-hist") == 0) {
 histswitch = 1;
 grayscalemode_start = 0;
      }else if(strcmp(argv[i], "-noboxx") == 0) {
 noboxx = 1;
      }else if(strcmp(argv[i], "-noboxy") == 0) {
 noboxy = 1;
      }else if(strcmp(argv[i], "-lp") == 0) {
 grayscalemode_start = 0;
      }else if(strcmp(argv[i], "-map") == 0) {
 grayscalemode_start = 1;
      }else if(strcmp(argv[i], "-showtop") == 0) {
 showtop = 1;
      }else if(strcmp(argv[i], "-showright") == 0) {
 showright = 1;
      }else if(strcmp(argv[i], "-showwedge") == 0) {
 showwedge = 1;
      }else if(strcmp(argv[i], "-interactive") == 0 || strcmp(argv[i], "-ia") == 0) {
 interactive_flag = 1;
      }else if(strcmp(argv[i], "-nopoly") == 0) {
 polymode_flag = 0;
      }else {
 if(argv[i][0] == '-') {
   printerror(application.verbose_state.debug, "pplot: Unknown option: %s\n\nRun pplot without command line arguments to show help", argv[i]);
   terminateApplication(&application);
   return 0;
 }else {
   if(applicationAddFilename(i, application.verbose_state) == 0)
     return 0;
 }
      }
    }
  }
  if(histswitch && grayscalemode_start) {
    printerror(application.verbose_state.debug, "ERROR pplot: Cannot use the -hist and -map options simultaneously.");
    return 0;
  }
  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pplot: No files specified");
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) > 1 && interactive_flag) {
    printerror(application.verbose_state.debug, "ERROR pplot: Cannot plot multiple files in interactive mode.");
    return 0;
  }
  if(application.macro_ptr != NULL) {
    interactive_flag = 1;
  }
  if(viewportOptionsSet == 0) {
    if(nrpanelsx > 1 || nrpanelsy > 1) {
      viewport_endx = 0.95;
      viewport_startx = 0.05;
      viewport_starty = 0.1;
      viewport_endy = 0.95;
    }
  }
  ppgopen(application.pgplotdevice);
  pgplot_setWindowsize(application.windowwidth, application.windowheight, -1);
  ppgask(0);
  ppgslw(1);
  zappedVectors = malloc(max_nr_zap*sizeof(int));
  zappedSubints = malloc(max_nr_zap*sizeof(int));
  stack_state = malloc(max_nr_stack*sizeof(plot_state_def));
  if(zappedVectors == NULL || zappedSubints == NULL || stack_state == NULL) {
    printerror(application.verbose_state.debug, "ERROR pplot: Memory allocation error.");
    return 0;
  }
  int curpanelnrx, curpanelnry;
  curpanelnrx = 0;
  curpanelnry = 0;
  int goingToPlotFirstPanel;
  goingToPlotFirstPanel = 1;
  int dataSubset_allocated;
  dataSubset_allocated = 0;
  int stack_poly_allocated;
  stack_poly_allocated = 0;
  int maxy_allocated;
  maxy_allocated = 0;
  while((inputfilename = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
    int data_read;
    data_read = 0;
    int nrpolarizations;
    nrpolarizations = 0;
    int didtranspose_orig_nrbin;
    didtranspose_orig_nrbin = 0;
    float baseline;
    float viewport_startxCurrentPanel;
    float viewport_endxCurrentPanel;
    float scale, scale2;
    long *maxy;
    float dxshift;
    int bin1, bin2;
    int zapmode;
    int nx_last_clicked_zap, ny_last_clicked_zap;
    nx_last_clicked_zap = 0;
    ny_last_clicked_zap = 0;
    float *dataSubset;
    float *xstack_poly, *ystack_poly;
    float dummyf1, dummyf2, dummyf3, dummyf4;
    int dummyi;
    char txt1[1000], txt2[1000], txt3[1000];
    maxy = NULL;
    dataSubset = NULL;
    xstack_poly = ystack_poly = NULL;
    dxshift = dxshift_start;
    zapmode = 0;
    if(subint_range_defined == 0) {
      subint_start = -1;
    }
    if(goingToPlotFirstPanel == 0) {
      if(curpanelnrx == 0 && curpanelnry == 0) {
 if(nokeypresses_flag == 0) {
   printf("Press a key to continue\n");
   fflush(stdout);
   pgetch_macro(&application, application.verbose_state);
 }
      }
    }
    goingToPlotFirstPanel = 0;
    current_stack_pos = 0;
    printf("Plotting: %s\n", inputfilename);
    stack_state[0].nrZappedVectors = 0;
    stack_state[0].nrZappedSubints = 0;
    cleanPSRData(&fin, application.verbose_state);
    bin1 = bin1_start;
    bin2 = bin2_start;
    scale = scale_start;
    scale2 = scale2_start;
    stack_state[0].grayscalemode = grayscalemode_start;
    stack_state[1].grayscalemode = grayscalemode_start;
    if(notitleset) {
      strcpy(title, inputfilename);
    }
    int firstPulsedataSubset;
    int ytitle_set_to_intensity, ytitle_set_to_pulsenumber, ytitle_set_to_subint, ytitle_set_to_freqchannel, ytitle_set_to_fluctuationFreq, ytitle_set_to_padist, ytitle_set_to_elldist, ytitle_set_to_pulselongitude, ytitle_set_to_p3fold, ytitle_set_to_spectral_power, ytitle_set_to_harmonic_number, ytitle_set_to_lag_number;
    ytitle_set_to_intensity = 0;
    ytitle_set_to_pulsenumber = 0;
    ytitle_set_to_subint = 0;
    ytitle_set_to_freqchannel = 0;
    ytitle_set_to_fluctuationFreq = 0;
    ytitle_set_to_padist = 0;
    ytitle_set_to_elldist = 0;
    ytitle_set_to_pulselongitude = 0;
    ytitle_set_to_p3fold = 0;
    ytitle_set_to_spectral_power = 0;
    ytitle_set_to_harmonic_number = 0;
    ytitle_set_to_lag_number = 0;
    firstPulsedataSubset = 0;
    do {
      float currentPanelScaling, currentPanelScalingx, currentPanelScalingy;
      float viewport_startyCurrentPanel, viewport_endyCurrentPanel;
      float xmin, xmax, ymin, ymax, xminshow, xmaxshow;
      float x, x2, y, y2;
      int redraw;
      char *newtext, ch;
      x = x2 = y = y2 = ch = 0;
      currentPanelScaling = 1;
      viewport_startxCurrentPanel = viewport_startx;
      if(nrpanelsx > 1) {
 if(appendframes_flag == 0)
   currentPanelScalingx = (viewport_endx-viewport_startx)/(viewport_endx-viewport_startx + nrpanelsx -1);
 else
   currentPanelScalingx = 1.0/nrpanelsx;
 currentPanelScaling = currentPanelScalingx;
 if(appendframes_flag == 0 || curpanelnrx == 0)
   viewport_startxCurrentPanel += curpanelnrx*currentPanelScalingx*(viewport_startx);
 viewport_startxCurrentPanel += curpanelnrx*currentPanelScalingx*(viewport_endx-viewport_startx);
 if(appendframes_flag == 0)
   viewport_startxCurrentPanel += curpanelnrx*currentPanelScalingx*(1-viewport_endx);
      }
      viewport_endxCurrentPanel = viewport_endx;
      if(nrpanelsx > 1) {
 viewport_endxCurrentPanel = viewport_startxCurrentPanel;
 viewport_endxCurrentPanel += currentPanelScalingx*(viewport_endx-viewport_startx);
      }
      viewport_startyCurrentPanel = viewport_starty;
      if(nrpanelsy > 1) {
 currentPanelScalingy = (viewport_endy-viewport_starty)/(viewport_endy-viewport_starty + nrpanelsy -1);
 if(currentPanelScalingy < currentPanelScaling)
   currentPanelScaling = currentPanelScalingy;
 viewport_startyCurrentPanel += (nrpanelsy-curpanelnry-1)*currentPanelScalingy*(viewport_starty);
 viewport_startyCurrentPanel += (nrpanelsy-curpanelnry-1)*currentPanelScalingy*(viewport_endy-viewport_starty);
 viewport_startyCurrentPanel += (nrpanelsy-curpanelnry-1)*currentPanelScalingy*(1-viewport_endy);
      }
      viewport_endyCurrentPanel = viewport_endy;
      if(nrpanelsy > 1) {
 viewport_endyCurrentPanel = viewport_startyCurrentPanel;
 viewport_endyCurrentPanel += currentPanelScalingy*(viewport_endy-viewport_starty);
      }
      if(curpanelnrx == 0 && curpanelnry == 0)
 ppgpage();
      if(data_read == 0 || data_read == 2) {
 int iformat;
 iformat = application.iformat;
 if(access(inputfilename, F_OK) != 0) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot open file %s for reading.", inputfilename);
   free(zappedVectors);
   free(zappedSubints);
   free(stack_state);
   closePSRData(&fin, 0, application.verbose_state);
   terminateApplication(&application);
   return 0;
 }
 if(iformat <= 0)
   iformat = guessPSRData_format(inputfilename, 0, application.verbose_state);
 if(isValidPSRDATA_format(iformat) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Please specify a valid input format with the -iformat option.\n");
   free(zappedVectors);
   free(zappedSubints);
   free(stack_state);
   closePSRData(&fin, 0, application.verbose_state);
   terminateApplication(&application);
   return 0;
 }
 closePSRData(&fin, 0, application.verbose_state);
 if(!openPSRData(&fin, inputfilename, iformat, 0, 1, 0, application.verbose_state)) {
   printerror(application.verbose_state.debug, "ERROR pplot: Error opening file.\n");
   return 0;
 }
 if(PSRDataHeader_parse_commandline(&fin, argc, argv, application.verbose_state) == 0)
   return 0;
 for(i = 1; i < argc; i++) {
   if(strcmp(argv[i], "-header") == 0) {
     printwarning(application.verbose_state.debug, "WARNING pplot: If using the -header option, be aware it applied BEFORE the preprocessing.");
     break;
   }
 }
 if(preprocessApplication(&application, &fin) == 0) {
   return 0;
 }
 if(fin.gentype == GENTYPE_PENERGY) {
   if(application.verbose_state.verbose) {
     printwarning(application.verbose_state.debug, "The file appears to be penergy output. The \"polarization channels\" correspond to peak intensity (on & off-pulse), integrated energy (on & off-pulse), rms (on & off-pulse) and the S/N of the on-pulse region. The \"subintegrations\" correspond to the different polarizations in the original file.\n");
   }
 }
 if((fin.poltype == POLTYPE_ILVPAdPA || fin.poltype == POLTYPE_ILVPAdPATEldEl) && fin.NrSubints == 1 && fin.NrFreqChan == 1) {
   if(application.verbose_state.verbose) {
     printwarning(application.verbose_state.debug, "The file appears to be ppol output. Plotting this data with ppol allows the profile and PA-swing to be shown in a single figure.\n");
   }
 }
 if(fin.isFolded && fin.foldMode == FOLDMODE_FIXEDPERIOD) {
   if(fin.fixedPeriod <= 0.0) {
     printwarning(application.verbose_state.debug, "WARNING pplot: Period does not appear to be set, assuming it is 1 sec.");
     fin.fixedPeriod = 1.0;
   }
 }
 if(fin.tsampMode == TSAMPMODE_LONGITUDELIST) {
   if(convert_to_fixed_tsamp(&fin, application.verbose_state) != 1) {
     printerror(application.verbose_state.debug, "ERROR pplot: Cannot convert sampling to a regular grid.");
     return 0;
   }
 }
 if(fin.isFolded && fin.tsampMode == TSAMPMODE_FIXEDTSAMP) {
   if(get_tsamp(fin, 0, application.verbose_state) <= 0.0) {
     printwarning(application.verbose_state.debug, "WARNING pplot: Assuming full period is stored.");
     double period;
     int ret;
     ret = get_period(fin, 0, &period, application.verbose_state);
     if(ret == 2) {
       printerror(application.verbose_state.debug, "ERROR pplot (%s): Cannot obtain period", fin.filename);
       return 0;
     }
     fin.fixedtsamp = period/(double)fin.NrBins;
   }
 }
 region_frac_to_int(&(application.onpulse), fin.NrBins, 0);
 if(data_read == 0) {
   if(setBaselineParams(fin, &baseline, &dxshift, &xUnitsSwitch, application.verbose_state) == 0)
     return 0;
   if(fin.yrangeset) {
     yUnitsSwitch = 1;
   }
   if(longitudeRangeSet_flag == 0) {
     longitude_start = 0 + dxshift;
     longitude_end = baseline + dxshift;
   }
 }
 nrpolarizations = fin.NrPols;
 if(fin.NrPols > 1) {
   datafile_definition clone;
   if(current_polnr >= fin.NrPols)
     current_polnr = 0;
   if(preprocess_polselect(fin, &clone, current_polnr, application.verbose_state) == 0) {
     printerror(application.verbose_state.debug, "ERROR pplot: Error selecting polarization channel %d.", current_polnr);
     return 0;
   }
   swap_orig_clone(&fin, &clone, application.verbose_state);
 }
 didtranspose_orig_nrbin = 0;
 ytitle_set_to_intensity = 0;
 ytitle_set_to_pulsenumber = 0;
 ytitle_set_to_subint = 0;
 ytitle_set_to_freqchannel = 0;
 ytitle_set_to_fluctuationFreq = 0;
 ytitle_set_to_padist = 0;
 ytitle_set_to_elldist = 0;
 ytitle_set_to_pulselongitude = 0;
 ytitle_set_to_p3fold = 0;
 ytitle_set_to_spectral_power = 0;
 ytitle_set_to_harmonic_number = 0;
 ytitle_set_to_lag_number = 0;
 if(fin.NrFreqChan > 1) {
   datafile_definition clone;
   didtranspose_orig_nrbin = fin.NrBins;
   if(preprocess_transposeRawFBdata(fin, &clone, application.verbose_state) == 0) {
     printerror(application.verbose_state.debug, "ERROR pplot: Error transposing data.");
     return 0;
   }
   if(longitudeRangeSet_flag == 0) {
     longitude_start = 0;
     longitude_end = clone.NrSubints;
     if(fin.NrSubints > 1)
       xUnitsSwitch = XUNIT_PHASE;
     else
       xUnitsSwitch = XUNIT_BINS;
   }
   clone.NrSubints = fin.NrFreqChan;
   clone.NrBins = fin.NrSubints * fin.NrBins;
   clone.NrFreqChan = 1;
   if(ytitleset == 0)
     ytitle_set_to_freqchannel = 1;
   if(yUnitsMHz)
     yUnitsSwitch = 1;
   if(xUnitsSwitch == XUNIT_BINS) {
     baseline = fin.NrSubints * fin.NrBins;
   }else if(xUnitsSwitch == XUNIT_DEG) {
     baseline = fin.NrSubints * 360;
   }else if(xUnitsSwitch == XUNIT_PHASE) {
     baseline = fin.NrSubints;
   }else if(xUnitsSwitch == XUNIT_TIME) {
     baseline = fin.NrSubints * fin.NrBins * get_tsamp(fin, 0, application.verbose_state);
   }
   swap_orig_clone(&fin, &clone, application.verbose_state);
 }else {
   if(ytitleset == 0) {
     if((fin.gentype == GENTYPE_2DFS || fin.gentype == GENTYPE_LRFS || fin.gentype == GENTYPE_S2DFSP3 || fin.gentype == GENTYPE_S2DFSP2) && (vectorMinMaxSpecified == 0 && fin.yrangeset == 1)) {
       ytitle_set_to_fluctuationFreq = 1;
     }else if(fin.gentype == GENTYPE_P3FOLD) {
       ytitle_set_to_p3fold = 1;
     }else if(fin.gentype == GENTYPE_HRFS_UNFOLDED) {
       ytitle_set_to_spectral_power = 1;
     }else if(fin.gentype == GENTYPE_HRFS) {
       ytitle_set_to_harmonic_number = 1;
     }else if(fin.gentype == GENTYPE_LRAC) {
       ytitle_set_to_lag_number = 1;
     }else if((fin.gentype == GENTYPE_LRCC) && (vectorMinMaxSpecified == 0 && fin.yrangeset == 1)) {
       ytitle_set_to_pulselongitude = 1;
     }else if(fin.gentype == GENTYPE_PADIST&& (vectorMinMaxSpecified == 0 && fin.yrangeset == 1)) {
       ytitle_set_to_padist = 1;
     }else if(fin.gentype == GENTYPE_ELLDIST&& (vectorMinMaxSpecified == 0 && fin.yrangeset == 1)) {
       ytitle_set_to_elldist = 1;
     }else if((fin.gentype == GENTYPE_RMMAP) && (vectorMinMaxSpecified == 0 && fin.yrangeset == 1)) {
       strcpy(ytitle, "RM (rad/m\\u2\\d)");
     }else if((fin.gentype == GENTYPE_PULSESTACK || fin.gentype == GENTYPE_PROFILE) && fin.NrSubints == 1) {
       ytitle_set_to_intensity = 1;
     }else if((fin.gentype == GENTYPE_SEARCHMODE) && fin.NrSubints == 1 && fin.NrFreqChan == 1) {
       ytitle_set_to_intensity = 1;
     }else if(fin.gentype == GENTYPE_PULSESTACK) {
       ytitle_set_to_pulsenumber = 1;
       if(yUnitsTobs && fin.NrSubints > 1)
  yUnitsSwitch = 1;
     }else {
       ytitle_set_to_subint = 1;
       if(yUnitsTobs && fin.NrSubints > 1)
  yUnitsSwitch = 1;
     }
   }
 }
 if(data_read == 0) {
   if(subint_start >= fin.NrSubints) {
     subint_start = fin.NrSubints -1;
     subint_end = fin.NrSubints -1;
   }
   if(subint_start < 0) {
     subint_start = 0;
     subint_end = fin.NrSubints -1;
   }
   if(subint_end >= fin.NrSubints) {
     subint_end = fin.NrSubints -1;
   }
   if(bin1 < 0) {
     bin1 = 0;
     bin2 = fin.NrBins -1;
   }
   if(xUnitsSwitch != XUNIT_BINS) {
     bin1 = fin.NrBins*(longitude_start-dxshift)/baseline;
     bin2 = fin.NrBins*(longitude_end-dxshift)/baseline - 1;
     if(bin2 >= fin.NrBins)
       bin2 = fin.NrBins-1;
     if(bin1 < 0)
       bin1 = 0;
   }
   if(application.verbose_state.verbose) printf("Specified binrange: %d - %d\n", bin1, bin2);
 }
 if(dataSubset_allocated)
   free(dataSubset);
 dataSubset = (float *)malloc((subint_end-subint_start+1)*fin.NrBins*sizeof(float));
 dataSubset_allocated = 1;
 if(stack_poly_allocated) {
   free(xstack_poly);
   free(ystack_poly);
 }
 if(maxy_allocated)
   free(maxy);
 if(polymode_flag) {
   if(histswitch) {
     xstack_poly = (float *)malloc(2*(fin.NrBins+2)*sizeof(float));
     ystack_poly = (float *)malloc(2*(fin.NrBins+2)*sizeof(float));
   }else {
     xstack_poly = (float *)malloc((fin.NrBins+2)*sizeof(float));
     ystack_poly = (float *)malloc((fin.NrBins+2)*sizeof(float));
   }
   stack_poly_allocated = 1;
 }else {
   maxy = (long *)malloc(fin.NrBins*sizeof(long));
   maxy_allocated = 1;
 }
 if(dataSubset == NULL
    || (maxy_allocated && maxy == NULL)
    || (stack_poly_allocated && (xstack_poly == NULL || ystack_poly == NULL))) {
   printerror(application.verbose_state.debug, "ERROR pplot: Memory allocation error.");
   if(dataSubset == NULL && application.verbose_state.debug) {
     printerror(application.verbose_state.debug, "Cannot allocate %ld bytes", (subint_end-subint_start+1)*fin.NrBins*sizeof(float));
     printerror(application.verbose_state.debug, "subint_end=%ld subint_start=%ld fin.NrBins=%ld", subint_end, subint_start, fin.NrBins);
   }
   return 0;
 }
 for(i = subint_start; i <= subint_end; i++) {
   if(!readPulsePSRData(&fin, i, 0, 0, 0, fin.NrBins, &dataSubset[(i-subint_start)*fin.NrBins], application.verbose_state)) {
     printerror(application.verbose_state.debug, "ERROR pplot: Cannot read data.");
     return 0;
   }
 }
 if(stack_state[current_stack_pos].grayscalemode == -1) {
   if(fin.NrSubints > 1)
     stack_state[current_stack_pos].grayscalemode = 1;
   else
     stack_state[current_stack_pos].grayscalemode = 0;
 }
 if(data_read == 0) {
   if(change_pulse_numbering_flag) {
     if(interactive_flag != 0 || stack_state[current_stack_pos].grayscalemode != 0) {
       printwarning(application.verbose_state.debug, "WARNING: -map mode and interactive mode are incompatable with the -0 and -1 options. Commands are ignored.");
     }else {
       subint_end -= subint_start - change_pulse_numbering_value;
       subint_start = change_pulse_numbering_value;
     }
   }
   firstPulsedataSubset = subint_start;
   stack_state[current_stack_pos].subint_start = subint_start;
   stack_state[current_stack_pos].subint_end = subint_end;
   stack_state[current_stack_pos].x1 = bin1;
   stack_state[current_stack_pos].x2 = bin2;
   stack_state[current_stack_pos].nrZappedVectors = 0;
   stack_state[current_stack_pos].nrZappedSubints = 0;
   current_stack_pos++;
 }
 if(current_stack_pos > 0) {
   if(stack_state[current_stack_pos-1].nrZappedSubints > 0 || stack_state[current_stack_pos-1].nrZappedVectors) {
     if(zapVectors(stack_state[current_stack_pos-1].nrZappedVectors, zappedVectors, subint_start, fin, dataSubset, -1, application.verbose_state) == 0)
       return 0;
     if(zapSubints(stack_state[current_stack_pos-1].nrZappedSubints, zappedSubints, subint_start, didtranspose_orig_nrbin, fin, dataSubset, -1, application.verbose_state) == 0)
       return 0;
   }
 }
      }
      if(setBaselineParams(fin, &baseline, &dxshift, &xUnitsSwitch, application.verbose_state) == 0)
 return 0;
      if(heading_string_index && curpanelnrx == 0 && curpanelnry == 0) {
 ppgsch(heading_font.characterheight);
 ppgslw(heading_font.linewidth);
 ppgscf(heading_font.font);
 ppgsvp(0.1, 0.9, 0.1, 0.9);
 ppgswin(0, 1, 0, 1);
 ppgptxt(0.5, 1.07, 0, 0.5, argv[heading_string_index]);
      }
      ppgsvp(viewport_startxCurrentPanel, viewport_endxCurrentPanel, viewport_startyCurrentPanel, viewport_endyCurrentPanel);
      xmin = stack_state[current_stack_pos-1].x1;
      xmax = stack_state[current_stack_pos-1].x2;
      if(histswitch || stack_state[current_stack_pos-1].grayscalemode) {
 xmin -= 0.5;
 xmax += 0.5;
      }
      if(xUnitsSwitch != XUNIT_BINS) {
 xmin *= baseline/(float)fin.NrBins;
 xmin += dxshift;
 xmax *= baseline/(float)fin.NrBins;
 xmax += dxshift;
      }
      y = stack_state[current_stack_pos-1].subint_start;
      y2 = stack_state[current_stack_pos-1].subint_end;
      if(stack_state[current_stack_pos-1].subint_start == stack_state[current_stack_pos-1].subint_end) {
 if(stack_state[current_stack_pos-1].grayscalemode == 0) {
   y = y2 = ypos(dataSubset, stack_state[current_stack_pos-1].x1, stack_state[current_stack_pos-1].subint_start, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
 }
      }
      if(fixverticalscale_flag) {
 y = 0;
 y2 = 1;
      }
      if(fixverticalscale_flag == 0 && stack_state[current_stack_pos-1].grayscalemode == 0) {
 for(i = stack_state[current_stack_pos-1].subint_start; i <= stack_state[current_stack_pos-1].subint_end; i++) {
   for(j = stack_state[current_stack_pos-1].x1; j <= stack_state[current_stack_pos-1].x2; j++) {
     if(ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) < y)
       y = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
     if(ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) > y2)
       y2 = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
   }
 }
      }else if(y == y2 && stack_state[current_stack_pos-1].grayscalemode == 0){
 y2 += 1;
      }
      ymin = y-ymargin_both;
      ymax = y2+ymargin_both+ymargin_top;
      xminshow = xmin;
      xmaxshow = xmax;
      pgplot_frame.swin_showtwice = showTwice_flag;
      pgplot_frame.swin_x1 = xminshow;
      pgplot_frame.swin_x2 = xmaxshow;
      pgplot_frame.swin_y1 = ymin;
      pgplot_frame.swin_y2 = ymax;
      if(stack_state[current_stack_pos-1].grayscalemode == 0) {
 ppgswin(xminshow,xmaxshow,ymin,ymax);
      }
      if(stack_state[current_stack_pos-1].grayscalemode) {
 float max, min, oldmin, oldmax;
 float xleft, xright, xleft2, xright2;
 min = max = dataSubset[stack_state[current_stack_pos-1].x1+ignorebins];
 for(i = stack_state[current_stack_pos-1].x1; i <= stack_state[current_stack_pos-1].x2; i++) {
   for(j = 0; j <= stack_state[current_stack_pos-1].subint_end-stack_state[current_stack_pos-1].subint_start; j++) {
     if(dataSubset[(j+stack_state[current_stack_pos-1].subint_start-firstPulsedataSubset)*fin.NrBins+i] > max)
       max = dataSubset[(j+stack_state[current_stack_pos-1].subint_start-firstPulsedataSubset)*fin.NrBins+i];
     if(dataSubset[(j+stack_state[current_stack_pos-1].subint_start-firstPulsedataSubset)*fin.NrBins+i] < min)
       min = dataSubset[(j+stack_state[current_stack_pos-1].subint_start-firstPulsedataSubset)*fin.NrBins+i];
   }
 }
 oldmin = min;
 oldmax = max;
 max = min + (max-min)/scale;
 min = oldmin + (oldmax-oldmin)*scale2;
 xleft = stack_state[current_stack_pos-1].x1+ignorebins;
 xright = stack_state[current_stack_pos-1].x2-ignorebins-ignorebins2;
 xleft2 = 0;
 xright2 = fin.NrBins-1;
 xleft -= 0.5;
 xright += 0.5;
 if(xUnitsSwitch != XUNIT_BINS) {
   xleft *= baseline/(float)(fin.NrBins);
   xright *= baseline/(float)(fin.NrBins);
   xleft += dxshift;
   xright += dxshift;
   xleft2 = 0+dxshift;
   xright2 = baseline*(fin.NrBins-1)/(float)(fin.NrBins)+dxshift;
 }
 printf("Plotting vectors (vertical axis): %ld - %ld\n", stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end);
 if(nonumside_flag) {
   if(showright)
     showright = 2;
   if(showtop)
     showtop = 2;
 }
 ppgsvp(viewport_startxCurrentPanel, viewport_endxCurrentPanel, viewport_startyCurrentPanel, viewport_endyCurrentPanel);
 dummyf1 = (0+change_pulse_numbering_value)*yUnitCmdLine+dyshift;
 dummyf2 = (fin.NrSubints-1+change_pulse_numbering_value)*yUnitCmdLine+dyshift;
 dummyf3 = (stack_state[current_stack_pos-1].subint_start-0.5)*yUnitCmdLine+dyshift;
 dummyf4 = (stack_state[current_stack_pos-1].subint_end+0.5)*yUnitCmdLine+dyshift;
 dummyi = showright;
 if((fin.yrangeset && yUnitsSwitch) || vectorMinMaxSpecified) {
   if(showright)
     dummyi += 2;
 }
 pgplot_options_definition pgplot_options;
 pgplot_clear_options(&pgplot_options);
 pgplot_options.viewport.windowwidth = application.windowwidth;
 pgplot_options.viewport.windowheight = application.windowheight;
 pgplot_options.viewport.dxplot = viewport_startxCurrentPanel-0.15;
 pgplot_options.viewport.xsize = (viewport_endxCurrentPanel-viewport_startxCurrentPanel)/(0.9-0.15);
 pgplot_options.viewport.dyplot = viewport_startyCurrentPanel-0.15;
 pgplot_options.viewport.ysize = (viewport_endyCurrentPanel-viewport_startyCurrentPanel)/(0.9-0.15);
 pgplot_options.viewport.noclear = 1;
 strcpy(pgplot_options.viewport.plotDevice, application.pgplotdevice);
 pgplot_options.viewport.dontopen = 1;
 pgplot_options.viewport.dontclose = 1;
 newtext = str_replace_header_params(fin, title, application.verbose_state);
 if(newtext == NULL) {
   fflush(stdout);
   printwarning(application.verbose_state.debug, "WARNING pplot: Cannot substitute keyword in title");
   pgplot_options.box.title[0] = 0;
 }else {
   strcpy(pgplot_options.box.title, newtext);
   free(newtext);
 }
 pgplot_options.box.title_ch = title_font.characterheight*currentPanelScaling;
 pgplot_options.box.title_lw = title_font.linewidth;
 pgplot_options.box.title_f = title_font.font;
 pgplot_options.box.label_f = title_font.font;
 pgplot_options.box.label_ch = label_font.characterheight*currentPanelScaling;
 pgplot_options.box.box_labelsize = box_font.characterheight*currentPanelScaling;
 pgplot_options.box.box_lw = box_font.linewidth;
 pgplot_options.box.label_lw = box_font.linewidth;
 if(wedgelabel_set)
   strcpy(pgplot_options.box.wedgelabel, wedgelabel);
 int levelset = 1;
 if(didtranspose_orig_nrbin)
   levelset = 0;
 if(scalerange) {
   levelset = 1;
   min = scalerange_min;
   max = scalerange_max;
 }
 if(pgplotMap(&pgplot_options, fin.data, fin.NrBins, fin.NrSubints, xleft2, xright2, xleft, xright, dummyf1, dummyf2, dummyf3, dummyf4, application.cmap, application.itf, 0, 0, NULL, 1, 0, 1, levelset, min, max, 1, 2, dummyi, 0, showtop, 0, plotlw, showwedge, !application.do_noplotsubset, showTwice_flag, application.verbose_state) == 0) {
   printerror(application.verbose_state.debug, "ERROR pplot: Cannot plot data.");
   return 0;
 }
      }else {
 printf("Plotting vectors (vertical axis): %ld - %ld\n", stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end);
 ppgbbuf();
 if(polymode_flag) {
   ppgslw(1);
   for(i = stack_state[current_stack_pos-1].subint_end - ignorelastpulses; i >= stack_state[current_stack_pos-1].subint_start; i--) {
     k = 0;
     x = 0;
     for(j = stack_state[current_stack_pos-1].x1+ignorebins; j <= stack_state[current_stack_pos-1].x2-ignorebins; j++) {
       if(!histswitch) {
  x = j;
       }else {
  x = j - 0.5;
       }
       if(xUnitsSwitch != XUNIT_BINS) {
  x *= baseline/(float)fin.NrBins;
  x += dxshift;
       }
       y = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
       if(k == 0) {
  xstack_poly[k] = x-100*0;
  ystack_poly[k] = ymin-100;
  k++;
       }
       int loopnr;
       for(loopnr = 0; loopnr < 2; loopnr++) {
  if(loopnr == 0) {
    xstack_poly[k] = x;
    ystack_poly[k] = y;
    k++;
  }else {
    float x2;
    x2 = j + 0.5;
    if(xUnitsSwitch != XUNIT_BINS) {
      x2 *= baseline/(float)fin.NrBins;
      x2 += dxshift;
    }
    xstack_poly[k] = x2;
    ystack_poly[k] = y;
    k++;
  }
  if(k == MaxNrPointsInPolygon) {
    xstack_poly[k] = x+100*0;
    ystack_poly[k] = ymin-100;
    k++;
    ppgsfs(1);
    ppgsci(0);
    ppgpoly(k, xstack_poly, ystack_poly);
    xstack_poly[0] = xstack_poly[k-3];
    ystack_poly[0] = ymin-100;
    xstack_poly[1] = xstack_poly[k-3];
    ystack_poly[1] = ystack_poly[k-3];
    xstack_poly[2] = xstack_poly[k-2];
    ystack_poly[2] = ystack_poly[k-2];
    k = 3;
  }
  if(!histswitch) {
    break;
  }
       }
     }
     xstack_poly[k] = x+100*0;
     ystack_poly[k] = ymin-100;
     k++;
     ppgsfs(1);
     ppgsci(0);
     ppgpoly(k, xstack_poly, ystack_poly);
     ppgslw(plotlw);
     ppgsci(1);
     for(j = stack_state[current_stack_pos-1].x1+ignorebins; j <= stack_state[current_stack_pos-1].x2-ignorebins; j++) {
       if(!histswitch) {
  x = j;
       }else {
  x = j - 0.5;
       }
       if(xUnitsSwitch != XUNIT_BINS) {
  x *= baseline/(float)fin.NrBins;
  x += dxshift;
       }
       y = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
       if(j == stack_state[current_stack_pos-1].x1+ignorebins) {
  ppgmove(x, y);
       }else {
  ppgdraw(x, y);
       }
       if(histswitch) {
  x = j + 0.5;
  if(xUnitsSwitch != XUNIT_BINS) {
    x *= baseline/(float)fin.NrBins;
    x += dxshift;
  }
  ppgdraw(x, y);
       }
     }
   }
 }else {
   int skip, didfirstmove;
   double xdouble, ydouble;
   skip = 0;
   ppgslw(plotlw);
   for(j = 0; j < fin.NrBins; j++)
     maxy[j] = stack_state[current_stack_pos-1].subint_start;
   for(i = stack_state[current_stack_pos-1].subint_start; i <= stack_state[current_stack_pos-1].subint_end - ignorelastpulses; i++) {
     x = 0;
     if(xUnitsSwitch != XUNIT_BINS) {
       x *= baseline/(float)fin.NrBins;
       x += dxshift;
     }
     ppgmove(x, ypos(dataSubset, 0,i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift));
     didfirstmove = 0;
     for(j = stack_state[current_stack_pos-1].x1+ignorebins; j <= stack_state[current_stack_pos-1].x2-ignorebins; j++) {
       y = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
       if(y > ypos(dataSubset, j, maxy[j], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) || i == maxy[j]) {
  if(skip == 0) {
    if(!histswitch) {
      x = j;
      if(xUnitsSwitch != XUNIT_BINS) {
        x *= baseline/(float)fin.NrBins;
        x += dxshift;
      }
      if((i != stack_state[current_stack_pos-1].subint_end - ignorelastpulses || x <= xmax_oscil) && (i != stack_state[current_stack_pos-1].subint_start || x >= xmin_oscil)) {
        if(didfirstmove == 0) {
   ppgmove(x, y);
   didfirstmove = 1;
        }else {
   ppgdraw(x, y);
        }
      }
    }else {
      x = j-0.5;
      y = ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
      if(xUnitsSwitch != XUNIT_BINS) {
        x *= baseline/(float)fin.NrBins;
        x += dxshift;
      }
      ppgmove(x, y);
      y = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
      if((i != stack_state[current_stack_pos-1].subint_end - ignorelastpulses || x <= xmax_oscil) && (i != stack_state[current_stack_pos-1].subint_start || x >= xmin_oscil))
        ppgdraw(x, y);
      x = j+0.5;
      if(xUnitsSwitch != XUNIT_BINS) {
        x *= baseline/(float)fin.NrBins;
        x += dxshift;
      }
      if((i != stack_state[current_stack_pos-1].subint_end - ignorelastpulses || x <= xmax_oscil) && (i != stack_state[current_stack_pos-1].subint_start || x >= xmin_oscil))
        ppgdraw(x, y);
    }
  }else {
    if(j != stack_state[current_stack_pos-1].x1) {
      xdouble = ((double)ypos(dataSubset, j-1, maxy[j-1], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) - (double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift))/(double)((double)ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) - (double)((double)ypos(dataSubset, j, maxy[j-1], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, maxy[j-1], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)));
      ydouble = (double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) + xdouble*((double)ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift));
      x = xdouble;
      y = ydouble;
      if(x >= 0 && x <= 1) {
        x += j-1;
        if(xUnitsSwitch != XUNIT_BINS) {
   x *= baseline/(float)fin.NrBins;
   x += dxshift;
        }
        ppgmove(x, y);
      }
    }
    y = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
    x = j;
    if(xUnitsSwitch != XUNIT_BINS) {
      x *= baseline/(float)fin.NrBins;
      x += dxshift;
    }
    if((i != stack_state[current_stack_pos-1].subint_end - ignorelastpulses || x <= xmax_oscil) && (i != stack_state[current_stack_pos-1].subint_start || x >= xmin_oscil))
      ppgdraw(x, y);
    skip = 0;
  }
       }else {
  if(skip == 1) {
    ppgmove(j, y);
  }else {
    skip = 1;
    xdouble = ((double)ypos(dataSubset, j-1, maxy[j], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) - (double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift))/(double)((double)ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) - ((double)ypos(dataSubset, j, maxy[j], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, maxy[j], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)));
    ydouble = (double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) + xdouble*((double)ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift));
    x = xdouble;
    y = ydouble;
    x += j-1;
    if(xUnitsSwitch != XUNIT_BINS) {
      x *= baseline/(float)fin.NrBins;
      x += dxshift;
    }
    if(x > 0 && x < 10000 && y > -stack_state[current_stack_pos-1].subint_end && y < 2*stack_state[current_stack_pos-1].subint_end) {
      if((i != stack_state[current_stack_pos-1].subint_end - ignorelastpulses || x <= xmax_oscil) && (i != stack_state[current_stack_pos-1].subint_start || x >= xmin_oscil)) {
        if(y < ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) || y < ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift))
   ppgdraw(x, y);
      }
    }
  }
       }
       if(y >= ypos(dataSubset, j, maxy[j], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift))
  maxy[j] = i;
     }
   }
 }
 ppgebuf();
      }
      ppgsci(1);
      ppgsch(box_font.characterheight);
      ppgslw(box_font.linewidth);
      sprintf(txt1, "bc");
      sprintf(txt2, "bc");
      if(noboxx) {
 txt1[0] = 0;
 sprintf(txt2, "b");
      }
      if(noboxy) {
 sprintf(txt1, "b");
 txt2[0] = 0;
      }
      pgplot_options_definition pgplot_options;
      pgplot_clear_options(&pgplot_options);
      pgplot_options.box.drawbox = 0;
      pgplot_options.box.drawtitle = 0;
      pgplot_options.box.drawlabels = 0;
      pgplot_options.box.title_ch = title_font.characterheight*currentPanelScaling;
      pgplot_options.box.title_lw = title_font.linewidth;
      pgplot_options.box.title_f = title_font.font;
      pgplot_options.box.label_f = title_font.font;
      pgplot_options.box.label_ch = label_font.characterheight*currentPanelScaling;
      pgplot_options.box.box_labelsize = box_font.characterheight*currentPanelScaling;
      pgplot_options.box.box_lw = box_font.linewidth;
      pgplot_options.box.box_f = box_font.font;
      pgplot_options.box.label_lw = box_font.linewidth;
      if(noboxx == 0 || noboxy == 0) {
 if(curpanelnrx == 0 || appendframes_flag == 0) {
   if(disable_x_numbers && disable_y_numbers ) {
     strcat(txt1, "st");
     strcat(txt2, "t");
   }else if(disable_x_numbers && !disable_y_numbers) {
     strcat(txt1, "st");
     strcat(txt2, "nt");
   }else if(!disable_x_numbers && !disable_y_numbers) {
     strcat(txt1, "nst");
     strcat(txt2, "nt");
   }else if(!disable_x_numbers && disable_y_numbers) {
     strcat(txt1, "nst");
     strcat(txt2, "t");
   }
 }else {
   if(disable_x_numbers) {
     strcat(txt1, "st");
     strcat(txt2, "t");
   }else {
     strcat(txt1, "nst");
     strcat(txt2, "t");
   }
 }
 pgplot_options.box.drawbox = 1;
 if((fin.yrangeset && yUnitsSwitch)
    || (vectorMinMaxSpecified && yUnitsSwitch)
    || (((fin.gentype == GENTYPE_SUBINTEGRATIONS || fin.gentype == GENTYPE_PULSESTACK) && fin.NrFreqChan == 1 && fin.NrSubints > 1) && yUnitsSwitch)
    || (yUnitsSwitch && didtranspose_orig_nrbin)
    ) {
   pgplot_frame.swin_showtwice = showTwice_flag;
   pgplot_frame.swin_x1 = xminshow;
   pgplot_frame.swin_x2 = xmaxshow;
   pgplot_frame.swin_y1 = ymin;
   pgplot_frame.swin_y2 = ymax;
   strcpy(pgplot_options.box.box_xopt, txt1);
   strcpy(pgplot_options.box.box_yopt, txt2);
   float vectorMin2, vectorMax2;
   yvec2unit(fin, yUnitsSwitch, vectorMinMaxSpecified, vectorMin, vectorMax, didtranspose_orig_nrbin, 1, ymin, &vectorMin2, stack_state[current_stack_pos-1].grayscalemode, 0, application.verbose_state);
   yvec2unit(fin, yUnitsSwitch, vectorMinMaxSpecified, vectorMin, vectorMax, didtranspose_orig_nrbin, 2, ymax, &vectorMax2, stack_state[current_stack_pos-1].grayscalemode, 0, application.verbose_state);
   if(showTwice_flag == 0)
     ppgswin(xminshow, xmaxshow, yUnitCmdLine*vectorMin2, yUnitCmdLine*vectorMax2);
   else
     ppgswin(xminshow, xmaxshow, yUnitCmdLine*vectorMin2, yUnitCmdLine*(vectorMax2+(vectorMax2-vectorMin2)));
   pgplot_drawbox(&pgplot_options.box);
   if(stack_state[current_stack_pos-1].grayscalemode == 0)
     ppgswin(xminshow, xmaxshow, ymin, ymax);
   else {
     ppgswin(xminshow, xmaxshow, ymin-0.5, ymax+0.5);
   }
 }else {
   pgplot_options.box.drawbox = 1;
   strcpy(pgplot_options.box.box_xopt, txt1);
   strcpy(pgplot_options.box.box_yopt, txt2);
   pgplot_drawbox(&pgplot_options.box);
 }
      }
      pgplot_options.box.drawbox = 0;
      ppgsch(label_font.characterheight);
      pgplot_options.box.label_f = label_font.font;
      if(curpanelnrx == 0 || appendframes_flag == 0) {
 if(xtitle_set == 0) {
   if(xUnitsSwitch != XUNIT_BINS) {
     if(xUnitsSwitch == XUNIT_TIME) {
       sprintf(xtitle, "Time (sec)");
     }else if(xUnitsSwitch == XUNIT_PHASE) {
       if(didtranspose_orig_nrbin == 0)
  sprintf(xtitle, "Pulse phase");
       else
  sprintf(xtitle, "Subint");
     }else if(xUnitsSwitch == XUNIT_DEG) {
       sprintf(xtitle, "Pulse longitude (deg)");
     }
     if((fin.gentype == GENTYPE_2DFS) && (longitudeRangeSet_flag == 0 && fin.xrangeset == 1)) {
       strcpy(xtitle, "fluctuation frequency (cycles/period)");
     }else if((fin.gentype == GENTYPE_HRFS_UNFOLDED || fin.gentype == GENTYPE_HRFS) && (longitudeRangeSet_flag == 0 && fin.xrangeset == 1)) {
       strcpy(xtitle, "fluctuation frequency (cycles/period)");
     }
   }else {
     if(fin.gentype == GENTYPE_PROFILE || fin.gentype == GENTYPE_PULSESTACK || fin.gentype == GENTYPE_SUBINTEGRATIONS || fin.gentype == GENTYPE_DYNAMICSPECTRUM || fin.gentype == GENTYPE_LRFS || fin.gentype == GENTYPE_P3FOLD || fin.gentype == GENTYPE_LRCC || fin.gentype == GENTYPE_PADIST || fin.gentype == GENTYPE_ELLDIST) {
       sprintf(xtitle, "Pulse longitude (bins)");
     }else if(fin.gentype == GENTYPE_S2DFSP3 || fin.gentype == GENTYPE_S2DFSP2) {
       sprintf(xtitle, "Block number (pulses)");
     }else {
       sprintf(xtitle, "Bin");
     }
   }
 }
 if(ytitleset == 0) {
   if(ytitle_set_to_intensity) {
     if(stack_state[current_stack_pos-1].grayscalemode != 1)
       strcpy(ytitle, "Intensity");
     else
       strcpy(ytitle, "Pulse number");
   }else if(ytitle_set_to_freqchannel && yUnitsSwitch == 0) {
     strcpy(ytitle, "Frequency channel");
   }else if(ytitle_set_to_pulsenumber && yUnitsSwitch == 0) {
     strcpy(ytitle, "Pulse number");
     if(stack_state[current_stack_pos-1].subint_start == stack_state[current_stack_pos-1].subint_end && stack_state[current_stack_pos-1].grayscalemode != 1) {
       strcpy(ytitle, "Intensity");
     }
   }else if(ytitle_set_to_subint && yUnitsSwitch == 0) {
     strcpy(ytitle, "Subint number");
     if(stack_state[current_stack_pos-1].subint_start == stack_state[current_stack_pos-1].subint_end && stack_state[current_stack_pos-1].grayscalemode != 1) {
       strcpy(ytitle, "Intensity");
     }
   }else if(ytitle_set_to_p3fold && yUnitsSwitch) {
     strcpy(ytitle, "Pulse number");
   }else if(ytitle_set_to_spectral_power) {
     strcpy(ytitle, "Spectral Power");
   }else if(ytitle_set_to_harmonic_number) {
     strcpy(ytitle, "Harmonic Number");
   }else if(ytitle_set_to_lag_number) {
     strcpy(ytitle, "Lag number");
   }else if((ytitle_set_to_pulsenumber || ytitle_set_to_subint) && yUnitsSwitch) {
     strcpy(ytitle, "Time (sec)");
   }else if(ytitle_set_to_fluctuationFreq && yUnitsSwitch) {
     strcpy(ytitle, "fluctuation frequency (cycles/period)");
   }else if(ytitle_set_to_fluctuationFreq && yUnitsSwitch == 0) {
     strcpy(ytitle, "fluctuation frequency (bin)");
   }else if(ytitle_set_to_padist && yUnitsSwitch) {
     strcpy(ytitle, "PA (deg)");
   }else if(ytitle_set_to_padist && yUnitsSwitch == 0) {
     strcpy(ytitle, "PA (bin)");
   }else if(ytitle_set_to_elldist && yUnitsSwitch) {
     strcpy(ytitle, "\\gx (deg)");
   }else if(ytitle_set_to_elldist && yUnitsSwitch == 0) {
     strcpy(ytitle, "\\gx (bin)");
   }else if(ytitle_set_to_pulselongitude && yUnitsSwitch) {
     strcpy(ytitle, "Pulse longitude (deg)");
   }else if(ytitle_set_to_pulselongitude && yUnitsSwitch == 0) {
     strcpy(ytitle, "Pulse longitude (bins)");
   }else if(yUnitsSwitch && didtranspose_orig_nrbin) {
     strcpy(ytitle, "Frequency (MHz)");
   }
 }
 pgplot_options.box.drawlabels = 1;
 strcpy(pgplot_options.box.xlabel, xtitle);
 strcpy(pgplot_options.box.ylabel, ytitle);
 pgplot_drawbox(&pgplot_options.box);
      }else {
 pgplot_options.box.drawlabels = 1;
 strcpy(pgplot_options.box.xlabel, xtitle);
 strcpy(pgplot_options.box.ylabel, "");
 pgplot_drawbox(&pgplot_options.box);
      }
      if(!stack_state[current_stack_pos-1].grayscalemode) {
 pgplot_options.box.drawlabels = 0;
 pgplot_options.box.drawtitle = 1;
 newtext = str_replace_header_params(fin, title, application.verbose_state);
 if(newtext == NULL) {
   printwarning(application.verbose_state.debug, "WARNING pplot: Cannot substitute keyword in title");
   pgplot_options.box.title[0] = 0;
 }else {
   strcpy(pgplot_options.box.title, newtext);
   free(newtext);
 }
 pgplot_drawbox(&pgplot_options.box);
 pgplot_options.box.drawtitle = 0;
      }
      redraw = 0;
      if(interactive_flag) {
 do {
   int key, key2;
   FILE *fout;
   if(data_read == 2)
     data_read = 1;
   if(data_read == 0) {
     key = '?';
     data_read = 1;
   }else {
     printf("Option: ");
     fflush(stdout);
     do {
       key = pgetch_macro(&application, application.verbose_state);
     }while(key == '\n' || key == '\r');
     printf("%c\n", key);
   }
   switch(key) {
   case '?':
     printf("a       Auto-scale, only works with line plot\n");
     printf("b       Toggle units of the x-axis between bins and other units\n");
     printf("I       Info (statistical), limited by current x-axis selection\n");
     printf("l       Specify left/right bounds, i.e. the x-range\n");
     printf("M       Change between colour map/line plot mode\n");
     printf("n       Step in vector range (vertical range)\n");
     printf("p       Pop stack\n");
     printf("P       Switch polarization channel\n");
     printf("q       Quit\n");
     printf("r       Redraw\n");
     printf("R       Read cursor position\n");
     printf("s       Set scaling of data (similar to -scale/-scale2 options).\n");
     printf("v       Set vector range (vertical range)\n");
     printf("y       Toggle units of the y-axis\n");
     if(zapmode == 0)
       printf("z       Toggle from zap vectors (vertical) to zap subints (horizontal) mode\n");
     else
       printf("z       Toggle from zap subints (horizontal) to zap vectors (vertical) mode\n");
     if(zapmode == 0)
       printf("Z       Zap vectors (vertical axis, change to subints with 'z'). This only affects the plot, not the input-data.\n");
     else
       printf("Z       Zap subints (horizontal axis, change to vectors with 'z'). This only affects the plot, not the input-data.\n");
     if(change_filename_extension(inputfilename, txt3, "zap", 999, application.verbose_state) == 1) {
       if(zapmode == 0)
  printf("W       Write out list of zapped vectors to %s (change to subints with 'z')\n", txt3);
       else
  printf("W       Write out list of zapped subints to %s (change to vectors with 'z')\n", txt3);
     }
     printf(",/. Move to the left/right (if zoomed in)\n");
     printf("?   This help\n");
     break;
   case 'a':
     if(stack_state[current_stack_pos-1].grayscalemode) {
       printf("Resetting scale to 1.\n");
       scale = 1;
       redraw = 1;
     }else {
       scale = 1;
       y = y2 = -1;
       float float_tmp;
       for(i = stack_state[current_stack_pos-1].subint_start; i <= stack_state[current_stack_pos-1].subint_end; i++) {
  for(j = stack_state[current_stack_pos-1].x1; j <= stack_state[current_stack_pos-1].x2; j++) {
    float_tmp = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
    if(stack_state[current_stack_pos-1].subint_start != stack_state[current_stack_pos-1].subint_end)
      float_tmp -= i;
    if(-float_tmp > y || (i == stack_state[current_stack_pos-1].subint_start && j == stack_state[current_stack_pos-1].x1)) {
      y = float_tmp;
    }
    if(float_tmp > y2 || (i == stack_state[current_stack_pos-1].subint_start && j == stack_state[current_stack_pos-1].x1)) {
      y2 = float_tmp;
    }
  }
       }
       if(y2 > 0) {
  scale = fabs(1.0/y2);
       }
       if(fabs(1.0/y) < scale && y > 0) {
  scale = fabs(1.0/y);
       }
       redraw = 1;
     }
     break;
   case 'M':
     printf("  m = colour map\n");
     printf("  l = line drawing\n");
     fflush(stdout);
     if(current_stack_pos >= max_nr_stack) {
       printwarning(application.verbose_state.debug, "WARNING: Stack is full");
       current_stack_pos--;
     }else {
       copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
     }
     do {
       key2 = pgetch_macro(&application, application.verbose_state);
     }while(key2 == '\n' || key2 == '\r');
     if(key2 == 'm') {
       printf("set colour map mode\n");
       stack_state[current_stack_pos].grayscalemode = 1;
       scale = 1;
     }else {
       printf("set line plot mode\n");
       stack_state[current_stack_pos].grayscalemode = 0;
     }
     current_stack_pos++;
     redraw = 1;
     break;
   case 'p':
     current_stack_pos--;
     if(current_stack_pos <= 0) {
       printwarning(application.verbose_state.debug, "WARNING: Stack is empty");
       current_stack_pos = 1;
     }else {
       if(stack_state[current_stack_pos].nrZappedVectors != stack_state[current_stack_pos-1].nrZappedVectors || stack_state[current_stack_pos].nrZappedSubints != stack_state[current_stack_pos-1].nrZappedSubints) {
  data_read = 2;
       }
     }
     redraw = 1;
     break;
   case 'P':
     current_polnr++;
     if(current_polnr >= nrpolarizations)
       current_polnr = 0;
     printf("Selected polarization channel %d\n", current_polnr);
     redraw = 1;
     data_read = 2;
     break;
   case 'q':
   case 3:
     interactive_flag = 0;
     redraw = 1;
     break;
   case 'b':
     xUnitsSwitch += 1;
     if(xUnitsSwitch == XUNIT_ENDOFLIST || fin.gentype == GENTYPE_S2DFSP3 || fin.gentype == GENTYPE_S2DFSP2)
       xUnitsSwitch = 0;
     redraw = 1;
     break;
   case 'y':
     yUnitsSwitch += 1;
     if(yUnitsSwitch == 2)
       yUnitsSwitch = 0;
     redraw = 1;
     break;
   case 'r':
     redraw = 1;
     break;
   case 's':
     if(stack_state[current_stack_pos-1].grayscalemode) {
       printf("Set scale (default=1, current value %f). A value > 1 results in clipping, making weak features clearer.\n", scale);
     }else {
       printf("Set scale (default=1, current value %f). All intensities are multiplied with this value.\n", scale);
     }
     printf("Scale = ");
     fflush(stdout);
     if(application.macro_ptr == NULL) {
       scanf("%f", &dummyf1);
     }else {
       int ret;
       ret = fscanf(application.macro_ptr, "%f", &dummyf1);
       if(ret != 1) {
  fclose(application.macro_ptr);
  application.macro_ptr = NULL;
  printf("\nReached end of macro, switching to keyboard input: type in two numbers\n");
  scanf("%f", &dummyf1);
       }else {
  printf("%f\n", dummyf1);
       }
     }
     scale = dummyf1;
     if(stack_state[current_stack_pos-1].grayscalemode) {
       printf("Set second scale (default=0, current value %f). A value closer to 1 results in clipping of low intesity values to emphasize the bright features.\n", scale2);
       printf("Scale = ");
       fflush(stdout);
       if(application.macro_ptr == NULL) {
  scanf("%f", &dummyf1);
       }else {
  int ret;
  ret = fscanf(application.macro_ptr, "%f", &dummyf1);
  if(ret != 1) {
    fclose(application.macro_ptr);
    application.macro_ptr = NULL;
    printf("\nReached end of macro, switching to keyboard input: type in two numbers\n");
    scanf("%f", &dummyf1);
  }else {
    printf("%f\n", dummyf1);
  }
       }
       scale2 = dummyf1;
     }
     redraw = 1;
     break;
   case 'R':
     printf("Left click to read positions, press other key to quit.\n");
     do {
       int nx, ny;
       ppgband(7, 0, 0, 0, &x, &y, &ch);
       if(ch == 65) {
  pgplotMapCoordinate(x, y, &nx, &ny);
  yvec2unit(fin, yUnitsSwitch, vectorMinMaxSpecified, vectorMin, vectorMax, didtranspose_orig_nrbin, 0, y, &dummyf1, stack_state[current_stack_pos-1].grayscalemode, 0, application.verbose_state);
  printf("%f %f (%d %d)\n", x, dummyf1, nx, ny);
       }
     }while(ch == 65);
     break;
   case 'z':
     if(zapmode == 0) {
       if(didtranspose_orig_nrbin) {
  printf("Swapping to zapping of subints\n");
  zapmode = 1;
       }else {
  printf("Swapping to zapping of subints only works if there is more than one subint and more than one frequency channel\n");
       }
     }else {
       printf("Swapping to zapping of vectors (vertical axis)\n");
       zapmode = 0;
     }
     break;
   case 'W':
     if(zapmode == 0) {
       if(stack_state[current_stack_pos-1].nrZappedVectors == 0) {
  printwarning(application.verbose_state.debug, "WARNING pplot: No vectors are zapped, so writing out of a zapfile is ignored.");
       }else {
  int ret;
  if(fin.NrFreqChan == 1 && didtranspose_orig_nrbin == 0) {
    ret = change_filename_extension(inputfilename, txt3, "subint.zap", 999, application.verbose_state);
  }else {
    ret = change_filename_extension(inputfilename, txt3, "freq.zap", 999, application.verbose_state);
  }
  if(ret == 1) {
    fout = fopen(txt3, "w");
    if(fout != NULL) {
      long largest, count;
      count = 0;
      largest = zappedVectors[0];
      for(i = 0; i < stack_state[current_stack_pos-1].nrZappedVectors; i++) {
        if(zappedVectors[i] > largest) {
   largest = zappedVectors[i];
        }
      }
      for(j = 0; j <= largest; j++) {
        for(i = 0; i < stack_state[current_stack_pos-1].nrZappedVectors; i++) {
   if(zappedVectors[i] == j) {
     fprintf(fout, "%d\n", zappedVectors[i]);
     count++;
     break;
   }
        }
      }
      printf("pplot: Written %ld vector numbers to %s\n", count, txt3);
      fclose(fout);
    }else {
      printerror(application.verbose_state.debug, "ERROR pplot: Cannot open zapfile %s.", txt3);
      return 0;
    }
  }else {
    printerror(application.verbose_state.debug, "ERROR pplot: Filename is too long.");
    return 0;
  }
       }
     }else {
       if(stack_state[current_stack_pos-1].nrZappedSubints == 0) {
  printwarning(application.verbose_state.debug, "WARNING pplot: No subints are zapped, so writing out of a zapfile is ignored.");
       }else {
  if(change_filename_extension(inputfilename, txt3, "subint.zap", 999, application.verbose_state) == 1) {
    fout = fopen(txt3, "w");
    if(fout != NULL) {
      long largest, count;
      count = 0;
      largest = zappedSubints[0];
      for(i = 0; i < stack_state[current_stack_pos-1].nrZappedSubints; i++) {
        if(zappedSubints[i] > largest) {
   largest = zappedSubints[i];
        }
      }
      for(j = 0; j <= largest; j++) {
        for(i = 0; i < stack_state[current_stack_pos-1].nrZappedSubints; i++) {
   if(zappedSubints[i] == j) {
     fprintf(fout, "%d\n", zappedSubints[i]);
     count++;
     break;
   }
        }
      }
      printf("pplot: Written %ld subint numbers to %s\n", count, txt3);
      fclose(fout);
    }else {
      printerror(application.verbose_state.debug, "ERROR pplot: Cannot open zapfile %s.", txt3);
      return 0;
    }
  }else {
    printerror(application.verbose_state.debug, "ERROR pplot: Filename is too long.");
    return 0;
  }
       }
     }
     break;
   case 'Z':
     if(current_stack_pos >= max_nr_stack) {
       printwarning(application.verbose_state.debug, "WARNING: Stack is full");
       current_stack_pos--;
     }else {
       copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
     }
     if(zapmode == 0)
       printf("Left click to zap vector (vertical axis), left click followed by a right click removes a range. Press other key to quit. The result only visible after you're done with clicking.\n");
     else
       printf("Left click to zap subint (horizontal axis), left click followed by a right click removes a range. Press other key to quit. The result only visible after you're done with clicking.\n");
     do {
       int nx, ny, nx_final, ny_final, ret, ch_value;
       if(application.macro_ptr == NULL) {
  if(zapmode == 0)
    ppgband(5, 0, 0, 0, &x, &y, &ch);
  else
    ppgband(6, 0, 0, 0, &x, &y, &ch);
       }else {
  if(application.verbose_state.verbose) {
    printf("In macro mode, add lines with \"keyvalue xvalue yvalue\"\n");
    printf("keyvalue = 65 = Left click\n");
    printf("keyvalue = 88 = Right click\n");
  }
  ret = fscanf(application.macro_ptr, "%d %f %f", &ch_value, &x, &y);
  if(ret != 3) {
    printf("\nReached end of macro, switching to keyboard input\n");
    fclose(application.macro_ptr);
    application.macro_ptr = NULL;
    ch_value = 27;
  }
  ch = ch_value;
       }
       if(ch == 65 || ch == 88) {
  pgplotMapCoordinate(x, y, &nx, &ny);
  if(zapmode == 0) {
    if(ch == 65) {
      printf("Zapping vector %d\n", ny);
      ny_last_clicked_zap = ny;
      ny_final = ny;
    }else {
      if(ny < ny_last_clicked_zap) {
        int oldint;
        oldint = ny;
        ny = ny_last_clicked_zap;
        ny_last_clicked_zap = oldint;
      }
      printf("Zapping vectors %d - %d\n", ny_last_clicked_zap, ny);
      ny_final = ny;
    }
    for(ny = ny_last_clicked_zap; ny <= ny_final; ny++) {
      if(stack_state[current_stack_pos].nrZappedVectors < max_nr_zap-1) {
        zappedVectors[stack_state[current_stack_pos].nrZappedVectors++] = ny;
      }else {
        printwarning(application.verbose_state.debug, "WARNING: Maximum number of interactively zapped vectors exceeded.");
      }
      if(zapVectors(stack_state[current_stack_pos].nrZappedVectors, zappedVectors, subint_start, fin, dataSubset, ny, application.verbose_state) == 0)
        return 0;
    }
  }else {
    nx /= didtranspose_orig_nrbin;
    if(ch == 65) {
      printf("Zapping subint %d\n", nx);
      nx_last_clicked_zap = nx;
      nx_final = nx;
    }else {
      if(nx < nx_last_clicked_zap) {
        int oldint;
        oldint = nx;
        nx = nx_last_clicked_zap;
        nx_last_clicked_zap = oldint;
      }
      printf("Zapping subints %d - %d\n", nx_last_clicked_zap, nx);
      nx_final = nx;
    }
    for(nx = nx_last_clicked_zap; nx <= nx_final; nx++) {
      if(stack_state[current_stack_pos].nrZappedSubints < max_nr_zap-1) {
        zappedSubints[stack_state[current_stack_pos].nrZappedSubints++] = nx;
      }else {
        printwarning(application.verbose_state.debug, "WARNING: Maximum number of interactively zapped subints exceeded.");
      }
      if(zapSubints(stack_state[current_stack_pos].nrZappedSubints, zappedSubints, subint_start, didtranspose_orig_nrbin, fin, dataSubset, nx, application.verbose_state) == 0)
        return 0;
    }
  }
  nx_final = nx;
  ny_final = ny;
       }else {
  printf("Got key: %d\n", ch);
       }
     }while(ch == 65);
     redraw = 1;
     current_stack_pos++;
     break;
   case 'v':
     printf("Specify vector (vertical) range (two numbers): ");
     fflush(stdout);
     if(current_stack_pos >= max_nr_stack) {
       printwarning(application.verbose_state.debug, "WARNING: Stack is full");
       current_stack_pos--;
     }else {
       copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
     }
     if(application.macro_ptr == NULL) {
       scanf("%f %f", &dummyf1, &dummyf2);
     }else {
       int ret;
       ret = fscanf(application.macro_ptr, "%f %f", &dummyf1, &dummyf2);
       if(ret != 2) {
  printf("\nReached end of macro, switching to keyboard input: type in two numbers\n");
  scanf("%f %f", &dummyf1, &dummyf2);
  fclose(application.macro_ptr);
  application.macro_ptr = NULL;
       }else {
  printf("%f %f\n", dummyf1, dummyf2);
       }
     }
     if((fin.yrangeset && yUnitsSwitch)
        || (vectorMinMaxSpecified && yUnitsSwitch)
        || (((fin.gentype == GENTYPE_SUBINTEGRATIONS || fin.gentype == GENTYPE_PULSESTACK) && fin.NrFreqChan == 1 && fin.NrSubints > 1) && yUnitsSwitch)
        || (yUnitsSwitch && didtranspose_orig_nrbin)
        ) {
       float rangemin, rangemax;
       yvec2unit(fin, yUnitsSwitch, vectorMinMaxSpecified, vectorMin, vectorMax, didtranspose_orig_nrbin, 0, dummyf1, &rangemin, stack_state[current_stack_pos-1].grayscalemode, 1, application.verbose_state);
       yvec2unit(fin, yUnitsSwitch, vectorMinMaxSpecified, vectorMin, vectorMax, didtranspose_orig_nrbin, 0, dummyf2, &rangemax, stack_state[current_stack_pos-1].grayscalemode, 1, application.verbose_state);
       stack_state[current_stack_pos].subint_start = round(rangemin);
       stack_state[current_stack_pos].subint_end = round(rangemax);
     }else {
       stack_state[current_stack_pos].subint_start = dummyf1/yUnitCmdLine;
       stack_state[current_stack_pos].subint_end = dummyf2/yUnitCmdLine;
     }
     if(stack_state[current_stack_pos].subint_start < subint_start) {
       printwarning(application.verbose_state.debug, "WARNING: set start vector to %ld", subint_start);
       stack_state[current_stack_pos].subint_start = subint_start;
     }
     if(stack_state[current_stack_pos].subint_start > subint_end) {
       printwarning(application.verbose_state.debug, "WARNING: set end vector to %ld", subint_end);
       stack_state[current_stack_pos].subint_start = subint_end;
     }
     if(stack_state[current_stack_pos].subint_end < subint_start) {
       printwarning(application.verbose_state.debug, "WARNING: set end vector to %ld", subint_start);
       stack_state[current_stack_pos].subint_end = subint_start;
     }
     if(stack_state[current_stack_pos].subint_end > subint_end) {
       printwarning(application.verbose_state.debug, "WARNING: set end vector to %ld", subint_end);
       stack_state[current_stack_pos].subint_end = subint_end;
     }
     current_stack_pos++;
     redraw = 1;
     break;
   case 'l':
     printf("Specify x-range (two numbers, or R to use mouse): ");
     fflush(stdout);
     if(current_stack_pos >= max_nr_stack) {
       printwarning(application.verbose_state.debug, "WARNING: Stack is full");
       current_stack_pos--;
     }else {
       copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
     }
     if(application.macro_ptr == NULL) {
       fgets(txt3, 990, stdin);
       i = sscanf(txt3, "%s %s", txt1, txt2);
       if(i == 1 && strcmp(txt1, "R") != 0) {
  fgets(txt3, 990, stdin);
  sscanf(txt3, "%s", txt2);
       }
       if(strcmp(txt1, "R") == 0) {
  printf("Click on left edge\n");
  ppgband(6, 0, 0, 0, &x, &y, &ch);
       }else {
  sscanf(txt1, "%f", &x);
       }
       if(strcmp(txt1, "R") == 0) {
  printf("Click on right edge\n");
  x2 = x;
  ppgband(4, 0, x, y, &x2, &y, &ch);
       }else {
  sscanf(txt2, "%f", &x2);
       }
     }else {
       int ret;
       ret = fscanf(application.macro_ptr, "%f %f", &x, &x2);
       if(ret != 2) {
  printf("\nReached end of macro, switching to keyboard input: type in two numbers\n");
  scanf("%f %f", &dummyf1, &dummyf2);
  fclose(application.macro_ptr);
  application.macro_ptr = NULL;
       }
     }
     printf("Select %f-%f\n", x, x2);
     if(xUnitsSwitch != XUNIT_BINS) {
       x -= dxshift;
       x /= baseline/(float)fin.NrBins;
       x2 -= dxshift;
       x2 /= baseline/(float)fin.NrBins;
     }
     stack_state[current_stack_pos].x1 = x;
     stack_state[current_stack_pos].x2 = x2;
     if(stack_state[current_stack_pos].x1 < 0) {
       printwarning(application.verbose_state.debug, "WARNING: set x1 to 0");
       stack_state[current_stack_pos].x1 = 0;
     }
     if(stack_state[current_stack_pos].x1 >= fin.NrBins) {
       printwarning(application.verbose_state.debug, "WARNING: set x1 to %ld", fin.NrBins-1);
       stack_state[current_stack_pos].x1 = fin.NrBins-1;
     }
     if(stack_state[current_stack_pos].x2 < 0) {
       printwarning(application.verbose_state.debug, "WARNING: set x2 to 0");
       stack_state[current_stack_pos].x2 = 0;
     }
     if(stack_state[current_stack_pos].x2 >= fin.NrBins) {
       printwarning(application.verbose_state.debug, "WARNING: set x2 to %ld", fin.NrBins-1);
       stack_state[current_stack_pos].x2 = fin.NrBins-1;
     }
     current_stack_pos++;
     redraw = 1;
     break;
   case '.':
     if(current_stack_pos >= max_nr_stack) {
       printwarning(application.verbose_state.debug, "WARNING: Stack is full");
       current_stack_pos--;
     }else {
       copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
     }
     y = stack_state[current_stack_pos].x2 - stack_state[current_stack_pos].x1;
     x = stack_state[current_stack_pos].x1 + y;
     x2 = stack_state[current_stack_pos].x2 + y;
     if(x < 0) {
       printwarning(application.verbose_state.debug, "WARNING: set x1 to 0");
       x = 0;
       x2 = y;
     }
     if(x2 >= fin.NrBins) {
       printwarning(application.verbose_state.debug, "WARNING: set x2 to %ld", fin.NrBins-1);
       x2 = fin.NrBins-1;
       x = x2 - y;
       if(x < 0)
  x = 0;
     }
     stack_state[current_stack_pos].x1 = x;
     stack_state[current_stack_pos].x2 = x2;
     current_stack_pos++;
     redraw = 1;
     break;
   case ',':
     if(current_stack_pos >= max_nr_stack) {
       printwarning(application.verbose_state.debug, "WARNING: Stack is full");
       current_stack_pos--;
     }else {
       copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
     }
     y = stack_state[current_stack_pos].x2 - stack_state[current_stack_pos].x1;
     x = stack_state[current_stack_pos].x1 - y;
     x2 = stack_state[current_stack_pos].x2 - y;
     if(x < 0) {
       printwarning(application.verbose_state.debug, "WARNING: set x1 to 0");
       x = 0;
       x2 = y;
     }
     if(x2 >= fin.NrBins) {
       printwarning(application.verbose_state.debug, "WARNING: set x2 to %ld", fin.NrBins-1);
       x2 = fin.NrBins-1;
       x = x2 - y;
       if(x < 0)
  x = 0;
     }
     stack_state[current_stack_pos].x1 = x;
     stack_state[current_stack_pos].x2 = x2;
     current_stack_pos++;
     redraw = 1;
     break;
   case 'n':
     i = stack_state[current_stack_pos-1].subint_end-stack_state[current_stack_pos-1].subint_start+1;
     j = stack_state[current_stack_pos-1].subint_start;
     k = stack_state[current_stack_pos-1].subint_end;
     if(current_stack_pos >= max_nr_stack) {
       printwarning(application.verbose_state.debug, "WARNING: Stack is full");
       current_stack_pos--;
     }else {
       copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
     }
     stack_state[current_stack_pos].subint_start = j+i;
     stack_state[current_stack_pos].subint_end = k+i;
     if(stack_state[current_stack_pos].subint_start < subint_start) {
       printwarning(application.verbose_state.debug, "WARNING: set vector to %ld", subint_start);
       stack_state[current_stack_pos].subint_start = subint_start;
     }
     if(stack_state[current_stack_pos].subint_start > subint_end) {
       printwarning(application.verbose_state.debug, "WARNING: set vector to %ld", subint_end);
       stack_state[current_stack_pos].subint_start = subint_end;
     }
     if(stack_state[current_stack_pos].subint_end > subint_end) {
       printwarning(application.verbose_state.debug, "WARNING: set vector to %ld", subint_end);
       stack_state[current_stack_pos].subint_end = subint_end;
     }
     current_stack_pos++;
     redraw = 1;
     break;
   case 'I':
     printf("Vector -   min            at binnr  max            at binnr  mean           rms            centroid bin\n");
     for(i = stack_state[current_stack_pos-1].subint_start; i <= stack_state[current_stack_pos-1].subint_end; i++) {
       double Imin, Imax, mean, rms, x_centroid;
       long bin_min, bin_max;
       mean = 0;
       rms = 0;
       Imin = Imax = 0;
       bin_min = bin_max = 0;
       x_centroid = 0;
       for(j = stack_state[current_stack_pos-1].x1; j <= stack_state[current_stack_pos-1].x2; j++) {
  double float_tmp;
  float_tmp = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
  if(stack_state[current_stack_pos-1].subint_start != stack_state[current_stack_pos-1].subint_end)
    float_tmp -= i;
  mean += float_tmp;
  x_centroid += j*float_tmp;
  if(j == stack_state[current_stack_pos-1].x1 || float_tmp > Imax) {
    Imax = float_tmp;
    bin_max = j;
  }
  if(j == stack_state[current_stack_pos-1].x1 || float_tmp < Imin) {
    Imin = float_tmp;
    bin_min = j;
  }
       }
       x_centroid /= (double)mean;
       mean /= (double)(stack_state[current_stack_pos-1].x2-stack_state[current_stack_pos-1].x1+1);
       for(j = stack_state[current_stack_pos-1].x1; j <= stack_state[current_stack_pos-1].x2; j++) {
  double float_tmp;
  float_tmp = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
  if(stack_state[current_stack_pos-1].subint_start != stack_state[current_stack_pos-1].subint_end)
    float_tmp -= i;
  rms += (float_tmp-mean)*(float_tmp-mean);
       }
       rms /= (float)(stack_state[current_stack_pos-1].x2-stack_state[current_stack_pos-1].x1+1);
       rms = sqrt(rms);
       printf("%6ld - %14e %9ld %14e %9ld %14e %14e %14e\n", i, Imin, bin_min, Imax, bin_max, mean, rms, x_centroid);
     }
     break;
   default:
     printf("Unknown option: '%c' (%d)\nPress '?' for help.\n", key, key);
   }
 }while(redraw == 0);
      }
      curpanelnrx++;
      if(curpanelnrx == nrpanelsx) {
 curpanelnrx = 0;
 curpanelnry++;
      }
      if(curpanelnry == nrpanelsy)
 curpanelnry = 0;
    }while(interactive_flag);
    if(dataSubset_allocated) {
      free(dataSubset);
      dataSubset_allocated = 0;
    }
    if(stack_poly_allocated) {
      free(xstack_poly);
      free(ystack_poly);
      stack_poly_allocated = 0;
    }
    if(maxy_allocated) {
      free(maxy);
      maxy_allocated = 0;
    }
    closePSRData(&fin, 0, application.verbose_state);
    printf("Finished plotting of: %s\n", inputfilename);
  }
  ppgend();
  free(zappedVectors);
  free(zappedSubints);
  free(stack_state);
  terminateApplication(&application);
  return 0;
}
float ypos(float *stackI, int bin, long PulseNr, float scale, int NrBins, long pulse_bottom, long pulse_top, long subint_start, float yUnitCmdLine, float dyshift)
{
  long i;
  i = (PulseNr-subint_start)*NrBins+bin;
  if(pulse_bottom != pulse_top)
    return stackI[i]*scale+PulseNr*yUnitCmdLine + dyshift;
  else
    return stackI[i]*scale + dyshift;
}
int setBaselineParams(datafile_definition fin, float *baseline, float *dxshift, int *xUnitsSwitch, verbose_definition verbose)
{
  double period;
  int ret;
  if(fin.isFolded) {
    ret = get_period(fin, 0, &period, verbose);
    if(ret == 2) {
      printerror(verbose.debug, "ERROR pplot (%s): Cannot obtain period", fin.filename);
      return 0;
    }
  }else {
    ret = 1;
    period = -1;
  }
  if(fin.xrangeset) {
    *baseline = (fin.xrange[1]-fin.xrange[0])*(fin.NrBins)/(float)(fin.NrBins-0.999);
    *dxshift = fin.xrange[0];
    *xUnitsSwitch = XUNIT_DEG;
  }else {
    if(fin.isFolded != 0 || fin.gentype != GENTYPE_SEARCHMODE) {
      if(period < 0.001) {
 fflush(stdout);
 printwarning(verbose.debug, "pplot: The period does not appear to be set in the header. Consider using the -header option.");
 if(xUnitsSwitch != XUNIT_BINS) {
   printerror(verbose.debug, "       Terminating.");
   return 0;
 }else {
   printwarning(verbose.debug, "       (warning only)");
 }
      }
    }
    if((get_tsamp(fin, 0, verbose) < 0.0000001 || get_tsamp(fin, 0, verbose) >= 100)) {
      fflush(stdout);
      printwarning(verbose.debug, "pplot: The sampling time does not appear to be set correctly in the header. Consider using the -header option");
      if(xUnitsSwitch != XUNIT_BINS) {
 printerror(verbose.debug, "       Terminating.");
 return 0;
      }else {
 printwarning(verbose.debug, "       (warning only)");
      }
    }
    if(fin.isFolded != 0 || fin.gentype != GENTYPE_SEARCHMODE)
      *baseline = 360.0*get_tsamp(fin, 0, verbose)*fin.NrBins/period;
    else
      *baseline = 360.0;
    if(*xUnitsSwitch == XUNIT_PHASE) {
      if(fin.isFolded != 0 || fin.gentype != GENTYPE_SEARCHMODE)
 *baseline = get_tsamp(fin, 0, verbose)*fin.NrBins/period;
      else
 *baseline = 1.0;
    }else if(*xUnitsSwitch == XUNIT_TIME) {
      *baseline = get_tsamp(fin, 0, verbose)*fin.NrBins;
    }
    if(verbose.verbose && xUnitsSwitch != XUNIT_BINS) {
      printf("Based on sampling time and pulse period the baseline appears to be %f.\n", *baseline);
    }
  }
  return 1;
}
int zapVectors(int nrZappedVectors, int *zappedVectors, long subint_start, datafile_definition fin, float *dataSubset, int onlysinglechannel, verbose_definition verbose)
{
  int i, nx, ny;
  float I;
  if(onlysinglechannel >= 0)
    nrZappedVectors = 1;
  for(i = 0; i < nrZappedVectors; i++) {
    I = 0;
    if(onlysinglechannel < 0)
      ny = zappedVectors[i];
    else
      ny = onlysinglechannel;
    for(nx = 0; nx < fin.NrBins; nx++) {
      dataSubset[(ny-subint_start)*fin.NrBins + nx] = 0;
      if(writePulsePSRData(&fin, ny, 0, 0, nx, 1, &I, verbose) != 1) {
 printerror(verbose.debug, "ERROR pplot: Error writing data while zapping");
 return 0;
      }
    }
  }
  return 1;
}
int zapSubints(int nrZappedSubints, int *zappedSubints, long subint_start, int didtranspose_orig_nrbin, datafile_definition fin, float *dataSubset, int onlysinglesubint, verbose_definition verbose)
{
  float I;
  long binnr;
  int i, nx, ny;
  if(didtranspose_orig_nrbin) {
    if(onlysinglesubint >= 0)
      nrZappedSubints = 1;
    for(i = 0; i < nrZappedSubints; i++) {
      I = 0;
      if(onlysinglesubint < 0)
 nx = zappedSubints[i];
      else
 nx = onlysinglesubint;
      for(binnr = nx*didtranspose_orig_nrbin; binnr < (nx+1)*didtranspose_orig_nrbin; binnr++) {
 for(ny = 0; ny < fin.NrSubints; ny++) {
   dataSubset[(ny-subint_start)*fin.NrBins + binnr] = 0;
   if(writePulsePSRData(&fin, ny, 0, 0, binnr, 1, &I, verbose) != 1) {
     printerror(verbose.debug, "ERROR pplot: Error writing data while zapping");
     return 0;
   }
 }
      }
    }
  }
  return 1;
}
void yvec2unit(datafile_definition fin, int yUnitsSwitch, int vectorMinMaxSpecified, float vectorMin, float vectorMax, int didtranspose_orig_nrbin, int type, float valuein, float *valueout, int mapmode, int inverse, verbose_definition verbose)
{
  if(vectorMinMaxSpecified == 0) {
    if(fin.yrangeset && yUnitsSwitch) {
      vectorMin = fin.yrange[0];
      vectorMax = fin.yrange[1];
      vectorMinMaxSpecified = 1;
    }else if(((fin.gentype == GENTYPE_SUBINTEGRATIONS || fin.gentype == GENTYPE_PULSESTACK) && fin.NrFreqChan == 1 && fin.NrSubints > 1) && yUnitsSwitch && didtranspose_orig_nrbin == 0) {
      vectorMin = 0;
      vectorMax = get_tobs(fin, verbose);
      vectorMinMaxSpecified = 1;
    }else if(yUnitsSwitch && didtranspose_orig_nrbin) {
      if(fin.freqMode != FREQMODE_UNIFORM) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING pplot: Frequency range is for first subint.");
      }
      vectorMin = get_nonweighted_channel_freq(fin, 0, verbose);
      vectorMax = get_nonweighted_channel_freq(fin, fin.NrSubints-1, verbose);
      vectorMinMaxSpecified = 1;
    }
  }
  if(vectorMinMaxSpecified) {
    float binsize;
    binsize = (vectorMax - vectorMin)/(float)(fin.NrSubints-1);
    if(inverse == 0) {
      *valueout = valuein*binsize + vectorMin;
      if(mapmode) {
 if(type == 1)
   *valueout -= 0.5*binsize;
 else if(type == 2)
   *valueout += 0.5*binsize;
      }
    }else {
      if(mapmode) {
 if(type == 1)
   valuein += 0.5*binsize;
 else if(type == 2)
   valuein -= 0.5*binsize;
      }
      *valueout = (valuein - vectorMin)/binsize;
    }
  }else {
    *valueout = valuein;
  }
}

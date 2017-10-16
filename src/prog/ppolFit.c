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
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "psrsalsa.h"
#include "cpgplot.h"

#define M_PI 3.14159265358979323846
#define MaxNrJumps 100


void PrintHelp();
void convertAlphaBeta(double *alpha, double *beta
        );
void calcBeamWidths(int nalpha, int nbeta, double alphastart, double alphaend, double betastart, double betaend,
      int calculate_beam_widths, int calculate_interpulse_widths, double pulse_width2, float *rhogrid, float *rhogrid2, int nocounters);

int internal_fit_pa_or_l0;

void PlotGrid(float *chigrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, double level, double suppress_fac, int GridDeviceID, double alpha, double beta, double lwbox, double labelcharheight, double boxlabelcharheight, int drawCross, int draw_title, int drawcontours, int nogray,
       double chimax, double chimin, int maptype, int showwedge, char *showwedge_label);

void DoFitting(double alpha0, double beta0, double pa0, double dpa0, double l0, double dl0, double dh0, double ddh0, double ftol, double *fit_pa0, double *fit_l0, double *fit_a, double *fit_b, double *fit_dh, double *chi, int *nfunk, int searchAll, int report, FILE *reportStream, int finderrors, double nrofsigmas, int *nfitparameters, int amoeba_algorithm, verbose_definition verbose);

void PlotPAswing(double alpha, double beta, double pa0, double l0, int PlotFit, double leftPulseLongitude, double rightPulseLongitude, double dh);

void PlotContours(float *rhogrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, int nrlevels, float *TR, int GridDeviceID, int contour_txt, int contourcolor, int fixedContours, int nruserContours, float *userContours, int lwbox, int dotted);

void calcIntersectionRhoAndBanana(float *rhogrid, float *chigrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, int nruserContours, float *userContours, double chimax, double chimin, verbose_definition verbose);

void print_steepness(double alpha, double beta, double l0, double pa0, int verbose, double dh, double *sina_b);

double funk(double x[]);

double internal_funk_gsl(double *x, void *params);

double dy_180(double y1, double y2);
double dy_90(double y1, double y2);




struct {

  double *data_l;
  float *data_pa, *data_dpa;
  int NrDataPoints;

  double jump_longitude[MaxNrJumps], jump_offset[MaxNrJumps];
  int nrJumps, autojump;

  double add_height_longitude;


  double l0_start, max_l0_diff;

  int force_set;
  double force_l, force_pa, force_dpa;
  double pulse_width, sigma_width, rho_bcw, sigma_rho;
}fitterinfo;


int main(int argc, char **argv)
{
  char dumpfile[1000], c, device1[100], device2[100], txt[MaxStringLength], *txtptr, *showwedge_label;
  char prefix[1000];
  int i, j, nfunk, loadresults, GridDeviceID, PADeviceID, PSDeviceID, macrofilename, ret;
  int calculate_beam_widths, nrcontourlevels, nrcontourlevels2, redraw;
  int fixedContours, printsmooth, nogray, drawcontours, boxlw, nrfitparams;
  int contour_plot, contour_txt, calculate_interpulse_widths, iformat;
  int alphaset, betaset, lset, paset, printnofit, printfit, GridSearch, nalpha, nbeta;
  int onlyshowbest, showgraphics, drawCross, title_txt;
  int enableerrors, nruserContours, doBruteForce;
  int showwedge, version;
  int gridDevice_resx, gridDevice_resy, paDevice_resx, paDevice_resy;
  int invertGrayscale, amoeba_algorithm, doMC, devicenores, fixseed;
  int suppress_greyscale, beamwidth_params_only_w;
  float *chigrid, *l0grid, *pa0grid, *dhgrid, *rhogrid, *rhogrid2, TR[6], userContours[500];
  double ftol;
  double l0, dl0, pa0, dpa0, dh0, ddh0, alpha0, beta0, chi, chi_old, chimax, chimin;
  double fit_pa0, fit_l0, fit_dh0, fit_alpha, fit_beta, bestalpha, bestbeta;
  double labelcharheight, boxlabelcharheight, suppress_fac;
  double alphastart, alphaend, betastart, betaend, add_longitude_shift;
  double leftPulseLongitude, rightPulseLongitude, nrofsigmas;
  double level, pulse_width2, bestl0, bestpa0, beshdh;
  double paErrorFac;
  FILE *fin, *macrofile, *beamwidth_params_fin;
  double sina_b;
  int contourcolor, contourcolorIP;
  int doContourRange, contcol_nr_edge_steps, contcol_nr_fid_plane_steps;
  double contcol_edge1, contcol_edge1plus, contcol_edge1min, contcol_edge2, contcol_edge2plus, contcol_edge2min, contcol_phi0, contcol_phi0plus, contcol_phi0min, contcol_l0, contcol_l0plus, contcol_l0min, contcol_p;
  datafile_definition datain;
  int usegsl_error_grid_search;
  psrsalsaApplication application;
  verbose_definition beamwidth_params_vebose;

  initApplication(&application, "ppolFit", "[options] paswing_file");

  application.switch_verbose = 1;
  application.switch_nocounters = 1;
  application.switch_debug = 1;
  application.switch_history_cmd_only = 1;


  suppress_greyscale = 6;
  fitterinfo.pulse_width = 0;
  nogray = 0;
  ftol = 1e-3;
  alphaset = betaset = lset = paset = printnofit = printfit = GridSearch = loadresults = 0;
  alphastart = 0;
  alphaend = 180;
  betastart = -90;
  betaend = 90;
  fitterinfo.max_l0_diff = 45;
  calculate_beam_widths = 0;
  calculate_interpulse_widths = 0;
  strcpy(dumpfile, "dump.dat");
  contour_txt = 1;
  nrcontourlevels = 19;
  nrcontourlevels2 = 10;
  fixedContours = 0;
  printsmooth = 0;
  onlyshowbest = 0;
  showgraphics = 1;
  iformat = -1;
  fitterinfo.nrJumps = 0;
  fitterinfo.autojump = 0;
  add_longitude_shift = 0;
  fitterinfo.add_height_longitude = 1e9;
  drawCross = 1;
  fitterinfo.force_set = 0;
  title_txt = 1;
  boxlw = 4;
  drawcontours = 3;
  labelcharheight = 1.7;
  boxlabelcharheight = 1;
  suppress_fac = 1;
  dh0 = -1e9;
  ddh0 = -1e9;
  enableerrors = 0;
  contour_plot = 0;
  doBruteForce = 0;
  sprintf(device1, "?");
  sprintf(device2, "?");
  nruserContours = 0;
  showwedge_label = "Reduced-\\gx\\u2";
  macrofilename = 0;
  beamwidth_params_fin = NULL;
  showwedge = 0;
  nrofsigmas = 3;
  contourcolor = 2;
  contourcolorIP = 3;
  gridDevice_resx = 680;
  gridDevice_resy = 680;
  paDevice_resx = 680;
  paDevice_resy = 340;
  amoeba_algorithm = 0;

  usegsl_error_grid_search = 1;
  devicenores = 0;
  fixseed = 0;
  cleanVerboseState(&beamwidth_params_vebose);

  fin = NULL;
  rhogrid = rhogrid2 = NULL;
  bestalpha = bestbeta = -1e10;
  chimin = chimax = -1e10;
  PADeviceID = GridDeviceID = 0;
  doMC = 0;
  paErrorFac = 1;
  sprintf(prefix, "chi2");
  doContourRange = 0;

  if(argc < 2) {
    printf("Program to fit the rotating vector model to a position angle swing as generated by ppol. Usage:\n\n");
    printApplicationHelp(&application);
    fprintf(stdout, "Define a grid search:\n");
    fprintf(stdout, "  -g \"a b\"      Do gridsearch of a alpha and b beta beta points\n");
    fprintf(stdout, "  -A \"min max\"  Specify search range in alpha [def=\"0 180\" degrees]\n");
    fprintf(stdout, "  -B \"min max\"  Specify search range in beta  [def=\"-90 90\" degrees]\n");
    fprintf(stdout, "  -l \"l0 dl\"    Set start value and stepsize for the longitude of magnetic axis\n");
    fprintf(stdout, "                [def=\"180 90\" degrees], same as -l0.\n");
    fprintf(stdout, "  -pa \"pa0 dpa\" Set start value and stepsize for the PA at the position of the\n");
    fprintf(stdout, "                magnetic axis [def=\"0 90\" degrees], same as -pa0.\n\n");
    fprintf(stdout, "Alternative: fixed alpha and beta fitting (only fit for pa_0 and l_0):\n");
    fprintf(stdout, "  -a            Fixed alpha value in deg (use together with -b)\n");
    fprintf(stdout, "  -b            Fixed beta value in deg (use together with -a)\n");
    fprintf(stdout, "  -printfit     Print out fit of pa-swing to the terminal\n");
    fprintf(stdout, "  -printnofit   Like -printfit, but pa_0 and l_0 are taken from the -l0 and -pa0\n");
    fprintf(stdout, "                options, so no fit will be done\n");
    fprintf(stdout, "  -printsmooth  Phase-wraps are tried to be avoided with -printfit\n");
    fprintf(stdout, "\nBeam-width (rho) contour options, which can be used in grid search mode.\n              (Consider using the -cont or 'r' options):\n");
    fprintf(stdout, "  -wmp        Define pulse width of main-pulse in degrees\n");
    fprintf(stdout, "  -wmp2       Alternative contours are derived using this width in degrees\n");
    fprintf(stdout, "  -wip        Same as -wmp, but now for the interpulse\n");
    fprintf(stdout, "\nOptions affecting grid-search operation:\n");
    fprintf(stdout, "  -best       Show best solution and quit program\n");
    fprintf(stdout, "  -brute      Do normal search, then create a lattice of 16 points in l0, pa0\n");
    fprintf(stdout, "              space with point separations given by dl0 and dpa0 centered at l0\n");
    fprintf(stdout, "              and pa0.  This has a better chance to find solutions, but it will\n");
    fprintf(stdout, "              take ~17 times as long to run.\n");
    fprintf(stdout, "  -macro      Execute commands in specified file rather than via manual input\n");
    fprintf(stdout, "\nGeneral fit options:\n");
    fprintf(stdout, "  -autoopm            Shift each pa point by 90 degrees to see if it fits better\n");
    fprintf(stdout, "  -fitdh \"l h dh\"     Fit for emission height difference (as fraction of light\n");
    fprintf(stdout, "                      cylinder) at this pulse longitude l (deg), initial guess h\n");
    fprintf(stdout, "                      and initial step size dh.\n");
    fprintf(stdout, "  -maxdl              Set maximum allowed deviation of the position of the\n");
    fprintf(stdout, "                      magnetic axis from start value [def=%.1f deg].\n", fitterinfo.max_l0_diff);
    fprintf(stdout, "  -forcepa \"l PA dPA\" Forces fit to go through PA with error bar dPA\n");
    fprintf(stdout, "                      at longitude l.\n");
    fprintf(stdout, "  -opm \"l o\"          Put an OPM at longitude l with offset o in degrees in\n");
    fprintf(stdout, "                      the model. You can use this option multiple times.\n");
    fprintf(stdout, "  -tol                Set tolerance for fitting (default is %f, and 1000x\n", ftol);
    fprintf(stdout, "                      better when not doing a gridsearch).\n");
    fprintf(stdout, "\nOptions operating on input file:\n");
    fprintf(stdout, "  -dh \"l h\"   Put at pulse longitude l (in degrees) a height difference h\n");
    fprintf(stdout, "              (as fraction of the light cylinder) in the model.\n");
    fprintf(stdout, "  -dl         Apply this longitude shift in degrees to data\n");
    fprintf(stdout, "  -mc         When reading in the PA-values, each PA is taken to be a value\n");
    fprintf(stdout, "              from a Gaussian distribution defined by its error-bar. The\n");
    fprintf(stdout, "              errorbar is set to a fixed value. This allows Monte-Carlo\n");
    fprintf(stdout, "              type of analysis.\n");
    fprintf(stdout, "  -paerrfac   Multiply all PA errorbars with this factor\n");
    fprintf(stdout, "\nOther functionality:\n");
    fprintf(stdout, "  -contcol    \"edg1 edg1+ edg1- edg2 edg2+ edg2- fi fi+ fi- l0 l0+ l0- P N1 N2\"\n");
    fprintf(stdout, "              Rather than doing a fit/showing a chi^2 grid, construct a\n");
    fprintf(stdout, "              collection of contours which can be used with the R option in\n");
    fprintf(stdout, "              interactive mode. The observed pulse starts at pulse longitude\n");
    fprintf(stdout, "              edg1 with pos/neg error edg1+/edg1- in deg. Similarly the end of\n");
    fprintf(stdout, "              the pulse is defined. The fiducial plane is at fi (in deg) with \n");
    fprintf(stdout, "              errorbars and the inflection point at l0 (in deg) with errorbars. \n");
    fprintf(stdout, "              The period of the pulsar is P seconds. All numbers are expected to\n");
    fprintf(stdout, "              be positive. N1 and N2 are the number of W and rho values generated.\n");
    fprintf(stdout, "              the overall range of pulse longitudes covered by the open field line\n");
    fprintf(stdout, "              region is based on the furthest edge from the fiducial plane position.\n");
    fprintf(stdout, "\nPlot options:\n");
    fprintf(stdout, "  -device1         Specify the pgplot device for the chi^2 grid\n");
    fprintf(stdout, "  -device2         Specify the pgplot device for the PA-swing\n");
    fprintf(stdout, "  -device1res      Change the resolution of device 1 [def=\"%d %d\"]\n", gridDevice_resx, gridDevice_resy);
    fprintf(stdout, "  -device2res      Change the resolution of device 2 [def=\"%d %d\"]\n", paDevice_resx, paDevice_resy);
    fprintf(stdout, "  -devicenores     Use the default resolution decided by pgplot\n");
    fprintf(stdout, "  -cont            Enables the plotting of rho contours, without waiting for\n");
    fprintf(stdout, "                   user input in interactive mode\n");
    fprintf(stdout, "  -showwedge       Show a colour wedge indicating the reduced chi^2 scale\n");
    fprintf(stdout, "  -showwedge_label Specify label to be shown next to the wedge\n");
    fprintf(stdout, "\nFile options:\n");
    fprintf(stdout, "  -load       Load the specified dump file containing the chi^2 grid and\n");
    fprintf(stdout, "              PA-swing. A dumpfile is automatically generated by ppolFit\n");
    fprintf(stdout, "              after a grid-search. The default output name is %s.\n", dumpfile);
    fprintf(stdout, "  -save       Change the default name of the output dump file\n");
    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about fitting position-angle swings and using the beam-width information can be found in:\n");
    printf(" - 	Rookyard et al. 2015, MNRAS, 446, 3367\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {
    int index;
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
 i = index;
      }else if(strcmp(argv[i], "-a") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &alpha0, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 alphaset = 1;
 i++;
      }else if(strcmp(argv[i], "-b") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &beta0, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 betaset = 1;
 i++;
      }else if(strcmp(argv[i], "-device1") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%100s", device1, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-device2") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%100s", device2, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-device1res") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &gridDevice_resx, &gridDevice_resy, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-device2res") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &paDevice_resx, &paDevice_resy, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-devicenores") == 0) {
 devicenores = 1;
      }else if(strcmp(argv[i], "-wmp") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &fitterinfo.pulse_width, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 calculate_beam_widths = 1;
 i++;
      }else if(strcmp(argv[i], "-macro") == 0) {
 macrofilename = i+1;
 i++;
      }else if(strcmp(argv[i], "-dl") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &add_longitude_shift, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-dh") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &fitterinfo.add_height_longitude, &dh0, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-fitdh") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %lf", &fitterinfo.add_height_longitude, &dh0, &ddh0, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-wip") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &pulse_width2, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 calculate_interpulse_widths = 1;
 i++;
      }else if(strcmp(argv[i], "-wmp2") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &pulse_width2, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 calculate_interpulse_widths = 2;
 i++;
      }else if(strcasecmp(argv[i], "-paerrfac") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &paErrorFac, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-contcol") == 0
        ) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", &contcol_edge1, &contcol_edge1plus, &contcol_edge1min, &contcol_edge2, &contcol_edge2plus, &contcol_edge2min, &contcol_phi0, &contcol_phi0plus, &contcol_phi0min, &contcol_l0, &contcol_l0plus, &contcol_l0min, &contcol_p, &contcol_nr_edge_steps, &contcol_nr_fid_plane_steps, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }





 doContourRange = 1;
 i++;
      }else if(strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "-l0") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &l0, &dl0, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 lset = 1;
 i++;
      }else if(strcmp(argv[i], "-maxdl") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &fitterinfo.max_l0_diff, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-brute") == 0) {
 doBruteForce = 1;
      }else if(strcmp(argv[i], "-opm") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &fitterinfo.jump_longitude[fitterinfo.nrJumps], &fitterinfo.jump_offset[fitterinfo.nrJumps], NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 fitterinfo.nrJumps++;
 i++;
      }else if(strcmp(argv[i], "-forcepa") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %lf", &fitterinfo.force_l, &fitterinfo.force_pa, &fitterinfo.force_dpa, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 fitterinfo.force_set = 1;
 i++;
      }else if(strcmp(argv[i], "-pa") == 0 || strcmp(argv[i], "-pa0") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &pa0, &dpa0, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 paset = 1;
 i++;
      }else if(strcmp(argv[i], "-A") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &alphastart, &alphaend, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-B") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &betastart, &betaend, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-tol") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &ftol, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "showwedge_label") == 0) {
 showwedge_label = argv[i+1];
 i++;
      }else if(strcasecmp(argv[i], "-mc") == 0) {
 doMC = 1;
      }else if(strcasecmp(argv[i], "-fixseed") == 0) {
 fixseed = 1;
      }else if(strcmp(argv[i], "-cont") == 0) {
 contour_plot = 1;


      }else if(strcmp(argv[i], "-best") == 0) {
 onlyshowbest = 1;
      }else if(strcmp(argv[i], "-autoopm") == 0) {
 fitterinfo.autojump = 1;
      }else if(strcmp(argv[i], "-showwedge") == 0) {
 showwedge = 1;
      }else if(strcmp(argv[i], "-showwedge_label") == 0) {
 showwedge_label = argv[i+1];
 i++;
      }else if(strcmp(argv[i], "-load") == 0) {
 loadresults = 1;
 strcpy(dumpfile, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-save") == 0) {
 strcpy(dumpfile, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-printnofit") == 0) {
 printnofit = 1;
      }else if(strcmp(argv[i], "-printfit") == 0) {
 printfit = 1;
      }else if(strcmp(argv[i], "-printsmooth") == 0) {
 printsmooth = 1;
      }else if(strcmp(argv[i], "-g") == 0) {
 GridSearch = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &nalpha, &nbeta, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else {

 if(argv[i][0] == '-') {
   printerror(application.verbose_state.debug, "ppolFit: Unknown option: %s\n\nRun ppolFit without command line arguments to show help", argv[i]);
   terminateApplication(&application);
   return 0;
 }else {
   if(applicationAddFilename(i, application.verbose_state) == 0)
     return 0;
 }
      }
    }
  }
  if(doContourRange) {
    double wopen, hem, theta_pc, rho_pc, shiftBCW;
    fprintf(stderr, "For most the most likely input values:\n");
    fprintf(stderr, "  Pulse starts/ends at:   %lf (+%f -%f) and %lf (+%f -%f) deg\n", contcol_edge1, contcol_edge1plus, contcol_edge1min, contcol_edge2, contcol_edge2plus, contcol_edge2min);
    fprintf(stderr, "  Fiducial plane at:      %lf (+%f -%f) deg\n", contcol_phi0, contcol_phi0plus, contcol_phi0min);
    wopen = fabs(contcol_edge2 - contcol_phi0);
    if(fabs(contcol_edge1 - contcol_phi0) > wopen)
      wopen = fabs(contcol_edge1 - contcol_phi0);
    wopen *= 2.0;
    fprintf(stderr, "  This make Wopen:        %lf deg\n", wopen);
    fprintf(stderr, "  Inflection point at:    %lf (+%f -%f) deg\n", contcol_l0, contcol_l0plus, contcol_l0min);
    shiftBCW = contcol_l0 - contcol_phi0;
    fprintf(stderr, "  BCW shift:              %lf deg\n", contcol_l0 - contcol_phi0);
    shiftBCW *= M_PI/180.0;
    if(contcol_edge1plus < 0 || contcol_edge1min < 0 || contcol_edge2plus < 0 || contcol_edge2min < 0 || contcol_l0min < 0 || contcol_l0plus < 0 || contcol_phi0min < 0 || contcol_phi0plus < 0) {
      printerror(application.verbose_state.debug, "ERROR ppolFit: All errors defined for the -contcol option should be positive.\n");
      return 0;
    }
    if(shiftBCW > 0) {
      hem = contcol_p*299792458.0*shiftBCW/(8.0*M_PI);
      fprintf(stderr, "  Emission height:        %lf km\n", 0.001*hem);
      theta_pc = asin(sqrt(2.0*M_PI*hem/(contcol_p*299792458.0)));
      rho_pc = theta_pc + atan(0.5*tan(theta_pc));
      fprintf(stderr, "  Opening angle beam:     %lf deg (half-opening angle)\n", rho_pc*180.0/M_PI);
    }else {
      fprintf(stderr, "  Emission height:        Inconsistent (BCW shift is negative)\n");
    }
    double rho_min, rho_max, shiftBCW_min, shiftBCW_max;
    shiftBCW = contcol_l0 + contcol_l0plus - (contcol_phi0 - contcol_phi0min);
    if(shiftBCW < 0) {
      printerror(application.verbose_state.debug, "ERROR ppolFit: No positive values for the BCW shift allowed, so rho is undefined\n");
      return 0;
    }
    shiftBCW_max = shiftBCW;
    shiftBCW *= M_PI/180.0;
    hem = contcol_p*299792458.0*shiftBCW/(8.0*M_PI);
    theta_pc = asin(sqrt(2.0*M_PI*hem/(contcol_p*299792458.0)));
    rho_pc = theta_pc + atan(0.5*tan(theta_pc));
    rho_pc *= 180.0/M_PI;
    rho_max = rho_pc;
    shiftBCW = contcol_l0 - contcol_l0min - (contcol_phi0 + contcol_phi0plus);
    if(shiftBCW < 0) {
      shiftBCW = 0;
    }
    shiftBCW_min = shiftBCW;
    shiftBCW *= M_PI/180.0;
    hem = contcol_p*299792458.0*shiftBCW/(8.0*M_PI);
    theta_pc = asin(sqrt(2.0*M_PI*hem/(contcol_p*299792458.0)));
    rho_pc = theta_pc + atan(0.5*tan(theta_pc));
    rho_pc *= 180.0/M_PI;
    rho_min = rho_pc;
    fprintf(stderr, "Predicted range in rho:    %lf ... %lf deg\n", rho_min, rho_max);
      long k, l;
      double phi0, edge1, edge2, wopen2, wopen_min, wopen_max;
      for(l = 0; l < contcol_nr_fid_plane_steps; l++) {
 shiftBCW = shiftBCW_min + (shiftBCW_max-shiftBCW_min)*(double)l/(double)(contcol_nr_fid_plane_steps-1);
 shiftBCW *= M_PI/180.0;
 hem = contcol_p*299792458.0*shiftBCW/(8.0*M_PI);
 theta_pc = asin(sqrt(2.0*M_PI*hem/(contcol_p*299792458.0)));
 rho_pc = theta_pc + atan(0.5*tan(theta_pc));
 rho_pc *= 180.0/M_PI;
 wopen_min = wopen_max = -1.0;
 for(k = 0; k < contcol_nr_fid_plane_steps*10; k++) {
   phi0 = contcol_phi0 - contcol_phi0min + (contcol_phi0plus+contcol_phi0min)*(double)k/(double)(contcol_nr_fid_plane_steps*10-1);
   l0 = phi0 + shiftBCW*180.0/M_PI;
   if(l0 >= contcol_l0 - contcol_l0min && l0 <= contcol_l0 + contcol_l0plus) {
     for(i = 0; i < contcol_nr_edge_steps; i++) {
       edge1 = contcol_edge1 - contcol_edge1min + (contcol_edge1plus+contcol_edge1min)*(double)i/(double)(contcol_nr_edge_steps-1);
       for(j = 0; j < contcol_nr_edge_steps; j++) {
  edge2 = contcol_edge2 - contcol_edge2min + (contcol_edge2plus+contcol_edge2min)*(double)j/(double)(contcol_nr_edge_steps-1);
  wopen = edge2 - phi0;
  if(wopen < 0)
    wopen *= -1.0;
  wopen2 = phi0 - edge1;
  if(wopen2 < 0)
    wopen2 *= -1.0;
  if(wopen2 > wopen)
    wopen = wopen2;
  wopen *= 2.0;
  if(wopen_min < 0 || wopen < wopen_min) {
    wopen_min = wopen;
  }
  if(wopen_max < 0 || wopen > wopen_max) {
    wopen_max = wopen;
  }
       }
     }
   }
 }
 if(wopen_min < 0 || wopen_max < 0) {
   printerror(application.verbose_state.debug, "ERROR ppolFit: Determining range in wopen failed. Possibly N2 is too small, which can especially be an issue if the error on l0 is small.\n");
   return 0;
 }
 for(i = 0; i < contcol_nr_edge_steps; i++) {
   wopen = wopen_min + (wopen_max-wopen_min)*(double)i/(double)(contcol_nr_edge_steps-1);
   printf("%e %e\n", wopen, rho_pc);
 }
      }
    fprintf(stderr, "\nThe output can be re-directed to a file and used with the R option of ppolFit in interactive mode.\n");
    terminateApplication(&application);
    return 0;
  }
  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) > 0 && loadresults) {
    printerror(application.verbose_state.debug, "ERROR ppolFit: Unknown option \"%s\".", getNextFilenameFromList(&application, argv, application.verbose_state));
    printerror(application.verbose_state.debug, "Note: when using the -load you're not allowed to specify the PA-swing file, as it is already stored in dump file.");
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) > 1) {
    printerror(application.verbose_state.debug, "ERROR ppolFit: Only one input file can be selected on the command line.");
    return 0;
  }
  if(lset == 0) {
    if(loadresults == 0) {
      printwarning(application.verbose_state.debug, "  Assumed -l '180 90'");
      l0 = 180;
      dl0 = 90;
    }
  }
  if(paset == 0) {
    if(loadresults == 0) {
      printwarning(application.verbose_state.debug, "  Assumed -pa '0 90'");
      pa0 = 0;
      dpa0 = 90;
    }
  }
  gsl_rng *rand_num_gen;
  const gsl_rng_type *rand_num_gen_type;
  long idnum;
  gsl_rng_env_setup();
  rand_num_gen_type = gsl_rng_default;
  rand_num_gen = gsl_rng_alloc(rand_num_gen_type);
  if(fixseed)
    idnum = 1;
  else
    randomize_idnum(&idnum);
  gsl_rng_set(rand_num_gen, idnum);
  fit_dh0 = dh0;
  if(printnofit == 1) {
    for(i = 0; i < 360; i++)
      printf("%d %lf\n", i, paswing_double(alpha0, beta0, i, pa0, l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, dh0));
    terminateApplication(&application);
    gsl_rng_free(rand_num_gen);
    return 0;
  }
  if(loadresults == 0) {
    iformat = guessPSRData_format(argv[argc-1], 0, application.verbose_state);
    if(iformat != PPOL_format && iformat != PPOL_SHORT_format) {
      printerror(application.verbose_state.debug, "ppolFit: Input file does not appear to be in a ppol format.");
      return 0;
    }
    if(!openPSRData(&datain, argv[argc-1], iformat, 0, 1, 0, application.verbose_state))
      return 0;
    fitterinfo.NrDataPoints = filterPApoints(&datain, application.verbose_state);
    if(fitterinfo.NrDataPoints == 0) {
      printerror(application.verbose_state.debug, "ppolFit: No significant pa points in data.");
      return 0;
    }
    if(paErrorFac != 1.0) {
      if(application.verbose_state.verbose)
 printf("ppolFit: Multiplying input PA errors with %lf.\n", paErrorFac);
      for(i = 0; i < datain.NrBins; i++) {
 if(datain.poltype == POLTYPE_ILVPAdPA || datain.poltype == POLTYPE_ILVPAdPATEldEl) {
   datain.data[i+4*datain.NrBins] *= paErrorFac;
 }else {
   datain.data[i+1*datain.NrBins] *= paErrorFac;
 }
      }
    }
    if(doMC) {
      if(application.verbose_state.verbose)
 printf("ppolFit: Randomizing input PAs and set errors to one.\n");
      for(i = 0; i < datain.NrBins; i++) {
 if(datain.poltype == POLTYPE_ILVPAdPA || datain.poltype == POLTYPE_ILVPAdPATEldEl) {
   datain.data[i+3*datain.NrBins] += datain.data[i+4*datain.NrBins]*gsl_ran_gaussian(rand_num_gen, 1.0);
   datain.data[i+4*datain.NrBins] = 1;
 }else {
   datain.data[i ] += datain.data[i+1*datain.NrBins]*gsl_ran_gaussian(rand_num_gen, 1.0);
   datain.data[i+1*datain.NrBins] = 1;
 }
      }
    }
    if(datain.tsampMode != TSAMPMODE_LONGITUDELIST || datain.tsamp_list == NULL) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR ppolFit: Pulse longitudes were expected to be defined.");
      return 0;
    }
    fitterinfo.data_l = datain.tsamp_list;
    if(datain.poltype == POLTYPE_ILVPAdPA || datain.poltype == POLTYPE_ILVPAdPATEldEl) {
      fitterinfo.data_pa = &(datain.data[3*datain.NrBins]);
      fitterinfo.data_dpa = &(datain.data[4*datain.NrBins]);
    }else {
      fitterinfo.data_pa = datain.data;
      fitterinfo.data_dpa = &(datain.data[datain.NrBins]);
    }
    if(add_longitude_shift != 0.0) {
      if(application.verbose_state.verbose)
 printf("ppolFit: Add longitude offset of %lf to PA points.\n", add_longitude_shift);
      for(i = 0; i < fitterinfo.NrDataPoints; i++)
 fitterinfo.data_l[i] += add_longitude_shift;
    }
  }else {
    if(paErrorFac != 1.0) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "The -paerrfac option cannot be used together with the -load option");
      return 0;
    }
    if(add_longitude_shift != 0.0) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "The -dl option cannot be used together with the -load option");
      return 0;
    }
    if(doMC) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "The -mc option cannot be used together with the -load option");
      return 0;
    }
  }
  if(loadresults) {
    fin = fopen(dumpfile, "rb");
    if(fin == NULL) {
      printerror(application.verbose_state.debug, "Cannot open %s", dumpfile);
      return 0;
    }
    if(fread(txt, sizeof(char), 12, fin) != 12) {printerror(application.verbose_state.debug, "Read error from %s.", dumpfile); return 0; }
    txt[13] = 0;
    if(strcmp(txt, "fit_paswing") == 0) {
      if(fread(&version, sizeof(int), 1, fin) != 1) { printerror(application.verbose_state.debug, "Read error from %s.", dumpfile); return 0; }
    }else {
      printwarning(application.verbose_state.debug, "WARNING Unrecognized file version in %s. I assume it is written by an obsolete version of ppolFit.", dumpfile);
      rewind(fin);
      version = 0;
    }
    if(application.verbose_state.verbose)
      printf("File version of %s: %d\n", dumpfile, version);
    if(version < 0 || version > 5) {
      printerror(application.verbose_state.debug, "ERROR: ppolFit does not recognize file version %d", version);
    }
    if(version >= 2) {
      fread(&alphastart, sizeof(double), 1, fin);
      fread(&alphaend, sizeof(double), 1, fin);
      fread(&betastart, sizeof(double), 1, fin);
      fread(&betaend, sizeof(double), 1, fin);
      fread(&nalpha, sizeof(int), 1, fin);
      fread(&nbeta, sizeof(int), 1, fin);
      fread(&l0, sizeof(double), 1, fin);
      fread(&dl0, sizeof(double), 1, fin);
      fread(&pa0, sizeof(double), 1, fin);
      fread(&dpa0, sizeof(double), 1, fin);
    }else {
      float dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      alphastart = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      alphaend = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      betastart = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      betaend = dummyf;
      fread(&nalpha, sizeof(int), 1, fin);
      fread(&nbeta, sizeof(int), 1, fin);
      fread(&dummyf, sizeof(float), 1, fin);
      l0 = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      dl0 = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      pa0 = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      dpa0 = dummyf;
    }
    if(version >= 2) {
 int dummy;
 fread(&dummy, sizeof(int), 1, fin);
 fread(&dummy, sizeof(int), 1, fin);
    }
    if(version >= 3) {
      if(fitterinfo.nrJumps != 0) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "You cannot use -opm together with the -load option as information is stored in the dump file");
 return 0;
      }
      fread(&fitterinfo.nrJumps, sizeof(int), 1, fin);
      for(i = 0; i < fitterinfo.nrJumps; i++) {
 fread(&fitterinfo.jump_longitude[i], sizeof(double), 1, fin);
 fread(&fitterinfo.jump_offset[i], sizeof(double), 1, fin);
      }
      if(fitterinfo.autojump != 0) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "You cannot use -autoopm together with the -load option as information is stored in the dump file");
 return 0;
      }
      fread(&fitterinfo.autojump, sizeof(int), 1, fin);
    }
    fread(&fitterinfo.NrDataPoints, sizeof(int), 1, fin);
    if(fitterinfo.NrDataPoints < 0) {
      printerror(application.verbose_state.debug, "ppolFit: Cannot load %d data points, something appears to be wrong with dumpfile\n", fitterinfo.NrDataPoints);
      return 0;
    }
    if(fitterinfo.NrDataPoints > 10000) {
      printwarning(application.verbose_state.debug, "WARNING ppolFit: Going to load %d PA data points, which seems a large number. Maybe something is wrong with dumpfile?\n", fitterinfo.NrDataPoints);
    }
    if(application.verbose_state.debug) {
      fprintf(stdout, "  Going to load %d data points\n", fitterinfo.NrDataPoints);
      fflush(stdout);
    }
    fitterinfo.data_l = (double *)malloc(fitterinfo.NrDataPoints*sizeof(double));
    fitterinfo.data_pa = (float *)malloc(fitterinfo.NrDataPoints*sizeof(float));
    fitterinfo.data_dpa = (float *)malloc(fitterinfo.NrDataPoints*sizeof(float));
    if(fitterinfo.data_l == NULL || fitterinfo.data_pa == NULL || fitterinfo.data_dpa == NULL) {
      printerror(application.verbose_state.debug, "Cannot allocate memory");
      return 0;
    }
    for(i = 0; i < fitterinfo.NrDataPoints; i++) {
      float dummy_float;
      if(version >= 4) {
 fread(&fitterinfo.data_l[i], sizeof(double), 1, fin);
      }else {
 fread(&dummy_float, sizeof(float), 1, fin);
 fitterinfo.data_l[i] = dummy_float;
      }
      fread(&fitterinfo.data_pa[i], sizeof(float), 1, fin);
      fread(&fitterinfo.data_dpa[i], sizeof(float), 1, fin);
    }
    if(version >= 5) {
      int dummyint;
      if(fread(&dummyint, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Read error."); return 0; }
      if(dummyint < 0) {
 printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot load a string of %d bytes, something appears to be wrong with dumpfile\n", dummyint);
 return 0;
      }
      if(dummyint > MaxStringLength) {
 printerror(application.verbose_state.debug, "ERROR ppolFit: Command line stored in dump file appears to be very long (%d bytes). Maybe something is wrong with dumpfile?\n", dummyint);
      }
      if(fread(txt, 1, dummyint, fin) != dummyint) {printerror(application.verbose_state.debug, "Read error."); return 0; }
      txt[dummyint] = 0;
      if(application.verbose_state.verbose) {
 printf("Dump file created with command: %s\n", txt);
      }
    }
    GridSearch = 1;
  }
  if(alphaset == 0 || betaset == 0) {
    if(GridSearch == 0) {
      printerror(application.verbose_state.debug, "Need to specify the -g option to do a grid search over alpha and beta, or use the -a and -b option to fix their values.");
      return 0;
    }
  }
  if(GridSearch == 1) {
    chigrid = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    l0grid = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    pa0grid = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    dhgrid = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    rhogrid = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    rhogrid2 = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    if(chigrid == NULL || l0grid == NULL || pa0grid == NULL || dhgrid == NULL || rhogrid == NULL || rhogrid2 == NULL) {
      printerror(application.verbose_state.debug, "Cannot allocate memory");
      return 0;
    }
    if(loadresults == 0) {
      for(i = 0; i < nalpha; i++) {
 for(j = 0; j < nbeta; j++) {
   chigrid[nalpha*j+i] = 1e10;
   l0grid[nalpha*j+i] = 1e10;
   pa0grid[nalpha*j+i] = 1e10;
   dhgrid[nalpha*j+i] = 1e10;
 }
      }
      for(i = 0; i < nalpha; i++) {
 for(j = 0; j < nbeta; j++) {
   alpha0 = i*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
   beta0 = j*(betaend-betastart)/(double)(nbeta-1)+betastart;
   convertAlphaBeta(&alpha0, &beta0
      );
   DoFitting(alpha0, beta0, pa0, dpa0, l0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 0, stdout, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
   chi_old = chi;
   bestl0 = fit_l0;
   bestpa0 = fit_pa0;
   beshdh = fit_dh0;
   if(doBruteForce) {
     int bruteForceSignPa0, bruteForceSignL0, pa0step, l0step;
     double newpa0, newl0;
     for(bruteForceSignPa0 = -1; bruteForceSignPa0 <= 1; bruteForceSignPa0 += 2) {
       for(bruteForceSignL0 = -1; bruteForceSignL0 <= 1; bruteForceSignL0 += 2) {
  for(pa0step = 1; pa0step <= 3; pa0step += 2) {
    for(l0step = 1; l0step <= 3; l0step += 2) {
      newpa0 = pa0+0.5*(double)bruteForceSignPa0*((double)pa0step*dpa0);
      newl0 = l0+0.5*(double)bruteForceSignL0*((double)l0step*dl0);
      DoFitting(fit_alpha, fit_beta, newpa0, dpa0, newl0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 0, stdout, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
      if(chi < chi_old) {
        chi_old = chi;
        bestl0 = fit_l0;
        bestpa0 = fit_pa0;
        beshdh = fit_dh0;
      }
    }
  }
       }
     }
   }
   chi = chi_old;
   chigrid[nalpha*j+i] = chi/(double)(fitterinfo.NrDataPoints-nrfitparams);
   l0grid[nalpha*j+i] = bestl0;
   pa0grid[nalpha*j+i] = bestpa0;
   dhgrid[nalpha*j+i] = fit_dh0;
 }
 if(application.verbose_state.nocounters == 0)
   fprintf(stderr, "%.1f%%     \r",(100.0*(i+1))/(double)(nalpha));
      }
      fin = fopen(dumpfile, "wb");
      if(fin == NULL) {
 printerror(application.verbose_state.debug, "Cannot open %s", dumpfile);
 return 0;
      }
      sprintf(txt, "fit_paswing");
      txt[12] = 0;
      if(fwrite(txt, sizeof(char), 12, fin) != 12) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      version = 5;
      if(fwrite(&version, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&alphastart, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&alphaend, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&betastart, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&betaend, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&nalpha, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&nbeta, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&l0, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&dl0, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&pa0, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&dpa0, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
 int dummy;
 dummy = 0;
 if(fwrite(&dummy, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
 if(fwrite(&dummy, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&fitterinfo.nrJumps, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      for(i = 0; i < fitterinfo.nrJumps; i++) {
 if(fwrite(&fitterinfo.jump_longitude[i], sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
 if(fwrite(&fitterinfo.jump_offset[i], sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      }
      if(fwrite(&fitterinfo.autojump, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&fitterinfo.NrDataPoints, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      for(i = 0; i < fitterinfo.NrDataPoints; i++) {
 if(fwrite(&fitterinfo.data_l[i], sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
 if(fwrite(&fitterinfo.data_pa[i], sizeof(float), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
 if(fwrite(&fitterinfo.data_dpa[i], sizeof(float), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      }
      int dummyint;
      constructCommandLineString(txt, MaxStringLength, argc, argv, application.verbose_state);
      dummyint = strlen(txt);
      if(fwrite(&dummyint, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(txt, 1, dummyint, fin) != dummyint) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(chigrid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(l0grid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(pa0grid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(dhgrid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      fprintf(stderr, "Dumped data to %s\n", dumpfile);
      if(fclose(fin) != 0) {printerror(application.verbose_state.debug, "File close error."); return 0; }
    }else {
      if(fin == NULL) {
 printerror(application.verbose_state.debug, "ppolFit: BUG!!!!!!!!!!!!");
 return 0;
      }
      if(fread(chigrid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Read error."); return 0; }
      if(version >= 1) {
 if(fread(l0grid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Read error."); return 0; }
 if(fread(pa0grid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Read error."); return 0; }
 if(fread(dhgrid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Read error."); return 0; }
      }
      fclose(fin);
      if(application.verbose_state.verbose) fprintf(stderr, "Loaded data from %s\n", dumpfile);
    }
    if(calculate_beam_widths != 0 || calculate_interpulse_widths != 0) {
      calcBeamWidths(nalpha, nbeta, alphastart, alphaend, betastart, betaend,
       calculate_beam_widths, calculate_interpulse_widths, pulse_width2, rhogrid, rhogrid2, application.verbose_state.nocounters);
    }
    for(i = 0; i < nalpha; i++) {
      for(j = 0; j < nbeta; j++) {
 if(i == 0 && j == 0) {
   chimax = chigrid[nalpha*j+i];
   chimin = chigrid[nalpha*j+i];
 }
 if(chigrid[nalpha*j+i] > chimax)
   chimax = chigrid[nalpha*j+i];
 if(chigrid[nalpha*j+i] < chimin) {
   chimin = chigrid[nalpha*j+i];
   bestalpha = i*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
   bestbeta = j*(betaend-betastart)/(double)(nbeta-1)+betastart;
 }
      }
    }
    for(i = 0; i < nalpha; i++) {
      for(j = 0; j < nbeta; j++) {
 chigrid[nalpha*j+i] = (chigrid[nalpha*j+i]-chimin)/(chimax-chimin);
 chigrid[nalpha*j+i] = 1-chigrid[nalpha*j+i];
      }
    }
    fprintf(stderr, "Found chi2 values between %e and %e\n", chimin, chimax);
    fprintf(stderr, "Best fit in grid: alpha = %f and beta = %f\n", bestalpha, bestbeta);
    alpha0 = bestalpha;
    beta0 = bestbeta;
    TR[0] = alphastart-0.5*(alphaend-alphastart)/(double)nalpha;
    TR[1] = (alphaend-alphastart)/(double)nalpha;
    TR[2] = 0;
    TR[3] = betastart-0.5*(betaend-betastart)/(double)nbeta;
    TR[4] = 0;
    TR[5] = (betaend-betastart)/(double)nbeta;
      leftPulseLongitude = 0;
      rightPulseLongitude = 360;
    if(showgraphics) {
      if(onlyshowbest == 0) {
 PADeviceID = ppgopen(device2);
 if(devicenores == 0)
   pgplot_setWindowsize(paDevice_resx, paDevice_resy, -1);
 PlotPAswing(alpha0, beta0, 0, 0, 0, leftPulseLongitude, rightPulseLongitude, fit_dh0);
 GridDeviceID = ppgopen(device1);
 if(devicenores == 0)
   pgplot_setWindowsize(gridDevice_resx, gridDevice_resy, -1);
      }
      level = suppress_greyscale*chimin;
      level = (chimax - level)/(chimax - chimin);
      PrintHelp();
    }
    if(macrofilename) {
      macrofile = fopen(argv[macrofilename], "r");
      if(macrofile == NULL) {
 printerror(application.verbose_state.debug, "Cannot open %s", argv[macrofilename]);
 return 0;
      }
    }
    redraw = 1;
    do {
      if(showgraphics) {
 if(redraw) {
   if(onlyshowbest == 2) {
     GridDeviceID = ppgopen(device1);
     if(devicenores == 0)
       pgplot_setWindowsize(gridDevice_resx, gridDevice_resy, -1);
   }
   if(redraw == 2) {
     if(onlyshowbest != 1) {
       cpgslct(GridDeviceID);
       PlotGrid(chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, level, suppress_fac, GridDeviceID, bestalpha, bestbeta, boxlw, labelcharheight, boxlabelcharheight, drawCross, title_txt, drawcontours, nogray,
         chimax, chimin, PPGPLOT_GRAYSCALE, showwedge, showwedge_label);
     }
   }else {
     if(onlyshowbest != 1) {
       cpgslct(GridDeviceID);
       PlotGrid(chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, level, suppress_fac, GridDeviceID, -1, -1, boxlw, labelcharheight, boxlabelcharheight, drawCross, title_txt, drawcontours, nogray,
         chimax, chimin, PPGPLOT_GRAYSCALE, showwedge, showwedge_label);
     }
   }
   if(calculate_beam_widths != 0 && contour_plot) {
     if(onlyshowbest != 1) {
       cpgslct(GridDeviceID);
       int dotted;
       if(beamwidth_params_fin != NULL && beamwidth_params_only_w == 0) {
  rewind(beamwidth_params_fin);
  dotted = 0;
  int ret;
  do {
    if(calculate_interpulse_widths == 0)
      ret = fscanf(beamwidth_params_fin, "%lf %f", &fitterinfo.pulse_width, &userContours[0]);
    else
      ret = fscanf(beamwidth_params_fin, "%lf %f", &pulse_width2, &userContours[0]);
    if(ret == 2) {
      nruserContours = 1;
      calculate_beam_widths = 1;
      if(calculate_interpulse_widths == 0)
        printf("Calculating contours with W=%f deg and rho=%f deg\n", fitterinfo.pulse_width, userContours[0]);
      else
        printf("Calculating contours for interpulse with W=%f deg and rho=%f deg\n", pulse_width2, userContours[0]);
      calcBeamWidths(nalpha, nbeta, alphastart, alphaend, betastart, betaend,
       calculate_beam_widths, calculate_interpulse_widths, pulse_width2, rhogrid, rhogrid2, application.verbose_state.nocounters);
      PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, GridDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, dotted);
      calcIntersectionRhoAndBanana(rhogrid, chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nruserContours, userContours, chimax, chimin, beamwidth_params_vebose);
    }
  }while(ret == 2);
       }else {
  if(onlyshowbest != 1) {
    dotted = 1;
    PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, GridDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, dotted);
  }
       }
     }
   }
   if(calculate_interpulse_widths != 0 && contour_plot) {
     if(onlyshowbest != 1) {
       cpgslct(GridDeviceID);
       int dotted;
       dotted = 1;
       if(beamwidth_params_fin != NULL && beamwidth_params_only_w == 0) {
  rewind(beamwidth_params_fin);
  dotted = 0;
  int ret;
  do {
    if(calculate_interpulse_widths == 0)
      ret = fscanf(beamwidth_params_fin, "%lf %f", &fitterinfo.pulse_width, &userContours[0]);
    else
      ret = fscanf(beamwidth_params_fin, "%lf %f", &pulse_width2, &userContours[0]);
    if(ret == 2) {
      nruserContours = 1;
      calculate_interpulse_widths = 1;
      if(calculate_interpulse_widths == 0)
        printf("Calculating contours with W=%f deg and rho=%f deg\n", fitterinfo.pulse_width, userContours[0]);
      else
        printf("Calculating contours for interpulse with W=%f deg and rho=%f deg\n", pulse_width2, userContours[0]);
      calcBeamWidths(nalpha, nbeta, alphastart, alphaend, betastart, betaend,
       calculate_beam_widths, calculate_interpulse_widths, pulse_width2, rhogrid, rhogrid2, application.verbose_state.nocounters);
      PlotContours(rhogrid2, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels2, TR, GridDeviceID, contour_txt, contourcolorIP, fixedContours, nruserContours, userContours, boxlw, dotted);
      calcIntersectionRhoAndBanana(rhogrid2, chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nruserContours, userContours, chimax, chimin, beamwidth_params_vebose);
    }
  }while(ret == 2);
       }else {
  dotted = 1;
  PlotContours(rhogrid2, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels2, TR, GridDeviceID, contour_txt, contourcolorIP, fixedContours, nruserContours, userContours, boxlw, dotted);
       }
     }
   }
   redraw = 0;
   if(onlyshowbest == 2) {
     ppgclos();
   }
 }
 if(onlyshowbest == 0)
    cpgslct(GridDeviceID);
      }
 if(onlyshowbest == 1) {
   c = 98;
   onlyshowbest = 2;
      }else if(onlyshowbest == 2) {
 c = 27;
      }else {
 if(showgraphics) {
   if(macrofilename == 0) {
     float dummyf1, dummyf2;
     cpgband(7, 0, 0, 0, &dummyf1, &dummyf2, &c);
     alpha0 = dummyf1;
     beta0 = dummyf2;
   }else {
     int ret;
     do {
       ret = fscanf(macrofile, "%c", &c);
       if(ret != 1) {
  printf("Reached EOF of macro\n");
  c = 27;
       }
     }while(c == 10);
   }
 }else {
   c = 27;
 }
      }
      switch(c) {
      case 27:
      case 'q': break;
      case 108:
 printf("Set new grayscale level to this number of sigma's (it is now %f): ", (chimax - level*(chimax - chimin))/chimin - 1);
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%lf", &level);
 }else {
   fscanf(macrofile, "%lf", &level);
   printf("%lf\n", level);
 }
 level = (level+1)*chimin;
 level = (chimax - level)/(chimax - chimin);
 printf("chi2 value corresponding to black = %f and white = %f\n", chimax - level*(chimax - chimin), chimin);
 redraw = 1;
 break;
      case 65:
 {
   convertAlphaBeta(&alpha0, &beta0
      );
   i = (nalpha-1)*(alpha0 - alphastart)/(double)(alphaend-alphastart);
   j = (nbeta-1)*(beta0 - betastart)/(double)(betaend-betastart);
   double optimum_l0, optimum_pa0;
   optimum_l0 = l0grid[nalpha*j+i];
   optimum_pa0 = pa0grid[nalpha*j+i];
   if(application.verbose_state.verbose) {
     printf("Clicked on grid position: %d %d\n", i, j);
     printf("  alpha: %lf deg\n", alpha0);
     printf("  beta:  %lf deg\n", beta0);
     printf("  l0:    %lf deg\n", optimum_l0);
     printf("  pa0:   %lf deg\n", optimum_pa0);
   }
   cpgslct(PADeviceID);
   DoFitting(alpha0, beta0, pa0, dpa0, l0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 1, stdout, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
   if(doBruteForce || loadresults) {
     PlotPAswing(alpha0, beta0, optimum_pa0, optimum_l0, 1, leftPulseLongitude, rightPulseLongitude, fit_dh0);
     if(doBruteForce) {
       printwarning(application.verbose_state.debug, "WARNING: Showing solution as found during the grid-search. The fit done after clicking might not be accurate as it did not considered the range of initial parameters.  Run ppolFit with -v to see the parameters used in the plotted RVM curve.");
     }else {
       printwarning(application.verbose_state.debug, "WARNING: Showing solution as found during the grid-search when using the -load option. The fit done after clicking might not be accurate as it did not considered the range of initial parameters if the -brute option was used to generate the grid. Run ppolFit with -v to see the parameters used in the plotted RVM curve.");
     }
   }else {
     PlotPAswing(alpha0, beta0, fit_pa0, fit_l0, 1, leftPulseLongitude, rightPulseLongitude, fit_dh0);
   }
 }
 break;
      case 122:
 printf("Set left edge pulse longitude window: ");
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%lf", &leftPulseLongitude);
 }else {
   fscanf(macrofile, "%lf", &leftPulseLongitude);
   printf("%lf\n", leftPulseLongitude);
 }
 printf("Set right edge pulse longitude window: ");
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%lf", &rightPulseLongitude);
 }else {
   fscanf(macrofile, "%lf", &rightPulseLongitude);
   printf("%lf\n", rightPulseLongitude);
 }
 cpgslct(PADeviceID);
 PlotPAswing(alpha0, beta0, 0, 0, 0, leftPulseLongitude, rightPulseLongitude, fit_dh0);
 break;
      case 98:
 printf("Plotting best solution\n");
 if(onlyshowbest == 0)
   cpgslct(PADeviceID);
 alpha0 = bestalpha;
 beta0 = bestbeta;
 convertAlphaBeta(&alpha0, &beta0
    );
 DoFitting(alpha0, beta0, pa0, dpa0, l0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 1, stdout, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
 DoFitting(alpha0, beta0, fit_pa0, 10, fit_l0, 10, dh0, ddh0, ftol*0.01, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 1, 1, stdout, enableerrors, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
 alpha0 = fit_alpha;
 beta0 = fit_beta;
 printf("     Could do -paswing '%lf %lf %lf %lf' ", alpha0, beta0, fit_pa0, fit_l0);
 int jumpnr;
 if(fitterinfo.nrJumps > 0) {
   for(jumpnr = 0; jumpnr < fitterinfo.nrJumps; jumpnr++) {
     printf("-opm '%f %f' ", fitterinfo.jump_longitude[jumpnr], fitterinfo.jump_offset[jumpnr]);
   }
 }
 printf("in ppolFig\n");
 print_steepness(alpha0, beta0, fit_l0, fit_pa0, application.verbose_state.verbose, fit_dh0, &sina_b);
 if(showgraphics) {
   if(onlyshowbest != 0) {
     PADeviceID = ppgopen(device2);
     if(devicenores == 0)
       pgplot_setWindowsize(paDevice_resx, paDevice_resy, -1);
   }
   PlotPAswing(alpha0, beta0, fit_pa0, fit_l0, 1, leftPulseLongitude, rightPulseLongitude, fit_dh0);
   if(onlyshowbest != 0) {
     ppgclos();
   }
 }
 redraw = 2;
 break;
      case 116:
 if(contour_txt == 0) {
   printf("Drawing of labels for countours is now switched on\n");
   contour_txt = 1;
 }else {
   printf("Drawing of labels for countours is now switched off\n");
   contour_txt = 0;
 }
 redraw = 1;
 break;
      case 101:
 if(enableerrors != 1) {
   printf("\nThe error estimation is done by starting at the best solution. Then one by one each fit parameter is changed until the chi2 becomes (nr of sigma+1) times bigger. In this process a downhill-simplex search is done for all parameters to take into account the covariance between the parameters.\n\n");
   printf("Error estimation on best fit is now switched on\n");
   enableerrors = 1;
 }else {
   printf("Error estimation on best fit is now switched off\n");
   enableerrors = 0;
 }
 redraw = 1;
 break;
      case 84:
 if(title_txt == 0) {
   printf("Drawing title is switched on\n");
   title_txt = 1;
 }else {
   printf("Drawing title is switched off\n");
   title_txt = 0;
 }
 redraw = 1;
 break;
      case 103:
 if(nogray == 0) {
   printf("Drawing grayscale plot switched off\n");
   nogray = 1;
 }else {
   printf("Drawing grayscale plot switched on\n");
   nogray = 0;
 }
 redraw = 1;
 break;
      case 115:
 printf("Suppress factor for grayscale (>= 1, currently is set to %lf): ", suppress_fac);
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%lf", &suppress_fac);
 }else {
   fscanf(macrofile, "%lf", &suppress_fac);
   printf("%f\n", suppress_fac);
 }
 redraw = 1;
 break;
      case 83:
 printf("Set the number of sigmas in error calculation to (it is now %.1f): ", nrofsigmas);
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%lf", &nrofsigmas);
 }else {
   fscanf(macrofile, "%lf", &nrofsigmas);
   printf("%f\n", nrofsigmas);
 }
 redraw = 1;
 break;
      case 71:
 {
 double minchi2, maxchi2, chi2, minbeta, maxbeta, beta, minalpha, maxalpha, alpha, minl0, maxl0, minpa0, maxpa0, mindh, maxdh, minslope, maxslope;
 double minl0_2, maxl0_2, minpa0_2, maxpa0_2, mindh_2, maxdh_2, optimum_l0, optimum_pa0, optimum_chi2;
 int chi2unset, betaunset, fullsearch, reducel0, reducepa0, method;
 chi2unset = 1;
 betaunset = 1;
 fullsearch = 0;
 optimum_l0 = 0;
 do {
   printf("What algorithm do you want to use: \n");
   printf("  1. Look for range parameters by exploring best solutions for each alpha/beta pair (quick).\n");
   printf("  2. For each alpha/beta pair, explore how far l0 and pa0 can be pushed while solving for the other parameter (slow, but more precise).\n");
   printf("  3. Similar to option 2, but using downhill-simplex rather than Brent method.\n\n");
   printf("Type number in terminal: ");
   fflush(stdout);
   if(macrofilename == 0) {
     scanf("%d", &method);
   }else {
     fscanf(macrofile, "%d", &method);
   }
 }while(method < 1 || method > 3);
 if(method == 1) {
   fullsearch = 0;
   printf("\nDoing a quick search\n");
 }else if(method == 2) {
   fullsearch = 1;
   usegsl_error_grid_search = 1;
   printf("\nDoing a slow search with Brent method\n");
 }else if(method == 3) {
   fullsearch = 1;
   usegsl_error_grid_search = 0;
   printf("\nDoing a slow search with downhill-simplex method\n");
 }
 reducel0 = 1;
 reducepa0 = 1;
 for(i = 0; i < nalpha; i++) {
   for(j = 0; j < nbeta; j++) {
     chi2 = 1.0 - chigrid[nalpha*j+i];
     chi2 = chi2*(chimax - chimin) + chimin;
     if((i == 0 && j == 0) || chi2 < optimum_chi2) {
       optimum_chi2 = chi2;
       optimum_l0 = l0grid[nalpha*j+i];
       optimum_pa0 = pa0grid[nalpha*j+i];
     }
   }
 }
 printf("Found angles for l0 will be de-wrapped to an angle close to %lf degrees\n", optimum_l0);
 printf("Found angles for pa0 will be de-wrapped to an angle close to %lf degrees for chi2=%f\n", optimum_pa0, optimum_chi2);
 for(i = 0; i < nalpha; i++) {
   for(j = 0; j < nbeta; j++) {
     chi2 = 1.0 - chigrid[nalpha*j+i];
     chi2 = chi2*(chimax - chimin) + chimin;
     if(chi2 <= (nrofsigmas+1)*chimin) {
       if(chi2unset) {
  maxchi2 = minchi2 = chi2;
  chi2unset = 0;
       }
       if(chi2 > maxchi2)
  maxchi2 = chi2;
       if(chi2 < minchi2)
  minchi2 = chi2;
       alpha = i*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
       beta = j*(betaend-betastart)/(double)(nbeta-1)+betastart;
       convertAlphaBeta(&alpha, &beta
          );
       double tmp_l0, tmp_pa0, tmp_slope;
       tmp_l0 = l0grid[nalpha*j+i];
       if(reducel0) {
  tmp_l0 = derotate_180_small_double(tmp_l0 - optimum_l0);
       }
       tmp_pa0 = pa0grid[nalpha*j+i];
       if(reducepa0)
  tmp_pa0 = derotate_180_small_double(tmp_pa0 - optimum_pa0);
       tmp_slope = sin(alpha*M_PI/180.0)/sin(beta*M_PI/180.0);
       if(betaunset) {
  minbeta = maxbeta = beta;
  minalpha = maxalpha = alpha;
  minl0 = maxl0 = tmp_l0;
  minpa0 = maxpa0 = tmp_pa0;
  mindh = maxdh = dhgrid[nalpha*j+i];
  minl0_2 = maxl0_2 = tmp_l0;
  minpa0_2 = maxpa0_2 = tmp_pa0;
  mindh_2 = maxdh_2 = dhgrid[nalpha*j+i];
  minslope = maxslope = tmp_slope;
  betaunset = 0;
       }
       if(beta < minbeta)
  minbeta = beta;
       if(beta > maxbeta)
  maxbeta = beta;
       if(alpha < minalpha)
  minalpha = alpha;
       if(alpha > maxalpha)
  maxalpha = alpha;
       if(tmp_l0 < minl0)
  minl0 = tmp_l0;
       if(tmp_l0 > maxl0)
  maxl0 = tmp_l0;
       if(tmp_pa0 < minpa0)
  minpa0 = tmp_pa0;
       if(tmp_pa0 > maxpa0)
  maxpa0 = tmp_pa0;
       if(dhgrid[nalpha*j+i] < mindh)
  mindh = dhgrid[nalpha*j+i];
       if(dhgrid[nalpha*j+i] > maxdh)
  maxdh = dhgrid[nalpha*j+i];
       if(tmp_slope > maxslope)
  maxslope = tmp_slope;
       if(tmp_slope < minslope)
  minslope = tmp_slope;
       if(fullsearch) {
  double dx[5], xfit[5], dplus, dmin;
  int fixed[5], param_i, param_j;
  for(param_i = 0; param_i < 3; param_i++) {
    for(param_j = 0; param_j < 5; param_j++) {
      dx[param_j] = 1;
      fixed[param_j] = 1;
    }
    xfit[0] = pa0grid[nalpha*j+i];
    xfit[1] = l0grid[nalpha*j+i];
    fitterinfo.l0_start = xfit[1];
    xfit[2] = alpha;
    xfit[3] = beta;
    xfit[4] = dhgrid[nalpha*j+i];
    fixed[0] = 0;
    fixed[1] = 0;
    fixed[4] = 0;
    nrfitparams = 3;
    if(dh0 < -2 ) {
      fixed[4] = 1;
      fixed[2] = 0;
      dx[2] = 0;
      nrfitparams -= 1;
      if(param_i == 2)
        break;
    }else {
      if(usegsl_error_grid_search) {
        printerror(application.verbose_state.debug, "Searching for errorbars when allowing emission height difference will not work without using downhill-simplex, as it is no longer a 1D fitting process.");
        exit(0);
      }
    }
    int debug_verbose;
    debug_verbose = 0;
    double chi2_notreduced;
    chi2_notreduced = chimin*(fitterinfo.NrDataPoints-nrfitparams);
    if(param_i == 0) {
      if(usegsl_error_grid_search == 0) {
        if(find_errors_amoeba_d(amoeba_algorithm, dx, fixed, xfit, chi2_notreduced, 5, funk, ftol, 0, &dplus, &dmin, nrofsigmas) != 0)
   return 0;
      }else {
        internal_fit_pa_or_l0 = 1;
        ret = find_1D_error(internal_funk_gsl, xfit, 0, 5, 1, -1, NULL, nrofsigmas, chi2_notreduced, 2000, 0.01, 0.0, &dplus, debug_verbose);
        if(ret != 0) {
   fflush(stdout);
   printwarning(application.verbose_state.debug, "WARNING: Finding error failed for this trial");
   if(ret == 1) {
     printwarning(application.verbose_state.debug, "WARNING: Maximum nr of itterations exceeded");
   }else if(ret == 2) {
     printwarning(application.verbose_state.debug, "WARNING: Did not found root");
   }else if(ret == 3) {
     printwarning(application.verbose_state.debug, "WARNING: Lower and upper limit do not bracket a root");
   }else {
     printwarning(application.verbose_state.debug, "WARNING: Unknown error in minimize_1D_double");
   }
        }
        ret = find_1D_error(internal_funk_gsl, xfit, 0, 5, 1, -1, NULL, -nrofsigmas, chi2_notreduced, 2000, 0.01, 0.0, &dmin, debug_verbose);
        dmin *= -1.0;
        if(ret != 0) {
   fflush(stdout);
   printwarning(application.verbose_state.debug, "WARNING: Finding error failed for this trial");
   if(ret == 1) {
     printwarning(application.verbose_state.debug, "WARNING: Maximum nr of itterations exceeded");
   }else if(ret == 2) {
     printwarning(application.verbose_state.debug, "WARNING: Did not found root");
   }else if(ret == 3) {
     printwarning(application.verbose_state.debug, "WARNING: Lower and upper limit do not bracket a root");
   }else {
     printwarning(application.verbose_state.debug, "WARNING: Unknown error in minimize_1D_double");
   }
        }
      }
      tmp_pa0 = xfit[0]+dplus;
      if(reducepa0)
        tmp_pa0 = derotate_180_small_double(tmp_pa0 - optimum_pa0);
      if(tmp_pa0 > maxpa0_2)
        maxpa0_2 = tmp_pa0;
      tmp_pa0 = xfit[0]+dmin;
      if(reducepa0)
        tmp_pa0 = derotate_180_small_double(tmp_pa0 - optimum_pa0);
      if(tmp_pa0 < minpa0_2)
        minpa0_2 = tmp_pa0;
    }else if(param_i == 1) {
      if(usegsl_error_grid_search == 0) {
        if(find_errors_amoeba_d(amoeba_algorithm, dx, fixed, xfit, chi2_notreduced, 5, funk, ftol, 1, &dplus, &dmin, nrofsigmas) != 0)
   return 0;
      }else {
        internal_fit_pa_or_l0 = 2;
        ret = find_1D_error(internal_funk_gsl, xfit, 1, 5, 1, fitterinfo.max_l0_diff, NULL, nrofsigmas, chi2_notreduced, 2000, 0.01, 0.0, &dplus, debug_verbose);
        if(ret != 0) {
   fflush(stdout);
   printwarning(application.verbose_state.debug, "WARNING: Finding error failed for this trial");
   if(ret == 1) {
     printwarning(application.verbose_state.debug, "WARNING: Maximum nr of itterations exceeded");
   }else if(ret == 2) {
     printwarning(application.verbose_state.debug, "WARNING: Did not found root");
   }else if(ret == 3) {
     printwarning(application.verbose_state.debug, "WARNING: Lower and upper limit do not bracket a root");
   }else {
     printwarning(application.verbose_state.debug, "WARNING: Unknown error in minimize_1D_double");
   }
        }
        ret = find_1D_error(internal_funk_gsl, xfit, 1, 5, 1, fitterinfo.max_l0_diff, NULL, -nrofsigmas, chi2_notreduced, 2000, 0.01, 0.0, &dmin, debug_verbose);
        dmin *= -1.0;
        if(ret != 0) {
   fflush(stdout);
   printwarning(application.verbose_state.debug, "WARNING: Finding error failed for this trial");
   if(ret == 1) {
     printwarning(application.verbose_state.debug, "WARNING: Maximum nr of itterations exceeded");
   }else if(ret == 2) {
     printwarning(application.verbose_state.debug, "WARNING: Did not found root");
   }else if(ret == 3) {
     printwarning(application.verbose_state.debug, "WARNING: Lower and upper limit do not bracket a root");
   }else {
     printwarning(application.verbose_state.debug, "WARNING: Unknown error in minimize_1D_double");
   }
        }
      }
      tmp_l0 = xfit[1]+dplus;
      if(reducel0)
        tmp_l0 = derotate_180_small_double(tmp_l0 - optimum_l0);
      if(tmp_l0 > maxl0_2)
        maxl0_2 = tmp_l0;
      tmp_l0 = xfit[1]+dmin;
      if(reducel0)
        tmp_l0 = derotate_180_small_double(tmp_l0 - optimum_l0);
      if(tmp_l0 < minl0_2)
        minl0_2 = tmp_l0;
    }else if(param_i == 2) {
      if(find_errors_amoeba_d(amoeba_algorithm, dx, fixed, xfit, chi2_notreduced, 5, funk, ftol, 4, &dplus, &dmin, nrofsigmas) != 0)
        return 0;
      if(xfit[4]+dplus > maxdh_2)
        maxdh_2 = xfit[4]+dplus;
      if(xfit[4]+dmin < mindh_2)
        mindh_2 = xfit[4]+dmin;
    }
  }
       }
     }
     if(application.verbose_state.nocounters == 0)
       fprintf(stderr, "%.1f%%     \r",(100.0*(i+1))/(double)(nalpha));
   }
 }
 fprintf(stderr, "\n");
 printf("After considering all the best solutions in the alpha-beta grid the following %.1f sigma ranges were found:\n", nrofsigmas);
 printf("The range in chi2                 = %f to %f", minchi2, maxchi2);
 printf("\n");
 printf("The range in alpha                = %f to %f deg\n", minalpha, maxalpha);
 printf("The range in beta                 = %f to %f deg\n", minbeta, maxbeta);
 printf("The range in l0                   = %f to %f deg", minl0, maxl0);
 if(reducel0) {
   printf(" offset with respect to %lf\n", optimum_l0);
 }else {
   printf("\n");
 }
 printf("The range in pa0                  = %f to %f", minpa0, maxpa0);
 if(reducepa0) {
   printf(" offset with respect to %lf\n", optimum_pa0);
 }else {
   printf("\n");
 }
 printf("The range in sin(alpha)/sin(beta) = %f to %f\n", minslope, maxslope);
 printf("The range in dh                   = %f to %f (fraction of light cylinder radius)\n", mindh, maxdh);
 if(fullsearch) {
   printf("\nEach of those solutions were taken as a start point of a calculation similar to the 'e' option. For each solution the allowed ranges for l0 (and pa0 and dh if used) were explored while fitting for pa0 and dh if used (or the other permutations in case the range in pa0 or dh are solved for) the following ranges were found:\n");
   printf("The range in l0    = %f to %f deg", minl0_2, maxl0_2);
   if(reducel0) {
     printf(" offset with respect to %lf\n", optimum_l0);
   }else {
     printf("\n");
   }
   printf("The range in pa0   = %f to %f deg", minpa0_2, maxpa0_2);
   if(reducepa0) {
     printf(" offset with respect to %lf\n", optimum_pa0);
   }else {
     printf("\n");
   }
   printf("The range in dh    = %f to %f deg\n", mindh_2, maxdh_2);
   redraw = 1;
 }
 }
 break;
      case 99:
 printf("Nr contour levels: ");
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%d", &drawcontours);
 }else {
   fscanf(macrofile, "%d", &drawcontours);
   printf("%d\n", drawcontours);
 }
 redraw = 1;
 break;
      case 111:
 printf("Change output plot parameters\n");
 printf("line width of box (now set to %d): ", boxlw);
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%d", &boxlw);
 }else {
   fscanf(macrofile, "%d", &boxlw);
   printf("%d\n", boxlw);
 }
 printf("character height labels (now set to %f): ", labelcharheight);
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%lf", &labelcharheight);
 }else {
   fscanf(macrofile, "%lf", &labelcharheight);
   printf("%f\n", labelcharheight);
 }
 printf("label numbers character height (now set to %f): ", boxlabelcharheight);
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%lf", &boxlabelcharheight);
 }else {
   fscanf(macrofile, "%lf", &boxlabelcharheight);
   printf("%f\n", boxlabelcharheight);
 }
 printf("Colour index for contours (now set to %d): ", contourcolor);
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%d", &contourcolor);
 }else {
   fscanf(macrofile, "%d", &contourcolor);
   printf("%d\n", contourcolor);
 }
 printf("Colour index for contours of interpulse (now set to %d): ", contourcolorIP);
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%d", &contourcolorIP);
 }else {
   fscanf(macrofile, "%d", &contourcolorIP);
   printf("%d\n", contourcolorIP);
 }
 redraw = 1;
 break;
      case 114:
 if(contour_plot == 0) {
   printf("Countour plot for beam radius is switched on\n");
   contour_plot = 1;
 }else {
   printf("Countour plot for beam radius is switched off\n");
   contour_plot = 0;
 }
 redraw = 1;
 break;
      case 110:
 fixedContours = 0;
 nruserContours = 0;
 printf("Nr contour levels: ");
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%d", &nrcontourlevels);
 }else {
   fscanf(macrofile, "%d", &nrcontourlevels);
   printf("%d\n", nrcontourlevels);
 }
 printf("nr contour levels interpulse: ");
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%d", &nrcontourlevels2);
 }else {
   fscanf(macrofile, "%d", &nrcontourlevels2);
   printf("%d\n", nrcontourlevels2);
 }
 redraw = 1;
 break;
      case 102:
 if(fixedContours) {
   fixedContours = 0;
   printf("Don't use fixed contours.\n");
 }else {
   fixedContours = 1;
   printf("Use fixed contours.\n");
 }
 redraw = 1;
 break;
      case 70:
 printf("Specify contours (e.g. 5,10,20, hit return without specifying anything to disable drawing contours):\n");
 if(macrofilename == 0) {
   fgets(txt, 990, stdin);
 }else {
   int index, ret;
   char character;
   index = 0;
   while((ret = fscanf(macrofile, "%c", &character))==1) {
     txt[index++] = character;
     if(index == 990 || character == '\n' || character == '\r') {
       txt[index] = 0;
       break;
     }
   }
 }
 if(strlen(txt) == 0) {
   printwarning(application.verbose_state.debug, "WARNING: empty string. If not correct, try F again.");
 }
 txtptr = strtok(txt, ",");
 nruserContours = 0;
 do {
   if(txtptr != NULL) {
     sscanf(txtptr, "%f", &userContours[nruserContours++]);
     printf("contour %d = %f\n", nruserContours, userContours[nruserContours-1]);
   }
   txtptr = strtok(NULL, ",");
 }while(txtptr != NULL);
 redraw = 1;
 break;
      case 82:
 {
   char beamwidth_params_file[1000];
   printf("Specify ascii file with pulse a width and rho contour (could use the -contcol option): ");
   fflush(stdout);
   if(macrofilename == 0) {
     scanf("%s", beamwidth_params_file);
   }else {
     int ret;
     ret = fscanf(macrofile, "%s", beamwidth_params_file);
     if(ret != 1) {
       printerror(application.verbose_state.debug, "ERROR ppolFit: failed to read in string from macrofile");
       return 0;
     }
   }
   printf("Will read parameters from %s\n", beamwidth_params_file);
   beamwidth_params_fin = fopen(beamwidth_params_file, "r");
   if(beamwidth_params_fin == NULL) {
     printerror(application.verbose_state.debug, "ERROR ppolFit: cannot open '%s'", beamwidth_params_file);
     return 0;
   }
 }
 printf("Specify MP or IP: ");
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%s", txt);
 }else {
   ret = fscanf(macrofile, "%s", txt);
   if(ret != 1) {
     printerror(application.verbose_state.debug, "ERROR ppolFit: failed to read in string from macrofile");
     return 0;
   }
 }
 if(strcasecmp(txt, "IP") == 0) {
   calculate_interpulse_widths = 1;
   calculate_beam_widths = 0;
   printf("Take interpulse\n");
 }else {
   calculate_interpulse_widths = 0;
   calculate_beam_widths = 1;
   printf("Take main pulse\n");
 }
 printf("Plot intersections between chi2 surface and contours (y/n): ");
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%s", txt);
 }else {
   ret = fscanf(macrofile, "%s", txt);
   if(ret != 1) {
     printerror(application.verbose_state.debug, "ERROR ppolFit: cannot read in string from macrofile");
     return 0;
   }
 }
 if(strcasecmp(txt, "y") == 0) {
   beamwidth_params_vebose.verbose = 1;
   printf("Do show intersections between chi2 surface and contours\n");
 }else {
   beamwidth_params_vebose.verbose = 0;
   printf("Do not show intersections between chi2 surface and contours\n");
 }
 printf("Draw contours ONLY with w option, i.e. not to current device (y/n): ");
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%s", txt);
 }else {
   ret = fscanf(macrofile, "%s", txt);
   if(ret != 1) {
     printerror(application.verbose_state.debug, "ERROR ppolFit: cannot read in string from macrofile");
     return 0;
   }
 }
 if(strcasecmp(txt, "y") == 0) {
   beamwidth_params_only_w = 1;
   printf("Drawing contours ONLY with w option\n");
 }else {
   beamwidth_params_only_w = 0;
   printf("Always drawing contours\n");
 }
 redraw = 1;
 break;
      case 67:
 if(drawCross == 2) {
   drawCross = 0;
   printf("Disabled drawing of the best-fit cross\n");
 }else if(drawCross == 1) {
   drawCross = 2;
   printf("Enabled drawing of the best-fit cross (Black)\n");
 }else if(drawCross == 0) {
   drawCross = 1;
   printf("Enabled drawing of the best-fit cross (Red)\n");
 }
 redraw = 1;
 break;
      case 119:
 invertGrayscale = 1;
 printf("Do you want to specify pgplot device for chi2 grid yourself (yes/no)? ");
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%s", txt);
 }else {
   ret = fscanf(macrofile, "%s", txt);
   if(ret != 1) {
     printerror(application.verbose_state.debug, "ERROR ppolFit: cannot read in string from macrofile");
     return 0;
   }
 }
 if(strcasecmp(txt, "yes") == 0) {
   printf("Specify pgplot device: ");
   fflush(stdout);
   if(macrofilename == 0) {
     scanf("%s", txt);
   }else {
     ret = fscanf(macrofile, "%s", txt);
     if(ret != 1) {
       printerror(application.verbose_state.debug, "ERROR ppolFit: cannot read in string from macrofile");
       return 0;
     }
   }
   printf("Writing output to device %s\n", txt);
   PSDeviceID = ppgopen(txt);
   if(devicenores == 0)
     pgplot_setWindowsize(gridDevice_resx, gridDevice_resy, -1);
   printf("Do you want to invert grayscale (yes/no)? ");
   fflush(stdout);
   if(macrofilename == 0) {
     scanf("%s", txt);
   }else {
     ret = fscanf(macrofile, "%s", txt);
     if(ret != 1) {
       printerror(application.verbose_state.debug, "ERROR ppolFit: cannot read in string from macrofile");
       return 0;
     }
   }
   if(strcasecmp(txt, "yes") == 0) {
     invertGrayscale = 1;
     printf("Inverting grayscale\n");
   }else {
     invertGrayscale = 0;
     printf("Not inverting grayscale\n");
   }
 }else {
   printf("Writing output to %s.ps\n", prefix);
   sprintf(txt, "%s.ps/cps", prefix);
   PSDeviceID = ppgopen(txt);
 }
 ppgask(0);
 cpgslct(PSDeviceID);
 if(invertGrayscale) {
   PlotGrid(chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, level, suppress_fac, PSDeviceID, bestalpha, bestbeta, boxlw, labelcharheight, boxlabelcharheight, drawCross, title_txt, drawcontours, nogray,
     chimax, chimin, PPGPLOT_INVERTED_GRAYSCALE, showwedge, showwedge_label);
 }else {
   PlotGrid(chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, level, suppress_fac, PSDeviceID, bestalpha, bestbeta, boxlw, labelcharheight, boxlabelcharheight, drawCross, title_txt, drawcontours, nogray,
     chimax, chimin, PPGPLOT_GRAYSCALE, showwedge, showwedge_label);
 }
 if(calculate_beam_widths != 0 && contour_plot) {
   int dotted;
   if(beamwidth_params_fin != NULL) {
     rewind(beamwidth_params_fin);
     dotted = 0;
     int ret;
     do {
       if(calculate_interpulse_widths == 0)
  ret = fscanf(beamwidth_params_fin, "%lf %f", &fitterinfo.pulse_width, &userContours[0]);
       else
  ret = fscanf(beamwidth_params_fin, "%lf %f", &pulse_width2, &userContours[0]);
       if(ret == 2) {
  nruserContours = 1;
  calculate_beam_widths = 1;
  if(calculate_interpulse_widths == 0)
    printf("Calculating contours with W=%lf deg and rho=%f deg\n", fitterinfo.pulse_width, userContours[0]);
  else
    printf("Calculating contours for ip with W=%f deg and rho=%f deg\n", pulse_width2, userContours[0]);
  calcBeamWidths(nalpha, nbeta, alphastart, alphaend, betastart, betaend,
          calculate_beam_widths, calculate_interpulse_widths, pulse_width2, rhogrid, rhogrid2, application.verbose_state.nocounters);
  PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, PSDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, dotted);
       }
     }while(ret == 2);
   }else {
     dotted = 1;
     PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, PSDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, dotted);
     PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, PSDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, dotted);
   }
 }
 if(calculate_interpulse_widths != 0 && contour_plot) {
   int dotted;
   if(beamwidth_params_fin != NULL) {
     rewind(beamwidth_params_fin);
     dotted = 0;
     int ret;
     do {
       if(calculate_interpulse_widths == 0)
  ret = fscanf(beamwidth_params_fin, "%lf %f", &fitterinfo.pulse_width, &userContours[0]);
       else
  ret = fscanf(beamwidth_params_fin, "%lf %f", &pulse_width2, &userContours[0]);
       if(ret == 2) {
  nruserContours = 1;
  calculate_interpulse_widths = 1;
  if(calculate_interpulse_widths == 0)
    printf("Calculating contours with W=%lf deg and rho=%f deg\n", fitterinfo.pulse_width, userContours[0]);
  else
    printf("Calculating contours for ip with W=%f deg and rho=%f deg\n", pulse_width2, userContours[0]);
  calcBeamWidths(nalpha, nbeta, alphastart, alphaend, betastart, betaend,
          calculate_beam_widths, calculate_interpulse_widths, pulse_width2, rhogrid, rhogrid2, application.verbose_state.nocounters);
  PlotContours(rhogrid2, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels2, TR, GridDeviceID, contour_txt, contourcolorIP, fixedContours, nruserContours, userContours, boxlw, dotted);
       }
     }while(ret == 2);
   }else {
     dotted = 1;
     PlotContours(rhogrid2, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels2, TR, GridDeviceID, contour_txt, contourcolorIP, fixedContours, nruserContours, userContours, boxlw, dotted);
   }
 }
 ppgclos();
 printf("Writing chi2 grid done. Making a PA-swing plot is next.\n");
 printf("alpha: ");
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%lf", &alpha0);
 }else {
   fscanf(macrofile, "%lf", &alpha0);
   printf("%f\n", alpha0);
 }
 printf("beta: ");
 fflush(stdout);
 if(macrofilename == 0) {
   scanf("%lf", &beta0);
 }else {
   fscanf(macrofile, "%lf", &beta0);
   printf("%f\n", beta0);
 }
 DoFitting(alpha0, beta0, pa0, dpa0, l0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 1, stdout, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
 printf("Writing output to pa.ps\n");
 ppgopen("pa.ps/cps");
 PlotPAswing(alpha0, beta0, fit_pa0, fit_l0, 1, leftPulseLongitude, rightPulseLongitude, fit_dh0);
 ppgclos();
 cpgslct(GridDeviceID);
 printf("Making plots done.\n");
 break;
      case 104:
 PrintHelp();
 break;
      default: printf("Unknown key: %d\n", c);
      }
    }while(c != 27 && c != 'q');
    if(showgraphics)
      ppgend();
    free(chigrid);
    free(l0grid);
    free(pa0grid);
    free(dhgrid);
    free(rhogrid);
    free(rhogrid2);
  }else {
    DoFitting(alpha0, beta0, pa0, dpa0, l0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 1, stderr, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
    alpha0 = fit_alpha;
    beta0 = fit_beta;
    if(application.verbose_state.verbose) printf("\n\n#Fitting procedure: ");
    for(i = 0; i < argc; i++)
      if(application.verbose_state.verbose) printf("%s ", argv[i]);
    if(application.verbose_state.verbose) printf("\n#\n#After %d steps the downhill-simplex found:\n#pa0 = %f and l0 = %f with chi2 = %e\n", nfunk, fit_pa0, fit_l0, chi);
    if(application.verbose_state.verbose) printf("#\n");
    if(printfit == 1) {
      dpa0 = paswing_double(alpha0, beta0, i, fit_pa0, fit_l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, fit_dh0);
      for(i = 0; i < 360; i++) {
 pa0 = paswing_double(alpha0, beta0, i, fit_pa0, fit_l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, fit_dh0);
 if(printsmooth) {
   if(pa0 - dpa0 > 90)
     pa0 -= 180;
   else if(dpa0 - pa0 > 90)
     pa0 += 180;
   if(pa0 - dpa0 > 90)
     pa0 -= 180;
   else if(dpa0 - pa0 > 90)
     pa0 += 180;
   dpa0 = pa0;
 }
 if(i == 0)
   printf("# BinNr PA[deg] alpha[deg] beta[deg] longitude_0[deg] PA_0[deg]\n");
 printf("%d %f %f %f %f %f\n", i, pa0, alpha0, beta0, fit_l0, fit_pa0);
      }
    }
  }
  if(loadresults == 0) {
    closePSRData(&datain, 0, application.verbose_state);
  }
  if(macrofilename)
    fclose(macrofile);
  terminateApplication(&application);
  gsl_rng_free(rand_num_gen);
  return 0;
}
double funk(double x[])
{
  double pa, y1, dy, chi2, alpha, beta, heightshift;
  int i;
  alpha = x[2];
  beta = x[3];
  heightshift = x[4];
  chi2 = 0;
  if(alpha == 0.0000)
    alpha = 0.00001;
  if(beta == 0.000)
    beta = 0.00001;
  if(fabs(fitterinfo.l0_start-x[1]) > fitterinfo.max_l0_diff) {
    return 1e10;
  }
  if(fitterinfo.force_set) {
    pa = paswing_double(alpha, beta, fitterinfo.force_l, x[0], x[1], fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, heightshift);
    if(fitterinfo.autojump)
      dy = dy_90(pa, fitterinfo.force_pa);
    else
      dy = dy_180(pa, fitterinfo.force_pa);
    if(dy > fitterinfo.force_dpa)
      return 1e10;
  }
  for(i = 0; i < fitterinfo.NrDataPoints; i++) {
    pa = paswing_double(alpha, beta, fitterinfo.data_l[i], x[0], x[1], fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, heightshift);
    y1 = fitterinfo.data_pa[i];
    if(fitterinfo.autojump)
      dy = dy_90(pa, y1);
    else
      dy = dy_180(pa, y1);
    chi2 += dy*dy/((double)(fitterinfo.data_dpa[i])*(double)(fitterinfo.data_dpa[i]));
  }
  return chi2;
}
double internal_PA0_funk(double pa0, void *params)
{
  double chi2;
  double *x;
  x = params;
  x[0] = pa0;
  chi2 = funk(x);
  return chi2;
}
double internal_L0_funk(double l0, void *params)
{
  double chi2;
  double *x;
  x = params;
  x[1] = l0;
  chi2 = funk(x);
  return chi2;
}
double internal_funk_gsl(double *x, void *params)
{
  int ret, nrpoints, debug_verbose, debug_verbose2, trial, higher_res_step;
  double chi2, chi2_before, pa0, l0, xnew[5], ftol;
  nrpoints = 200;
  ftol = 0.00001;
  debug_verbose = 0;
  debug_verbose2 = 0;
  for(higher_res_step = 0; higher_res_step < 2; higher_res_step++) {
    memcpy(xnew, x, 5*sizeof(double));
    if(higher_res_step != 0) {
      fflush(stdout);
      printwarning(0, "Failed to find global minimum, trying to use a higher resolution grid");
      nrpoints *= 100;
    }
    chi2_before = funk(xnew);
    for(trial = 0; trial < 2; trial++) {
      if(internal_fit_pa_or_l0 == 1)
 ret = minimize_1D_double(0, internal_L0_funk, x, fitterinfo.l0_start-fitterinfo.max_l0_diff, fitterinfo.l0_start+fitterinfo.max_l0_diff, nrpoints, 1, 2, &l0, 2000, ftol, 0.0, debug_verbose, debug_verbose2);
      else
 ret = minimize_1D_double(0, internal_PA0_funk, x, 0, 180, nrpoints, 0, 2, &pa0, 2000, ftol, 0.0, debug_verbose, debug_verbose2);
      if(ret == 1) {
 fflush(stdout);
 ftol *= 100;
 printwarning(debug_verbose, "WARNING: Maximum nr of itterations exceeded in internal_funk_gsl (alpha=%f, beta=%f)\nTrying to lower tolerance to %f", x[2], x[3], ftol);
      }else {
 break;
      }
    }
    if(ret != 0) {
      if(ret == 1) {
 fflush(stdout);
 printwarning(debug_verbose, "WARNING: Maximum nr of itterations exceeded in internal_funk_gsl (alpha=%f, beta=%f)\nTrying to lower tolerance", x[2], x[3]);
      }else if(ret == 2) {
 fflush(stdout);
 printwarning(debug_verbose, "WARNING: Did not found root");
      }else if(ret == 3) {
 fflush(stdout);
 printwarning(debug_verbose, "WARNING: Lower and upper limit do not bracket a root in internal_funk_gsl");
      }else {
 fflush(stdout);
 printwarning(debug_verbose, "WARNING: Unknown error in minimize_1D_double in internal_funk_gsl");
 exit(0);
      }
    }
    if(internal_fit_pa_or_l0 == 1) {
      xnew[1] = l0;
    }else {
      xnew[0] = pa0;
    }
    chi2 = funk(xnew);
    if(1.05*chi2_before < chi2) {
      fflush(stdout);
      printwarning(0, "Optimization failed: %f > %f for alpha=%f beta=%f", chi2, chi2_before, x[2], x[3]);
      if(higher_res_step > 0)
 exit(0);
    }else {
      break;
    }
  }
  return chi2;
}
void DoFitting(double alpha0, double beta0, double pa0, double dpa0, double l0, double dl0, double dh0, double ddh0, double ftol, double *fit_pa0, double *fit_l0, double *fit_a, double *fit_b, double *fit_dh, double *chi, int *nfunk, int searchAll, int report, FILE *reportStream, int finderrors, double nrofsigmas, int *nfitparameters, int amoeba_algorithm, verbose_definition verbose)
{
  double xstart[5], dx[5], xfit[5], dplus[5], dmin[5], chi_d;
  int fixed[5];
  fitterinfo.l0_start = l0;
  xstart[0] = pa0;
  xstart[1] = l0;
  xstart[2] = alpha0;
  xstart[3] = beta0;
  xstart[4] = dh0;
  dx[0] = dpa0;
  dx[1] = dl0;
  dx[2] = 10;
  dx[3] = 10;
  dx[4] = ddh0;
  fixed[0] = 0;
  fixed[1] = 0;
  fixed[2] = 0;
  fixed[3] = 0;
  fixed[4] = 0;
  if(searchAll) {
    *nfitparameters = 4;
  }else {
    fixed[2] = 1;
    fixed[3] = 1;
    *nfitparameters = 2;
  }
  if(ddh0 >= -2) {
    (*nfitparameters) += 1;
  }else {
    fixed[4] = 1;
    xstart[4] = dh0;
  }
  if(finderrors && *nfitparameters <= 2) {
    printwarning(verbose.debug, "WARNING: Need at least three fit parameters to do error estimation.");
    finderrors = 0;
  }
  do {
    if(finderrors) {
      printf("param0 = pa0, param1 = l0, param2 = alpha, param3 = beta, param4=dh\n");
    }
    if(doAmoeba_d(amoeba_algorithm, xstart, dx, fixed, xfit, &chi_d, 5, funk, ftol, nfunk, 0, finderrors, nrofsigmas, dplus, dmin) == 1) {
      printwarning(verbose.debug, "WARNING: Adjusting downhill-simplex tollerance to try to converge.");
      ftol *= 10;
    }else {
      *chi = chi_d;
      break;
    }
  }while(ftol < 0.01);
  *fit_pa0 = xfit[0];
  *fit_l0 = xfit[1];
  *fit_a = xfit[2];
  *fit_b = xfit[3];
  *fit_dh = xfit[4];
  if(fixed[4]) {
    if(fitterinfo.add_height_longitude <= 360)
      *fit_dh = dh0;
    else
      *fit_dh = 0;
  }
  if(report) {
    if(searchAll == 0) {
      fprintf(reportStream, "Fit using specified alpha and beta value:");
    }else {
      fprintf(reportStream, "Fit keeping alpha and beta as a free parameter:");
    }
    fprintf(reportStream, "\n     alpha = %15f deg", *fit_a);
    if(finderrors) fprintf(reportStream, " (%.1f sigma error: %+15f %15f -> range = %f to %f)", nrofsigmas, dplus[2], dmin[2], *fit_a+dmin[2], *fit_a+dplus[2]);
    fprintf(reportStream, "\n     beta  = %15f deg", *fit_b);
    if(finderrors) fprintf(reportStream, " (%.1f sigma error: %+15f %15f -> range = %f to %f)", nrofsigmas, dplus[3], dmin[3], *fit_b+dmin[3], *fit_b+dplus[3]);
    fprintf(reportStream, "\n           = %15f for the other pole", *fit_b+2.0*(*fit_a)-180.0);
    fprintf(reportStream, "\n     l0    = %15f deg", *fit_l0);
    if(finderrors) fprintf(reportStream, " (%.1f sigma error: %+15f %15f -> range = %f to %f)", nrofsigmas, dplus[1], dmin[1], *fit_l0+dmin[1], *fit_l0+dplus[1]);
    fprintf(reportStream, "\n     pa0   = %15f deg", *fit_pa0);
    if(finderrors) fprintf(reportStream, " (%.1f sigma error: %+15f %15f -> range = %f to %f)", nrofsigmas, dplus[0], dmin[0], *fit_pa0+dmin[0], *fit_pa0+dplus[0]);
    fprintf(reportStream, "\n     dh    = %15f", *fit_dh);
    if(finderrors) fprintf(reportStream, " (%.1f sigma error: %+15f %15f -> range = %f to %f)", nrofsigmas, dplus[4], dmin[4], *fit_dh+dmin[4], *fit_dh+dplus[4]);
    fprintf(reportStream, "\n     reduced chi^2=%f (tot=%f) %d params and %d points\n", *chi/(double)(fitterinfo.NrDataPoints-(*nfitparameters)), *chi, *nfitparameters, fitterinfo.NrDataPoints);
  }
}
void PlotGrid(float *chigrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, double level, double suppress_fac, int GridDeviceID, double alpha, double beta, double lwbox, double labelcharheight, double boxlabelcharheight, int drawCross, int draw_title, int drawcontours, int nogray,
       double chimax, double chimin, int maptype, int showwedge, char *showwedge_label)
{
  int i;
  char txt1[100], txt2[100], txt3[100];
  float contours[200];
  verbose_definition noverbose;
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;
    for(i = 0; i < drawcontours; i++) {
      contours[i] = (chimax - (i+2)*chimin)/(chimax - chimin);
    }
  if(drawcontours > 200) {
    printwarning(noverbose.debug, "Maximum allowed number contours is 200");
    drawcontours = 200;
  }
  if(draw_title) {
    sprintf(txt1, "\\(2148)\\u2\\d grid");
  }else {
    txt1[0] = 0;
  }
  sprintf(txt2, "\\(2128) [deg]");
  sprintf(txt3, "\\(2127) [deg]");
  ppgask(0);
  ppgpage();
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);
  pgplot_options.viewport.dontopen = 1;
  pgplot_options.viewport.dontclose = 1;
  strcpy(pgplot_options.box.xlabel, txt3);
  strcpy(pgplot_options.box.ylabel, txt2);
  strcpy(pgplot_options.box.title, txt1);
  pgplot_options.box.label_ch = labelcharheight;
  pgplot_options.box.box_labelsize = boxlabelcharheight;
  pgplot_options.box.box_lw = lwbox;
  pgplot_options.box.label_lw = lwbox;
  pgplotMap(&pgplot_options, chigrid, nalpha, nbeta, alphastart, alphaend, alphastart, alphaend, betastart, betaend, betastart, betaend, maptype, 0, nogray, drawcontours, contours, lwbox, 0, 1, 1, level, suppress_fac, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, noverbose);
  if(showwedge) {
    ppgsch(boxlabelcharheight*0.5);
    ppgslw(lwbox);
    ppgwedg("RI", 0, 6, (1.0-suppress_fac)*(chimax-chimin)+chimin, (1.0-level)*(chimax-chimin)+chimin, showwedge_label);
    ppgslw(1);
    ppgsch(1);
  }
  if(alpha >= 0 && drawCross) {
    if(drawCross == 1)
      ppgsci(2);
    else
      ppgsci(0);
    ppgmove(alpha-2, beta);
    ppgdraw(alpha+2, beta);
    ppgmove(alpha, beta-2);
    ppgdraw(alpha, beta+2);
    ppgsci(1);
  }
  ppgslw(1);
}
int define_contours_from_list_points(float *alpha, float *beta, int *contnr, int npoints, double maxsep)
{
  int i, j, curcont, found, findnextcontour;
  double sep;
  if(npoints < 1) {
    return 0;
  }
  for(i = 0; i < npoints; i++) {
    contnr[i] = 0;
  }
  curcont = 0;
  do {
    curcont++;
    found = 0;
    findnextcontour = 0;
    for(i = 0; i < npoints; i++) {
      if(contnr[i] == 0) {
 contnr[i] = curcont;
 found = 1;
 break;
      }
    }
    if(found == 0)
      break;
    else
      findnextcontour = 1;
    while(found) {
      found = 0;
      for(i = 0; i < npoints; i++) {
 if(contnr[i] == 0) {
   for(j = 0; j < npoints; j++) {
     if(contnr[j] == curcont) {
       sep = sqrt((alpha[j]-alpha[i])*(alpha[j]-alpha[i])+(beta[j]-beta[i])*(beta[j]-beta[i]));
       if(sep <= maxsep) {
  contnr[i] = curcont;
  found = 1;
  break;
       }
     }
   }
 }
      }
    }
    found = 1;
  }while(findnextcontour);
  return curcont-1;
}
void calcIntersectionRhoAndBanana(float *rhogrid, float *chigrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, int nruserContours, float *userContours, double chimax, double chimin, verbose_definition verbose)
{
  int i, j, k, ialpha, ibeta, ialpha2, ibeta2, nfound, nfoundbest, *contnr, nrcontoursfound;
  float *alpha, *beta, *chi, chibest, maxsep;
  alpha = malloc(2*nalpha*nbeta*sizeof(float));
  beta = malloc(2*nalpha*nbeta*sizeof(float));
  chi = malloc(2*nalpha*nbeta*sizeof(float));
  contnr = malloc(2*nalpha*nbeta*sizeof(int));
  if(alpha == NULL || beta == NULL || chi == NULL || contnr == NULL) {
    printerror(verbose.debug, "ERROR ppolFit: Cannot allocate memory");
    return;
  }
  if(nruserContours > 0) {
    for(i = 0; i < nruserContours; i++) {
      nfound = 0;
      for(ialpha = 0; ialpha < nalpha; ialpha++) {
 for(ibeta = 0; ibeta < nbeta; ibeta++) {
   ialpha2 = ialpha + 1;
   ibeta2 = ibeta;
   if(ialpha2 >= nalpha)
     ialpha2 = ialpha;
   if((rhogrid[nalpha*ibeta+ialpha] < userContours[i] && rhogrid[nalpha*ibeta2+ialpha2] > userContours[i]) || (rhogrid[nalpha*ibeta+ialpha] > userContours[i] && rhogrid[nalpha*ibeta2+ialpha2] < userContours[i])) {
     alpha[nfound] = 0.5*(ialpha+ialpha2)*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
     beta[nfound] = 0.5*(ibeta+ibeta2)*(betaend-betastart)/(double)(nbeta-1)+betastart;
     chi[nfound] = 0.5*(chigrid[nalpha*ibeta+ialpha] + chigrid[nalpha*ibeta2+ialpha2]);
     if(chi[nfound] > (chimax - 4*chimin)/(chimax - chimin))
        nfound++;
   }
   ialpha2 = ialpha;
   ibeta2 = ibeta+1;
   if(ibeta2 >= nbeta)
     ibeta2 = ibeta;
   if((rhogrid[nalpha*ibeta+ialpha] < userContours[i] && rhogrid[nalpha*ibeta2+ialpha2] > userContours[i]) || (rhogrid[nalpha*ibeta+ialpha] > userContours[i] && rhogrid[nalpha*ibeta2+ialpha2] < userContours[i])) {
     alpha[nfound] = 0.5*(ialpha+ialpha2)*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
     beta[nfound] = 0.5*(ibeta+ibeta2)*(betaend-betastart)/(double)(nbeta-1)+betastart;
     chi[nfound] = 0.5*(chigrid[nalpha*ibeta+ialpha] + chigrid[nalpha*ibeta2+ialpha2]);
     if(chi[nfound] > (chimax - 4*chimin)/(chimax - chimin))
        nfound++;
   }
 }
      }
      if(nfound <= 0) {
 printf("No intersections found between 3-sigma banana and contour\n");
      }else {
 maxsep = sqrt(2.0*(alphaend-alphastart)*(alphaend-alphastart)/(double)((nalpha)*(nalpha)) + 2.0*(betaend-betastart)*(betaend-betastart)/(double)((nbeta)*(nbeta)));
 nrcontoursfound = define_contours_from_list_points(alpha, beta, contnr, nfound, maxsep);
 printf("Found %d cross points between contour %d and 3-sigma chi2 contour\n", nrcontoursfound, i+1);
 if(verbose.verbose) {
   for(j = 0; j < nfound; j++) {
     ppgsci(contnr[j]+4);
     ppgpt1(alpha[j], beta[j], -1);
   }
 }
 for(k = 1; k <= nrcontoursfound; k++) {
   nfoundbest = 0;
   chibest = NAN;
   for(j = 0; j < nfound; j++) {
     if(contnr[j] == k) {
       if(chi[j] > chibest || isnan(chibest)) {
  chibest = chi[j];
  nfoundbest = j;
       }
     }
   }
   if(verbose.verbose) {
     ppgsci(8);
     ppgslw(5);
     ppgpt1(alpha[nfoundbest], beta[nfoundbest], 2);
   }
   chibest = 1.0-chibest;
   chibest = chibest*(chimax-chimin) + chimin;
   printf("Cross point %d: alpha=%f beta=%f chi=%f\n", k, alpha[nfoundbest], beta[nfoundbest], chibest);
 }
      }
    }
  }
  free(alpha);
  free(beta);
  free(chi);
  free(contnr);
  ppgsci(1);
  ppgslw(1);
}
void PlotContours(float *rhogrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, int nrlevels, float *TR, int GridDeviceID, int contour_txt, int contourcolor, int fixedContours, int nruserContours, float *userContours, int lwbox, int dotted)
{
  float C[500];
  char txt[100];
  int i;
  if(fixedContours == 0) {
    for(i = 0; i < nrlevels; i++)
      C[i] = 180.0*i/(double)(nrlevels-1.0);
  }else {
    C[0] = 1;
    C[1] = 2;
    C[2] = 3;
    C[3] = 4;
    C[4] = 5;
    C[5] = 10;
    C[6] = 20;
    C[7] = 30;
    C[8] = 40;
    C[9] = 50;
    C[10] = 60;
    C[11] = 70;
    C[12] = 80;
    C[13] = 90;
    C[14] = 100;
    C[15] = 110;
    C[16] = 120;
    C[17] = 130;
    C[18] = 140;
    C[19] = 150;
    C[20] = 160;
    C[21] = 170;
    C[22] = 180;
    nrlevels = 23;
  }
  if(nruserContours > 0) {
    for(i = 0; i < nruserContours; i++)
      C[i] = userContours[i];
    nrlevels = nruserContours;
  }else if( nruserContours == -1){
    C[0] = 5;
    C[1] = 10;
    C[2] = 15;
    C[3] = 20;
    C[4] = 25;
    C[5] = 30;
    nrlevels = 6;
  }
  ppgsci(contourcolor);
  ppgslw(lwbox);
  if(dotted)
    ppgsls(4);
  ppgcont(rhogrid, nalpha, nbeta, 1, nalpha, 1, nbeta, C, -nrlevels, TR);
  ppgsls(1);
  ppgslw(1);
  if(contour_txt) {
    for(i = 0; i < nrlevels; i++) {
      sprintf(txt, "%.0f", C[i]);
      ppgconl(rhogrid, nalpha, nbeta, 1, nalpha, 1, nbeta, C[i], TR, txt, nalpha, 0.01*nalpha);
    }
  }
  ppgsci(1);
}
double dy_180(double y1, double y2)
{
  double dy;
  y1 = derotate_180_double(y1);
  y2 = derotate_180_double(y2);
  dy = fabs(y1-y2);
  if(fabs(y1-y2-180) < dy) {
    y2 += 180.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-180) < dy) {
    y2 += 180.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-180) < dy) {
    y2 += 180.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+180) < dy) {
    y2 -= 180.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+180) < dy) {
    y2 -= 180.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+180) < dy) {
    y2 -= 180.0;
    dy = fabs(y1-y2);
  }
  return dy;
}
double dy_90(double y1, double y2)
{
  double dy;
  y1 = derotate_180_double(y1);
  y2 = derotate_180_double(y2);
  dy = fabs(y1-y2);
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  return dy;
}
void PlotPAswing(double alpha, double beta, double pa0, double l0, int PlotFit, double leftPulseLongitude, double rightPulseLongitude, double dh)
{
  int i;
  double oldpa, newpa;
  ppgpage();
  ppgask(0);
  ppgslw(1);
  ppgsvp(0.1, 0.9, 0.1, 0.9);
  ppgswin(leftPulseLongitude, rightPulseLongitude, 0, 180);
  ppglab("Pulse longitude", "PA", "PA-fit");
  ppgbox("bcnsti",0.0,0,"bcnmsti",0.0,0);
  for(i = 0; i < fitterinfo.NrDataPoints; i++) {
    ppgerr1(6, fitterinfo.data_l[i], derotate_180(fitterinfo.data_pa[i]), fitterinfo.data_dpa[i], 3);
  }
  if(PlotFit) {
    ppgsci(2);
    oldpa = paswing_double(alpha, beta, 0, pa0, l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, dh);
    ppgmove(0, oldpa);
    for(i = 1; i < 3600; i++) {
      newpa = paswing_double(alpha, beta, 0.1*i, pa0, l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, dh);
      if(fabs(newpa-oldpa) < 100)
 ppgdraw(0.1*i, newpa);
      else
 ppgmove(0.1*i, newpa);
      oldpa = newpa;
    }
    if(l0 > 360)
      l0 -= 360;
    if(l0 < 0)
      l0 += 360;
    ppgsci(3);
    ppgmove(l0, 0);
    ppgdraw(l0, 180);
    l0 += 180;
    ppgsci(4);
    if(l0 > 360)
      l0 -= 360;
    ppgmove(l0, 0);
    ppgdraw(l0, 180);
    ppgsci(1);
  }
}
void PrintHelp()
{
  printf("General options:\n");
  printf("  ESC (or q): Exit\n");
  printf("  w:   Write out postscript files or other file\n");
  printf("Fit options:\n");
  printf("  b:   Find and plot best solution\n");
  printf("  e:   Toggle error estimation on best fit\n");
  printf("  S:   Set the number of sigma's used for error calculation\n");
  printf("  G:   Find allowed range of fit parameters by considering the alpha-beta grid search\n");
  printf("Plot refinement options:\n");
  printf("  C:   Toggle showing of the cross of the best fit\n");
  printf("  g:   Toggle drawing chi2 surface in grayscale\n");
  printf("  l:   Set upper chi^2 level (for grayscale)\n");
  printf("  o:   Change output plot parameters (such as line widths)\n");
  printf("  s:   Set suppress factor of grayscale\n");
  printf("  T:   Toggle showing of the title in the chi^2 plot\n");
  printf("  z:   Zoom in on pulse longitude\n");
  printf("Contour options:\n");
  printf("  c:   Plot chi^2 contours\n");
  printf("  f:   Use build-in non-equally spaced specification of contour levels for beam radius\n");
  printf("  F:   Manually specify shown contour levels for beam radius\n");
  printf("  n:   Set nr of equally spaced contour levels for beam radius\n");
  printf("  r:   Toggle plotting contours of the beam radius\n");
  printf("  R:   Input pulse widths and rho values from ascii file which might help to get an idea of the allowed spread in the contours\n");
  printf("  t:   Switch drawing of labels for countours\n");
  printwarning(0, "WARNING: -maxdl is set to %f. Changing this parameter might change the results considerably.\n\n", fitterinfo.max_l0_diff);
}
void print_steepness(double alpha, double beta, double l0, double pa0, int verbose, double dh, double *sina_b)
{
  int i;
  double l, l1, l2, psi, psi_old, dpsi, dpsidphi, max_dpsidphi, resolution;
  *sina_b = sin(alpha*M_PI/180.0)/sin(beta*M_PI/180.0);
  printf("\nsin(a)/sin(b)   = %f\n", sin(alpha*M_PI/180.0)/sin(beta*M_PI/180.0));
  max_dpsidphi = 0;
  l1 = l2 = -1;
  resolution = fitterinfo.data_l[1] - fitterinfo.data_l[0];
  for(i = 1; i < fitterinfo.NrDataPoints; i++) {
    if(fitterinfo.data_l[i] - fitterinfo.data_l[i-1] < resolution)
      resolution = fitterinfo.data_l[i] - fitterinfo.data_l[i-1];
  }
  if(verbose) printf("Resolution = %f degrees\n", resolution);
  for(i = 0; i < fitterinfo.NrDataPoints; i++) {
    l = fitterinfo.data_l[i];
    psi = paswing_double(alpha, beta, l, pa0, l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, dh);
    if(i > 0) {
      dpsi = psi-psi_old;
      dpsidphi = dpsi/(fitterinfo.data_l[i]-fitterinfo.data_l[i-1]);
      if(fabs(dpsidphi) > fabs(max_dpsidphi)) {
 if(fabs(dpsi) > fabs(dpsi-180)) {
   dpsi -= 180;
 }
 if(fabs(dpsi) > fabs(dpsi+180)) {
   dpsi += 180;
 }
 dpsidphi = dpsi/(fitterinfo.data_l[i]-fitterinfo.data_l[i-1]);
 if(fabs(dpsidphi) > fabs(max_dpsidphi) && (fitterinfo.data_l[i] - fitterinfo.data_l[i-1] < 2.5*resolution)) {
   max_dpsidphi = dpsidphi;
   l1 = fitterinfo.data_l[i-1];
   l2 = fitterinfo.data_l[i];
 }
      }
    }
    psi_old = psi;
  }
  printf("measured steepest part = %f (between longitude %f - %f)\n", max_dpsidphi, l1, l2);
}
void convertAlphaBeta(double *alpha, double *beta
        )
{
  double alpha0, beta0;
  alpha0 = *alpha;
  beta0 = *beta;
  *alpha = alpha0;
  *beta = beta0;
}
void calcBeamWidths(int nalpha, int nbeta, double alphastart, double alphaend, double betastart, double betaend,
      int calculate_beam_widths, int calculate_interpulse_widths, double pulse_width2, float *rhogrid, float *rhogrid2, int nocounters)
{
  int i, j;
  double dummyf, chi;
  double alpha0, beta0;
  for(i = 0; i < nalpha; i++) {
    for(j = 0; j < nbeta; j++) {
      alpha0 = i*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
      beta0 = j*(betaend-betastart)/(double)(nbeta-1)+betastart;
      convertAlphaBeta(&alpha0, &beta0
         );
      if(calculate_beam_widths) {
 dummyf = cos(alpha0*M_PI/180.0)*cos((alpha0+beta0)*M_PI/180.0) + sin(alpha0*M_PI/180.0)*sin((alpha0+beta0)*M_PI/180.0)*cos(0.5*fitterinfo.pulse_width*M_PI/180.0);
 chi = acos(dummyf)*180.0/M_PI;
 rhogrid[nalpha*j+i] = chi;
      }
      if(calculate_interpulse_widths == 1) {
 dummyf = -cos(alpha0*M_PI/180.0)*cos((alpha0+beta0)*M_PI/180.0) + sin(alpha0*M_PI/180.0)*sin((alpha0+beta0)*M_PI/180.0)*cos(0.5*pulse_width2*M_PI/180.0);
 chi = acos(dummyf)*180.0/M_PI;
 rhogrid2[nalpha*j+i] = chi;
      }
      if(calculate_interpulse_widths == 2) {
 beta0 = -2.0*alpha0-beta0;
 dummyf = cos(alpha0*M_PI/180.0)*cos((alpha0+beta0)*M_PI/180.0) + sin(alpha0*M_PI/180.0)*sin((alpha0+beta0)*M_PI/180.0)*cos(0.5*pulse_width2*M_PI/180.0);
 chi = acos(dummyf)*180.0/M_PI;
 rhogrid2[nalpha*j+i] = chi;
      }
    }
    if(nocounters == 0)
      fprintf(stderr, "%.1f%%     \r",(100.0*(i+1))/(double)(nalpha));
  }
}

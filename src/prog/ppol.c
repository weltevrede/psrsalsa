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
#include <stdlib.h>
#include "psrsalsa.h"
#define MaxNrJumps 100
int main(int argc, char **argv)
{
  int index, firstfiletoopen, nokeypresses, compute_PA, extendedPol, writeout, writeoutFilename;
  int manualOnpulseSelection, deviceID, normalize, correctLbias, correctPbias, nxsub, boxlinewidth;
  int titleset, titlelw, dashed, datalinewidth, noylabel, overlayPA, overlayFine, nrJumps;
  int outline, outline_color, extended, onlysignificantPA, twoprofiles, sigmaI;
  int show_pa_dist, pa_dist_normalise, elldist, doprojection, spstat;
  int subtract, subtractFile, decompose_method;
  int projection_binnr, projection_weighting, proj_type, projection_nolabels, projection_nrx, projection_nry;
  int sigmaset, padist_weighttype, onpulseonly;
  long i, j, f;
  float *profileI, loffset, correctQV, correctV, sigma_limit;
  float xsize, ysize, ysizepa, xtick, plotI1, plotI2, titlech;
  float plotl1, plotl2, plotp1, plotp2, PAoffset_data;
  float overlayalpha, overlaybeta, overlaypa0, overlayl0;
  float pamask_value, projection_threshold;
  float jump_longitude[MaxNrJumps], jump_offset[MaxNrJumps];
  float *mapProjection, proj_rot_long, proj_rot_lat;
  char *filename_ptr, title[1000], txt[100], PlotDevice[100], PlotDevice2[100], *extension;
  char ofilename[1000];
  datafile_definition datain, dataout, subtract_fin, padist_data;
  psrsalsaApplication application;
  pgplot_options_definition pgplot_options;
  initApplication(&application, "ppol", "[options] inputfile(s)");
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.switch_formatlist = 1;
  application.switch_header = 1;
  application.switch_headerlist = 1;
  application.switch_onpulse = 1;
  application.switch_onpulsef = 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_rebin = 1;
  application.switch_TSCR = 1;
  application.switch_FSCR = 1;
  application.switch_debase = 1;
  application.switch_rotateStokes = 1;
  application.switch_nocounters = 1;
  application.switch_stokes = 1;
  application.switch_deparang = 1;
  application.switch_deFaraday = 1;
  application.switch_filelist = 1;
  application.switch_dedisperse = 1;
  application.switch_deFaraday = 1;
  application.switch_changeRefFreq = 1;
  application.switch_fchan = 1;
  application.switch_history_cmd_only = 1;
  application.switch_nskip = 1;
  application.switch_nread = 1;
  application.switch_forceUniformFreqLabelling = 1;
  application.switch_cmaplist = 1;
  application.switch_cmap = 1;
  application.cmap = -1;
  application.oformat = PPOL_format;
  nokeypresses = 0;
  extendedPol = 0;
  writeout = 0;
  writeoutFilename = -1;
  strcpy(PlotDevice, "?");
  strcpy(PlotDevice2, "?");
  manualOnpulseSelection = 0;
  extension = NULL;
  normalize = 1;
  correctLbias = 1;
  correctPbias = 0;
  correctQV = 1;
  correctV = 1;
  loffset = 0;
  xsize = 0.8;
  ysize = 1.0;
  ysizepa = 1.0;
  xtick = 0.0;
  nxsub = 0;
  sigma_limit = 3.0;
  boxlinewidth = 1;
  titleset = 0;
  title[0] = 0;
  titlech = 1;
  titlelw = 1;
  plotl1 = 0;
  plotl2 = 360;
  plotp1 = -180;
  plotp2 = 180;
  plotI1 = plotI2 = 0;
  PAoffset_data = 0;
  datalinewidth = 1;
  dashed = 0;
  noylabel = 0;
  overlayPA = 0;
  overlayFine = 1;
  outline = -1;
  nrJumps = 0;
  twoprofiles = 0;
  outline_color = 0;
  overlayalpha = 0;
  overlaybeta = 0;
  overlaypa0 = 0;
  overlayl0 = 0;
  sigmaI = 0;
  doprojection = 0;
  spstat = 0;
  subtract = 0;
  subtractFile = 0;
  decompose_method = 0;
  show_pa_dist = 0;
  elldist = 0;
  pa_dist_normalise = 0;
  pamask_value = sqrt(-1);
  projection_nolabels = 0;
  projection_threshold = -1.0;
  sigmaset = 0;
  padist_weighttype = 0;
  onpulseonly = 0;
  if(argc < 2) {
    printf("Program to convert Stokes parameters in other polarization products such as\n");
    printf("postition angle and linear intensity. The results can be written out and plotted\n");
    printf("using ppolFig or fitted using ppolFit.\n\n");
    printApplicationHelp(&application);
    fprintf(stdout, "General input/output options:\n");
    fprintf(stdout, "-ext               Write out polarised profile to file with this extension.\n");
    fprintf(stdout, "-stdout            Like -ofile, but now write to stdout.\n");
    fprintf(stdout, "-ofile             Like -stdout, but now write to the specified file.\n");
    fprintf(stdout, "\nOptions related to generating polarimetric data:\n");
    fprintf(stdout, "-addfile           Add the PA-swing from this datafile to the data.\n");
    fprintf(stdout, "                   This is the opposite from -subtractfile.\n");
    fprintf(stdout, "-loffset           Shift longitudes by this amount.\n");
    fprintf(stdout, "-medianLdebias     Naively subtract the median L of offpulse region rather than\n");
    fprintf(stdout, "                   applying the better de-bias method of Wardle & Kronberg.\n");
    fprintf(stdout, "-noLdebias         Simply use L^2=Q^2+U^2, which is biased.\n");
    fprintf(stdout, "-nonorm            Do not normalize the output profile.\n");
    fprintf(stdout, "-paoffset          Add this angle to PA (in degrees).\n");
    fprintf(stdout, "-selectonpulse     Enables manual graphical selection of more on-pulse regions\n");
    fprintf(stdout, "                   in addition to any provided on the command line.\n");
    fprintf(stdout, "-sigma             Set sigma limit on L required for PA calculation [def=%.1f].\n", sigma_limit);
    fprintf(stdout, "                   Here sigma is defined as the RMS of the not de-biased off-pulse L.\n");
    fprintf(stdout, "-subtractfile      Subtract the PA-swing from this datafile from the data.\n");
    fprintf(stdout, "                   This is the opposite from -addfile.\n");
    fprintf(stdout, "-2                 Write out two pulse periods (so there is duplicated data).\n");
    fprintf(stdout, "-extendedpol       Also generate the total amount of polarization\n");
    fprintf(stdout, "                   and ellipticity angles.\n");
    fprintf(stdout, "-padist N      Calculates the PA distribution with N PA bins.\n");
    fprintf(stdout, "-elldist       When specified combined with -padist, those options\n");
    fprintf(stdout, "               now use the ellipticity angle rather than the PA.\n");
    fprintf(stdout, "-projection    Shows a projection of the Poincare sphere. Provide\n");
    fprintf(stdout, "               \"binnr hsize weighting type dlong dlat\". If binnr is set, only\n");
    fprintf(stdout, "               that pulse longitude bin is considered. If -1 it only uses the\n");
    fprintf(stdout, "               onpulse region and if -2 it will use all bins. Set the width in\n");
    fprintf(stdout, "               pixels of the map with hsize. If weighting is 0, the number of\n");
    fprintf(stdout, "               samples which point in a given direction are shown. If weighting\n");
    fprintf(stdout, "               is 1, pixel value is the sum of sqrt(Q^2+U^2+V^2). If weighting\n");
    fprintf(stdout, "               is 2, this sum is normalised by the sum of Stokes I. If weighting\n");
    fprintf(stdout, "               is 3, then the samples are weighted by Stokes I.\n");
    fprintf(stdout, "               Type 1 is the (equal area) Hammer-Aitoff, type 2 is a normal\n");
    fprintf(stdout, "               spherical projection, 3 is a long/lat map. The orientation of\n");
    fprintf(stdout, "               the projection can be set by dlong and dlat. No S/N limit is\n");
    fprintf(stdout, "               imposed when the input is Stokes data. ppol -extendedpol can\n");
    fprintf(stdout, "               be used first to calculate PA's and ellipticity, which allows\n");
    fprintf(stdout, "               ppol -sigma to be used.\n");
    fprintf(stdout, "\nOptions affecting the plotting:\n");
    fprintf(stdout, "-1             Only plot PA-swing once (equivalent to -yrange \"0 180\").\n");
    fprintf(stdout, "-device        Specify plotting device for onpulse region selection.\n");
    fprintf(stdout, "-device2       Specify plotting device final plot.\n");
    fprintf(stdout, "-projection_nobox  Do not show a rectangular box with numerical labels around the\n");
    fprintf(stdout, "                   plot produced by the -projection option.\n");
    fprintf(stdout, "-xrange        Specify longitude range covered in the plot in degrees.\n");
    fprintf(stdout, "-xrange_phase  Specify longitude range covered in the plot in phase.\n");
    fprintf(stdout, "-yrange        Specify PA-range covered in the plot.\n");
    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about fitting position-angle swings and using the beam-width information can be found in:\n");
    printf(" - 	Rookyard et al. 2015, MNRAS, 446, 3367\n");
    printf("More information about polarization distributions Poincare sphere projections (-projection) can be found in:\n");
    printf(" - 	Ilie et al. 2020, MNRAS, 491, 3385\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else if(argc >= 2) {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
 i = index;
      }else if(strcmp(argv[i], "-device") == 0) {
 strcpy(PlotDevice,argv[i+1]);
        i++;
      }else if(strcmp(argv[i], "-device2") == 0) {
 strcpy(PlotDevice2,argv[i+1]);
        i++;
      }else if(strcmp(argv[i], "-stdout") == 0) {
 writeout = 1;
      }else if(strcmp(argv[i], "-nonorm") == 0) {
 normalize = 0;
      }else if(strcmp(argv[i], "-ofile") == 0) {
 writeout = 1;
 writeoutFilename = i+1;
 i++;
      }else if(strcmp(argv[i], "-ext") == 0) {
 writeout = 1;
 extension = argv[i+1];
        i++;
      }else if(strcmp(argv[i], "-paoffset") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &PAoffset_data, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-loffset") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &loffset, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-nokeypress") == 0) {
 nokeypresses = 1;
      }else if(strcmp(argv[i], "-extendedpol") == 0) {
 extendedPol = 1;
      }else if(strcmp(argv[i], "-projection") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d %d %d %f %f", &projection_binnr, &projection_nrx, &projection_weighting, &proj_type, &proj_rot_long, &proj_rot_lat, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 if(proj_type < 1 || proj_type > 3) {
   printerror(application.verbose_state.debug, "ERROR ppol: In the '%s' option an invalid type was specified.", argv[i]);
   return 0;
 }
 if(projection_weighting < 0 || projection_weighting > 3) {
   printerror(application.verbose_state.debug, "ERROR ppol: In the '%s' option an invalid weighting was specified.", argv[i]);
   return 0;
 }
 projection_nry = 0.5*projection_nrx;
 doprojection = 1;
        i++;
      }else if(strcmp(argv[i], "-projection_nobox") == 0) {
 projection_nolabels = 1;
      }else if(strcmp(argv[i], "-sigma") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &sigma_limit, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 sigmaset = 1;
        i++;
      }else if(strcmp(argv[i], "-padist") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &show_pa_dist, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-elldist") == 0) {
 elldist = 1;
 extendedPol = 1;
      }else if(strcmp(argv[i], "-xrange") == 0 || strcmp(argv[i], "-xrange_phase") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &plotl1, &plotl2, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 if(strcmp(argv[i], "-xrange_phase") == 0) {
   plotl1 *= 360.0;
   plotl2 *= 360.0;
 }
        i++;
      }else if(strcmp(argv[i], "-yrange") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &plotp1, &plotp2, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-yrange2") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &plotI1, &plotI2, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-subtractfile") == 0) {
 subtract = 1;
 subtractFile = i+1;
 i++;
      }else if(strcmp(argv[i], "-addfile") == 0) {
 subtract = 2;
 subtractFile = i+1;
 i++;
      }else if(strcmp(argv[i], "-1") == 0) {
 plotp1 = 0.0;
 plotp2 = 180.0;
      }else if(strcmp(argv[i], "-selectonpulse") == 0) {
 manualOnpulseSelection = 1;
      }else if(strcmp(argv[i], "-2") == 0) {
 twoprofiles = 1;
      }else if(strcmp(argv[i], "-medianLdebias") == 0) {
 correctLbias = 0;
      }else if(strcmp(argv[i], "-noLdebias") == 0) {
 correctLbias = -1;
      }else if(strcmp(argv[i], "-medianPdebias") == 0) {
 correctPbias = 0;
      }else {
 if(argv[i][0] == '-') {
   printerror(application.verbose_state.debug, "ppol: Unknown option: %s", argv[i]);
   printerror(application.verbose_state.debug, "\nRun ppol without command line arguments to show help");
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
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR ppol: No files specified");
    return 0;
  }
  if(doprojection) {
    extendedPol = 0;
    writeout = 0;
  }
  firstfiletoopen = 1;
  while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
    if(firstfiletoopen == 0) {
      if(nokeypresses == 0) {
 i = pgplot_device_type(PlotDevice2, application.verbose_state);
 if(i < 3 || i > 10) {
   printf("Press a key in the terminal to continue\n");
   fflush(stdout);
   pgetch();
 }
      }
    }
    cleanPSRData(&dataout, application.verbose_state);
    if(!openPSRData(&datain, filename_ptr, application.iformat, 0, 1, 0, application.verbose_state))
      return 0;
    if(application.verbose_state.verbose) {
      fflush(stdout);
      printwarning(application.verbose_state.debug, "Input data contains %ld bins, %ld pulses, %ld polarizations and %ld frequencies.", (datain.NrBins), datain.NrSubints, (datain.NrPols), datain.NrFreqChan);
    }
    if(PSRDataHeader_parse_commandline(&datain, argc, argv, application.verbose_state) == 0)
      return 0;
    if(datain.poltype == POLTYPE_ILVPAdPA || datain.poltype == POLTYPE_PAdPA || datain.poltype == POLTYPE_ILVPAdPATEldEl) {
      if(doprojection == 0) {
 printerror(application.verbose_state.debug, "ERROR ppol: File already contains PA data. You can use ppolFig to plot this data.");
 return 0;
      }else {
 if(datain.poltype != POLTYPE_ILVPAdPATEldEl) {
   printerror(application.verbose_state.debug, "ERROR ppol: The -projection option expects either IQUV data, or polarization information including a PA and ellipticity (i.e. generated with ppol -extendedpol).");
   return 0;
 }
      }
    }
    compute_PA = 1;
    if(datain.isFolded && datain.foldMode == FOLDMODE_FIXEDPERIOD) {
      if(datain.fixedPeriod <= 0.0) {
 printwarning(application.verbose_state.debug, "WARNING ppol: Period does not appear to be set, assuming it is 1 sec.");
 datain.fixedPeriod = 1.0;
      }
    }
    if(datain.isFolded && datain.tsampMode == TSAMPMODE_FIXEDTSAMP) {
      if(get_tsamp(datain, 0, application.verbose_state) <= 0.0) {
 printwarning(application.verbose_state.debug, "WARNING ppol: Assuming full period is stored.");
 double period;
 int ret;
 ret = get_period(datain, 0, &period, application.verbose_state);
 if(ret == 2) {
   printerror(application.verbose_state.debug, "ERROR ppol (%s): Cannot obtain period", datain.filename);
   return 0;
 }
 datain.fixedtsamp = period/(double)datain.NrBins;
      }
    }
    if(get_tsamp(datain, 0, application.verbose_state) < 0.0000001 || get_tsamp(datain, 0, application.verbose_state) >= 10) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR ppol: The sampling time does not appear to be set correctly in the header. Consider using the -header option.");
      return 0;
    }
    for(i = 1; i < argc; i++) {
      if(strcmp(argv[i], "-header") == 0) {
 fflush(stdout);
 printwarning(application.verbose_state.debug, "WARNING ppol: If using the -header option, be aware it applied BEFORE the preprocessing.");
 break;
      }
    }
    if(preprocessApplication(&application, &datain) == 0) {
      fflush(stdout);
      printwarning(application.verbose_state.debug, "WARNING ppol: Applying preprocess options failed.");
      return 0;
    }
    if(!(doprojection && datain.poltype == POLTYPE_ILVPAdPATEldEl)) {
      if(datain.poltype == POLTYPE_COHERENCY) {
 if(application.verbose_state.verbose) {
   printf("Data is written as coherency parameters. Going to convert them into Stokes parameters first.\n");
 }
 if(preprocess_stokes(&datain, application.verbose_state) == 0) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Conversion into Stokes parameters failed.");
   return 0;
 }
      }else if(datain.poltype == POLTYPE_UNKNOWN) {
 fflush(stdout);
 printwarning(application.verbose_state.debug, "WARNING ppol: It is assumed that the data are Stokes parameters.");
      }else if(datain.poltype != POLTYPE_STOKES) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "ERROR ppol: The polarization type of the data cannot be handled by ppol.");
 return 0;
      }
      if(datain.NrPols != 4) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "ERROR ppol: Need Stokes IQUV data, but the number of polarization channels != 4!");
 return 0;
      }
    }
    if(writeout) {
      if(datain.NrSubints > 1 && spstat == 0) {
 if(application.oformat == PPOL_format) {
   printf("Data contains multiple subints. Data will be written out in FITS format.\n");
   application.oformat = FITS_format;
 }
      }
      if(datain.NrFreqChan > 1 && spstat == 0) {
 if(application.oformat == PPOL_format) {
   printf("Data contains multiple frequency channels. Data will be written out in FITS format.\n");
   application.oformat = FITS_format;
 }
      }
    }
      region_frac_to_int(&(application.onpulse), datain.NrBins, 0);
    long total_nr_offpulse_file_bins;
    if(doprojection == 0 || (doprojection == 1 && projection_binnr == -1)) {
      strcpy(txt, "Select on-pulse region ");
      strcat(txt, datain.psrname);
 total_nr_offpulse_file_bins = datain.NrBins;
      profileI = (float *)malloc(total_nr_offpulse_file_bins*sizeof(float));
      if(profileI == NULL) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "ERROR ppol: Memory allocation error.");
 return 0;
      }
      int ret;
 ret = read_profilePSRData(datain, profileI, NULL, 0, application.verbose_state);
      if(ret != 1) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "ERROR ppol: Cannot form pulse profile");
 return 0;
      }
      if(manualOnpulseSelection == 1 || (application.onpulse.nrRegions == 0 && manualOnpulseSelection != -1)) {
 pgplot_clear_options(&pgplot_options);
 strcpy(pgplot_options.viewport.plotDevice, PlotDevice);
 strcpy(pgplot_options.box.xlabel, "Bin");
 strcpy(pgplot_options.box.ylabel, "Intensity");
 strcpy(pgplot_options.box.title, txt);
 selectRegions(profileI, total_nr_offpulse_file_bins, &pgplot_options, 0, 0, 0, &(application.onpulse), application.verbose_state);
 if(firstfiletoopen == 0) {
   ppgslct(deviceID);
 }
      }else {
 pgplot_clear_options(&pgplot_options);
 strcpy(pgplot_options.box.xlabel, "Bin");
 strcpy(pgplot_options.box.ylabel, "Intensity");
 strcpy(pgplot_options.box.title, txt);
 strcpy(pgplot_options.viewport.plotDevice, PlotDevice);
 pgplotGraph1(&pgplot_options, profileI, NULL, NULL, total_nr_offpulse_file_bins, 0, total_nr_offpulse_file_bins, 0, 0, total_nr_offpulse_file_bins, 0, 0, 0, 1, 0, 1, 0, 1, 1, &(application.onpulse), -1, application.verbose_state);
 if(firstfiletoopen == 0) {
   ppgslct(deviceID);
 }
      }
      region_int_to_frac(&(application.onpulse), 1.0/(float)total_nr_offpulse_file_bins, 0);
      regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);
    }
    if(decompose_method == 0) {
      if(subtractFile) {
 if(!openPSRData(&subtract_fin, argv[subtractFile], 0, 0, 1, 0, application.verbose_state))
   return 0;
 if(subtract_fin.NrBins != datain.NrBins) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: The reference PA-swing has a different number of pulse phase bins compared to the data file.");
   return 0;
 }
 if(subtract_fin.NrFreqChan != 1) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: The reference PA-swing should have a single frequency channel.");
   return 0;
 }
      }
    }
    if(doprojection == 0) {
      verbose_definition verbose2;
      copyVerboseState(application.verbose_state, &verbose2);
      if(application.verbose_state.verbose && datain.NrSubints == 1 && datain.NrFreqChan == 1)
 verbose2.verbose = 1;
      else
 verbose2.verbose = 0;
      int nolongitudes;
      if((datain.NrFreqChan > 1 || datain.NrSubints > 1) && spstat == 0) {
 nolongitudes = 1;
      }else {
 nolongitudes = 0;
      }
 if(make_paswing_fromIQUV(&datain, extendedPol, spstat, sigma_limit, sigmaI, application.onpulse, normalize, correctLbias, correctPbias, correctQV, correctV, nolongitudes, loffset, PAoffset_data, NULL, 1.0, onpulseonly, verbose2) == 0) {
   return 0;
 }
    }else {
      mapProjection = malloc(projection_nrx*projection_nry*sizeof(float));
      if(mapProjection == NULL) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "ERROR ppol: Cannot allocate memory");
 return 0;
      }
      if(sigmaset) {
 if(datain.poltype == POLTYPE_ILVPAdPATEldEl) {
   printwarning(application.verbose_state.debug, "WARNING ppol: A sigma limit is set by the user, but since the provided data has already a sigma limit applied when calculating the PA and ellipticity, it has no effect.");
 }else {
   printwarning(application.verbose_state.debug, "WARNING ppol: A sigma limit is set by the user, but isn't used when producing a Poincare sphere from Stokes parameters. Consider generating PA's and ellipticities first with ppol -extendedpol after assigning a sigma limit. The generated output file can then be read in with a separate ppol -projection command to generated a Poincare sphere.");
 }
      }
      if(subtractFile) {
 if(subtract == 2) {
   printerror(application.verbose_state.debug, "ERROR ppol: The -addfile option cannot be used together with the -projection option.");
   return 0;
 }
 if(make_polarization_projection_map(datain, mapProjection, projection_nrx, projection_nry, 0, projection_binnr, application.onpulse, projection_weighting, projection_threshold, proj_type, proj_rot_long+2.0*(PAoffset_data), proj_rot_lat, -1.0, &subtract_fin, application.verbose_state) == 0)
   return 0;
      }else {
 if(make_polarization_projection_map(datain, mapProjection, projection_nrx, projection_nry, 0, projection_binnr, application.onpulse, projection_weighting, projection_threshold, proj_type, proj_rot_long+2.0*(PAoffset_data), proj_rot_lat, -1.0, NULL, application.verbose_state) == 0) {
   return 0;
 }
      }
    }
    if(doprojection == 0 && decompose_method == 0) {
      if(subtract && subtractFile) {
 if(subtract == 2) {
   if(paswing_remove_observed_PA_swing(&datain, subtract_fin, 1, application.verbose_state) == 0) {
     printerror(application.verbose_state.debug, "ERROR ppol: Subtraction of PA-swing failed.");
     return 0;
   }
 }else {
   if(paswing_remove_observed_PA_swing(&datain, subtract_fin, 0, application.verbose_state) == 0) {
     printerror(application.verbose_state.debug, "ERROR ppol: Subtraction of PA-swing failed.");
     return 0;
   }
 }
      }
      if(datain.NrSubints > 1 || datain.NrFreqChan > 1) {
 fflush(stdout);
 printwarning(application.verbose_state.debug, "WARNING ppol: Only showing polarization of the first subint/frequency channel");
      }
      pgplot_clear_options(&pgplot_options);
      pgplot_options.viewport.dxplot = 0;
      pgplot_options.viewport.dyplot = 0;
      pgplot_options.viewport.xsize = xsize/0.8;
      pgplot_options.viewport.ysize = ysize;
      ppgpage();
      strcpy(pgplot_options.viewport.plotDevice, PlotDevice2);
      pgplot_options.viewport.noclear = 1;
      pgplot_options.viewport.dontclose = 1;
      pgplot_options.viewport.dontopen = !firstfiletoopen;
      pgplot_options.box.box_lw = boxlinewidth;
      pgplot_options.box.label_lw = boxlinewidth;
      pgplot_options.box.box_xtick = xtick;
      pgplot_options.box.box_nxsub = nxsub;
      if(titleset)
 strcpy(pgplot_options.box.title, title);
      else
 strcpy(pgplot_options.box.title, datain.filename);
      pgplot_options.box.title_ch = titlech;
      pgplot_options.box.title_lw = titlelw;
      pgplotPAplot(datain, extendedPol, 0, extendedPol, &pgplot_options, "Pulse longitude [deg]", "I,Linear,V", "PA [deg]", "\\gx [deg]", "Fraction", plotl1, plotl2, 0, 0.0, plotI1, plotI2, plotp1, plotp2, 0.0, sigma_limit, datalinewidth, ysizepa, dashed, noylabel, "-text", "-text_pa", 0.0, 0, "-text_ell", 0.0, 0, "-text_padist", "-text_elldist", "-text_spstat", "-herrorbar", "-herrorbarpa", "-herrorbar2", NULL, NULL, NULL, "-verrorbar", "-verrorbarpa", "-verrorbar2", NULL, argc, argv, outline, outline, outline_color, overlayPA, overlayalpha, overlaybeta, overlaypa0, overlayl0, overlayFine, nrJumps, jump_longitude, jump_offset, NULL, -90, 90, 1.0, 0, 0, NULL, 1.0, 0, NULL, application.verbose_state);
      if(firstfiletoopen) {
 ppgqid(&deviceID);
      }
    }else if(decompose_method) {
    }else {
      pgplot_clear_options(&pgplot_options);
      strcpy(pgplot_options.viewport.plotDevice, PlotDevice2);
      pgplot_options.viewport.dontclose = 1;
      pgplot_options.viewport.dontopen = 1;
      ppgopen(pgplot_options.viewport.plotDevice);
      int cmaptype;
      if(application.cmap < 0) {
 cmaptype = pgplot_device_type(NULL, application.verbose_state);
 if(cmaptype <= 2)
   cmaptype = PPGPLOT_HEAT;
 else
   cmaptype = PPGPLOT_INVERTED_HEAT4;
      }else {
 cmaptype = application.cmap;
      }
      if(projection_nolabels) {
 pgplot_options.box.drawbox = 0;
      }else {
 pgplot_options.box.drawbox = 1;
      }
      if(titleset) {
 strcpy(pgplot_options.box.title, title);
      }else {
 strcpy(pgplot_options.box.title, datain.filename);
      }
      pgplot_options.box.title_ch = titlech;
      pgplot_options.box.title_lw = titlelw;
      pgplot_options.box.box_lw = boxlinewidth;
      pgplot_options.box.label_lw = boxlinewidth;
      pgplot_options.box.box_xtick = xtick;
      pgplot_options.box.box_nxsub = nxsub;
      if(proj_type == 3) {
 if(pgplotMap(&pgplot_options, mapProjection, projection_nrx, projection_nry, -180, 180, -180, 180, -90, 90, -90, 90, cmaptype, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, application.verbose_state) != 1) {
   return 0;
 }
      }else {
 if(pgplotMap(&pgplot_options, mapProjection, projection_nrx, projection_nry, -2.25, 2.25, -2.25, 2.25, -1.125, 1.125, -1.125, 1.125, cmaptype, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, application.verbose_state) != 1) {
   return 0;
 }
      }
      drawSphericalGrid(30, 30, proj_rot_long, proj_rot_lat, 3, proj_type);
      int outline_txt = 0;
      if(outline > 0)
 outline_txt = 1;
      if(pgplot_process_text_option("-text", outline_txt, outline, outline_color, argc, argv, 0, datain, application.verbose_state) == 0) {
 printerror(application.verbose_state.debug, "ERROR ppol: Error processing %s option.", "-text");
 return 0;
      }
    }
    if(show_pa_dist && doprojection == 0) {
      if(nokeypresses == 0) {
 printf("Press a key in the terminal key to continue\n");
 pgetch();
      }
      pgplot_clear_options(&pgplot_options);
      pgplot_options.viewport.dontopen = 1;
      pgplot_options.viewport.noclear = 0;
      strcpy(pgplot_options.box.xlabel, "Pulse longitude bin");
      if(elldist == 0)
 strcpy(pgplot_options.box.ylabel, "PA [deg]");
      else
 strcpy(pgplot_options.box.ylabel, "\\gx [deg]");
      if(titleset)
 strcpy(pgplot_options.box.title, title);
      else
 strcpy(pgplot_options.box.title, datain.filename);
      datafile_definition *pamask_ptr;
      pamask_ptr = NULL;
      if(make_pa_distribution(datain, &padist_data, show_pa_dist, pa_dist_normalise, padist_weighttype, pamask_ptr, pamask_value, elldist, application.verbose_state) == 0)
 return 0;
      int ret;
      int cmaptype;
      if(application.cmap < 0) {
 cmaptype = PPGPLOT_GRAYSCALE;
      }else {
 cmaptype = application.cmap;
      }
      if(elldist) {
 ret = pgplotMap(&pgplot_options, padist_data.data, padist_data.NrBins, padist_data.NrSubints, 0.5, padist_data.NrBins-0.5, padist_data.NrBins*plotl1/360.0, padist_data.NrBins*plotl2/360.0, -45+0.5*(90.0/(float)show_pa_dist), 45-0.5*(90.0/(float)show_pa_dist), -45, 45, cmaptype, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, application.verbose_state);
      }else {
 ret = pgplotMap(&pgplot_options, padist_data.data, padist_data.NrBins, padist_data.NrSubints, 0.5, padist_data.NrBins-0.5, padist_data.NrBins*plotl1/360.0, padist_data.NrBins*plotl2/360.0, -90+0.5*(180.0/(float)show_pa_dist), 90-0.5*(180.0/(float)show_pa_dist), -90, 90, cmaptype, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, application.verbose_state);
      }
      if(ret == 0) {
 return 0;
      }
      if(writeout
  ) {
 if(extension != NULL) {
   if(change_filename_extension(filename_ptr, ofilename, extension, MaxFilenameLength, application.verbose_state) == 0) {
     fflush(stdout);
     printerror(application.verbose_state.debug, "ERROR ppol: Changing filename failed");
     return 0;
   }
 }
 if(writeoutFilename != -1) {
   strcpy(ofilename, argv[writeoutFilename]);
 }
 if(writeoutFilename != -1 || extension != NULL) {
   int oformat;
   oformat = application.oformat;
   if(oformat == PPOL_format) {
     oformat = FITS_format;
   }
   if(!openPSRData(&padist_data, ofilename, oformat, 1, 0, 0, application.verbose_state))
     return 0;
   if(!writeHeaderPSRData(&padist_data, argc, argv, application.history_cmd_only, NULL, application.verbose_state))
     return 0;
   if(writePSRData(&padist_data, padist_data.data, application.verbose_state) != 1)
     return 0;
 }
      }
      closePSRData(&padist_data, 0, 0, application.verbose_state);
    }
    if(subtractFile) {
      closePSRData(&subtract_fin, 0, 0, application.verbose_state);
    }
    if(writeout && (show_pa_dist == 0
      ) && doprojection == 0) {
      copy_params_PSRData(datain, &dataout, application.verbose_state);
      if(extension != NULL) {
 if(change_filename_extension(filename_ptr, ofilename, extension, MaxFilenameLength, application.verbose_state) == 0) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Changing filename failed");
   return 0;
 }
      }
      if(writeoutFilename != -1) {
 strcpy(ofilename, argv[writeoutFilename]);
      }
      if(writeoutFilename != -1 || extension != NULL) {
 if(!openPSRData(&dataout, ofilename, application.oformat, 1, 0, 0, application.verbose_state))
   return 0;
      }else {
 if(datain.NrSubints > 1 || datain.NrFreqChan > 1) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot write data to stdout when there are multiple subints and/or frequency channels");
   return 0;
 }
 dataout.fptr = stdout;
 dataout.fptr_hdr = stdout;
 dataout.format = application.oformat;
      }
      if(!writeHeaderPSRData(&dataout, argc, argv, application.history_cmd_only, NULL, application.verbose_state))
 return 0;
      if(application.oformat == PPOL_format) {
 extended = 1;
 onlysignificantPA = 0;
      }else {
 extended = 0;
 onlysignificantPA = 1;
      }
      if(dataout.format == PPOL_format || dataout.format == PPOL_SHORT_format) {
 dataout.offpulse_rms = datain.offpulse_rms;
 if(writePPOLfile(dataout, datain.data, extended, onlysignificantPA, twoprofiles, 0.0, application.verbose_state) == 0) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot write data");
   return 0;
 }
      }else {
 for(j = 0; j < datain.NrPols; j++) {
   for(i = 0; i < datain.NrSubints; i++) {
     for(f = 0; f < datain.NrFreqChan; f++) {
       if(writePulsePSRData(&dataout, i, j, f, 0, datain.NrBins, &datain.data[datain.NrBins*(j+datain.NrPols*(f+i*datain.NrFreqChan))], application.verbose_state) != 1) {
  fflush(stdout);
  printerror(application.verbose_state.debug, "ERROR ppol: Cannot write data");
  return 0;
       }
     }
   }
 }
      }
    }
    closePSRData(&datain, 0, 0, application.verbose_state);
    if(dataout.format == PPOL_format || dataout.format == PPOL_SHORT_format) {
      dataout.offpulse_rms = NULL;
    }
    closePSRData(&dataout, 0, 0, application.verbose_state);
    firstfiletoopen = 0;
    if(compute_PA) {
      free(profileI);
    }
    if(doprojection) {
      free(mapProjection);
    }
  }
  ppgend();
  terminateApplication(&application);
  return 0;
}

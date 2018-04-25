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
  int manualOnpulseSelection, deviceID, normalize, correctLbias, nxsub, boxlinewidth;
  int titleset, titlelw, dashed, datalinewidth, noylabel, overlayPA, overlayFine, nrJumps;
  int outline, outline_color, extended, onlysignificantPA, twoprofiles;
  long i, j, f;
  float *profileI, loffset, correctQV, correctV, sigma_limit;
  float xsize, ysize, ysizepa, xtick, plotI1, plotI2, titlech;
  float plotl1, plotl2, plotp1, plotp2, PAoffset_data;
  float overlayalpha, overlaybeta, overlaypa0, overlayl0;
  float jump_longitude[MaxNrJumps], jump_offset[MaxNrJumps];
  char *filename_ptr, title[1000], txt[100], PlotDevice[100], PlotDevice2[100], *extension;
  char ofilename[1000];
  datafile_definition datain, dataout;
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
  if(argc < 2) {
    printf("Program to convert Stokes parameters in other polarization products such as\n");
    printf("postition angle and linear intensity. The results can be written out and plotted\n");
    printf("using ppolFig or fitted using ppolFit.\n\n");
    printApplicationHelp(&application);
    fprintf(stdout, "General input/output options:\n");
    fprintf(stdout, "-ext               Write out polarised profile to file with this extension.\n");
    fprintf(stdout, "-stdout            Like -ofile, but now write to stdout.\n");
    fprintf(stdout, "-ofile             Like -stdout, but now write to the specified file.\n");
    fprintf(stdout, "\nOptions related to generating polimetric data:\n");
    fprintf(stdout, "-loffset           Shift longitudes by this amount.\n");
    fprintf(stdout, "-medianLdebias     Naively subtract the median L of offpulse region rather than\n");
    fprintf(stdout, "                   applying the better de-bias method of Wardle & Kronberg.\n");
    fprintf(stdout, "-noLdebias         Simply use L^2=Q^2+U^2, which is biased.\n");
    fprintf(stdout, "-nonorm            Do not normalize the output profile.\n");
    fprintf(stdout, "-paoffset          Add this angle to PA (in degrees).\n");
    fprintf(stdout, "-selectonpulse     Enables manual graphical selection of more on-pulse regions\n");
    fprintf(stdout, "                   in addition to any provided on the command line.\n");
    fprintf(stdout, "-sigma             Set sigma limit on L required for PA calculation [def=%.1f].\n", sigma_limit);
    fprintf(stdout, "-2                 Write out two pulse periods (so there is duplicated data).\n");
    fprintf(stdout, "-extendedpol       Also generate the total amount of polarization\n");
    fprintf(stdout, "\nOptions affecting the plotting:\n");
    fprintf(stdout, "-1             Only plot PA-swing once (equivalent to -yrange \"0 180\").\n");
    fprintf(stdout, "-device        Specify plotting device for onpulse region selection.\n");
    fprintf(stdout, "-device2       Specify plotting device final plot.\n");
    fprintf(stdout, "-xrange        Specify longitude range covered in the plot in degrees.\n");
    fprintf(stdout, "-xrange_phase  Specify longitude range covered in the plot in phase.\n");
    fprintf(stdout, "-yrange        Specify PA-range covered in the plot.\n");
    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about fitting position-angle swings and using the beam-width information can be found in:\n");
    printf(" - 	Rookyard et al. 2015, MNRAS, 446, 3367\n\n");
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
      }else if(strcmp(argv[i], "-extendedpol") == 0) {
 extendedPol = 1;
      }else if(strcmp(argv[i], "-sigma") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &sigma_limit, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
        i++;
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
      printerror(application.verbose_state.debug, "ERROR ppol: File already contains PA data. You can use ppolFig to plot this data.");
      return 0;
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
    if(datain.NrSubints > 1) {
      if(application.oformat == PPOL_format) {
 printf("Data contains multiple subints. Data will be written out in FITS format.\n");
 application.oformat = FITS_format;
      }
    }
    if(datain.NrFreqChan > 1) {
      if(application.oformat == PPOL_format) {
 printf("Data contains multiple frequency channels. Data will be written out in FITS format.\n");
 application.oformat = FITS_format;
      }
    }
    if(datain.NrPols != 4) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR ppol: Need Stokes IQUV data, but the number of polarization channels != 4!");
      return 0;
    }
      region_frac_to_int(&(application.onpulse), datain.NrBins, 0);
    long total_nr_offpulse_file_bins;
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
 pgplotGraph1(&pgplot_options, profileI, NULL, NULL, total_nr_offpulse_file_bins, 0, total_nr_offpulse_file_bins, 0, 0, total_nr_offpulse_file_bins, 0, 0, 0, 1, 0, 0, 1, 1, &(application.onpulse), application.verbose_state);
 if(firstfiletoopen == 0) {
   ppgslct(deviceID);
 }
      }
      region_int_to_frac(&(application.onpulse), 1.0/(float)total_nr_offpulse_file_bins, 0);
      regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);
      verbose_definition verbose2;
      copyVerboseState(application.verbose_state, &verbose2);
      if(application.verbose_state.verbose && datain.NrSubints == 1 && datain.NrFreqChan == 1)
 verbose2.verbose = 1;
      else
 verbose2.verbose = 0;
      int nolongitudes;
      if(datain.NrFreqChan > 1 || datain.NrSubints > 1) {
 nolongitudes = 1;
      }else {
 nolongitudes = 0;
      }
 if(make_paswing_fromIQUV(&datain, extendedPol, application.onpulse, normalize, correctLbias, correctQV, correctV, nolongitudes, loffset, PAoffset_data, NULL, 1.0, verbose2) == 0)
   return 0;
      if(sigma_limit > 0
  ) {
 int pachannel;
 if(datain.poltype == POLTYPE_ILVPAdPATEldEl) {
   pachannel = 3;
 }else {
   pachannel = datain.NrPols-2;
 }
 for(i = 0; i < datain.NrSubints; i++) {
   for(f = 0; f < datain.NrFreqChan; f++) {
     for(j = 0; j < datain.NrBins; j++) {
       if(sigma_limit > 0) {
  if(datain.offpulse_rms == NULL) {
    if(i == 0 && f == 0 && j == 0) {
      fflush(stdout);
      printwarning(application.verbose_state.debug, "WARNING ppol: Cannot apply sigma limit on PA-points for this type of input file. The used sigma limit will be the same as when the input file was generated.");
    }
  }else {
    if(datain.data[j+datain.NrBins*(1+datain.NrPols*(f+i*datain.NrFreqChan))] < sigma_limit*datain.offpulse_rms[1+datain.NrPols*(f + datain.NrFreqChan*i)]) {
      datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] = 0;
      datain.data[j+datain.NrBins*(pachannel+1 + datain.NrPols*(f+datain.NrFreqChan*i))] = -1;
    }
  }
       }
     }
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
      pgplotPAplot(datain, extendedPol, 0, 0
     , &pgplot_options, "Pulse longitude [deg]", "I,Linear,V", "PA [deg]", "\\gx [deg]", plotl1, plotl2, plotI1, plotI2, plotp1, plotp2, 0.0, sigma_limit, datalinewidth, ysizepa, dashed, noylabel, "-text", "-herrorbar", "-herrorbar2", "-verrorbar", "-verrorbar2", argc, argv, outline, outline, outline_color, overlayPA, overlayalpha, overlaybeta, overlaypa0, overlayl0, overlayFine, nrJumps, jump_longitude, jump_offset, NULL, NULL, application.verbose_state);
      if(firstfiletoopen) {
 ppgqid(&deviceID);
      }
    if(writeout
) {
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
      if(!writeHeaderPSRData(&dataout, argc, argv, application.history_cmd_only, application.verbose_state))
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
    closePSRData(&datain, 0, application.verbose_state);
    if(dataout.format == PPOL_format || dataout.format == PPOL_SHORT_format) {
      dataout.offpulse_rms = NULL;
    }
    closePSRData(&dataout, 0, application.verbose_state);
    firstfiletoopen = 0;
    if(compute_PA)
      free(profileI);
  }
  ppgend();
  terminateApplication(&application);
  return 0;
}

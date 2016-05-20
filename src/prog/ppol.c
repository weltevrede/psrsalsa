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
  int index, firstfiletoopen, curpanelnrx, curpanelnry, nokeypresses, compute_PA, showtotPol, writeout, writeoutFilename;
  int manualOnpulseSelection, deviceID, normalize, correctLbias, nrpanelsx, nrpanelsy, nxsub, boxlinewidth;
  int titleset, titlelw, dashed, datalinewidth, noylabel, overlayPA, overlayFine, nrJumps;
  int outline, outline_color, extended, onlysignificantPA, twoprofiles;
  long i, j, f;
  float *profileI, *Ppulse, loffset, correctQV, correctV, sigma_limit;
  float xsize, ysize, ysizepa, xtick, plotI1, plotI2, titlech;
  float plotl1, plotl2, plotp1, plotp2, PAoffset;
  float overlayalpha, overlaybeta, overlaypa0, overlayl0;
  float jump_longitude[MaxNrJumps], jump_offset[MaxNrJumps];
  char *filename_ptr, title[1000], txt[100], PlotDevice[100], PlotDevice2[100], *extension;
  char ofilename[MaxOutputNameLength];
  pgplot_viewport_def viewport;
  pgplot_box_def pgplotbox;
  datafile_definition datain, dataout;
  psrsalsaApplication application;
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
  application.switch_rotateQU = 1;
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
  application.oformat = PPOL_format;
  curpanelnrx = 0;
  curpanelnry = 0;
  nokeypresses = 0;
  showtotPol = 0;
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
  nrpanelsx = 1;
  nrpanelsy = 1;
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
  PAoffset = 0;
  datalinewidth = 1;
  dashed = 0;
  noylabel = 0;
  overlayPA = 0;
  overlayFine = 1;
  outline = -1;
  nrJumps = 0;
  twoprofiles = 0;
  if(argc < 2) {
    printf("Program to convert Stokes parameters in postition angle and linear intensity.\n");
    printf("The result can be plotted and/or written out.\n\n");
    printApplicationHelp(application);
    fprintf(stdout, "Other options affecting the output file:\n");

    fprintf(stdout, "-ext               Write out polarised profile to file with this extension.\n");
    fprintf(stdout, "-stdout            Like -ofile, but now write to stdout.\n");
    fprintf(stdout, "-ofile             Like -stdout, but now write to the specified file.\n");
    fprintf(stdout, "-loffset           Shift longitudes by this amount.\n");
    fprintf(stdout, "-medianLdebias     Naively subtract the median L of offpulse region rather than\n");
    fprintf(stdout, "                   applying the better de-bias method of Wardle & Kronberg.\n");
    fprintf(stdout, "-noLdebias         Simply use L^2=Q^2+U^2, which is biased.\n");
    fprintf(stdout, "-nonorm            Do not normalize the output profile.\n");
    fprintf(stdout, "-paoffset          Add this angle to PA.\n");
    fprintf(stdout, "-selectonpulse     Enables manual graphical selection of more on-pulse regions\n");
    fprintf(stdout, "                   in addition to any provided on the command line.\n");
    fprintf(stdout, "-sigma             Set sigma limit on L required for PA calculation [def=%.1f].\n", sigma_limit);
    fprintf(stdout, "-2                 Write out two pulse periods (so there is duplicated data).\n");
    fprintf(stdout, "\nOther options affecting the plotting:\n");
    fprintf(stdout, "-1             Only plot PA-swing once (equivalent to -yrange \"0 180\").\n");
    fprintf(stdout, "-boxlw         Set linewidth of the box [def=%d].\n", boxlinewidth);
    fprintf(stdout, "-dash          Make L and V profiles dashed/dotted.\n");
    fprintf(stdout, "-device        Specify plotting device for onpulse region selection.\n");
    fprintf(stdout, "-device2       Specify plotting device final plot.\n");
    fprintf(stdout, "-herrorbar     \"xleft xcentre xright y SizeOfMarkers lineWidth colourIndex\"\n");
    fprintf(stdout, "               Draw a horizontal errorbar. The position y is normalised\n");
    fprintf(stdout, "               between 0 and 1.\n");
    fprintf(stdout, "-herrorbar2    Identical to to -herrorbar option, but now for the PA panel,\n");
    fprintf(stdout, "               except that y is in degrees rather than being normalized.\n");
    fprintf(stdout, "-lw            Set linewidth used to plot data [def=%d].\n", datalinewidth);
    fprintf(stdout, "-N             \"nrx nry\" Create nrx by nry panels, rather than a single plot\n");
    fprintf(stdout, "               per page\n");
    fprintf(stdout, "-nokeypress    Do not wait for key presses to go to next page in output plot\n");
    fprintf(stdout, "               (by default this shouldn't happen if plotting to a file)\n");
    fprintf(stdout, "-noylabel      Disable labeling and numbers along y-axis.\n");
    fprintf(stdout, "-opm           Put OPM at this longitude and with this amount of degrees when\n");
    fprintf(stdout, "               using -paswing. You can use this option multiple times.\n");
    fprintf(stdout, "-outline       \"lw color\" Makes text of -text option appear in outline\n");
    fprintf(stdout, "-paswing       \"alpha beta pa0 l0\" Overlay PA-swing with these parameters.\n");
    fprintf(stdout, "-text          \"x y ch lw font color\" \"text\" plots text at this position in\n");
    fprintf(stdout, "               the top panel. See -textkeywords for possible keywords to use.\n");
    fprintf(stdout, "-textkeywords  Lists to keywords you can use in the -text and -title options.\n");
    fprintf(stdout, "-title \"...\"   Set title (default is file name). The title supports keywords\n");
    fprintf(stdout, "               listed by the -textkeywords option.\n");
    fprintf(stdout, "-title_fmt     \"characterheight line_width\" (default being \"%.1f %d\").\n", titlech, titlelw);
    fprintf(stdout, "-totpol        Also shows the total amount of polarization.\n");
    fprintf(stdout, "-xrange        Specify longitude range covered in the plot in degrees.\n");
    fprintf(stdout, "-xrange_phase  Specify longitude range covered in the plot in phase.\n");
    fprintf(stdout, "-yrange        Specify PA-range covered in the plot.\n");
    fprintf(stdout, "-xsize         Set xsize of the plot [def=%.1f].\n", xsize);
    fprintf(stdout, "-ysize         Set ysize of the plot [def=%.1f].\n", ysize);
    fprintf(stdout, "-ysizepa       Set relative ysize of PA plot [def=%.1f].\n", ysizepa);
    fprintf(stdout, "-xticks        \"XTICK NXSUB\"  Adjust PGPLOT tickmarks [def=\"%.0f %d\"].\n", xtick, nxsub);

    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about fitting position-angle swings and using the beam-width information can be found in:\n");
    printf(" - 	Rookyard et al. 2015, MNRAS, 446, 3367\n\n");
    printCitationInfo();
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
 PAoffset = atof(argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-loffset") == 0) {
 loffset = atof(argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-nokeypress") == 0) {
 nokeypresses = 1;
      }else if(strcmp(argv[i], "-totPol") == 0 || strcmp(argv[i], "-totpol") == 0) {
 showtotPol = 1;
      }else if(strcmp(argv[i], "-dash") == 0) {
 dashed = 1;
      }else if(strcmp(argv[i], "-noylabel") == 0) {
 noylabel = 1;
      }else if(strcmp(argv[i], "-N") == 0) {
 j = sscanf(argv[i+1], "%d %d", &nrpanelsx, &nrpanelsy);
 if(j != 2) {
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing %s option", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-xsize") == 0) {
 j = sscanf(argv[i+1], "%f", &xsize);
 if(j != 1) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -xsize option");
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-ysize") == 0) {
 j = sscanf(argv[i+1], "%f", &ysize);
 if(j != 1) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -ysize option");
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-sigma") == 0) {
 j = sscanf(argv[i+1], "%f", &sigma_limit);
 if(j != 1) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -sigma option");
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-ysizepa") == 0) {
 j = sscanf(argv[i+1], "%f", &ysizepa);
 if(j != 1) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -ysizepa option");
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-title") == 0) {
 strcpy(title, argv[i+1]);
 titleset = 1;
        i++;
      }else if(strcmp(argv[i], "-title_fmt") == 0) {
 j = sscanf(argv[i+1], "%f %d", &titlech, &titlelw);
 if(j != 2) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing %s option", argv[i]);
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-boxlw") == 0) {
 j = sscanf(argv[i+1], "%d", &boxlinewidth);
 if(j != 1) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -boxlw option");
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-lw") == 0) {
 j = sscanf(argv[i+1], "%d", &datalinewidth);
 if(j != 1) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -lw option");
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-xrange") == 0 || strcmp(argv[i], "-xrange_phase") == 0) {
 j = sscanf(argv[i+1], "%f %f", &plotl1, &plotl2);
 if(j != 2) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing %s option", argv[i]);
   return 0;
 }
 if(strcmp(argv[i], "-xrange_phase") == 0) {
   plotl1 *= 360.0;
   plotl2 *= 360.0;
 }
        i++;
      }else if(strcmp(argv[i], "-xticks") == 0) {
 j = sscanf(argv[i+1], "%f %d", &xtick, &nxsub);
 if(j != 2) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -xticks option");
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-yrange") == 0) {
 j = sscanf(argv[i+1], "%f %f", &plotp1, &plotp2);
 if(j != 2) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -yrange option");
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-yrange2") == 0) {
 j = sscanf(argv[i+1], "%f %f", &plotI1, &plotI2);
 if(j != 2) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -yrange2 option");
   return 0;
 }
        i++;
      }else if(strcmp(argv[i], "-paswing") == 0) {
 j = sscanf(argv[i+1], "%f %f %f %f", &overlayalpha, &overlaybeta, &overlaypa0, &overlayl0);
 if(j != 4) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -paswing option");
   return 0;
 }
 overlayPA = 1;
        i++;
      }else if(strcmp(argv[i], "-1") == 0) {
 plotp1 = 0.0;
 plotp2 = 180.0;
      }else if(strcmp(argv[i], "-text") == 0) {
 i += 2;
      }else if(strcmp(argv[i], "-herrorbar") == 0) {
 i += 1;
      }else if(strcmp(argv[i], "-herrorbar2") == 0) {
 i += 1;
      }else if(strcmp(argv[i], "-outline") == 0) {
 j = sscanf(argv[i+1], "%d %d", &outline, &outline_color);
 if(j != 2) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -outline option");
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-selectonpulse") == 0) {
 manualOnpulseSelection = 1;
      }else if(strcmp(argv[i], "-opm") == 0) {
 j = sscanf(argv[i+1], "%f %f", &jump_longitude[nrJumps], &jump_offset[nrJumps]);
 if(j != 2) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse -opm option, need two values.");
   return 0;
 }
 nrJumps++;
 if(nrJumps == MaxNrJumps) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Maximum number of -opm options exceeded.");
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-textkeywords") == 0) {
 str_list_replace_keys(0);
 return 0;
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
  if(numberInApplicationFilenameList(application, argv, application.verbose_state) == 0) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR ppol: No files specified");
    return 0;
  }


  firstfiletoopen = 1;
  while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {

    if(firstfiletoopen == 0) {
      if(curpanelnrx == 0 && curpanelnry == 0) {
 if(nokeypresses == 0) {
   i = pgplot_device_type(PlotDevice2, application.verbose_state);
   if(i < 3 || i > 10) {
     printf("Press a key in the terminal to continue\n");
     fflush(stdout);
     pgetch();
   }
 }
      }
    }

    cleanPSRData(&datain, application.verbose_state);
    cleanPSRData(&dataout, application.verbose_state);
    if(!openPSRData(&datain, filename_ptr, application.iformat, 0, 1, 0, application.verbose_state))
      return 0;

    if(application.verbose_state.verbose) {
      fflush(stdout);
      printwarning(application.verbose_state.debug, "Input data contains %ld bins, %ld pulses, %ld polarizations and %ld frequencies.", (datain.NrBins), datain.NrSubints, (datain.NrPols), datain.NrFreqChan);
    }



    if(PSRDataHeader_parse_commandline(&datain, argc, argv, application.verbose_state) == 0)
      return 0;

    if(datain.poltype == POLTYPE_ILVPAdPA || datain.poltype == POLTYPE_PAdPA) {
      compute_PA = 0;
      fflush(stdout);
      printwarning(application.verbose_state.debug, "WARNING ppol: File already contains PA data. ppol will show data, but not regenerate the PA points.");
    }else {
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
   datain.fixedtsamp = get_period(datain, 0, application.verbose_state)/(double)datain.NrBins;
 }
      }
      if(get_tsamp(datain, 0, application.verbose_state) < 0.0000001 || get_tsamp(datain, 0, application.verbose_state) >= 10) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "ERROR ppol: The sampling time does not appear to be set correctly in the header. Consider using the -header option.");
 return 0;
      }
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

    if(compute_PA) {

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
    }




    region_frac_to_int(&(application.onpulse), datain.NrBins, 0);

    if(compute_PA) {
      if(datain.NrSubints > 1) {
 printf("Data contains multiple subints. Data will be written out in FITS format.\n");
 application.oformat = FITS_format;
 if(showtotPol) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: The -totpol option is not supported when multiple subints are present");
   return 0;
 }
      }
      if(datain.NrFreqChan > 1) {
 printf("Data contains multiple frequency channels. Data will be written out in FITS format.\n");
 application.oformat = FITS_format;
 if(showtotPol) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: The -totpol option is not supported when multiple frequency channels are present");
   return 0;
 }
      }
      if(datain.NrPols != 4) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "ERROR ppol: Need Stokes IQUV data, but the number of polarization channels != 4!");
 return 0;
      }
    }else {
      if(datain.NrPols != 5) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "ERROR ppol: When plotting previously calculated PA data, the number of polarization channels is expected to be 5, but it is %d", datain.NrPols);
 if(datain.NrPols == 2) {
   printerror(application.verbose_state.debug, "Note that the PPOLSHORT format can not be displayed with ppol.");
 }
 return 0;
      }
    }

    if(compute_PA
) {





      if(showtotPol) {
 Ppulse = (float *)malloc((datain.NrBins)*sizeof(float));
 if(Ppulse == NULL) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Memory allocation error.");
   return 0;
 }
      }else {
 Ppulse = NULL;
      }
    }else {
      showtotPol = 0;
      writeout = 0;
      Ppulse = NULL;
    }

    if((compute_PA
 )
       ) {
      strcpy(txt, "Select on-pulse region ");
      strcat(txt, datain.psrname);
      profileI = (float *)malloc(datain.NrBins*sizeof(float));
      if(profileI == NULL) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "ERROR ppol: Memory allocation error.");
 return 0;
      }
      if(read_profilePSRData(datain, profileI, NULL, 0, application.verbose_state) != 1) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "ERROR ppol: Cannot form pulse profile");
 return 0;
      }
      if(manualOnpulseSelection == 1 || (application.onpulse.nrRegions == 0 && manualOnpulseSelection != -1)) {
 pgplot_clear_viewport_def(&viewport);
 strcpy(viewport.plotDevice, PlotDevice);
 pgplot_box_def pgplotbox;
 clear_pgplot_box(&pgplotbox);
 strcpy(pgplotbox.xlabel, "Bin");
 strcpy(pgplotbox.ylabel, "Intensity");
 strcpy(pgplotbox.title, txt);
 selectRegions(profileI, datain.NrBins, viewport, pgplotbox, 0, 0, 0, &(application.onpulse), application.verbose_state);
 if(firstfiletoopen == 0) {
   ppgslct(deviceID);
 }
      }else {
 pgplot_clear_viewport_def(&viewport);
 clear_pgplot_box(&pgplotbox);
 strcpy(pgplotbox.xlabel, "Bin");
 strcpy(pgplotbox.ylabel, "Intensity");
 strcpy(pgplotbox.title, txt);
 strcpy(viewport.plotDevice, PlotDevice);
 pgplotGraph1(viewport, profileI, NULL, NULL, datain.NrBins, 0, datain.NrBins, 0, 0, datain.NrBins, 0, 0, 0, pgplotbox, 1, 0, 0, 1, 1, &(application.onpulse), application.verbose_state);
 if(firstfiletoopen == 0) {
   ppgslct(deviceID);
 }
      }
      region_int_to_frac(&(application.onpulse), 1.0/(float)datain.NrBins, 0);
      regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);

 verbose_definition verbose2;
 copyVerboseState(application.verbose_state, &verbose2);
 if(application.verbose_state.verbose && datain.NrSubints == 1 && datain.NrFreqChan == 1)
   verbose2.verbose = 1;
 else
   verbose2.verbose = 0;

 if(make_paswing_fromIQUV(&datain, Ppulse, application.onpulse, normalize, correctLbias, correctQV, correctV, loffset, verbose2) == 0)
 return 0;
    }


    float currentPanelScaling, currentPanelScalingx, currentPanelScalingy, viewport_left, viewport_right, viewport_top, viewport_bottom;
    currentPanelScaling = 1;
    viewport_left = 0.1;
    if(nrpanelsx > 1) {
      currentPanelScalingx = xsize/((nrpanelsx-1)*0.2 + nrpanelsx*xsize);
      currentPanelScaling = currentPanelScalingx;
      viewport_left += curpanelnrx*currentPanelScalingx*0.1;
      viewport_left += curpanelnrx*currentPanelScalingx*xsize;
      viewport_left += curpanelnrx*currentPanelScalingx*0.1;
    }
    viewport_right = viewport_left+xsize;
    if(nrpanelsx > 1) {
      viewport_right = viewport_left;
      viewport_right += currentPanelScalingx*xsize;
    }
    viewport_bottom = 0.35-0.25*ysize*ysizepa;
    if(nrpanelsy > 1) {
      currentPanelScalingy = (0.55*ysize + 0.25*ysize*ysizepa)/((nrpanelsy-1)*(0.35-0.25*ysize*ysizepa) + nrpanelsy*(0.55*ysize + 0.25*ysize*ysizepa) + (nrpanelsy-1)*(1-0.35-0.55*ysize));
      if(currentPanelScalingy < currentPanelScaling)
 currentPanelScaling = currentPanelScalingy;
      viewport_bottom += (nrpanelsy-curpanelnry-1)*currentPanelScalingy*(0.35-0.25*ysize*ysizepa);
      viewport_bottom += (nrpanelsy-curpanelnry-1)*currentPanelScalingy*(0.55*ysize + 0.25*ysize*ysizepa);
      viewport_bottom += (nrpanelsy-curpanelnry-1)*currentPanelScalingy*(1-0.35-0.55*ysize);
    }
    viewport_top = viewport_bottom + (0.55*ysize + 0.25*ysize*ysizepa);
    if(nrpanelsy > 1) {
      viewport_top = viewport_bottom;
      viewport_top += currentPanelScalingy*(0.55*ysize + 0.25*ysize*ysizepa);
    }
    if(curpanelnrx == 0 && curpanelnry == 0)
      ppgpage();
    curpanelnrx++;
    if(curpanelnrx == nrpanelsx) {
      curpanelnrx = 0;
      curpanelnry++;
    }
    if(curpanelnry == nrpanelsy)
      curpanelnry = 0;
      if(sigma_limit > 0
) {
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
      datain.data[j+datain.NrBins*(datain.NrPols-2 + datain.NrPols*(f+datain.NrFreqChan*i))] = 0;
      datain.data[j+datain.NrBins*(datain.NrPols-1 + datain.NrPols*(f+datain.NrFreqChan*i))] = -1;
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
      clear_pgplot_box(&pgplotbox);
      pgplot_clear_viewport_def(&viewport);
      viewport.dxplot = viewport_left - 0.1;
      viewport.xsize = (viewport_right - 0.1 - viewport.dxplot)/0.8;
      viewport.ysize = (viewport_top - viewport_bottom)/(0.55 + 0.25*ysizepa);
      viewport.dyplot = viewport_top - 0.35 - 0.55*(viewport_top - viewport_bottom)/(0.55 + 0.25*ysizepa);
      strcpy(viewport.plotDevice, PlotDevice2);
      viewport.noclear = 1;
      viewport.dontclose = 1;
      viewport.dontopen = !firstfiletoopen;
      pgplotbox.box_lw = boxlinewidth;
      pgplotbox.box_xtick = xtick;
      pgplotbox.box_nxsub = nxsub;
      if(titleset)
 strcpy(pgplotbox.title, title);
      else
 strcpy(pgplotbox.title, datain.filename);
      pgplotbox.title_ch = titlech*currentPanelScaling;
      pgplotbox.title_lw = titlelw;
      pgplotbox.label_ch *= currentPanelScaling;
      pgplotbox.box_labelsize *= currentPanelScaling;
      pgplotPAplot(datain, Ppulse, viewport, pgplotbox, "Pulse longitude [deg]", "I,Linear,V", "PA [deg]", plotl1, plotl2, plotI1, plotI2, plotp1, plotp2, PAoffset, sigma_limit, datalinewidth, ysizepa, dashed, noylabel, "-text", "-herrorbar", "-herrorbar2", argc, argv, outline, outline, outline_color, overlayPA, overlayalpha, overlaybeta, overlaypa0, overlayl0, overlayFine, nrJumps, jump_longitude, jump_offset, application.verbose_state);
      if(firstfiletoopen) {
 ppgqid(&deviceID);
      }
    if(writeout
) {
      copy_params_PSRData(datain, &dataout, application.verbose_state);
      if(extension != NULL) {
 if(change_filename_extension(filename_ptr, ofilename, extension, MaxOutputNameLength, application.verbose_state) == 0) {
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
      dataout.poltype = POLTYPE_ILVPAdPA;
      if(!writeHeaderPSRData(&dataout, argc, argv, application.history_cmd_only, application.verbose_state))
 return 0;
      if(application.oformat == PPOL_format) {
 extended = 1;
 onlysignificantPA = 0;
      }else {
 extended = 0;
 onlysignificantPA = 1;
      }
      if(dataout.format == FITS_format) {
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
      }else {
 dataout.offpulse_rms = datain.offpulse_rms;
 if(writePPOLfile(dataout, datain.data, extended, onlysignificantPA, twoprofiles, PAoffset, application.verbose_state) == 0) {
   fflush(stdout);
   printerror(application.verbose_state.debug, "ERROR ppol: Cannot write data");
   return 0;
 }
      }
    }
    closePSRData(&datain, application.verbose_state);
    closePSRData(&dataout, application.verbose_state);
    firstfiletoopen = 0;
    if(compute_PA)
      free(profileI);
  }
  ppgend();
  return 0;
}

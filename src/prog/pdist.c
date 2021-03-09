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
#include <gsl/gsl_sort.h>
#include "psrsalsa.h"
int main(int argc, char **argv)
{
  int read_log, centered_at_zero, file_column1, file_column2, colspecified, polspecified, filename, truncate, cdf, select, long_resolved, lrnorm;
  long ndata, i, j, k, nrbins, nrbinsy, *distr, nrdatacolumns;
  float *cmap, *distr_float;
  double **data, min_x_data, max_x_data, min_y_data, max_y_data, average_x, average_y, normalise;
  double min_x, max_x, min_y, max_y, dx, dy, x, y, s, rangex_min, rangex_max, rangey_min, rangey_max, extra_phase, select1, select2;
  int nrbins_specified, dx_specified, dy_specified, output_sigma, output_fraction, nrbins_specifiedy, twoDmode, showGraphics, rangex_set, rangey_set, plotlog, file_column2_defined, xlabelset, linewidth;
  char title[500], *filename_ptr, output_suffix[100], *oname;
  FILE *ofile;
  double min, max;
  psrsalsaApplication application;
  datafile_definition datain, dataout;
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);
  initApplication(&application, "pdist", "[options] inputfile(s)");
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_iformat = 1;
  application.switch_formatlist = 1;
  application.switch_device = 1;
  application.switch_cmap = 1;
  application.switch_cmaplist = 1;
  application.switch_nocounters = 1;
  application.switch_history_cmd_only = 1;
  application.oformat = FITS_format;
  nrbins_specified = 0;
  nrbins_specifiedy = 0;
  dx_specified = 0;
  dy_specified = 0;
  output_sigma = 0;
  output_fraction = 0;
  read_log = 0;
  twoDmode = 0;
  centered_at_zero = 1;
  extra_phase = 0;
  showGraphics = 0;
  xlabelset = 0;
  sprintf(title, "Distribution");
  rangex_set = 0;
  rangey_set = 0;
  file_column1 = 1;
  file_column2 = 2;
  colspecified = 0;
  polspecified = 0;
  plotlog = 0;
  file_column2_defined = 1;
  application.cmap = PPGPLOT_HEAT;
  strcpy(output_suffix, "hist");
  filename = 0;
  truncate = 0;
  cdf = 0;
  select = 0;
  linewidth = 1;
  long_resolved = 0;
  lrnorm = 0;
  normalise = 0.0;
  if(argc < 2) {
    printf("Program to generate or plot a histogram by binning data. Also a cumulative\n");
    printf("distribution can be generated. Usage:\n\n");
    printApplicationHelp(&application);
    printf("Input options:\n");
    printf("-2               Turn on 2D mode: Read in two columns of data, resulting in the\n");
    printf("                 distribution N(x,y) rathern than N(x).\n");
    printf("-col \"c1 c2\"     Specify two column numbers (counting from 1) to be used.\n");
    printf("                 Without the -2 option only c1 needs to be specified.\n");
    printf("                 This option implies the input files are ascii files with each\n");
    printf("                 line having an equal number of columns. Lines starting with a #\n");
    printf("                 will be ignored.\n");
    printf("-pol \"p1 p2\"     This option implies that the input file is recognized as a\n");
    printf("                 pulsar format. Polarization p1 (and p2 in 2D mode), counting\n");
    printf("                 from zero, are used (all bins, freqs and subints). The default\n");
    printf("                 default is \"0 1\", or the integrated onpulse pulse energy for\n");
    printf("                 penergy output files in mode 1.\n");
    printf("\nOutput options:\n");
    printf("-ext             Specify suffix, default is '%s'\n", output_suffix);
    printf("-output filename Write output to this file [def=.%s extension].\n", output_suffix);
    printf("\nDistribution range and binning option:\n");
    printf("-dx value        Specify bin width (can also use -n).\n");
    printf("-dy value        Specify bin width used for the second column if -2 is used.\n");
    printf("-n nr            The range of input values is divided in nr bins.\n");
    printf("-nx nr           Same as -n\n");
    printf("-ny nr           Similar to -nx, but now for the second column if -2 is used.\n");
    printf("-rangex \"x1 x2\"  The range of bins produced in output will cover x1..x2, which\n");
    printf("                 should include all input values (see also -trunc). It therefore\n");
    printf("                 allows more bins to be generated than strictly necessary.\n");
    printf("-rangey          As -rangex, but now for second input column.\n");
    printf("-trunc           Modifies behaviour of -rangex and -rangey. Values outside the\n");
    printf("                 specified range are set to the extremes of the range.\n");
    printf("-zero            By default one of the bins will be centred at zero. This option\n");
    printf("                 results in the edge of a bin to coincides with zero. In both\n");
    printf("                 cases the reported locations of the bins are their centres.\n");
    printf("                 This means that in gnuplot you want to use \"with histeps\".\n");
    printf("-zeroshift off   Put in an extra offset off to where zero would fall with\n");
    printf("                 respect to a bin. -zero is equivalent to -zeroshift 0.5. \n");
    printf("\nOther distribution options:\n");
    printf("-cdf             The cumulative distribution is generated. Use no -dx or -n.\n");
    printf("-frac            Output fraction of counts per bin rather than counts.\n");
    printf("-log             The base-10 log of the input values is used.\n");
    printf("-norm            Normalise such that the average of the distribution is 1.\n");
    printf("-sigma           Generate extra column with the sqrt of the nr of counts.\n");
    printf("\nPlot options:\n");
    printf("-plot            Plot the distribution rather than writing output to file.\n");
    printf("-plotlog         The log10 of the counts is plotted (implies -plot and bins with\n");
    printf("                 zero counts are set to 1).\n");
    printf("-overplot        Like plot, but multiple input files are overplotted. The ranges\n");
    printf("                 are determined by first input file.\n");
    printf("-title  'title'  Set the title of the plot.\n");
    printf("-xlabel 'label'  Set label on horizontal axis.\n");
    printf("-labels          \"label_ch label_lw label_f box_ch box_lw box_f\"\n");
    printf("                  default is \"%.1f %d %d %.1f %d %d\"\n", pgplot_options.box.label_ch, pgplot_options.box.label_lw, pgplot_options.box.label_f, pgplot_options.box.box_labelsize, pgplot_options.box.box_lw, pgplot_options.box.box_f);
    printf("-lw              Set line width of the histogram.\n");
    printf("-vp              \"left right bottom top\"\n");
    printf("                 Modify the dimensions of the plot. Default is \"%.2f %.2f %.2f %.2f\"\n", pgplot_options.viewport.dxplot, pgplot_options.viewport.xsize, pgplot_options.viewport.dyplot, pgplot_options.viewport.ysize);
    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about fitting distributions (in the context of pulse energies) can be found in:\n");
    printf(" - Weltevrede et al. 2006, A&A, 458, 269\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      int index;
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
 i = index;
      }else if(strcmp(argv[i], "-2") == 0) {
 twoDmode = 1;
      }else if(strcmp(argv[i], "-plot") == 0) {
 showGraphics = 1;
      }else if(strcmp(argv[i], "-overplot") == 0) {
 showGraphics = 2;
      }else if(strcmp(argv[i], "-plotlog") == 0) {
 plotlog = 1;
 if(showGraphics == 0)
   showGraphics = 1;
      }else if(strcmp(argv[i], "-xlabel") == 0) {
 xlabelset = i+1;
 i++;
      }else if(strcmp(argv[i], "-title") == 0) {
 strcpy(title, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-output") == 0) {
 filename = i+1;
 i++;
      }else if(strcmp(argv[i], "-ext") == 0) {
 strcpy(output_suffix, argv[i+1]);
 i++;
      }else if(strcmp(argv[i], "-labels") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %d %d %f %d %d", &(pgplot_options.box.label_ch), &(pgplot_options.box.label_lw), &(pgplot_options.box.label_f), &(pgplot_options.box.box_labelsize), &(pgplot_options.box.box_lw), &(pgplot_options.box.box_f), NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-lw") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &linewidth, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-vp") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f %f", &(pgplot_options.viewport.dxplot), &(pgplot_options.viewport.xsize), &(pgplot_options.viewport.dyplot), &(pgplot_options.viewport.ysize), NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-col") == 0 || strcmp(argv[i], "-pol") == 0) {
 if(strcmp(argv[i], "-pol") == 0)
   polspecified = 1;
 else
   colspecified = 1;
 int ret;
 ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%d %d", &file_column1, &file_column2, NULL);
 file_column2_defined = 0;
 if(ret == 2) {
   file_column2_defined = 1;
 }else if(ret == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse %s option, need 1 or 2 integer values.\n", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "-nx") == 0) {
 if(dx_specified) {
   printerror(application.verbose_state.debug, "Cannot use -n and -dx option simultaniously.\n");
   return 0;
 }
 nrbins_specified = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &nrbins, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-ny") == 0) {
 if(dy_specified) {
   printerror(application.verbose_state.debug, "Cannot use -ny and -dy option simultaniously.\n");
   return 0;
 }
 nrbins_specifiedy = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &nrbinsy, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-dx") == 0) {
 if(nrbins_specified) {
   printerror(application.verbose_state.debug, "Cannot use -n and -dx option simultaniously.\n");
   return 0;
 }
 dx_specified = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &dx, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-dy") == 0) {
 if(nrbins_specifiedy) {
   printerror(application.verbose_state.debug, "Cannot use -ny and -dy option simultaniously.\n");
   return 0;
 }
 dy_specified = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &dy, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-rangex") == 0) {
 rangex_set = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &rangex_min, &rangex_max, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-rangey") == 0) {
 rangey_set = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &rangey_min, &rangey_max, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-select") == 0) {
 select = 1;
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &select1, &select2, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-trunc") == 0) {
 truncate = 1;
      }else if(strcmp(argv[i], "-sigma") == 0) {
 output_sigma = 1;
      }else if(strcmp(argv[i], "-log") == 0) {
 read_log = 1;
      }else if(strcmp(argv[i], "-zero") == 0) {
 centered_at_zero = 0;
      }else if(strcmp(argv[i], "-zeroshift") == 0) {
 if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &extra_phase, NULL) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: Cannot parse '%s' option.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-norm") == 0) {
 normalise = 1.0;
      }else if(strcmp(argv[i], "-frac") == 0) {
 output_fraction = 1;
      }else if(strcmp(argv[i], "-cdf") == 0) {
 cdf = 1;
      }else {
 if(argv[i][0] == '-') {
   printerror(application.verbose_state.debug, "pdist: Unknown option: %s\n\nRun pdist without command line arguments to show help", argv[i]);
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
    printerror(application.verbose_state.debug, "ERROR pdist: No files specified");
    return 0;
  }
  if(select) {
    if(twoDmode) {
      printerror(application.verbose_state.debug, "ERROR pdist: The -2 option cannot be used with -select.\n");
      return 0;
    }
    if(cdf) {
      printerror(application.verbose_state.debug, "ERROR pdist: The -cdf option cannot be used with -select.\n");
      return 0;
    }
    if(showGraphics) {
      printerror(application.verbose_state.debug, "ERROR pdist: No plots can be generated when together using -select.\n");
      return 0;
    }
    if(output_fraction) {
      printerror(application.verbose_state.debug, "ERROR pdist: The -frac option cannot be used with -select.\n");
      return 0;
    }
    if(output_sigma) {
      printerror(application.verbose_state.debug, "ERROR pdist: The -sigma option cannot be used together with -select.\n");
      return 0;
    }
  }
  if(cdf) {
    if(nrbins_specified || dx_specified) {
      printerror(application.verbose_state.debug, "ERROR pdist: For a cdf -n and -dx should not be used.\n");
      return 0;
    }
    if(twoDmode) {
      printerror(application.verbose_state.debug, "ERROR pdist: For a cdf only one input column is expected.\n");
      return 0;
    }
    if(output_sigma) {
      printerror(application.verbose_state.debug, "ERROR pdist: For a cdf -sigma is not supported.\n");
      return 0;
    }
    if(rangex_set || rangey_set || truncate) {
      printerror(application.verbose_state.debug, "ERROR pdist: For a cdf the -rangex and -rangey options are not supported.\n");
      return 0;
    }
  }else {
    if(select == 0) {
      if(nrbins_specified == 0 && dx_specified == 0) {
 printerror(application.verbose_state.debug, "ERROR pdist: Specify resolution with the -n or -dx option.\n");
 return 0;
      }
      if(twoDmode) {
 if(nrbins_specifiedy == 0 && dy_specified == 0) {
 printerror(application.verbose_state.debug, "ERROR pdist: Specify resolution with the -ny or -dy option.\n");
 return 0;
 }
      }
    }
  }
  int initial_iformat, currentfile;
  currentfile = 1;
  initial_iformat = application.iformat;
  while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
    application.iformat = initial_iformat;
    if(application.iformat <= 0 && colspecified == 0) {
      application.iformat = guessPSRData_format(filename_ptr, 1, application.verbose_state);
      if(application.iformat == -2 || application.iformat == -3)
 return 0;
    }
    if(colspecified)
      application.iformat = -1;
    if(polspecified && application.iformat <= 0) {
      printerror(application.verbose_state.debug, "ERROR pdist: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized. Data is not recognized as pulsar data, but -pol was used. Use -col for ascii files.\n", filename_ptr);
      return 0;
    }
    if(application.iformat > 0) {
      verbose_definition verbose;
      cleanVerboseState(&verbose);
      copyVerboseState(application.verbose_state, &verbose);
      if(application.iformat == PSRCHIVE_ASCII_format) {
 if(verbose.debug == 0)
   verbose.verbose = 0;
      }
      i = openPSRData(&datain, filename_ptr, application.iformat, 0, 1, 1, verbose);
      if(i == 0) {
 printerror(application.verbose_state.debug, "ERROR pdist: Error opening data");
 return 0;
      }
      if(polspecified == 0) {
 file_column1 = 0;
 file_column2 = 1;
      }
      if(application.iformat == PSRCHIVE_ASCII_format && datain.gentype == GENTYPE_PENERGY) {
 if(application.verbose_state.verbose)
   printf("The file is a penergy file generated in mode 1.\n");
 if(polspecified == 0) {
   printwarning(application.verbose_state.debug, "WARNING pdist: -pol not specified, default in the integrated on-pulse energy (-pol 2)");
   file_column1 = 2;
 }
      }
      if(application.verbose_state.verbose)
 printf("Loading %ld points", datain.NrSubints*datain.NrFreqChan*datain.NrBins);
      if(twoDmode)
 printf(" times two input polarizations");
      printf("\n");
      if(long_resolved == 0) {
 if(twoDmode || cdf || select) {
   nrdatacolumns = 2;
 }else {
   nrdatacolumns = 1;
 }
      }else {
 nrdatacolumns = datain.NrBins;
      }
      data = malloc(nrdatacolumns*sizeof(double *));
      if(data == NULL) {
 printerror(application.verbose_state.debug, "ERROR pdist: Memory allocation error.\n");
 return 0;
      }
      for(i = 0; i < nrdatacolumns; i++) {
 if(long_resolved == 0) {
   data[i] = malloc(datain.NrSubints*datain.NrFreqChan*datain.NrBins*sizeof(double));
 }else {
   data[i] = malloc(datain.NrSubints*datain.NrFreqChan*sizeof(double));
 }
 if(data[i] == NULL) {
   printerror(application.verbose_state.debug, "ERROR pdist: Memory allocation error.\n");
   return 0;
 }
      }
      long subintnr, freqnr, binnr;
      long double total_x, total_y;
      ndata = 0;
      total_x = 0;
      total_y = 0;
      min_x_data = max_x_data = 0;
      for(subintnr = 0; subintnr < datain.NrSubints; subintnr++) {
 for(freqnr = 0; freqnr < datain.NrFreqChan; freqnr++) {
   for(binnr = 0; binnr < datain.NrBins; binnr++) {
     int polnr;
     int nrpols = 1;
     if(twoDmode) {
       nrpols = 2;
     }
     for(polnr = 0; polnr < nrpols; polnr++) {
       int file_column;
       if(polnr == 0) {
  file_column = file_column1;
       }else {
  file_column = file_column2;
       }
       float dummy_float;
       if(readPulsePSRData(&datain, subintnr, file_column, freqnr, binnr, 1, &dummy_float, application.verbose_state) == 0) {
  printerror(application.verbose_state.debug, "ERROR pdist: Read error, shouldn't happen.\n");
  return 0;
       }
       if(read_log) {
  if(dummy_float <= 0) {
    printerror(application.verbose_state.debug, "ERROR pdist: Cannot take logarithm of a value <= 0.\n");
    return 0;
  }
  dummy_float = log10(dummy_float);
       }
       if(long_resolved == 0) {
  data[polnr][ndata] = dummy_float;
       }else {
  data[binnr][ndata] = dummy_float;
       }
       if(polnr == 0) {
  if(dummy_float > max_x_data || ndata == 0) {
    max_x_data = dummy_float;
  }
  if(dummy_float < min_x_data || ndata == 0) {
    min_x_data = dummy_float;
  }
  total_x += dummy_float;
       }else {
  if(twoDmode) {
    if(ndata == 0) {
      max_y_data = dummy_float;
      min_y_data = dummy_float;
    }
    if(dummy_float > max_y_data) {
      max_y_data = dummy_float;
    }
    if(dummy_float < min_y_data) {
      min_y_data = dummy_float;
    }
    total_y += dummy_float;
  }
       }
     }
     if(long_resolved == 0) {
       ndata++;
     }
   }
   if(long_resolved) {
     ndata++;
   }
 }
      }
      if(long_resolved) {
 average_x = total_x/(long double)(ndata*datain.NrBins);
 if(twoDmode) {
   average_y = total_y/(long double)(ndata*datain.NrBins);
 }
      }else {
 average_x = total_x/(long double)ndata;
 if(twoDmode) {
   average_y = total_y/(long double)ndata;
 }
      }
      if(application.verbose_state.verbose) {
 printf("xrange = %e ... %e\n", min_x_data, max_x_data);
 printf("average = %e\n", average_x);
 if(twoDmode) {
   printf("yrange = %e ... %e\n", min_y_data, max_y_data);
   printf("average = %e\n", average_y);
 }
      }
      if(normalise != 0.0) {
 if(long_resolved == 0) {
   for(i = 0; i < ndata; i++) {
     data[0][i] *= normalise/average_x;
     if(twoDmode) {
       data[1][i] *= normalise/average_y;
     }
   }
 }else {
   for(j = 0; j < nrdatacolumns; j++) {
     for(i = 0; i < ndata; i++) {
       data[j][i] *= normalise/average_x;
     }
   }
 }
 min_x_data *= normalise/average_x;
 max_x_data *= normalise/average_x;
 average_x *= normalise/average_x;
 if(twoDmode) {
   min_y_data *= normalise/average_y;
   max_y_data *= normalise/average_y;
   average_y *= normalise/average_y;
 }
 if(application.verbose_state.verbose) {
   printf("new xrange after applying normalisation = %e ... %e\n", min_x_data, max_x_data);
 }
      }
      if(long_resolved && lrnorm == 1) {
 double total, scale;
 for(j = 0; j < nrdatacolumns; j++) {
   total = 0;
   for(i = 0; i < ndata; i++) {
     total += data[j][i];
   }
   scale = ndata/total;
   for(i = 0; i < ndata; i++) {
     data[j][i] *= scale;
     if(data[j][i] > max_x_data || (i == 0 && j == 0)) {
       max_x_data = data[j][i];
     }
     if(data[j][i] < min_x_data || (i == 0 && j == 0)) {
       min_x_data = data[j][i];
     }
   }
 }
 if(application.verbose_state.verbose) {
   printf("new xrange after applying normalisation = %e ... %e\n", min_x_data, max_x_data);
 }
      }
      if(long_resolved && lrnorm == 2) {
 double total, scale;
 for(j = 0; j < nrdatacolumns; j++) {
   total = 0;
   for(i = 0; i < ndata; i++) {
     total += data[j][i];
   }
   if(j == 0 || total > scale) {
     scale = total;
   }
 }
 scale = ndata/scale;
 for(j = 0; j < nrdatacolumns; j++) {
   for(i = 0; i < ndata; i++) {
     data[j][i] *= scale;
     if(data[j][i] > max_x_data || (i == 0 && j == 0)) {
       max_x_data = data[j][i];
     }
     if(data[j][i] < min_x_data || (i == 0 && j == 0)) {
       min_x_data = data[j][i];
     }
   }
 }
 if(application.verbose_state.verbose) {
   printf("new xrange after applying normalisation = %e ... %e\n", min_x_data, max_x_data);
 }
      }
      if(long_resolved) {
 cleanPSRData(&dataout, application.verbose_state);
 copy_params_PSRData(datain, &dataout, application.verbose_state);
      }
      closePSRData(&datain, 0, application.verbose_state);
    }else {
      if(twoDmode || cdf || select) {
 nrdatacolumns = 2;
      }else {
 nrdatacolumns = 1;
      }
      data = malloc(nrdatacolumns*sizeof(double *));
      if(data == NULL) {
 printerror(application.verbose_state.debug, "ERROR pdist: Memory allocation error.\n");
 return 0;
      }
      int skiplines = 0;
      if(twoDmode) {
 if(application.verbose_state.verbose)
   fprintf(stdout, "Loading x values from ascii file\n");
 if(read_ascii_column_double(filename_ptr, skiplines, '#', -1, 1, &ndata, 0, file_column1, 1.0, normalise, read_log, &(data[0]), &min_x_data, &max_x_data, &average_x, application.verbose_state, 1) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: cannot load file.\n");
   if(colspecified) {
     printwarning(application.verbose_state.debug, "WARNING pdist: Using the -col option implies the input file is a simple ascii file. For penergy output (in mode 1), or an other recognized pulsar format, use the -pol option instead.\n");
   }
   return 0;
 }
 if(application.verbose_state.verbose)
   fprintf(stdout, "Loading y values from ascii file\n");
 if(file_column2_defined == 0) {
   printerror(application.verbose_state.debug, "In 2D mode, two input columns should be specified with the -col option.\n");
   return 0;
 }
 if(read_ascii_column_double(filename_ptr, skiplines, '#', -1, 1, &ndata, 0, file_column2, 1.0, normalise, read_log, &(data[1]), &min_y_data, &max_y_data, &average_y, application.verbose_state, 1) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: cannot load file.\n");
   if(colspecified) {
     printwarning(application.verbose_state.debug, "WARNING pdist: Using the -col option implies the input file is a simple ascii file. For penergy output (in mode 1), or an other recognized pulsar format, use the -pol option instead.\n");
   }
   return 0;
 }
      }else {
 if(application.verbose_state.verbose)
   fprintf(stdout, "Loading values from ascii file\n");
 if(read_ascii_column_double(filename_ptr, skiplines, '#', -1, 1, &ndata, 0, file_column1, 1.0, normalise, read_log, &(data[0]), &min_x_data, &max_x_data, &average_x, application.verbose_state, 1) == 0) {
   printerror(application.verbose_state.debug, "ERROR pdist: cannot load file.\n");
   if(colspecified) {
     printwarning(application.verbose_state.debug, "WARNING pdist: Using the -col option implies the input file is a simple ascii file. For penergy output (in mode 1), or an other recognized pulsar format, use the -pol option instead.\n");
   }
   return 0;
 }
 if(cdf || select) {
   data[1] = malloc(ndata*sizeof(double));
   if(data[1] == NULL) {
     printerror(application.verbose_state.debug, "ERROR pdist: Memory allocation error.\n");
     return 0;
   }
 }
      }
    }
    if(cdf == 0 && select == 0) {
      switch(set_binning_histogram(min_x_data, max_x_data, rangex_set, rangex_min, rangex_max, nrbins_specified, nrbins, centered_at_zero, extra_phase, &min_x, &max_x, &dx, application.verbose_state)) {
      case 0: break;
      case 1:
      case 2: return 0;
      default: {
 printerror(application.verbose_state.debug, "ERROR pdist: Unknown return value for call to set_binning_histogram.\n");
 return 0;
      }
      }
      if(rangex_set) {
 if(min_x > min_x_data || max_x < max_x_data) {
   if(truncate == 0) {
     printerror(application.verbose_state.debug, "ERROR pdist: Values are outside specified -rangex option. You may want to use -trunc?\n");
     return 0;
   }
 }
      }
      if(twoDmode) {
 switch(set_binning_histogram(min_y_data, max_y_data, rangey_set, rangey_min, rangey_max, nrbins_specifiedy, nrbinsy, centered_at_zero, extra_phase, &min_y, &max_y, &dy, application.verbose_state)) {
 case 0: break;
 case 1:
 case 2: return 0;
 default: {
   printerror(application.verbose_state.debug, "ERROR pdist: Unknown return value for call to set_binning_histogram.\n");
   return 0;
 }
 }
 if(application.verbose_state.verbose) {
   fprintf(stdout, "Going to use binsize dx=%e and dy=%e.\n", dx, dy);
 }
      }
      if(rangey_set) {
 if(min_y > min_y_data || max_y < max_y_data) {
   if(truncate == 0) {
     printerror(application.verbose_state.debug, "ERROR pdist: Values are outside specified -rangey option. You may want to use -trunc?\n");
     return 0;
   }
 }
      }
    }
    if(cdf == 0 && truncate) {
      if(long_resolved == 0) {
 if(truncate) {
   if(rangex_set) {
     for(i = 0; i < ndata; i++) {
       if(data[0][i] < min_x)
  data[0][i] = min_x;
       if(data[0][i] > max_x)
  data[0][i] = max_x;
     }
   }
   if(rangey_set) {
     for(i = 0; i < ndata; i++) {
       if(data[1][i] < min_y)
  data[1][i] = min_y;
       if(data[1][i] > max_y)
  data[1][i] = max_y;
     }
   }
 }
      }else {
 if(truncate) {
   if(rangex_set) {
     for(j = 0; j < nrdatacolumns; j++) {
       for(i = 0; i < ndata; i++) {
  if(data[j][i] < min_x)
    data[j][i] = min_x;
  if(data[j][i] > max_x)
    data[j][i] = max_x;
       }
     }
   }
 }
      }
    }
    if(cdf == 0 && select == 0) {
      nrbins = calculate_bin_number(max_x, dx, min_x, centered_at_zero, extra_phase) + 1;
    }else {
      nrbins = ndata;
    }
    if(select == 0)
      fprintf(stdout, "Distribution will have %ld bins.\n", nrbins);
    if(twoDmode) {
      nrbinsy = calculate_bin_number(max_y, dy, min_y, centered_at_zero, extra_phase)+ 1;
      fprintf(stdout, "Distribution will have %ld y-bins.\n", nrbinsy);
    }else {
      nrbinsy = 1;
    }
    long binnr, nr_distributions_to_compute = 1;
    if(long_resolved) {
      nr_distributions_to_compute = nrdatacolumns;
    }
    for(binnr = 0; binnr < nr_distributions_to_compute; binnr++) {
      distr = (long *)malloc(nrbins*nrbinsy*sizeof(long));
      if(distr == NULL) {
 printerror(application.verbose_state.debug, "ERROR pdist: Cannot allocate memory.\n");
 return 0;
      }
      for(i = 0; i < nrbins*nrbinsy; i++)
 distr[i] = 0;
      if(binnr == 0) {
 if(showGraphics) {
   pgplot_options.box.drawtitle = 1;
   strcpy(pgplot_options.box.title, title);
   if(xlabelset) {
     strcpy(pgplot_options.box.xlabel, argv[xlabelset]);
   }else {
     if(read_log)
       strcpy(pgplot_options.box.xlabel, "log X");
     else
       strcpy(pgplot_options.box.xlabel, "X");
   }
   if(twoDmode) {
     if(read_log)
       strcpy(pgplot_options.box.ylabel, "log Y");
     else
       strcpy(pgplot_options.box.ylabel, "Y");
   }else {
     if(plotlog)
       strcpy(pgplot_options.box.ylabel, "log N");
     else
       strcpy(pgplot_options.box.ylabel, "N");
   }
   strcpy(pgplot_options.viewport.plotDevice, application.pgplotdevice);
   if(showGraphics == 2) {
     pgplot_options.viewport.dontclose = 1;
   }
   if(currentfile > 1 && showGraphics == 2) {
     pgplot_options.viewport.noclear = 1;
     pgplot_options.viewport.dontopen = 1;
     pgplot_options.box.drawbox = 0;
     pgplot_options.box.drawtitle = 0;
     pgplot_options.box.drawlabels = 0;
   }
 }else {
   if(filename != 0) {
     oname = (char *) calloc(strlen(argv[filename])+1, 1);
   }else {
     oname = (char *) calloc(strlen(filename_ptr)+strlen(output_suffix)+2, 1);
   }
   if(oname == NULL) {
     printerror(application.verbose_state.debug, "ERROR pdist: Memory allocation error.\n");
     return 0;
   }
   if(filename != 0) {
     strcpy(oname, argv[filename]);
   }else {
     sprintf(oname,"%s.%s", filename_ptr, output_suffix);
   }
   if(long_resolved == 0) {
     printf("Output Ascii to %s\n", oname);
     ofile = fopen(oname,"w+");
     if(ofile == NULL) {
       printerror(application.verbose_state.debug, "ERROR pdist: Cannot open '%s'", oname);
       perror("");
       return 0;
     }
   }
 }
      }
      if(twoDmode) {
 for(i = 0; i < ndata; i++) {
   j = calculate_bin_number(data[0][i], dx, min_x, centered_at_zero, extra_phase);
   k = calculate_bin_number(data[1][i], dy, min_y, centered_at_zero, extra_phase);
   if(j < 0 || j >= nrbins) {
     printerror(application.verbose_state.debug, "BUG!\n");
     return 0;
   }
   if(k < 0 || k >= nrbinsy) {
     printerror(application.verbose_state.debug, "BUG!\n");
     return 0;
   }
   distr[k*nrbins+j] += 1;
 }
 for(i = 0; i < nrbins; i++) {
   for(j = 0; j < nrbinsy; j++) {
     x = calculate_bin_location(i, dx, min_x, centered_at_zero, extra_phase);
     y = calculate_bin_location(j, dy, min_y, centered_at_zero, extra_phase);
     if(!showGraphics) {
  fprintf(ofile, "%e %e", x, y);
       if(!output_fraction)
  fprintf(ofile, " %ld", distr[j*nrbins+i]);
       else
  fprintf(ofile, " %e", distr[j*nrbins+i]/(double)ndata);
       if(output_sigma) {
  if(distr[j*nrbins+i] == 0)
    s = 1;
  else
    s = sqrt(distr[j*nrbins+i]);
  if(output_fraction)
    s /= (double)ndata;
  fprintf(ofile, " %e", s);
       }
     }
     if(!showGraphics) {
       fprintf(ofile, "\n");
     }
   }
 }
 if(showGraphics) {
   cmap = (float *)malloc(nrbins*nrbinsy*sizeof(float));
   if(cmap == NULL) {
     printerror(application.verbose_state.debug, "Cannot allocate memory.\n");
     return 0;
   }
   if(plotlog) {
     if(distr[0] > 0) {
       min = log10(distr[0]);
       max = log10(distr[0]);
     }else {
       min = max = 0;
     }
   }else {
     min = distr[0];
     max = distr[0];
   }
   for(i = 0; i < nrbins; i++) {
     for(j = 0; j <nrbinsy; j++) {
       if(plotlog) {
  if(distr[i*nrbinsy+j] > 0)
    cmap[i*nrbinsy+j] = log10(distr[i*nrbinsy+j]);
  else
    cmap[i*nrbinsy+j] = 0;
       }
       else
  cmap[i*nrbinsy+j] = distr[i*nrbinsy+j];
       if(cmap[i*nrbinsy+j] > max)
  max = cmap[i*nrbinsy+j];
       if(cmap[i*nrbinsy+j] < min)
  min = cmap[i*nrbinsy+j];
     }
   }
   printf("Count range: %f to %f\n", min, max);
   int showwedge = 0;
   if(pgplotMap(&pgplot_options, cmap, nrbins, nrbinsy,
         calculate_bin_location(0, dx, min_x, centered_at_zero, extra_phase),
         calculate_bin_location(0, dx, max_x, centered_at_zero, extra_phase),
         calculate_bin_location(0, dx, min_x, centered_at_zero, extra_phase)-0.5*dx,
         calculate_bin_location(0, dx, max_x, centered_at_zero, extra_phase)+0.5*dx,
         calculate_bin_location(0, dy, min_y, centered_at_zero, extra_phase),
         calculate_bin_location(0, dy, max_y, centered_at_zero, extra_phase),
         calculate_bin_location(0, dy, min_y, centered_at_zero, extra_phase)-0.5*dy,
         calculate_bin_location(0, dy, max_y, centered_at_zero, extra_phase)+0.5*dy,
         application.cmap, 0, 0, 0, NULL, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, showwedge, 0, 0, 0, application.verbose_state) == 0) {
     printerror(application.verbose_state.debug, "ERROR pdist: Plotting failed.\n");
     return 0;
   }
   free(cmap);
 }
      }else {
 if(select) {
   nrbins = 0;
   for(i = 0; i < ndata; i++) {
     if(data[0][i] >= select1 && data[0][i] <= select2) {
       nrbins++;
     }
   }
   j = 0;
   for(i = 0; i < ndata; i++) {
     if(data[0][i] >= select1 && data[0][i] <= select2) {
       distr[j] = i;
       data[1][j] = data[0][i];
       j++;
     }
   }
 }else {
   if(cdf == 0) {
     for(i = 0; i < ndata; i++) {
       j = calculate_bin_number(data[binnr][i], dx, min_x, centered_at_zero, extra_phase);
       if(j < 0 || j >= nrbins) {
  printerror(application.verbose_state.debug, "BUG! bin number %ld outside range.\n", j);
  exit(0);
       }
       distr[j] += 1;
     }
   }else {
     gsl_sort(data[binnr], 1, ndata);
     if(long_resolved == 0) {
       for(i = 0; i < ndata; i++) {
  data[1][i] = (i+1)/(double)(ndata);
       }
     }
   }
 }
 if(!showGraphics) {
   if(long_resolved) {
   }else {
     for(i = 0; i < nrbins; i++) {
       x = i*dx;
  if(cdf == 0) {
    if(select == 0) {
      fprintf(ofile, "%e", calculate_bin_location(i, dx, min_x, centered_at_zero, extra_phase));
    }else {
      fprintf(ofile, "%ld", distr[i]);
    }
  }else {
    fprintf(ofile, "%e", data[0][i]);
  }
       if(!output_fraction) {
    if(cdf == 0) {
      if(select == 0) {
        fprintf(ofile, " %ld", distr[i]);
      }else {
        fprintf(ofile, " %e", data[1][i]);
      }
    }else {
      fprintf(ofile, " %e", data[1][i]);
    }
       }else {
    if(cdf == 0) {
      fprintf(ofile, " %e", distr[i]/(double)ndata);
    }else {
      fprintf(ofile, " %e", data[1][i]);
    }
       }
       if(output_sigma) {
  if(distr[i] == 0)
    s = 1;
  else
    s = sqrt(distr[i]);
  if(output_fraction)
    s /= (double)ndata;
  fprintf(ofile, " %e", s);
       }
       fprintf(ofile, "\n");
     }
   }
 }else {
   float *distr_x_float;
   distr_float = malloc(nrbins*sizeof(float));
   if(distr_float == NULL) {
     printerror(application.verbose_state.debug, "ERROR pdist: Cannot allocate memory\n");
     return 0;
   }
   if(cdf == 0) {
     for(i = 0; i < nrbins; i++) {
       distr_float[i] = distr[i];
       if(output_fraction)
  distr_float[i] /= (double)ndata;
       if(plotlog) {
  if(distr[i] > 0) {
    distr_float[i] = log10(distr_float[i]);
  }else {
    if(output_fraction) {
      distr_float[i] = log10(1.0/(float)ndata);
    }else {
      distr_float[i] = 0;
    }
  }
       }
     }
   }else {
     distr_x_float = malloc(nrbins*sizeof(float));
     if(distr_x_float == NULL) {
       printerror(application.verbose_state.debug, "ERROR pdist: Cannot allocate memory\n");
       return 0;
     }
     for(i = 0; i < nrbins; i++) {
       distr_float[i] = data[1][i];
       distr_x_float[i] = data[0][i];
     }
   }
   int forceMinZero, dontsetranges, colour;
   forceMinZero = 1;
   if(plotlog)
     forceMinZero = 0;
   dontsetranges = 0;
   colour = 1;
   if(showGraphics == 2) {
     if(currentfile > 1)
       dontsetranges = 1;
     colour = currentfile;
   }
   if(cdf == 0) {
     if(pgplotGraph1(&pgplot_options, distr_float, NULL, NULL, nrbins,
       calculate_bin_location(0, dx, min_x, centered_at_zero, extra_phase),
       calculate_bin_location(nrbins-1, dx, min_x, centered_at_zero, extra_phase),
       dontsetranges,
       calculate_bin_location(0, dx, min_x, centered_at_zero, extra_phase) - dx,
       calculate_bin_location(nrbins-1, dx, min_x, centered_at_zero, extra_phase) + dx,
       0, 0, forceMinZero, 1, 0, linewidth, 0, colour, 1, NULL, -1, application.verbose_state) == 0) {
       printerror(application.verbose_state.debug, "ERROR pdist: Cannot plot graph\n");
       return 0;
     }
   }else {
     if(pgplotGraph1(&pgplot_options, distr_float, distr_x_float, NULL, nrbins,
       0,
       0,
       dontsetranges,
       0,
       0,
       0, 0, forceMinZero, 1*0, 0, linewidth, 0, colour, 1, NULL, -1, application.verbose_state) == 0) {
       printerror(application.verbose_state.debug, "ERROR pdist: Cannot plot graph\n");
       return 0;
     }
     free(distr_x_float);
   }
   free(distr_float);
 }
      }
      free(distr);
    }
    for(i = 0; i < nrdatacolumns; i++) {
      free(data[i]);
    }
    free(data);
    if(showGraphics == 0) {
      if(long_resolved == 0) {
 fclose(ofile);
      }else {
 closePSRData(&dataout, 0, application.verbose_state);
      }
      free(oname);
    }
    if(showGraphics == 1) {
      if(currentfile != numberInApplicationFilenameList(&application, argv, application.verbose_state)) {
 i = pgplot_device_type(application.pgplotdevice, application.verbose_state);
 if(i < 3 || i > 10) {
   printf("Press a key to continue\n");
   fflush(stdout);
   pgetch();
 }
      }
    }
    currentfile += 1;
  }
  if(showGraphics == 2) {
    ppgend();
  }
  terminateApplication(&application);
  return 0;
}

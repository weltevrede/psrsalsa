/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "psrsalsa.h"

int internal_application_cmdline_FilenameList[MaxNrApplicationFilenames];
int internal_application_cmdline_Nrfilenames;
int internal_application_cmdline_CurFilename;
int internal_application_filelist_CurFilename;
long internal_application_filelist_Nrfilenames;
int internal_application_filelist_fileopen;
FILE *internal_application_filelist_fptr;
char *internal_application_filelist_filename;

void initApplication(psrsalsaApplication *application, char *name, char *genusage)
{
  internal_application_cmdline_Nrfilenames = 0;
  internal_application_cmdline_CurFilename = 0;
  internal_application_filelist_CurFilename = 0;
  internal_application_filelist_Nrfilenames = 0;
  internal_application_filelist_fileopen = 0;
  strcpy(application->progname, name);
  application->genusage = malloc(strlen(genusage)+1);
  if(application->genusage == NULL) {
    printerror(0, "ERROR initApplication: Cannot allocate memory.");
    exit(0);
  }
  strcpy(application->genusage, genusage);
  application->switch_verbose = 0;
  application->switch_debug = 0;
  application->switch_nocounters = 0;
  cleanVerboseState(&(application->verbose_state));
  application->switch_formatlist = 0;
  application->switch_iformat = 0;
  application->iformat = -1;
  application->switch_oformat = 0;
  application->oformat = -1;
  application->switch_headerlist = 0;
  application->switch_header = 0;
  application->switch_onpulse = 0;
  application->switch_onpulsef = 0;
  clearRegion(&(application->onpulse));
  application->switch_polselect = 0;
  application->polselectnr = -1;
  application->switch_itf = 0;
  application->itf = 0;
  application->switch_rebin = 0;
  application->dorebin = 0;
  application->switch_nread = 0;
  application->nread = -1;
  application->switch_nskip = 0;
  application->nskip = 0;
  application->switch_conshift = 0;
  application->doconshift = 0;
  application->switch_circshift = 0;
  application->docircshift = 0;
  application->switch_rot = 0;
  application->switch_rotdeg = 0;
  application->doshiftphase = 0;
  application->shiftPhase_cmdline = 0;
  application->shiftPhase = 0;
  application->switch_filelist = 0;
  application->filelist = 0;
  application->switch_device = 0;
  strcpy(application->pgplotdevice, "?");
  application->switch_tscr = 0;
  application->dotscr = 0;
  application->switch_tscr_complete = 0;
  application->tscr_complete = 0;
  application->switch_TSCR = 0;
  application->doTSCR = 0;
  application->switch_fscr = 0;
  application->dofscr = 0;
  application->switch_FSCR = 0;
  application->doFSCR = 0;
  application->switch_dedisperse = 0;
  application->do_dedisperse = 0;
  application->switch_deFaraday = 0;
  application->do_deFaraday = 0;
  application->switch_changeRefFreq = 0;
  application->newRefFreq = -100;
  application->switch_stokes = 0;
  application->dostokes = 0;
  application->switch_coherence = 0;
  application->docoherence = 0;
  application->switch_noweights = 0;
  application->noweights = 0;
  application->switch_scale = 0;
  application->doscale = 0;
  application->switch_debase = 0;
  application->dodebase = 0;
  application->switch_onpulsegr = 0;
  application->doonpulsegr = 0;
  application->switch_size = 0;
  application->windowwidth = -1;
  application->windowheight = -1;
  application->switch_macro = 0;
  application->macro_ptr = NULL;
  application->switch_noplotsubset = 0;
  application->do_noplotsubset = 0;
  application->switch_cmaplist = 0;
  application->switch_cmap = 0;
  application->cmap = PPGPLOT_HEAT;
  application->switch_deparang = 0;
  application->switch_insertparang = 0;
  application->do_parang_corr = 0;
  application->switch_norm = 0;
  application->normvalue = 1;
  application->do_norm = 0;
  application->switch_normglobal = 0;
  application->do_normglobal = 0;
  application->switch_clip = 0;
  application->do_clip = 0;
  application->clipvalue = 0;
  application->switch_fchan = 0;
  application->fchan_select = -1;
  application->switch_history_cmd_only = 0;
  application->history_cmd_only = 0;
  application->switch_fixseed = 0;
  application->fixseed = 0;
  application->switch_template = 0;
  application->template_specified = 0;
  application->switch_templatedata = 0;
  application->template_data_index = 0;
  cleanPSRData(&(application->template_file), application->verbose_state);
  application->switch_align = 0;
  application->doalign = 0;
  application->switch_blocksize = 0;
  application->blocksize = 0;
  application->switch_ext = 0;
  application->extension = NULL;
  application->switch_output = 0;
  sprintf(application->outputname, "output.dat");
  application->switch_shuffle = 0;
  application->doshuffle = 0;
  application->switch_rotateStokes = 0;
  application->nr_rotateStokes = 0;

  application->fzapMask = NULL;
  application->doautot = 0;
}

void printCitationInfo()
{
  printf("If you make use of PSRSALSA, please cite \"Weltevrede 2016, A&A, 590, A109\" and refer to the following website: https://github.com/weltevrede/psrsalsa\n");
}


int parse_command_string(verbose_definition verbose, int argc, char **argv, int argv_index, char *format, ...)
{
  va_list args;
  long i, curargnr, nrarguments;
  void *ptr;

  va_start(args, format);

  if(argv_index >= argc) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR parse_command_string: the command line option to be parsed does not appear to be followed by input values/text.", format);
    return 1;
  }

  if(verbose.debug) {
    printf("parse_command_string: start parsing the \"%s\" option\n", argv[argv_index-1]);
    printf("parse_command_string: parsing \"%s\" as format \"%s\"\n", argv[argv_index], format);
  }
  nrarguments = 0;
  for(i = 0; i < strlen(format); i++) {
    if(format[i] == '%') {
      nrarguments++;
    }
  }
  if(verbose.debug) {
    printf("parse_command_string: expected number of arguments = %ld\n", nrarguments);
  }
  if(nrarguments == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR parse_command_string: format string '%s' suggest no variables need to be parsed.", format);
    return 1;
  }
  for(i = 0; i < nrarguments; i++) {
    ptr = va_arg(args, void *);

    if(ptr == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR parse_command_string: %ld pointers are expected to be passed to this function, but at least one of them appears to be the NULL pointer, suggestive of a bug in the programme.", nrarguments);
      return 1;
    }
  }
  ptr = va_arg(args, void *);

  if(ptr != NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR parse_command_string: %ld pointers are expected to be passed to this function followed by the NULL pointer. However, the NULL pointer is not detected, suggestive of a bug in the programme.", nrarguments);
    return 1;
  }
  va_end(args);

  int nrwords;
  ptr = pickWordFromString(argv[argv_index], 1, &nrwords, 1, ' ', verbose);
  if(nrwords != nrarguments) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR parse_command_string: string \"%s\" appears to contain %d words, while the program expects %ld options to be provided.", argv[argv_index], nrwords, nrarguments);
    return 1;
  }

  va_start(args, format);


  curargnr = 0;
  int expecttype;
  expecttype = 0;
  for(i = 0; i < strlen(format); i++) {
    if(expecttype) {
      char *word_ptr, *word;
      word_ptr = pickWordFromString(argv[argv_index], curargnr, &nrwords, 1, ' ', verbose);
      word = malloc(strlen(word_ptr));
      if(word == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR parse_command_string: Memory allocation error");
 return 1;
      }
      sscanf(word_ptr, "%s", word);
      if(format[i] == 'c') {
 if(verbose.debug) {
   printf("parse_command_string: Parsing \"%s\" as a character\n", word);
 }
 if(strlen(word) != 1) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR parse_command_string: Cannot parse '%s' as a character, since it is not 1 character in length", word);
   return 1;
 }
 char *char_ptr;
 char_ptr = va_arg(args, char *);
 char_ptr[0] = word[0];
 if(verbose.debug) {
   printf("parse_command_string: Parsing \"%s\" as a character '%c'\n", word, char_ptr[0]);
 }
      }else if(format[i] == 'f') {
 if(verbose.debug) {
   printf("parse_command_string: Parsing \"%s\" as a floating point\n", word);
 }

 float *float_ptr;
 char *endptr;
 float_ptr = va_arg(args, float *);
 *float_ptr = strtof(word, &endptr);
 if(endptr == word || endptr == NULL || endptr != word+strlen(word)) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR parse_command_string: Cannot parse '%s' as a floating point", word);
 }
 if(verbose.debug) {
   printf("parse_command_string: Parsing \"%s\" as a %f\n", word, *float_ptr);
 }
      }else {
 fflush(stdout);
 printerror(verbose.debug, "ERROR parse_command_string: format string \"%s\" contains a unrecognized data type '%c'.", format, format[i]);
 free(word);
 return 1;
      }
      expecttype = 0;
      free(word);
    }else {
      if(format[i] == '%') {
 curargnr++;
 expecttype = 1;
      }else if(format[i] == ' ' || format[i] == '\t') {
      }else {
 fflush(stdout);
 printerror(verbose.debug, "ERROR parse_command_string: format string \"%s\" does not appear of the expected format, indicating a bug.", format);
 return 1;
      }
    }
  }

  va_end(args);
  return 0;
}


void printApplicationHelp(psrsalsaApplication application)
{
  fprintf(stdout, "%s %s\n", application.progname, application.genusage);
  if(application.switch_iformat || application.switch_oformat || application.switch_formatlist || application.switch_headerlist || application.switch_header || application.switch_filelist || application.switch_noweights || application.switch_history_cmd_only || application.switch_ext || application.switch_output
) {
    fprintf(stdout, "\nGeneral Input/Output options:\n");
    if(application.switch_filelist) {
      fprintf(stdout, "  -filelist file    Specify file with input file names: only one file name\n");
      fprintf(stdout, "                    per line, and they can be commented out with a #.\n");
    }
    if(application.switch_ext)
      fprintf(stdout, "  -ext ext          Specify output extension ext\n");
    if(application.switch_output) {
      fprintf(stdout, "  -output file      ");
      if(application.switch_ext)
 fprintf(stdout, "instead ");
      fprintf(stdout, "specify output filename, default is %s\n", application.outputname);
    }
    if(application.switch_formatlist)
      fprintf(stdout, "  -formatlist       Show supported file formats\n");
    if(application.switch_iformat)
      fprintf(stdout, "  -iformat id       Specify input format (e.g.  -iformat PSRFITS)\n");
    if(application.switch_oformat)
      fprintf(stdout, "  -oformat id       Specify output format (e.g. -oformat PSRFITS)\n");
    if(application.switch_headerlist)
      fprintf(stdout, "  -headerlist       Show available parameters for the -header option\n");
    if(application.switch_header)
      fprintf(stdout, "  -header           Change header parameter, e.g. -header 'name J0123+4567'\n");
    if(application.switch_noweights)
      fprintf(stdout, "  -noweights        Ignore the weights in the PSRFITS input file\n");
    if(application.switch_history_cmd_only) {
      fprintf(stdout, "  -history_cmd_only Write the history without timestamp, hence re-running the\n");
      fprintf(stdout, "                    same command might result in identical files.\n");
    }
  }
  if(application.switch_templatedata || application.switch_align || application.switch_template
     ) {
    fprintf(stdout, "\nGeneral template options:\n");
    if(application.switch_align) {
      fprintf(stdout, "  -align             Rotate all subints with the same amount to align data\n");
      fprintf(stdout, "                     with template\n");
    }
    if(application.switch_template)
      fprintf(stdout, "  -template file     Use mathematical template (von-Mises ascii file)\n");
    if(application.switch_templatedata)
      fprintf(stdout, "  -templatedata file Use this data file as a template profile\n");
  }
  if(application.switch_polselect || application.switch_rebin || application.switch_nread || application.switch_nskip || application.switch_conshift || application.switch_circshift || application.switch_rot || application.switch_rotdeg || application.switch_tscr || application.switch_TSCR || application.switch_tscr_complete || application.switch_fscr || application.switch_FSCR || application.switch_dedisperse || application.switch_deFaraday || application.switch_stokes || application.switch_coherence || application.switch_changeRefFreq || application.switch_scale || application.switch_debase || application.switch_deparang || application.switch_insertparang || application.switch_norm || application.switch_normglobal || application.switch_fchan || application.switch_blocksize || application.switch_shuffle || application.switch_clip || application.switch_rotateStokes
) {
    fprintf(stdout, "\nGeneral preprocess options:\n");
    if(application.switch_blocksize)
      fprintf(stdout, "  -blocksize b    Only read in a multiple of b subintegrations.\n");
    if(application.switch_clip) {
      fprintf(stdout, "  -clip           Specify threshhold value above which the value clipped\n");
      fprintf(stdout, "                  Clipping happens if sample > threshold or < -threshold.\n");
    }
    if(application.switch_dedisperse)
      fprintf(stdout, "  -dedisp         Dedisperse the data.\n");
    if(application.switch_deFaraday) {
      fprintf(stdout, "  -defarad        De-Faraday rotate the data.\n");
    }
    if(application.switch_deparang)
      fprintf(stdout, "  -deparang       Remove parallactic angle effect from data\n");
    if(application.switch_coherence)
      fprintf(stdout, "  -coherence      Convert Stokes to coherency parameters\n");
    if(application.switch_debase)
      fprintf(stdout, "  -debase         Subtract baseline from data (use with -onpulse)\n");
    if(application.switch_fchan)
      fprintf(stdout, "  -fchan f        Use frequency channel f only\n");
    if(application.switch_norm) {
      fprintf(stdout, "  -norm           Normalize peak value of each subint/channel independently. If\n");
      fprintf(stdout, "                  an onpulse region is set, the peak of that region is used. The\n");
      fprintf(stdout, "                  normalization is determined by the first polarization channel.\n");
    }
    if(application.switch_normglobal) {
      fprintf(stdout, "  -norm_global    As -norm, but use global factor for all subints/channels.\n");
    }
    if(application.switch_nskip)
      fprintf(stdout, "  -nskip n        Skip n subintegrations in input data (default is 0)\n");
    if(application.switch_nread)
      fprintf(stdout, "  -nread n        Use n subintegrationss in input data (default is all)\n");
    if(application.switch_insertparang)
      fprintf(stdout, "  -parang         Insert parallactic angle effect (-deparang corrects data)\n");
    if(application.switch_polselect)
      fprintf(stdout, "  -polselect p    Use polarization channel p (start counting from 0)\n");
    if(application.switch_rebin)
      fprintf(stdout, "  -rebin n        Rebin to n bins (not necessarily a power of two)\n");
    if(application.switch_conshift) {
      fprintf(stdout, "  -conshift       When applying rotations (-rot or -rotdeg), instead of rotating\n");
      fprintf(stdout, "                  subints independent of each other, rotate the end of one\n");
      fprintf(stdout, "                  subint to the start of another. One subint will be lost.\n");
    }
    if(application.switch_circshift) {
      fprintf(stdout, "  -circshift      idem, but makes the last subint spill over in the first.\n");
      fprintf(stdout, "                  No subints will be lost.\n");
    }
    if(application.switch_changeRefFreq) {
      fprintf(stdout, "  -reffreq f      Change the reference frequency for dedispersion/de-Faraday\n");
      fprintf(stdout, "                  rotation to specified value f (in MHz, -1=inf frequency).\n");
      fprintf(stdout, "                  Data will be re-dedispersed/de-Farady rotated if required.\n");
    }
    if(application.switch_rot)
      fprintf(stdout, "  -rot ph         Rotate each individual subint by ph pulse phase\n");
    if(application.switch_rotdeg)
      fprintf(stdout, "  -rotdeg ph      ditto, but in degrees\n");
    if(application.switch_rotateStokes) {
      fprintf(stdout, "  -rotateStokes   \"S1 S2 ph\" Rotate the polarization vectors in Stokes space\n");
      fprintf(stdout, "                  (Q,U,V) by ph degrees. S1 and S2 specify the direction of the\n");
      fprintf(stdout, "                  rotation, such that for instance \"Q U 10\" would be a\n");
      fprintf(stdout, "                  rotation about the Stokes V axis such that a vector in the\n");
      fprintf(stdout, "                  Q direction gets rotated towards the U axis. This option can\n");
      fprintf(stdout, "                  be used multiple times, and they are executed in the order\n");
      fprintf(stdout, "                  as specified on the command line.\n");
    }
    if(application.switch_scale)
      fprintf(stdout, "  -scale          \"scale offset\". output = scale*(input+offset)\n");
    if(application.switch_shuffle)
      fprintf(stdout, "  -shuffle        Shuffle the subints in a random order\n");
    if(application.switch_stokes)
      fprintf(stdout, "  -stokes         Convert to Stokes parameters\n");
    if(application.switch_tscr) {
      fprintf(stdout, "  -tscr t         Add t successive subints together.\n");
      fprintf(stdout, "                  A negative number duplicates subints.\n");
    }
    if(application.switch_TSCR)
      fprintf(stdout, "  -TSCR           Add all successive subints together\n");
    if(application.switch_tscr_complete) {
      fprintf(stdout, "  -tscr_complete  Use in combination with -tscr. Throw away last subints to\n");
      fprintf(stdout, "                  ensure each subint is the sum of the same number of.\n");
      fprintf(stdout, "                  input subints.\n");
    }
    if(application.switch_fscr) {
      fprintf(stdout, "  -fscr f         Add f successive frequency channels together.\n");
      fprintf(stdout, "                  Data will be dedispersed and de-Faraday rotated if required.\n");
      fprintf(stdout, "                  A negative number duplicates frequency channels.\n");
    }
    if(application.switch_FSCR) {
      fprintf(stdout, "  -FSCR           Add all successive frequency channels together.\n");
      fprintf(stdout, "                  Data will be dedispersed and de-Faraday rotated if required.\n");
    }
  }
  if(application.switch_itf || application.switch_device || application.switch_size || application.switch_noplotsubset || application.switch_cmaplist || application.switch_cmap
     ) {
    fprintf(stdout, "\nGeneral plotting options:\n");
    if(application.switch_cmap)
      fprintf(stdout, "  -cmap type    Select color map type\n");
    if(application.switch_cmaplist)
      fprintf(stdout, "  -cmaplist     List available color map types\n");
    if(application.switch_device)
      fprintf(stdout, "  -device dev   (or -dev) Specify PGPLOT plotting device dev\n");
    if(application.switch_noplotsubset)
      fprintf(stdout, "  -noplotsubset Disable plotting subsets of data (this could increase file size)\n");
    if(application.switch_itf) {
      fprintf(stdout, "  -itf n        Set image transfer function for colour map plots\n");
      fprintf(stdout, "                (0=linear (default), 1=logarithmic, 2=square-root)\n");
    }
    if(application.switch_size)
      fprintf(stdout, "  -size         \"width height\". Specify resolution of plot device (in pixels).\n");
  }
  if(application.switch_onpulse || application.switch_onpulsef || application.switch_onpulsegr
) {
    fprintf(stdout, "\nGeneral data selection options:\n");
   if(application.switch_onpulse)
     fprintf(stdout, "  -onpulse      \"left right\" manually select on-pulse regions (in bins)\n");
   if(application.switch_onpulsef)
     fprintf(stdout, "  -onpulsef     \"left right\" manually select on-pulse regions (in phase)\n");
   if(application.switch_onpulsegr)
      fprintf(stdout, "  -onpulsegr    Graphically select (additional) onpulse regions\n");
  }
  if(application.switch_verbose || application.switch_debug || application.switch_nocounters || application.switch_macro || application.switch_fixseed
     ) {
    fprintf(stdout, "\nOther general options:\n");
    if(application.switch_verbose)
      fprintf(stdout, "  -v            Verbose mode (to get a better idea what is happening)\n");
    if(application.switch_debug)
      fprintf(stdout, "  -debug        Enable more output (where implemented)\n");
    if(application.switch_fixseed) {
      fprintf(stdout, "  -fixseed      Do not randomy initialise seed of random number generator\n");
      fprintf(stdout, "                thereby making results reproducable.\n");
    }
    if(application.switch_nocounters)
      fprintf(stdout, "  -nocounters   Don't show counters etc (useful when generating log files)\n");
    if(application.switch_macro) {
      fprintf(stdout, "  -macro        Instead of taking commands from keyboard, read them from\n");
      fprintf(stdout, "                this macro file (put a ^ in front of symbol for the ctrl key)\n");
    }
  }
  fprintf(stdout, "\n");
}


int processCommandLine(psrsalsaApplication *application, int argc, char **argv, int *index)
{
  int j;
  if(strcmp(argv[*index], "-v") == 0 && application->switch_verbose) {
    application->verbose_state.verbose = 1;
    return 1;
  }else if(strcmp(argv[*index], "-debug") == 0 && application->switch_debug) {
    application->verbose_state.debug = 1;
    if(application->switch_verbose)
      application->verbose_state.verbose = 1;
    return 1;
  }else if(strcmp(argv[*index], "-iformat") == 0 && application->switch_iformat) {
    application->iformat = parsePSRDataFormats(argv[++(*index)]);
    if(application->iformat == 0)
      exit(0);
    return 1;
  }else if(strcmp(argv[*index], "-oformat") == 0 && application->switch_oformat) {
    application->oformat = parsePSRDataFormats(argv[++(*index)]);
    if(application->oformat == 0)
      exit(0);
    return 1;
  }else if(strcmp(argv[*index], "-formatlist") == 0 && application->switch_formatlist) {
    fprintf(stdout, "Supported file formats are:\n");
    printPSRDataFormats(stdout, 2);
    exit(0);
  }else if(strcmp(argv[*index], "-headerlist") == 0 && application->switch_headerlist) {
    printHeaderCommandlineOptions(stdout);
    exit(0);
  }else if(strcmp(argv[*index], "-gentypelist") == 0 && application->switch_headerlist) {
    printHeaderGentypeOptions(stdout);
    exit(0);
  }else if(strcmp(argv[*index], "-header") == 0 && application->switch_header) {
    (*index)++;
    return 1;
  }else if(strcmp(argv[*index], "-macro") == 0 && application->switch_macro) {
    printf("Opening macro '%s'\n", argv[++(*index)]);
    application->macro_ptr = fopen(argv[*index], "r");
    if(application->macro_ptr == NULL) {
      printerror(application->verbose_state.debug, "  Opening macro '%s' failed", argv[*index]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-cmaplist") == 0 && application->switch_cmaplist) {
    printCMAPCommandlineOptions(stdout);
    exit(0);
  }else if(strcmp(argv[*index], "-filelist") == 0 && application->switch_filelist) {
    (*index) += 1;
    application->filelist = *index;
    return 1;
  }else if(strcmp(argv[*index], "-nread") == 0 && application->switch_nread) {
    j = sscanf(argv[++(*index)], "%ld", &application->nread);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one integer.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-polselect") == 0 && application->switch_polselect) {
    j = sscanf(argv[++(*index)], "%d", &application->polselectnr);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one integer", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-nskip") == 0 && application->switch_nskip) {
    j = sscanf(argv[++(*index)], "%ld", &application->nskip);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one integer.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-template") == 0 && application->switch_template) {
    application->template_specified = ++(*index);
    if(readVonMisesModel(argv[application->template_specified], &(application->vonMises_components), application->verbose_state) == 0) {
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-templatedata") == 0 && application->switch_templatedata) {
    application->template_data_index = ++(*index);

    cleanPSRData(&(application->template_file), application->verbose_state);
    if(openPSRData(&(application->template_file), argv[application->template_data_index], 0, 0, 1, 0, application->verbose_state) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot open template file.");
      exit(0);
    }
    datafile_definition clone;
    if(application->template_file.NrSubints > 1) {
      if(!preprocess_addsuccessivepulses(application->template_file, &clone, application->template_file.NrSubints, 1, application->verbose_state))
 return 0;
      swap_orig_clone(&(application->template_file), &clone, application->verbose_state);
    }
    if(application->template_file.NrFreqChan > 1) {
      if(!preprocess_dedisperse(&(application->template_file), 0, 0, application->verbose_state))
 return 0;
      if(!preprocess_addsuccessiveFreqChans(application->template_file, &clone, application->template_file.NrFreqChan, NULL, application->verbose_state))
 return 0;
      swap_orig_clone(&(application->template_file), &clone, application->verbose_state);
    }
    if(application->template_file.NrPols > 1) {
      printwarning(application->verbose_state.debug, "WARNING processCommandLine: Polarization channel 0 is used as a template.");
      if(preprocess_polselect(application->template_file, &clone, 0, application->verbose_state) == 0)
 return 0;
      swap_orig_clone(&(application->template_file), &clone, application->verbose_state);
    }
    return 1;
  }else if(strcmp(argv[*index], "-onpulse") == 0 && application->switch_onpulse) {
    j = sscanf(argv[++(*index)], "%d %d", &(application->onpulse.left_bin[application->onpulse.nrRegions]), &(application->onpulse.right_bin[application->onpulse.nrRegions]));
    if(j != 2) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need two integers.", argv[(*index)-1]);
      exit(0);
    }
    application->onpulse.bins_defined[application->onpulse.nrRegions] = 1;
    (application->onpulse.nrRegions)++;
    if(application->onpulse.nrRegions == maxNrRegions) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "processCommandLine: To many regions selected.");
      exit(-1);
    }
    return 1;
  }else if(strcmp(argv[*index], "-onpulsef") == 0 && application->switch_onpulsef) {
    j = sscanf(argv[++(*index)], "%f %f", &(application->onpulse.left_frac[application->onpulse.nrRegions]), &(application->onpulse.right_frac[application->onpulse.nrRegions]));
    if(j != 2) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need two integers.", argv[(*index)-1]);
      exit(0);
    }
    application->onpulse.frac_defined[application->onpulse.nrRegions] = 1;
    (application->onpulse.nrRegions)++;
    return 1;
  }else if(strcmp(argv[*index], "-ext") == 0 && application->switch_ext) {
    application->extension = argv[++(*index)];
    return 1;
  }else if(strcmp(argv[*index], "-output") == 0 && application->switch_output) {
    strcpy(application->outputname, argv[++(*index)]);
    return 1;
  }else if((strcmp(argv[*index], "-device") == 0 || strcmp(argv[*index], "-dev") == 0) && application->switch_device) {
    strcpy(application->pgplotdevice, argv[++(*index)]);
    return 1;
  }else if(strcmp(argv[*index], "-fchan") == 0 && application->switch_fchan) {
    j = sscanf(argv[++(*index)], "%d", &application->fchan_select);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one integer.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-tscr_complete") == 0 && application->switch_tscr_complete) {
    application->tscr_complete = 1;
    return 1;
  }else if(strcmp(argv[*index], "-tscr") == 0 && application->switch_tscr) {
    j = sscanf(argv[++(*index)], "%ld", &application->dotscr);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one integer.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-fscr") == 0 && application->switch_fscr) {
    j = sscanf(argv[++(*index)], "%ld", &application->dofscr);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one integer.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-cmap") == 0 && application->cmap) {
    application->cmap = cmap_parse_commandline(argc, argv, application->verbose_state.debug);
    (*index)++;
    return 1;
  }else if(strcmp(argv[*index], "-rebin") == 0 && application->switch_rebin) {
    j = sscanf(argv[++(*index)], "%ld", &application->rebin);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one integer.", argv[(*index)-1]);
      exit(0);
    }
    application->dorebin = 1;
    return 1;
  }else if(strcmp(argv[*index], "-rot") == 0 && application->switch_rot) {
    j = sscanf(argv[++(*index)], "%f", &application->shiftPhase);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one float.", argv[(*index)-1]);
      exit(0);
    }
    application->shiftPhase_cmdline = application->shiftPhase;
    application->doshiftphase = 1;
    return 1;
  }else if(strcmp(argv[*index], "-rotdeg") == 0 && application->switch_rotdeg) {
    j = sscanf(argv[++(*index)], "%f", &application->shiftPhase);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one float.", argv[(*index)-1]);
      exit(0);
    }
    application->doshiftphase = 1;
    application->shiftPhase /= 360.0;
    application->shiftPhase_cmdline = application->shiftPhase;
    return 1;
  }else if(strcmp(argv[*index], "-itf") == 0 && application->switch_itf) {
    j = sscanf(argv[++(*index)], "%d", &application->itf);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one integer.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-size") == 0 && application->switch_size) {
    j = sscanf(argv[++(*index)], "%d %d", &application->windowwidth, &application->windowheight);
    if(j != 2) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one integer.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-noplotsubset") == 0 && application->switch_noplotsubset) {
    application->do_noplotsubset = 1;
    return 1;
  }else if(strcmp(argv[*index], "-nocounters") == 0 && application->switch_nocounters) {
    application->verbose_state.nocounters = 1;
    return 1;
  }else if(strcmp(argv[*index], "-debase") == 0 && application->switch_debase) {
    application->dodebase = 1;
    return 1;
  }else if(strcmp(argv[*index], "-norm") == 0 && application->switch_norm) {
    application->do_norm = 1;
    return 1;
  }else if(strcmp(argv[*index], "-norm_global") == 0 && application->switch_normglobal) {
    application->do_normglobal = 1;
    return 1;
  }else if(strcmp(argv[*index], "-TSCR") == 0 && application->switch_TSCR) {
    application->doTSCR = 1;
    return 1;
  }else if(strcmp(argv[*index], "-FSCR") == 0 && application->switch_FSCR) {
    application->doFSCR = 1;
    return 1;
  }else if(strcmp(argv[*index], "-clip") == 0 && application->switch_clip) {
    j = sscanf(argv[++(*index)], "%f", &application->clipvalue);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one floating point value.", argv[(*index)-1]);
      exit(0);
    }
    application->do_clip = 1;
    return 1;
  }else if(strcmp(argv[*index], "-reffreq") == 0 && application->switch_changeRefFreq) {
    j = sscanf(argv[++(*index)], "%lf", &application->newRefFreq);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one floating point value.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-dedisp") == 0 && application->switch_dedisperse) {
    application->do_dedisperse = 1;
    return 1;
  }else if(strcmp(argv[*index], "-defarad") == 0 && application->switch_deFaraday) {
    application->do_deFaraday = 1;
    return 1;
  }else if(strcmp(argv[*index], "-blocksize") == 0 && application->switch_blocksize) {
    j = sscanf(argv[++(*index)], "%d", &application->blocksize);
    if(j != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one integer.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-align") == 0 && application->switch_align) {
    application->doalign = 1;
    return 1;
  }else if(strcmp(argv[*index], "-conshift") == 0 && application->switch_conshift) {
    application->doconshift = 1;
    return 1;
  }else if(strcmp(argv[*index], "-circshift") == 0 && application->switch_circshift) {
    application->doconshift = 1;
    application->docircshift = 1;
    return 1;
  }else if(strcmp(argv[*index], "-onpulsegr") == 0 && application->switch_onpulsegr) {
    application->doonpulsegr = 1;
    return 1;
  }else if(strcmp(argv[*index], "-stokes") == 0 && application->switch_stokes) {
    application->dostokes = 1;
    return 1;
  }else if(strcmp(argv[*index], "-coherence") == 0 && application->switch_coherence) {
    application->docoherence = 1;
    return 1;
  }else if(strcmp(argv[*index], "-shuffle") == 0 && application->switch_shuffle) {
    application->doshuffle = 1;
    return 1;
  }else if(strcmp(argv[*index], "-fixseed") == 0 && application->switch_fixseed) {
    application->fixseed = 1;
    return 1;
  }else if(strcmp(argv[*index], "-scale") == 0 && application->switch_scale) {
    application->doscale = 1;
    j = sscanf(argv[++(*index)], "%f %f", &application->scale_scale, &application->scale_offset);
    if(j != 2) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need two floats.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcasecmp(argv[*index], "-rotateStokes") == 0 && application->switch_rotateStokes) {
    if(application->nr_rotateStokes == maxNrRotateStokes) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Maximum number of uses of the '%s' option is exceeded.", argv[(*index)-1]);
      exit(0);
    }


    char s1, s2;
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), "%c %c %f", &s1, &s2, &(application->rotateStokesAngle[application->nr_rotateStokes]), NULL)) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option, need one float.", argv[(*index)-1]);
      exit(0);
    }

    if(s1 == 'i' || s1 == 'I') {
      application->rotateStokes1[application->nr_rotateStokes] = 0;
    }else if(s1 == 'q' || s1 == 'Q') {
      application->rotateStokes1[application->nr_rotateStokes] = 1;
    }else if(s1 == 'u' || s1 == 'U') {
      application->rotateStokes1[application->nr_rotateStokes] = 2;
    }else if(s1 == 'v' || s1 == 'V') {
      application->rotateStokes1[application->nr_rotateStokes] = 3;
    }else {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option: cannot interpret '%c' as a Stokes parameter.", argv[(*index)-1], s1);
      exit(0);
    }
    if(s2 == 'i' || s2 == 'I') {
      application->rotateStokes2[application->nr_rotateStokes] = 0;
    }else if(s2 == 'q' || s2 == 'Q') {
      application->rotateStokes2[application->nr_rotateStokes] = 1;
    }else if(s2 == 'u' || s2 == 'U') {
      application->rotateStokes2[application->nr_rotateStokes] = 2;
    }else if(s2 == 'v' || s2 == 'V') {
      application->rotateStokes2[application->nr_rotateStokes] = 3;
    }else {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option: cannot interpret '%c' as a Stokes parameter.", argv[(*index)-1], s2);
      exit(0);
    }
    if(application->rotateStokes1[application->nr_rotateStokes] == application->rotateStokes2[application->nr_rotateStokes]) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option: Cannot rotate a Stokes parameter towards the same Stokes parameter.", argv[(*index)-1]);
      exit(0);
    }
    (application->nr_rotateStokes)++;
    return 1;
  }else if(application->switch_deparang && strcmp(argv[*index], "-deparang") == 0) {
    application->do_parang_corr = 1;
    return 1;
  }else if(application->switch_insertparang && strcmp(argv[*index], "-parang") == 0) {
    application->do_parang_corr = 2;
    return 1;
  }else if(strcmp(argv[*index], "-noweights") == 0 && application->switch_noweights) {
    application->noweights = 1;
    psrfits_set_noweights(1);
    return 1;
  }else if(strcmp(argv[*index], "-history_cmd_only") == 0 && application->switch_history_cmd_only) {
    application->history_cmd_only = 1;
    return 1;
  }else {
    return 0;
  }
  fflush(stdout);
  printwarning(application->verbose_state.debug, "WARNING processCommandLine: This line shouldn't be executed!");
  return 0;
}
int getOutputName(psrsalsaApplication application, char *filename, char *outputname, verbose_definition verbose)
{
  if(application.extension != NULL) {
    if(change_filename_extension(filename, outputname, application.extension, MaxOutputNameLength, verbose) == 0) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "getOutputName: Changing filename failed");
      return 0;
    }
  }else {
    strcpy(outputname, application.outputname);
  }
  return 1;
}
int preprocessApplication(psrsalsaApplication *application, datafile_definition *psrdata)
{
  datafile_definition clone;
  int device, original_gentype, original_poltype, original_isDeDisp, original_isDeFarad, original_isDePar, original_isDebase;
  double original_freq_ref;
  float x;
  verbose_definition verbose1, verbose2;
  long i;
  pgplot_viewport_def viewport;
  original_gentype = psrdata->gentype;
  original_poltype = psrdata->poltype;
  original_isDeDisp = psrdata->isDeDisp;
  original_isDeFarad = psrdata->isDeFarad;
  original_isDePar = psrdata->isDePar;
  original_isDebase = psrdata->isDebase;
  original_freq_ref = psrdata->freq_ref;
  if(application->verbose_state.verbose) {
    printf("\nApplying preprocess options\n");
  }
  copyVerboseState(application->verbose_state, &verbose1);
  copyVerboseState(application->verbose_state, &verbose2);
  verbose1.indent = application->verbose_state.indent + 2;
  verbose2.indent = application->verbose_state.indent + 4;
  if(application->nskip != 0 || application->nread > 0) {
    if(application->nread <= 0)
      application->nread = psrdata->NrSubints-application->nskip;
    if(preprocess_pulsesselect(*psrdata, &clone, application->nskip, application->nread, verbose1) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Error selecting pulses.");
      return 0;
    }
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(application->dostokes) {
    if(preprocess_stokes(psrdata, verbose1) == 0)
      return 0;
  }
  if(application->docoherence) {
    if(preprocess_coherency(psrdata, verbose1) == 0)
      return 0;
  }
  if(application->nr_rotateStokes > 0) {
    for(i = 0; i < application->nr_rotateStokes; i++) {
      if(preprocess_rotateStokes(psrdata, &clone, 1, -1, application->rotateStokesAngle[i], application->rotateStokes1[i], application->rotateStokes2[i], verbose1) == 0)
 return 0;
    }
  }
  if(application->do_parang_corr > 0) {
    if(application->do_parang_corr == 2) {
      if(preprocess_corrParAng(psrdata, NULL, 1, verbose1) == 0)
 return 0;
    }else {
      if(preprocess_corrParAng(psrdata, NULL, 0, verbose1) == 0)
 return 0;
    }
  }
  if(application->blocksize > 0) {
    if(preprocess_blocksize(*psrdata, &clone, application->blocksize, verbose1) == 0)
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(application->fchan_select != -1) {
    if(preprocess_channelselect(*psrdata, &clone, application->fchan_select, verbose1) == 0)
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(application->polselectnr >= 0) {
    if(preprocess_polselect(*psrdata, &clone, application->polselectnr, verbose1) == 0)
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(application->newRefFreq > -2) {
    if(preprocess_changeRefFreq(psrdata, application->newRefFreq, verbose1) == 0) {
      printerror(application->verbose_state.debug, "preprocessApplication: Error changing reference frequency.");
      return 0;
    }
  }
  if(application->doFSCR) {
    application->dofscr = psrdata->NrFreqChan;
    if(psrdata->NrFreqChan <= 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: -FSCR expect a positive number of channels. Maybe header parameters not loaded by main program?");
      return 0;
    }
  }
  if(application->do_dedisperse || application->dofscr) {
    if(!preprocess_dedisperse(psrdata, 0, 0, verbose1))
      return 0;
  }
  if(application->do_deFaraday || application->dofscr) {
    int skip;
    skip = 0;
    if(psrdata->NrPols != 4) {
      if(application->do_deFaraday == 0) {
 skip = 1;
      }
    }
    if(skip == 0) {
      if(application->do_deFaraday == 2 && application->dofscr == 0) {
 if(!preprocess_deFaraday(psrdata, 1, 0, 0, verbose1))
   return 0;
      }else {
 if(application->do_deFaraday == 2) {
   fflush(stdout);
   printerror(application->verbose_state.debug, "preprocessApplication: You cannot frequency scrunch and use the -farad option at the same time.");
   return 0;
  }
 if(!preprocess_deFaraday(psrdata, 0, 0, 0, verbose1))
   return 0;
      }
    }
  }
  if(application->dofscr) {
    if(!preprocess_addsuccessiveFreqChans(*psrdata, &clone, application->dofscr, application->fzapMask, verbose1))
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(application->doTSCR) {
    application->dotscr = psrdata->NrSubints;
    if(psrdata->NrSubints <= 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: -TSCR expect a positive number of subints. Maybe header parameters not loaded by main program?");
      return 0;
    }
  }
  if(application->dotscr) {
    if(!preprocess_addsuccessivepulses(*psrdata, &clone, application->dotscr, application->tscr_complete, verbose1))
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(application->doalign) {
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)
 printf(" ");
      printf("Aligning data using template\n");
    }
    if(application->template_specified == 0 && application->template_data_index == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Can only use -align option together with a specified template on the command line.");
      return 0;
    }
    if(!preprocess_addsuccessivepulses(*psrdata, &clone, psrdata->NrSubints, application->tscr_complete, verbose2))
      return 0;
    if(!preprocess_dedisperse(&clone, 0, 0, verbose2))
      return 0;
    if(psrdata->NrPols == 4) {
      if(!preprocess_deFaraday(&clone, 0, 0, 0, verbose2))
 return 0;
    }
    datafile_definition clone2;
    if(!preprocess_addsuccessiveFreqChans(clone, &clone2, clone.NrFreqChan, NULL, verbose2))
      return 0;
    if(application->template_specified) {
      x = correlateVonMisesFunction(application->vonMises_components, clone2.NrBins, clone2.data, verbose2);
    }else {
      if(clone2.NrBins != application->template_file.NrBins) {
 fflush(stdout);
 printerror(application->verbose_state.debug, "preprocessApplication: The template and the data file have a different amount of bins (%ld != %ld).", clone2.NrBins, application->template_file.NrBins);
 return 0;
      }
      int lag;
      float correl_max;
      if(find_peak_correlation(clone2.data, application->template_file.data, clone2.NrBins, 0, 0, 1, &lag, &correl_max, verbose2) == 0) {
 return 0;
      }
      x = lag/(double)clone2.NrBins;
    }
    closePSRData(&clone2, verbose2);
    closePSRData(&clone, verbose2);
    if(application->doshiftphase) {
      application->shiftPhase -= x;
    }else {
      application->doshiftphase = 1;
      application->shiftPhase = -x;
    }
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)
 printf(" ");
      printf("  done       \n");
    }
  }
  if(application->doshiftphase) {
    if(application->doconshift) {
      i = application->shiftPhase*psrdata->NrBins;
     if(i >= psrdata->NrBins)
       i -= psrdata->NrBins;
     if(i < 0)
       i += psrdata->NrBins;
     if(verbose1.verbose) {
       int j;
       for(j = 0; j < verbose1.indent; j++)
  printf(" ");
       printf("Rotating data by %ld bins\n", i);
     }
     if(continuous_shift(*psrdata, &clone, i, application->docircshift, "preprocessApplication", MEMORY_format, 0, NULL, verbose2, 0) == 0)
 return 0;
      swap_orig_clone(psrdata, &clone, application->verbose_state);
    }else {
      if(preprocess_fftshift(*psrdata, application->shiftPhase, 0, 0, verbose2) == 0)
 return 0;
    }
  }
  if(application->dorebin) {
    if(!preprocess_rebin(*psrdata, &clone, application->rebin, verbose1))
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(application->doonpulsegr) {
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)
 printf(" ");
      printf("Select more onpulse regions\n");
    }
    if(!preprocess_addsuccessivepulses(*psrdata, &clone, psrdata->NrSubints, application->tscr_complete, verbose2))
      return 0;
    ppgqid(&device);
    printf("Device = %d\n", device);
    pgplot_clear_viewport_def(&viewport);
    pgplot_box_def pgplotbox;
    clear_pgplot_box(&pgplotbox);
    strcpy(pgplotbox.xlabel, "Bin");
    strcpy(pgplotbox.ylabel, "Intensity");
    strcpy(pgplotbox.title, "Select on-pulse region of ");
    strcat(pgplotbox.title, psrdata->psrname);
    selectRegions(clone.data, clone.NrBins, viewport, pgplotbox, 0, 0, 0, &(application->onpulse), verbose2);
    if(device)
      ppgslct(device);
    closePSRData(&clone, verbose2);
    regionShowNextTimeUse(application->onpulse, "-onpulse", "-onpulsef", stdout);
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)
 printf(" ");
      printf("  done\n");
    }
  }
  if(application->dodebase) {
    if(!preprocess_debase(psrdata, application->onpulse, verbose1))
      return 0;
  }
  if(application->do_norm) {
    if(preprocess_norm(*psrdata, application->normvalue, &(application->onpulse), 0, verbose1) == 0)
      return 0;
  }
  if(application->do_normglobal) {
    if(preprocess_norm(*psrdata, application->normvalue, &(application->onpulse), 1, verbose1) == 0)
      return 0;
  }
  if(application->do_clip) {
    if(preprocess_clip(*psrdata, application->clipvalue, verbose1) == 0)
      return 0;
  }
  if(application->doscale) {
    if(preprocess_scale(*psrdata, application->scale_scale, application->scale_offset, verbose1) == 0)
      return 0;
  }
  if(application->doshuffle) {
    if(preprocess_shuffle(*psrdata, &clone, application->fixseed, verbose1) == 0)
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(original_gentype != psrdata->gentype && application->verbose_state.debug) {
    for(i = 0; i < verbose1.indent; i++)
      printf(" ");
    printf("Gentype of data changed from %s into %s\n", returnGenType_str(original_gentype), returnGenType_str(psrdata->gentype));
  }
  if(original_poltype != psrdata->poltype || original_isDeDisp != psrdata->isDeDisp || original_isDeFarad != psrdata->isDeFarad || original_isDePar != psrdata->isDePar || original_freq_ref != psrdata->freq_ref) {
    fflush(stdout);
    char *txt, *txt2;
    txt = malloc(10000);
    txt2 = malloc(10000);
    if(txt == NULL || txt2 == NULL) {
      printerror(application->verbose_state.debug, "ERROR preprocessApplication: Memory allocation error");
      return 0;
    }
    sprintf(txt, "WARNING: Note that after preprocessing the data are now:");
    if(original_poltype != psrdata->poltype) {
      if(psrdata->poltype == POLTYPE_STOKES) {
 strcat(txt, " Stokes parameters");
      }else if(psrdata->poltype == POLTYPE_COHERENCY) {
 strcat(txt, " Coherency parameters");
      }else if(psrdata->poltype == POLTYPE_ILVPAdPA) {
 strcat(txt, " I, L, V, PA and its error");
      }else if(psrdata->poltype == POLTYPE_PAdPA) {
 strcat(txt, " PA and its error");
      }else {
 strcat(txt, " unknown polarisation state");
      }
    }
    if(original_isDeDisp != psrdata->isDeDisp) {
      if(psrdata->isDeDisp) {
 strcat(txt, " dedispersed");
      }else {
 strcat(txt, " not dedispersed");
      }
    }
    if(original_isDeFarad != psrdata->isDeFarad) {
      if(psrdata->isDeFarad) {
 strcat(txt, " de-Faraday rotated");
      }else {
 strcat(txt, " not de-Faraday rotated");
      }
    }
    if(original_isDePar != psrdata->isDePar) {
      if(psrdata->isDePar) {
 strcat(txt, " de-parallactic angle rotated");
      }else {
 strcat(txt, " not de-parallactic angle rotated");
      }
    }
    if(original_isDebase != psrdata->isDebase) {
      if(psrdata->isDebase) {
 strcat(txt, " with baseline removed");
      }else {
 strcat(txt, " without removed baseline");
      }
    }
    if(original_freq_ref != psrdata->freq_ref) {
      if(psrdata->freq_ref > 0)
 sprintf(txt2, " with reference freq=%lf MHz", psrdata->freq_ref);
      else if(psrdata->freq_ref > -1.01 && psrdata->freq_ref < -0.99)
 sprintf(txt2, " with reference freq=infinity");
      else
 sprintf(txt2, " with unknown reference freq");
      strcat(txt, txt2);
    }
    printwarning(application->verbose_state.debug, "%s", txt);
    free(txt);
    free(txt2);
  }
  if(verbose1.verbose) {
    printf("Preprocessing done\n\n");
  }
  if(application->verbose_state.debug) {
    printHeaderPSRData(*psrdata, 1, application->verbose_state);
  }
  return 1;
}
int applicationAddFilename(int argi, verbose_definition verbose)
{
  if(internal_application_cmdline_Nrfilenames < MaxNrApplicationFilenames) {
    internal_application_cmdline_FilenameList[internal_application_cmdline_Nrfilenames++] = argi;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "applicationAddFilename: Too many file names on command-line (hard-coded limit is currently %d). Consider using -filelist.", MaxNrApplicationFilenames);
    return 0;
  }
  return 1;
}
int applicationFilenameList_checkConsecutive(char **argv, verbose_definition verbose)
{
  int i, j, indx_last;
  char txt[1000], txt2[1000];
  if(internal_application_cmdline_Nrfilenames < 2)
    return 1;
  indx_last = internal_application_cmdline_FilenameList[0];
  for(i = 1; i < internal_application_cmdline_Nrfilenames; i++) {
    if(internal_application_cmdline_FilenameList[i] != indx_last+1) {
      if(argv != NULL) {
 fflush(stdout);
 sprintf(txt, "Cannot parse command line. The following list of arguments (those not recognized as command line options) was interpreted as a list of file names:\n\n");
 for(j = 0; j < internal_application_cmdline_Nrfilenames; j++) {
   fflush(stdout);
   sprintf(txt2, "'%s' ", argv[internal_application_cmdline_FilenameList[j]]);
   strcat(txt, txt2);
 }
 fflush(stdout);
 strcat(txt, "\n");
 strcat(txt, "\nThese should however be consecutive, so something appears to be wrong.");
 printerror(verbose.debug, "%s", txt);
      }
      return 0;
    }
    indx_last++;
  }
  return 1;
}
int internal_open_filelist(psrsalsaApplication application, char **argv, verbose_definition verbose)
{
  int nrColumns;
  internal_application_filelist_fptr = fopen(argv[application.filelist], "r");
  if(internal_application_filelist_fptr == NULL) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "internal_open_filelist: Cannot open %s", argv[application.filelist]);
    return 0;
  }
  internal_application_filelist_fileopen = 1;
  nrColumns = -1;
  if(ascii_file_stats(internal_application_filelist_fptr, '#', &internal_application_filelist_Nrfilenames, 2048, 0, &nrColumns, verbose) == 0) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "internal_open_filelist: Determining the number of lines in %s failed", argv[application.filelist]);
    return 0;
  }
  rewind(internal_application_filelist_fptr);
  if(application.verbose_state.verbose) {
    printf("Opened input file name file '%s' with %ld file names\n", argv[application.filelist], internal_application_filelist_Nrfilenames);
  }
  return 1;
}
int numberInApplicationFilenameList(psrsalsaApplication application, char **argv, verbose_definition verbose)
{
  long total;
  total = internal_application_cmdline_Nrfilenames;
  if(application.filelist) {
    if(internal_application_filelist_fileopen == 0) {
      if(internal_open_filelist(application, argv, verbose) == 0) {
 fflush(stdout);
 printerror(application.verbose_state.debug, "numberInApplicationFilenameList: Opening list with file names failed");
 return 0;
      }
    }
    total += internal_application_filelist_Nrfilenames;
  }
  return total;
}
void rewindFilenameList(psrsalsaApplication application) {
  internal_application_cmdline_CurFilename = 0;
  internal_application_filelist_CurFilename = 0;
  if(application.filelist)
    rewind(internal_application_filelist_fptr);
}
char *getNextFilenameFromList(psrsalsaApplication *application, char **argv, verbose_definition verbose)
{
  int i;
  if(application->doautot != 0 && internal_application_cmdline_CurFilename != 0) {
    application->template_specified = 0;
  }
  application->shiftPhase = application->shiftPhase_cmdline;
  if(internal_application_cmdline_CurFilename < internal_application_cmdline_Nrfilenames) {
    return argv[internal_application_cmdline_FilenameList[internal_application_cmdline_CurFilename++]];
  }
  if(application->filelist != 0) {
    if(internal_application_filelist_CurFilename < internal_application_filelist_Nrfilenames) {
      if(internal_application_filelist_CurFilename == 0)
 internal_application_filelist_filename = malloc(2048);
      if(ascii_file_get_next_line(internal_application_filelist_fptr, internal_application_filelist_filename, 2048, '#', verbose) == 0) {
 fflush(stdout);
 printerror(application->verbose_state.debug, "ERROR getNextFilenameFromList: reading from %s failed", argv[application->filelist]);
 return 0;
      }
      for(i = 0; i < strlen(internal_application_filelist_filename); i++) {
 if(internal_application_filelist_filename[i] == '\n') {
   internal_application_filelist_filename[i] = 0;
   break;
 }else if(internal_application_filelist_filename[i] == '\r') {
   internal_application_filelist_filename[i] = 0;
   break;
 }else if(internal_application_filelist_filename[i] == ' ') {
   internal_application_filelist_filename[i] = 0;
   break;
 }
      }
      internal_application_filelist_CurFilename++;
      return internal_application_filelist_filename;
    }
  }
  return NULL;
}

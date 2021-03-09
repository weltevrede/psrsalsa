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
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <pwd.h>
#include <unistd.h>
#include <netdb.h>
#include <termios.h>
#include <sys/ioctl.h>
#include "psrsalsa.h"
#include <fcntl.h>
#include <signal.h>
#include <stdarg.h>
#include <errno.h>
#include <math.h>
int pgetch(void)
{
  struct termios oldattr, newattr;
  int ch;
  tcgetattr(STDIN_FILENO, &oldattr);
  newattr = oldattr;
  newattr.c_lflag &= ~( ICANON | ECHO | ISIG);
  tcsetattr(STDIN_FILENO, TCSANOW, &newattr);
  int gotkey = 0;
  int n;
  while(gotkey == 0) {
    if(ioctl(fileno(stdin), FIONREAD, &n) == 0 && n > 0) {
      ch = getchar();
      gotkey = 1;
    }
    if(gotkey == 0) {
      usleep(10000);
    }
  }
  tcsetattr(STDIN_FILENO, TCSANOW, &oldattr);
  if(ch == 3) {
    fflush(stdout);
    fprintf(stderr, "pgetch: caught control-c\n");
    fprintf(stderr, "Terminating program\n");
    exit(0);
  }
  if(ch == 26) {
    fflush(stdout);
    fprintf(stderr, "pgetch: caught control-z\n");
    fprintf(stderr, "Sending suspend signal. Program is not terminated.\n");
    raise(SIGTSTP);
  }
  return ch;
}
int pgetch_macro(psrsalsaApplication *application, verbose_definition verbose)
{
  int key;
  int ctrlstate;
  ctrlstate = 0;
  do {
    if(application->macro_ptr == NULL) {
      return pgetch();
    }else {
      key = fgetc(application->macro_ptr);
      if(key == EOF) {
 printf("\nReached end of macro, switching to keyboard input\n");
 fclose(application->macro_ptr);
 application->macro_ptr = NULL;
 ctrlstate = 0;
      }
    }
    if(key == '^') {
      if(ctrlstate == 0) {
 ctrlstate = 1;
      }else {
 fflush(stdout);
 printerror(verbose.debug, "ERROR pgetch_macro: key sequence ^^ is not allowed");
 exit(0);
      }
    }
  }while(key == '\n' || key == '\r' || key == '^' || key == EOF);
  if(ctrlstate) {
    if(key == 'c' || key == 'C') {
      fflush(stdout);
      fprintf(stderr, "pgetch_macro: caught control-c\n");
      exit(0);
    }
    if(key >= 'a' && key <= 'z') {
      return key - 96;
    }else if(key >= 'A' && key <= 'Z') {
      return key - 64;
    }else if(key == 32) {
      return 0;
    }
  }
  return key;
}
int getUsername(char **username, verbose_definition verbose)
{
  struct passwd pw, *result;
  char *tmpmem;
  tmpmem = malloc(MaxStringLength);
  if(tmpmem == NULL) {
    printwarning(verbose.debug, "ERROR getUsername: Memory allocation error.");
    *username = malloc(8);
    if(*username != NULL) {
      sprintf(*username, "Unknown");
    }
    return 0;
  }
  getpwuid_r(geteuid(), &pw, tmpmem, MaxStringLength, &result);
  if(result != NULL) {
    *username = malloc(strlen(pw.pw_name)+1);
    if(*username == NULL) {
      printwarning(verbose.debug, "ERROR getUsername: Memory allocation error.");
      return 0;
    }
    strcpy(*username, pw.pw_name);
    if(verbose.debug) {
      printf("Username: %s\n", *username);
    }
    free(tmpmem);
    return 1;
  }
  free(tmpmem);
  char *user;
  if((user = getenv ("USER")) == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR: getUsername failed to retrieve user name from USER environment variable.");
    *username = malloc(8);
    if(*username != NULL) {
      sprintf(*username, "Unknown");
    }
    return 0;
  }else {
    *username = malloc(strlen(user)+1);
    if(*username == NULL) {
      printwarning(verbose.debug, "ERROR getUsername: Memory allocation error.");
      return 0;
    }
    strcpy(*username, user);
    if(verbose.debug) {
      printf("Username: %s\n", *username);
    }
  }
  return 1;
}
void constructCommandLineString(char *txt, int length, int argc, char **argv, verbose_definition verbose)
{
  int i;
  txt[0] = 0;
  for(i = 0; i < argc; i++) {
    if(strlen(txt) + strlen(argv[i]) + 4 > length-1) {
      printwarning(verbose.debug, "WARNING constructCommandLineString: Truncating command line which is too long.");
      break;
    }
    if(strchr(argv[i], ' ') == NULL) {
      strcat(txt, argv[i]);
    }else {
      strcat(txt, "\"");
      strcat(txt, argv[i]);
      strcat(txt, "\"");
    }
    if(i != argc-1)
      strcat(txt, " ");
  }
}
int getMachinename(char *hostname, int size, verbose_definition verbose)
{
  hostname[size-1] = '\0';
  if(gethostname(hostname, size-1) == 0) {
    if(verbose.debug) {
      printf("Hostname: %s\n", hostname);
    }
    return 1;
  }else {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING getMachinename failed to find the machine name.");
    sprintf(hostname, "Unknown");
    return 0;
  }
}
char * pickWordFromString(char *string, int n, int *nrwords, int replacetabs, char separator, verbose_definition verbose)
{
  char *ptr, *ret, *field, *string_mod, *format;
  int nrchars;
  if(n <= 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pickWordFromString: requested word number (%d) is not valid", n);
    if(verbose.debug) {
      printerror(verbose.debug, "ERROR pickWordFromString: string=%s", string);
    }
    exit(0);
  }
  string_mod = malloc(strlen(string)+1);
  field = malloc(MaxPickWordFromString_WordLength+1);
  format = malloc(20);
  if(string_mod == NULL || field == NULL || format == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pickWordFromString: Memory allocation error");
    exit(0);
  }
  strcpy(string_mod, string);
  if(replacetabs) {
    for(nrchars = 0; nrchars < strlen(string_mod); nrchars++) {
      if(string_mod[nrchars] == '\t')
 string_mod[nrchars] = ' ';
    }
  }
  for(nrchars = strlen(string_mod)-1; nrchars >= 0; nrchars--) {
    if(string_mod[nrchars] == '\n')
      string_mod[nrchars] = 0;
    else if(string_mod[nrchars] == '\r')
      string_mod[nrchars] = 0;
    else if(string_mod[nrchars] == ' ')
      string_mod[nrchars] = 0;
    else
      break;
  }
  ptr = string_mod;
  *nrwords = 0;
  ret = NULL;
  do {
    while(*ptr == ' ')
      ++ptr;
    if(*ptr == 0)
      break;
    nrchars = 0;
    sprintf(format, "%%%d[^%c]%%n", MaxPickWordFromString_WordLength, separator);
    sscanf(ptr, format, field, &nrchars);
    (*nrwords) ++;
    if(*nrwords == n)
      ret = ptr;
    ptr += nrchars;
    if ( *ptr != separator ) {
      break;
    }
    ++ptr;
  }while(1);
  if(ret != NULL)
    ret = ret-string_mod+string;
  free(string_mod);
  free(field);
  free(format);
  return ret;
}
int change_filename_extension(char *inputname, char *outputname, char *extension, int outputnamelength, verbose_definition verbose)
{
  int i;
  for(i = strlen(inputname); i >= 0; i--) {
    if(inputname[i] == '.') {
      break;
    }
  }
  if(i <= 0) {
    fflush(stdout);
    printerror(verbose.debug, "change_filename_extension: no extension in '%s'?", inputname);
    return 0;
  }
  if(strlen(extension) + i+1 >= outputnamelength) {
    fflush(stdout);
    printerror(verbose.debug, "change_filename_extension: outputnamelength too long");
    return 0;
  }
  memcpy(outputname, inputname, i+1);
  outputname[i+1] = 0;
  strcat(outputname, extension);
  return 1;
}
int skipLinesInFile(FILE *fptr, int skiplines, verbose_definition verbose)
{
  long i, j;
  if(skiplines > 0) {
    if(verbose.verbose)
      printf("Skipping %d lines\n", skiplines);
    for(i = 0; i < skiplines; i++) {
      do {
 j = fgetc(fptr);
 if(j == EOF) {
   return 0;
 }
      }while(j !='\n');
    }
  }
  return 1;
}
int ascii_file_get_next_line(FILE *fin, char *txt, int maxlinelength, int skipChar, verbose_definition verbose)
{
  char *ret_ptr;
  int ret;
  ret = 0;
  do {
    ret_ptr = fgets(txt, maxlinelength, fin);
    if(ret_ptr == NULL) {
      return 0;
    }else {
      ret++;
      if(txt[0] != skipChar) {
 return ret;
      }
    }
  }while(txt[0] == skipChar);
  printerror(verbose.debug, "ascii_file_get_next_line: BUG!!!");
  exit(0);
}
int ascii_file_stats(FILE *fin, char skipChar, long *nrlines, int maxlinelength, int autoNrColumns, int *nrColumns, verbose_definition verbose)
{
  char *txt;
  int nrwords, ret;
  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ascii_file_stats: Cannot allocate temporary memory");
    return 0;
  }
  if(verbose.debug) {
    printf("ascii_file_stats: skipChar=%c\n", skipChar);
  }
  *nrlines = 0;
  do {
    ret = ascii_file_get_next_line(fin, txt, maxlinelength, skipChar, verbose);
    if(verbose.debug) {
      printf("ascii_file_stats: line %ld=%s\n", *nrlines, txt);
    }
    if(ret != 0) {
      pickWordFromString(txt, 1, &nrwords, 1, ' ', verbose);
      if(autoNrColumns != 0 && *nrlines == 0) {
 *nrColumns = nrwords;
      }else {
 if(nrColumns != NULL) {
   if(*nrColumns >= 0) {
     if(*nrColumns != nrwords) {
       fflush(stdout);
       printerror(verbose.debug, "ascii_file_stats: Nr of columns on line %ld is not the expected %d. Possibly the number of columns is changing from line to line?", (*nrlines)+1, *nrColumns);
       free(txt);
       return 0;
     }
   }else {
     if(nrwords < -(*nrColumns)) {
       fflush(stdout);
       printerror(verbose.debug, "ascii_file_stats: Nr of columns on line %ld is smaller than the expected %d. Possibly the number of columns is changing from line to line?", (*nrlines)+1, *nrColumns);
       free(txt);
       return 0;
     }
   }
 }
      }
      (*nrlines)++;
    }
  }while(ret != 0);
  free(txt);
  return 1;
}
int read_ascii_column(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, long nrdatapoints_max, int colnum, double scale, double norm, int read_log, float **data, float *mindata, float *maxdata, float *avdata, verbose_definition verbose, int verbose_stderr)
{
  FILE *fin, *verbose_stream;
  long i, j, n, fpos, maxlinelength, linenr;
  int nrwords, ret;
  char *txt, *word_ptr;
  double minx, maxx, sumx;
  maxlinelength = 10240+1;
  if(verbose_stderr) {
    fflush(stdout);
    verbose_stream = stderr;
  }else {
    verbose_stream = stdout;
  }
  fin = fopen(fname, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column: Cannot open %s", fname);
    return 0;
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "Opened file '%s'\n", fname);
    }
  }
  if(skipLinesInFile(fin, skiplines, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column: Reached EOF while skipping first %d lines", skiplines);
    return 0;
  }
  fpos = ftell(fin);
  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column: Cannot allocate temporary memory");
    return 0;
  }
  if(ascii_file_stats(fin, skipChar, nrdatapoints, maxlinelength, autoNrColumns, &nrColumns, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column: Error in determining the nr of lines");
    free(txt);
    return 0;
  }
  if(nrdatapoints_max > 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "  Capping the %ld datapoints to %ld\n", *nrdatapoints, nrdatapoints_max);
    }
    if(*nrdatapoints > nrdatapoints_max) {
      *nrdatapoints = nrdatapoints_max;
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  There are %ld datapoints\n", *nrdatapoints);
    if(autoNrColumns) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fprintf(verbose_stream, "  There are %d columns\n", nrColumns);
    }
  }
  *data = (float *)malloc((*nrdatapoints)*sizeof(float));
  if(*data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column: Memory allocation error");
    free(txt);
    return 0;
  }
  fseek(fin, fpos, SEEK_SET);
  minx = maxx = NAN;
  sumx = 0;
  n = 0;
  linenr = 0;
  do {
    ret = ascii_file_get_next_line(fin, txt, maxlinelength, skipChar, verbose);
    if(ret != 0) {
      linenr += ret;
      word_ptr = pickWordFromString(txt, colnum, &nrwords, 1, ' ', verbose);
      if(n >= *nrdatapoints) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column: Nr of lines in file changed?????");
 free(txt);
 free(*data);
 return 0;
      }
      if(word_ptr == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column: Cannot find column %d on line %ld", colnum, linenr);
 free(txt);
 free(*data);
 return 0;
      }
      j = sscanf(word_ptr, "%f", &((*data)[n]));
      if(j != 1) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column: Cannot interpret column %d on line %ld as a float", colnum, linenr);
 free(txt);
 free(*data);
 return 0;
      }
      (*data)[n] *= scale;
      if(read_log) {
 if((*data)[n] <= 0) {
   printerror(verbose.debug, "read_ascii_column: Cannot take logarithm of a value <= 0");
   return 0;
 }
 (*data)[n] = log10((*data)[n]);
      }
      if((*data)[n] < minx || n == 0) {
 minx = (*data)[n];
      }
      if((*data)[n] > maxx || n == 0) {
 maxx = (*data)[n];
      }
      sumx += (*data)[n];
      n++;
    }
    if(nrdatapoints_max > 0) {
      if(n == nrdatapoints_max) {
 break;
      }
    }
  }while(ret != 0);
  free(txt);
  fclose(fin);
  if(verbose.verbose) {
    fflush(stdout);
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  %ld points loaded from %s with values between %lf and %lf\n", *nrdatapoints, fname, minx, maxx);
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    fprintf(verbose_stream, "  Average value = %lf\n", sumx/(double)(*nrdatapoints));
  }
  if(norm != 0.0) {
    scale = 1.0/(sumx/(double)(*nrdatapoints));
    for(n = 0; n < *nrdatapoints; n++) {
      (*data)[n] *= scale;
    }
    minx *= scale;
    maxx *= scale;
    sumx *= scale;
  }
  if(mindata != NULL)
    *mindata = minx;
  if(maxdata != NULL)
    *maxdata = maxx;
  if(avdata != NULL)
    *avdata = sumx/(double)(*nrdatapoints);
  return 1;
}
int read_ascii_column_double(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, long nrdatapoints_max, int colnum, double scale, double norm, int read_log, double **data, double *mindata, double *maxdata, double *avdata, verbose_definition verbose, int verbose_stderr)
{
  FILE *fin, *verbose_stream;
  long i, j, n, fpos, maxlinelength, linenr;
  int nrwords, ret;
  char *txt, *word_ptr;
  double minx, maxx, sumx;
  maxlinelength = 10240+1;
  if(verbose_stderr) {
    fflush(stdout);
    verbose_stream = stderr;
  }else {
    verbose_stream = stdout;
  }
  fin = fopen(fname, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_double: Cannot open %s", fname);
    return 0;
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "Opened file '%s'\n", fname);
    }
  }
  if(skipLinesInFile(fin, skiplines, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_double: Reached EOF while skipping first %d lines", skiplines);
    return 0;
  }
  fpos = ftell(fin);
  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_double: Cannot allocate temporary memory");
    return 0;
  }
  if(ascii_file_stats(fin, skipChar, nrdatapoints, maxlinelength, autoNrColumns, &nrColumns, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_double: Error in determining the nr of lines");
    free(txt);
    return 0;
  }
  if(nrdatapoints_max > 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "  Capping the %ld datapoints to %ld\n", *nrdatapoints, nrdatapoints_max);
    }
    if(*nrdatapoints > nrdatapoints_max) {
      *nrdatapoints = nrdatapoints_max;
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  There are %ld datapoints\n", *nrdatapoints);
    if(autoNrColumns) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fprintf(verbose_stream, "  There are %d columns\n", nrColumns);
    }
  }
  *data = (double *)malloc((*nrdatapoints)*sizeof(double));
  if(*data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_double: Memory allocation error");
    free(txt);
    return 0;
  }
  fseek(fin, fpos, SEEK_SET);
  minx = maxx = NAN;
  sumx = 0;
  n = 0;
  linenr = 0;
  do {
    ret = ascii_file_get_next_line(fin, txt, maxlinelength, skipChar, verbose);
    if(ret != 0) {
      linenr += ret;
      word_ptr = pickWordFromString(txt, colnum, &nrwords, 1, ' ', verbose);
      if(n >= *nrdatapoints) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column_double: Nr of lines in file changed?????");
 free(txt);
 free(*data);
 return 0;
      }
      if(word_ptr == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column_double: Cannot find column %d on line %ld", colnum, linenr);
 free(txt);
 free(*data);
 return 0;
      }
      j = sscanf(word_ptr, "%lf", &((*data)[n]));
      if(j != 1) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column_double: Cannot interpret column %d on line %ld as a double", colnum, linenr);
 free(txt);
 free(*data);
 return 0;
      }
      (*data)[n] *= scale;
      if(read_log) {
 if((*data)[n] <= 0) {
   printerror(verbose.debug, "read_ascii_column: Cannot take logarithm of a value <= 0");
   return 0;
 }
 (*data)[n] = log10((*data)[n]);
      }
      if((*data)[n] < minx || n == 0) {
 minx = (*data)[n];
      }
      if((*data)[n] > maxx || n == 0) {
 maxx = (*data)[n];
      }
      sumx += (*data)[n];
      n++;
    }
    if(nrdatapoints_max > 0) {
      if(n == nrdatapoints_max) {
 break;
      }
    }
  }while(ret != 0);
  free(txt);
  fclose(fin);
  if(verbose.verbose) {
    fflush(stdout);
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  %ld points loaded from %s with values between %lf and %lf\n", *nrdatapoints, fname, minx, maxx);
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    fprintf(verbose_stream, "  Average value = %lf\n", sumx/(double)(*nrdatapoints));
  }
  if(norm != 0.0) {
    scale = 1.0/(sumx/(double)(*nrdatapoints));
    for(n = 0; n < *nrdatapoints; n++) {
      (*data)[n] *= scale;
    }
    minx *= scale;
    maxx *= scale;
    sumx *= scale;
  }
  if(mindata != NULL)
    *mindata = minx;
  if(maxdata != NULL)
    *maxdata = maxx;
  if(avdata != NULL)
    *avdata = sumx/(double)(*nrdatapoints);
  return 1;
}
int read_ascii_column_int(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, long nrdatapoints_max, int colnum, int **data, int *mindata, int *maxdata, double *avdata, verbose_definition verbose, int verbose_stderr)
{
  FILE *fin, *verbose_stream;
  long i, j, n, fpos, maxlinelength, linenr;
  int nrwords, ret;
  char *txt, *word_ptr;
  long minx, maxx, sumx;
  maxlinelength = 10240+1;
  if(verbose_stderr) {
    fflush(stdout);
    verbose_stream = stderr;
  }else {
    verbose_stream = stdout;
  }
  fin = fopen(fname, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_int: Cannot open %s", fname);
    return 0;
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "Opened file '%s'\n", fname);
    }
  }
  if(skipLinesInFile(fin, skiplines, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_int: Reached EOF while skipping first %d lines", skiplines);
    return 0;
  }
  fpos = ftell(fin);
  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_int: Cannot allocate temporary memory");
    return 0;
  }
  if(ascii_file_stats(fin, skipChar, nrdatapoints, maxlinelength, autoNrColumns, &nrColumns, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_int: Error in determining the nr of lines");
    free(txt);
    return 0;
  }
  if(nrdatapoints_max > 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "  Capping the %ld datapoints to %ld\n", *nrdatapoints, nrdatapoints_max);
    }
    if(*nrdatapoints > nrdatapoints_max) {
      *nrdatapoints = nrdatapoints_max;
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  There are %ld datapoints\n", *nrdatapoints);
    if(autoNrColumns) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fprintf(verbose_stream, "  There are %d columns\n", nrColumns);
    }
  }
  *data = (int *)malloc((*nrdatapoints)*sizeof(int));
  if(*data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_int: Memory allocation error");
    free(txt);
    return 0;
  }
  fseek(fin, fpos, SEEK_SET);
  minx = maxx = 0;
  sumx = 0;
  n = 0;
  linenr = 0;
  do {
    ret = ascii_file_get_next_line(fin, txt, maxlinelength, skipChar, verbose);
    if(ret != 0) {
      linenr += ret;
      word_ptr = pickWordFromString(txt, colnum, &nrwords, 1, ' ', verbose);
      if(n >= *nrdatapoints) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column_int: Nr of lines in file changed?????");
 free(txt);
 free(*data);
 return 0;
      }
      if(word_ptr == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column_int: Cannot find column %d on line %ld", colnum, linenr);
 free(txt);
 free(*data);
 return 0;
      }
      j = sscanf(word_ptr, "%d", &((*data)[n]));
      if(j != 1) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column_int: Cannot interpret column %d on line %ld as a int", colnum, linenr);
 free(txt);
 free(*data);
 return 0;
      }
      if((*data)[n] < minx || n == 0) {
 minx = (*data)[n];
      }
      if((*data)[n] > maxx || n == 0) {
 maxx = (*data)[n];
      }
      sumx += (*data)[n];
      n++;
    }
    if(nrdatapoints_max > 0) {
      if(n == nrdatapoints_max) {
 break;
      }
    }
  }while(ret != 0);
  free(txt);
  fclose(fin);
  if(verbose.verbose) {
    fflush(stdout);
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  %ld points loaded from %s with values between %ld and %ld\n", *nrdatapoints, fname, minx, maxx);
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    fprintf(verbose_stream, "  Average value = %lf\n", sumx/(double)(*nrdatapoints));
  }
  if(mindata != NULL)
    *mindata = minx;
  if(maxdata != NULL)
    *maxdata = maxx;
  if(avdata != NULL)
    *avdata = sumx/(double)(*nrdatapoints);
  return 1;
}
int read_ascii_column_str(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, long nrdatapoints_max, int colnum, char ***data, verbose_definition verbose, int verbose_stderr)
{
  FILE *fin, *verbose_stream;
  long i, n, fpos, maxlinelength, linenr;
  int nrwords, ret;
  char *txt, *word, *word_ptr;
  maxlinelength = 10240+1;
  if(verbose_stderr) {
    fflush(stdout);
    verbose_stream = stderr;
  }else {
    verbose_stream = stdout;
  }
  fin = fopen(fname, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_str: Cannot open %s", fname);
    return 0;
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "Opened file '%s'\n", fname);
    }
  }
  if(skipLinesInFile(fin, skiplines, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_str: Reached EOF while skipping first %d lines", skiplines);
    return 0;
  }
  fpos = ftell(fin);
  txt = malloc(maxlinelength);
  word = malloc(maxlinelength);
  if(txt == NULL || word == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_str: Cannot allocate temporary memory");
    return 0;
  }
  if(ascii_file_stats(fin, skipChar, nrdatapoints, maxlinelength, autoNrColumns, &nrColumns, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_str: Error in determining the nr of lines");
    free(txt);
    free(word);
    return 0;
  }
  if(nrdatapoints_max > 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "  Capping the %ld datapoints to %ld\n", *nrdatapoints, nrdatapoints_max);
    }
    if(*nrdatapoints > nrdatapoints_max) {
      *nrdatapoints = nrdatapoints_max;
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  There are %ld datapoints\n", *nrdatapoints);
    if(autoNrColumns) {
      for(i = 0; i < verbose.indent; i++)
 printf(" ");
      fprintf(verbose_stream, "  There are %d columns\n", nrColumns);
    }
  }
  *data = (char **)malloc((*nrdatapoints)*sizeof(char *));
  if(*data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_str: Memory allocation error");
    free(word);
    free(txt);
    return 0;
  }
  fseek(fin, fpos, SEEK_SET);
  n = 0;
  linenr = 0;
  do {
    ret = ascii_file_get_next_line(fin, txt, maxlinelength, skipChar, verbose);
    if(ret != 0) {
      linenr += ret;
      word_ptr = pickWordFromString(txt, colnum, &nrwords, 1, ' ', verbose);
      if(n >= *nrdatapoints) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column_str: Nr of lines in file changed?????");
 free(word);
 free(txt);
 free(*data);
 return 0;
      }
      if(word_ptr == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column_str: Cannot find column %d on line %ld", colnum, linenr);
 free(word);
 free(txt);
 free(*data);
 return 0;
      }
      sscanf(word_ptr, "%s", word);
      (*data)[n] = (char *)malloc(strlen(word)+1);
      if((*data)[n] == NULL) {
 fflush(stdout);
 printerror(verbose.debug, "read_ascii_column_str: Memory allocation error");
 free(word);
 free(txt);
 return 0;
      }
      strcpy((*data)[n], word);
      n++;
    }
    if(nrdatapoints_max > 0) {
      if(n == nrdatapoints_max) {
 break;
      }
    }
  }while(ret != 0);
  free(word);
  free(txt);
  fclose(fin);
  if(verbose.verbose) {
    fflush(stdout);
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  %ld words loaded from %s\n", *nrdatapoints, fname);
  }
  return 1;
}
char *str_replace(char *orig, char *rep, char *with, verbose_definition verbose)
{
  char *result;
  char *ins;
  char *tmp;
  int len_rep;
  int len_with;
  int len_front;
  int count;
  int new_string_length;
  if(orig == NULL)
    return NULL;
  if(rep == NULL)
    rep = "";
  len_rep = strlen(rep);
  if(with == NULL)
    with = "";
  len_with = strlen(with);
  ins = orig;
  for(count = 0; (tmp = strstr(ins, rep)) != NULL; ++count) {
    ins = tmp + len_rep;
  }
  new_string_length = strlen(orig) + (len_with - len_rep) * count + 1;
  result = malloc(new_string_length);
  tmp = result;
  if(result == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR str_replace: Cannot allocate memory");
    return NULL;
  }
  while(count--) {
    ins = strstr(orig, rep);
    len_front = ins - orig;
    tmp = strncpy(tmp, orig, len_front) + len_front;
    tmp = strcpy(tmp, with) + len_with;
    orig += len_front + len_rep;
  }
  strcpy(tmp, orig);
  return result;
}
void fprintf_color(FILE *destination, int color, const char *format, ...)
{
  va_list args;
  if(isatty(fileno(destination))) {
    switch(color) {
    case 2: fprintf(destination, "\x1B[31m"); break;
    case 3: fprintf(destination, "\x1B[32m"); break;
    case 4: fprintf(destination, "\x1B[33m"); break;
    case 5: fprintf(destination, "\x1B[34m"); break;
    case 6: fprintf(destination, "\x1B[35m"); break;
    case 7: fprintf(destination, "\x1B[36m"); break;
    case 8: fprintf(destination, "\x1B[37m"); break;
    }
  }
  va_start(args, format);
  vfprintf(destination, format, args);
  if(isatty(fileno(destination))) {
    fprintf(destination, "\x1B[0m");
  }
}

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
#include <gsl/gsl_statistics_double.h>
#include "psrsalsa.h"
#define KSTEST 1
#define KSFLAT 2
#define KSSIN 3
#define CORREL 20
#define PEARSON 21
#define MOMENTS 30
#define MEDIAN 31
#define CHI2TEST_HIST 50
#define CHI2TEST_CDF 51
int main(int argc, char **argv)
{
  psrsalsaApplication application;
  long i, j;
  int file1_column1, file1_column2, file1_column3, file2_column1, file2_column2, file2_column3, typetest, read_log, output_idx;
  double threshold1, threshold2, threshold3;
  initApplication(&application, "pstat", "[options] inputfile(s)");
  application.switch_libversions = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  file1_column1 = 0;
  file1_column2 = 0;
  file1_column3 = 0;
  file2_column1 = 0;
  file2_column2 = 0;
  file2_column3 = 0;
  typetest = 0;
  read_log = 0;
  threshold1 = 0;
  threshold2 = 0;
  threshold3 = 0;
  output_idx = 0;
  if(argc < 2) {
    printf("Program to perform various statistical tests on input data. One or two input\n");
    printf("files are required depending on the statistical test to be performed. The input\n");
    printf("files should be ascii files with each line having an equal number of columns.\n");
    printf("Lines starting with a # will be ignored. Usage:\n\n");
    printApplicationHelp(&application);
    printf("Input/Output options:\n");
    printf("-col \"c1 c2 c3\"  Specify one, two or three column numbers (counting from 1) to\n");
    printf("                 be read in from the first input file. Some statistical tests\n");
    printf("                 require only one input column, others two or more.\n");
    printf("-col1            Equivalent to -col\n");
    printf("-col2            Like -col1, but for second input file.\n");
    printf("-log             The base 10 log of the input values is used.\n");
    printf("-output          Specify output filename to use rather than the stdout.\n");
    printf("\nStatistical tests:\n");
    printf("-chi2hist \"t1 t2 t3\" Input should be two files, each being a histogram and\n");
    printf("                     having two (or three, see below) columns: bin location and\n");
    printf("                 the height. The histograms should have equal bin widths,\n");
    printf("                 overlap at least partially and have aligned bins. Only bins\n");
    printf("                 in the 1st input file with a height >= t1 are considered (so if\n");
    printf("                 negative, everything should be included. t2 is the theshold for\n");
    printf("                 the 2nd file, while t3 is that for the sum of heights of the\n");
    printf("                 two histograms. Specifying t2 and t3 is optional. The optional\n");
    printf("                 3rd column is the std. dev. on each bin, which is set to zero\n");
    printf("                 for bins outside the covered range of the histograms, and the\n");
    printf("                 height is set to zero as well. With only two columns, the\n");
    printf("                 chi-square test is unweighted, i.e. each bin has an equal\n");
    printf("                 weight (set to 1) regardless of the height of the bin(s).\n");
    printf("                 Examples: \n");
    printf("                 pstat -chi2hist \"-1 -1 0.1\" -col1 \"1 2\" -col2 \"1 2\" file1 file2\n");
    printf("                 pstat -chi2hist -1 -col1 \"1 2 3\" -col2 \"1 2 3\" file1 file2\n");
    printf("-chi2cdf         This test is similar to -chi2hist, except that the inputs are\n");
    printf("                 not histograms, but the distribution of points itself (i.e. one\n");
    printf("                 column). This test therefore does not depend on binning. The\n");
    printf("                 vertical difference in their CDF is quantified with an non-\n");
    printf("                 weighted chi-square test (weights are set to 1). This is done\n");
    printf("                 for each point in the CDF of the shortest distribution and\n");
    printf("                 using a linear interpolation of the other. This test can be\n");
    printf("                 useful to decide which model distribution fits an observed\n");
    printf("                 distribution best, but the chi square itself is not immediately\n");
    printf("                 meaningful. Example: pstat -chi2cdf file1 file2\n");
    printf("-correl          Produces the cross correlation function (as function of lag)\n");
    printf("                 calculated in the Fourier domain assuming equal sampling.\n");
    printf("                 Examples: pstat -correl -col1 1 -col2 1 file1 file2\n");
    printf("                 or:       pstat -correl -col \"1 2\" file1\n");
    printf("-ks              Kolmogorov-Smirnov test: Assess difference between two\n");
    printf("                 distributions of points. The made approximations are similar to\n");
    printf("                 the implementation in Numerical Recipes in C, 2nd edition and.\n");
    printf("                 does not depend on any binning.\n");
    printf("                 Examples: pstat -ks -col1 1 -col2 1 file1 file2\n");
    printf("                 or:       pstat -ks -col \"1 2\" file1\n");
    printf("-ksflat          Like -ks, but compare single input distribution with a flat\n");
    printf("                 distribution covering the range of input values.\n");
    printf("                 Example:  pstat -ksflat -col 1 file1\n");
    printf("-kssin           Like -ks, but compare single input distribution with a\n");
    printf("                 sinusoidal distribution between 0 and 90 deg.\n");
    printf("                 Example:  pstat -kssin -col 1 file1\n");
    printf("-median          Compute the median of the distribution.\n");
    printf("                 Example:  pstat -median -col 1 file1\n");
    printf("-moments         Compute different moments of the distribution (mean, variance\n");
    printf("                 etc.)\n");
    printf("                 Example:  pstat -moments -col 1 file1\n");
    printf("-pearson         Pearson product-moment correlation coefficient:\n");
    printf("                 Measures the linear correlation between two variables.\n");
    printf("                 Examples: pstat -pearson -col1 1 -col2 1 file1 file2\n");
    printf("                 or:       pstat -pearson -col \"1 2\" file1\n");
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
      }else if(strcmp(argv[i], "-col") == 0 || strcmp(argv[i], "-col1") == 0) {
 int ret;
 ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%d %d %d", &file1_column1, &file1_column2, &file1_column3, NULL);
 if(ret == 1) {
   file1_column2 = 0;
   file1_column3 = 0;
 }else if(ret == 2) {
   file1_column3 = 0;
 }else if(ret != 3) {
   printerror(application.verbose_state.debug, "ERROR pstat: Cannot parse '%s' option, 1, 2 or 3 values.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-col2") == 0) {
 int ret;
 ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%d %d %d", &file2_column1, &file2_column2, &file2_column3, NULL);
 if(ret == 1) {
   file2_column2 = 0;
   file2_column3 = 0;
 }else if(ret == 2) {
   file2_column3 = 0;
 }else if(ret != 3) {
   printerror(application.verbose_state.debug, "ERROR pstat: Cannot parse '%s' option, 1, 2 or 3 values.", argv[i]);
   return 0;
 }
 i++;
      }else if(strcmp(argv[i], "-output") == 0) {
 output_idx = i+1;
 i++;
      }else if(strcmp(argv[i], "-correl") == 0) {
 if(typetest != 0) {
   printerror(application.verbose_state.debug, "pstat: Cannot specify more than one type of statistical test at the time");
   return 0;
 }
 typetest = CORREL;
      }else if(strcasecmp(argv[i], "-pearson") == 0) {
 if(typetest != 0) {
   printerror(application.verbose_state.debug, "pstat: Cannot specify more than one type of statistical test at the time");
   return 0;
 }
 typetest = PEARSON;
      }else if(strcmp(argv[i], "-ks") == 0) {
 if(typetest != 0) {
   printerror(application.verbose_state.debug, "pstat: Cannot specify more than one type of statistical test at the time");
   return 0;
 }
 typetest = KSTEST;
      }else if(strcmp(argv[i], "-ksflat") == 0) {
 if(typetest != 0) {
   printerror(application.verbose_state.debug, "pstat: Cannot specify more than one type of statistical test at the time");
   return 0;
 }
 typetest = KSFLAT;
      }else if(strcmp(argv[i], "-moments") == 0) {
 if(typetest != 0) {
   printerror(application.verbose_state.debug, "pstat: Cannot specify more than one type of statistical test at the time");
   return 0;
 }
 typetest = MOMENTS;
      }else if(strcmp(argv[i], "-median") == 0) {
 if(typetest != 0) {
   printerror(application.verbose_state.debug, "pstat: Cannot specify more than one type of statistical test at the time");
   return 0;
 }
 typetest = MEDIAN;
      }else if(strcmp(argv[i], "-kssin") == 0) {
 if(typetest != 0) {
   printerror(application.verbose_state.debug, "pstat: Cannot specify more than one type of statistical test at the time");
   return 0;
 }
 typetest = KSSIN;
      }else if(strcmp(argv[i], "-chi2hist") == 0) {
 int ret;
 ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%lf %lf %lf", &threshold1, &threshold2, &threshold3, NULL);
 if(ret == 1) {
   threshold2 = threshold3 = -1;
 }else if(ret == 2) {
   threshold3 = -1;
 }else if(ret != 3) {
   printerror(application.verbose_state.debug, "ERROR pstat: Cannot parse '%s' option, 1, 2 or 3 values.", argv[i]);
   return 0;
 }
 if(typetest != 0) {
   printerror(application.verbose_state.debug, "pstat: Cannot specify more than one type of statistical test at the time");
   return 0;
 }
 typetest = CHI2TEST_HIST;
 i++;
      }else if(strcmp(argv[i], "-chi2cdf") == 0) {
 if(typetest != 0) {
   printerror(application.verbose_state.debug, "pstat: Cannot specify more than one type of statistical test at the time");
   return 0;
 }
 typetest = CHI2TEST_CDF;
      }else if(strcmp(argv[i], "-log") == 0) {
 read_log = 1;
      }else {
 if(argv[i][0] == '-') {
   printerror(application.verbose_state.debug, "pstat: Unknown option: %s\n\nRun pstat without command line arguments to show help", argv[i]);
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
  {
    int nrInputFiles;
    int nrInputColumns;
    nrInputFiles = numberInApplicationFilenameList(&application, argv, application.verbose_state);
    if(nrInputFiles < 1) {
      printerror(application.verbose_state.debug, "ERROR pstat: No files specified");
      return 0;
    }
    if(nrInputFiles > 2) {
      printerror(application.verbose_state.debug, "ERROR pstat: Cannot specify more than two input files. %d input files are currently specified.", nrInputFiles);
      if(application.verbose_state.debug) {
 printerror(0, "  These are:");
 char *filename_ptr;
 for(i = 0; i < nrInputFiles; i++) {
   filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state);
   if(filename_ptr != NULL) {
     printerror(0, "    %s", filename_ptr);
   }
 }
      }
      return 0;
    }
    if(file1_column1 == 0) {
      {
 printwarning(application.verbose_state.debug, "WARNING pstat: No -col1 option. Using the default (-col1 1).");
 file1_column1 = 1;
 file1_column2 = 0;
 file1_column3 = 0;
      }
    }
    if(nrInputFiles == 2 && file2_column1 == 0) {
      {
 printwarning(application.verbose_state.debug, "WARNING pstat: No -col2 option, while there are two input files. Default is -col2 1.");
 file2_column1 = 1;
 file2_column2 = 0;
 file2_column3 = 0;
      }
    }
    if(file2_column1 && nrInputFiles != 2) {
      printerror(application.verbose_state.debug, "ERROR pstat: -col2 option suggest that there should be two input files, but there is only %d defined on command line.", nrInputFiles);
      return 0;
    }
    nrInputColumns = 0;
    if(file1_column1)
      nrInputColumns++;
    if(file1_column2)
      nrInputColumns++;
    if(file1_column3)
      nrInputColumns++;
    if(file2_column1)
      nrInputColumns++;
    if(file2_column2)
      nrInputColumns++;
    if(file2_column3)
      nrInputColumns++;
    if(typetest == KSTEST) {
      if(nrInputColumns != 2) {
 printerror(application.verbose_state.debug, "ERROR pstat: KS-test requires two columns of data to be read in. Example: pstat -ks -col1 1 -col2 1 file1 file2 or pstat -ks -col \"1 2\" file1");
 return 0;
      }
    }else if(typetest == CORREL) {
      if(nrInputColumns != 2) {
 printerror(application.verbose_state.debug, "ERROR pstat: Correlation requires two columns of data to be read in. Example: pstat -correl -col1 1 -col2 1 file1 file2 or pstat -correl -col \"1 2\" file1");
 return 0;
      }
    }else if(typetest == PEARSON) {
      if(nrInputColumns != 2) {
 printerror(application.verbose_state.debug, "ERROR pstat: Correlation requires two columns of data to be read in. Example: pstat -pearson -col1 1 -col2 1 file1 file2 or pstat -pearson -col \"1 2\" file1");
 return 0;
      }
    }else if(typetest == MOMENTS) {
      if(nrInputColumns != 1) {
 printerror(application.verbose_state.debug, "ERROR pstat: Computation of the moments of a distribution requires one columns of data to be read in. Example: pstat -moments -col 1 file1");
 return 0;
      }
    }else if(typetest == MEDIAN) {
      if(nrInputColumns != 1) {
 printerror(application.verbose_state.debug, "ERROR pstat: Computation of the median of a distribution requires one columns of data to be read in. Example: pstat -median -col 1 file1");
 return 0;
      }
    }else if(typetest == KSFLAT) {
      if(nrInputColumns != 1) {
 printerror(application.verbose_state.debug, "ERROR pstat: KS-test comparison with a flat distribution requires one columns of data to be read in. Example: pstat -ksflat -col 1 file1");
 return 0;
      }
    }else if(typetest == KSSIN) {
      if(nrInputColumns != 1) {
 printerror(application.verbose_state.debug, "ERROR pstat: KS-test comparison with a sinusoidal distribution requires one columns of data to be read in. Example: pstat -kssin -col 1 file1");
 return 0;
      }
    }else if(typetest == CHI2TEST_HIST) {
      if(nrInputColumns != 4 && nrInputColumns != 6) {
 printerror(application.verbose_state.debug, "ERROR pstat: The chi-square histogram test requires four or six columns of data to be read in. Examples: pstat -chi2hist -1 -col1 \"1 2\" -col2 \"1 2\" file1 file2 or pstat -chi2hist -1 -col1 \"1 2 3\" -col2 \"1 2 3\" file1 file2");
 return 0;
      }
    }else if(typetest == CHI2TEST_CDF) {
      if(nrInputColumns != 2) {
 printerror(application.verbose_state.debug, "ERROR pstat: The chi-square CDF test requires two columns of data to be read in. Examples: pstat -chi2cdf -col1 1 -col2 1 file1 file2 or pstat -chi2cdf -col1 \"1 2\" file1");
 return 0;
      }
    }else {
      printerror(application.verbose_state.debug, "ERROR pstat: No statistical test has been specified, nothing to do.");
      return 0;
    }
  }
  char *filename_ptr;
  double *input_array[6];
  long number_values[6];
  int number_input_arrays;
  number_input_arrays = 0;
  number_values[0] = 0;
  filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state);
  if(filename_ptr == NULL) {
    printerror(application.verbose_state.debug, "ERROR pstat: Bug!");
    return 0;
  }
  if(file1_column1) {
    double min_x, max_x, avrg;
    if(read_ascii_column_double(filename_ptr, 0, '#', -1, 1, &number_values[number_input_arrays], file1_column1, 1.0, read_log, &input_array[number_input_arrays], &min_x, &max_x, &avrg, application.verbose_state, 0) == 0) {
      printerror(application.verbose_state.debug, "ERROR pstat: cannot load file.\n");
      return 0;
    }
    number_input_arrays++;
  }
  if(file1_column2) {
    double min_x, max_x, avrg;
    if(read_ascii_column_double(filename_ptr, 0, '#', -1, 1, &number_values[number_input_arrays], file1_column2, 1.0, read_log, &input_array[number_input_arrays], &min_x, &max_x, &avrg, application.verbose_state, 0) == 0) {
      printerror(application.verbose_state.debug, "ERROR pstat: cannot load file.\n");
      return 0;
    }
    number_input_arrays++;
  }
  if(file1_column3) {
    double min_x, max_x, avrg;
    if(read_ascii_column_double(filename_ptr, 0, '#', -1, 1, &number_values[number_input_arrays], file1_column3, 1.0, read_log, &input_array[number_input_arrays], &min_x, &max_x, &avrg, application.verbose_state, 0) == 0) {
      printerror(application.verbose_state.debug, "ERROR pstat: cannot load file.\n");
      return 0;
    }
    number_input_arrays++;
  }
  if(file2_column1 || file2_column2 || file2_column3) {
    filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state);
    if(filename_ptr == NULL) {
      printerror(application.verbose_state.debug, "ERROR pstat: Bug!");
      return 0;
    }
  }
  if(file2_column1) {
    double min_x, max_x, avrg;
    if(read_ascii_column_double(filename_ptr, 0, '#', -1, 1, &number_values[number_input_arrays], file2_column1, 1.0, read_log, &input_array[number_input_arrays], &min_x, &max_x, &avrg, application.verbose_state, 0) == 0) {
      printerror(application.verbose_state.debug, "ERROR pstat: cannot load file.\n");
      return 0;
    }
    number_input_arrays++;
  }
  if(file2_column2) {
    double min_x, max_x, avrg;
    if(read_ascii_column_double(filename_ptr, 0, '#', -1, 1, &number_values[number_input_arrays], file2_column2, 1.0, read_log, &input_array[number_input_arrays], &min_x, &max_x, &avrg, application.verbose_state, 0) == 0) {
      printerror(application.verbose_state.debug, "ERROR pstat: cannot load file.\n");
      return 0;
    }
    number_input_arrays++;
  }
  if(file2_column3) {
    double min_x, max_x, avrg;
    if(read_ascii_column_double(filename_ptr, 0, '#', -1, 1, &number_values[number_input_arrays], file2_column3, 1.0, read_log, &input_array[number_input_arrays], &min_x, &max_x, &avrg, application.verbose_state, 0) == 0) {
      printerror(application.verbose_state.debug, "ERROR pstat: cannot load file.\n");
      return 0;
    }
    number_input_arrays++;
  }
  FILE *fout;
  if(output_idx) {
    fout = fopen(argv[output_idx], "w");
    if(fout == NULL) {
      printerror(application.verbose_state.debug, "ERROR pstat: Cannot open %s", argv[output_idx]);
      return 0;
    }
  }else {
    fout = stdout;
  }
  if(typetest == KSTEST) {
    double max_diff, prob;
    if(number_input_arrays != 2) {
      printerror(application.verbose_state.debug, "ERROR pstat: KS-test requires two input arrays to be specified");
      return 0;
    }
    kstest(input_array[0], number_values[0], input_array[1], number_values[1], 0, NULL, &max_diff, &prob, application.verbose_state);
    if(application.verbose_state.verbose == 0) {
      fprintf(fout, "probability = %e\n", prob);
    }
  }else if(typetest == KSFLAT) {
    double max_diff, prob;
    if(number_input_arrays != 1) {
      printerror(application.verbose_state.debug, "ERROR pstat: KS-test comparison with a flat distribution requires one input array to be specified");
      return 0;
    }
    kstest(input_array[0], number_values[0], NULL, 0, 1, NULL, &max_diff, &prob, application.verbose_state);
    if(application.verbose_state.verbose == 0) {
      fprintf(fout, "probability = %e\n", prob);
    }
  }else if(typetest == KSSIN) {
    double max_diff, prob;
    if(number_input_arrays != 1) {
      printerror(application.verbose_state.debug, "ERROR pstat: KS-test comparison with a sinusoidal distribution requires one input array to be specified");
      return 0;
    }
    kstest(input_array[0], number_values[0], NULL, 0, 2, NULL, &max_diff, &prob, application.verbose_state);
    if(application.verbose_state.verbose == 0) {
      fprintf(fout, "probability = %e\n", prob);
    }
  }else if(typetest == MOMENTS) {
    double mean, variance, skew, kurt;
    if(number_input_arrays != 1) {
      printerror(application.verbose_state.debug, "ERROR pstat: Computation of the moments of a distribution requires one columns of data to be read in.");
      return 0;
    }
    mean = gsl_stats_mean(input_array[0], 1, number_values[0]);
    fprintf(fout, "Mean               = %e\n", mean);
    variance = gsl_stats_variance_m(input_array[0], 1, number_values[0], mean);
    fprintf(fout, "Variance           = %e (unbiased, normalised by 1/(N-1))\n", variance);
    fprintf(fout, "Standard deviation = %e (unbiased, normalised by 1/(N-1))\n", sqrt(variance));
    variance *= (number_values[0]-1.0)/(double)number_values[0];
    fprintf(fout, "Variance           = %e (normalised by 1/N)\n", variance);
    fprintf(fout, "Standard deviation = %e (normalised by 1/N)\n", sqrt(variance));
    skew = gsl_stats_skew(input_array[0], 1, number_values[0]);
    fprintf(fout, "Skewness           = %e\n", skew);
    kurt = gsl_stats_kurtosis(input_array[0], 1, number_values[0]);
    fprintf(fout, "Kurtosis           = %e\n", kurt);
  }else if(typetest == MEDIAN) {
    double median;
    if(number_input_arrays != 1) {
      printerror(application.verbose_state.debug, "ERROR pstat: Computation of the moments of a distribution requires one columns of data to be read in.");
      return 0;
    }
    gsl_sort(input_array[0], 1, number_values[0]);
    median = gsl_stats_median_from_sorted_data(input_array[0], 1, number_values[0]);
    fprintf(fout, "Median = %e\n", median);
  }else if(typetest == PEARSON) {
    double cc;
    if(number_input_arrays != 2) {
      printerror(application.verbose_state.debug, "ERROR pstat: Correlation requires two input arrays to be specified");
      return 0;
    }
    if(number_values[0] != number_values[1]) {
      printerror(application.verbose_state.debug, "ERROR pstat: Correlation requires two input arrays of equal length");
      return 0;
    }
#if GSL_VERSION_NUMBER >= 110
    cc = gsl_stats_correlation(input_array[0], 1, input_array[1], 1, number_values[0]);
    fprintf(fout, "Pearson correlation coefficient = %e\n", cc);
#else
    printerror(application.verbose_state.debug, "ERROR pstat: Pearson correlation coefficient cannot be calculated if GSL < 1.10");
    return 0;
#endif
  }else if(typetest == CORREL) {
    if(number_input_arrays != 2) {
      printerror(application.verbose_state.debug, "ERROR pstat: Correlation requires two input arrays to be specified");
      return 0;
    }
    if(number_values[0] != number_values[1]) {
      printerror(application.verbose_state.debug, "ERROR pstat: Correlation requires two input arrays of equal length");
      return 0;
    }
    if(number_values[0] > 2147483646) {
      printerror(application.verbose_state.debug, "ERROR pstat: Input array too long.");
      return 0;
    }
    float *x1, *y1;
    x1 = malloc(number_values[0]*sizeof(float));
    y1 = malloc(number_values[1]*sizeof(float));
    if(x1 == NULL || y1 == NULL) {
      printerror(application.verbose_state.debug, "ERROR pstat: Memory allocation error");
      return 0;
    }
    for(i = 0; i < number_values[0]; i++)
      x1[i] = (input_array[0])[i];
    for(i = 0; i < number_values[1]; i++)
      y1[i] = (input_array[1])[i];
    int extrazeropad, cc_length;
    float *ans;
    extrazeropad = 0;
    if(crosscorrelation_fft_padding(x1, y1, number_values[0], extrazeropad, &ans, &cc_length, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR pstat: Cross correlation failed.");
      return 0;
    }
    free(x1);
    free(y1);
    int lagoutput;
    lagoutput = 1;
    if(lagoutput == 0) {
      for(i = 0; i < cc_length; i++)
 fprintf(fout, "%ld %e\n", i, ans[i]);
    }else {
      for(i = cc_length/2; i <= cc_length-1; i++)
 fprintf(fout, "%ld %e\n", i-cc_length, ans[i]);
      for(i = 0; i < cc_length/2; i++)
 fprintf(fout, "%ld %e\n", i, ans[i]);
    }
    free(ans);
  }else if(typetest == CHI2TEST_HIST) {
    double binwidth, ratio, offset, chi2;
    long offset_binnr, file2_x_col, dof, i2, nr_overlapping_bins;
    if(number_input_arrays != 4 && number_input_arrays != 6) {
      printerror(application.verbose_state.debug, "ERROR pstat: The chi-square histogram test requires four or six columns of data to be specified.");
      return 0;
    }
    if(number_input_arrays == 4) {
      file2_x_col = 2;
      printwarning(application.verbose_state.debug, "WARNING pstat: Since no column numbers with error-bars are provided, uniform weighting of the different bins is assumed with sigma=1. This is unlikely to be correct.");
    }else {
      file2_x_col = 3;
    }
    binwidth = (input_array[0])[1] - (input_array[0])[0];
    ratio = binwidth/((input_array[file2_x_col])[1] - (input_array[file2_x_col])[0]);
    if(application.verbose_state.verbose) {
      printf("Ratio bin widths of two histograms is %lf (should be very close to 1)\n", binwidth);
    }
    if(ratio < 0.999 || ratio > 1.001) {
      printerror(application.verbose_state.debug, "ERROR pstat: The binwidths of the two histograms appear to be different (%e != %e).", (input_array[0])[1] - (input_array[0])[0], (input_array[file2_x_col])[1] - (input_array[file2_x_col])[0]);
      return 0;
    }
    offset = (input_array[file2_x_col])[0];
    offset -= (input_array[0])[0];
    offset /= binwidth;
    if(application.verbose_state.verbose) {
      printf("Offset between two histograms is %lf bins (should be very close to an integer value)\n", offset);
    }
    offset_binnr = round(offset);
    offset = fabs(offset - offset_binnr);
    if(offset > 0.001) {
      printerror(application.verbose_state.debug, "ERROR pstat: The bins of the two histograms do not appear to be aligned, but have an offset of %lf.", offset);
      return 0;
    }
    double height1, height2;
    height1 = 0;
    height2 = 0;
    for(i = 0; i < number_values[0]; i++)
      height1 += (input_array[1])[i];
    for(i = 0; i < number_values[file2_x_col]; i++)
      height2 += (input_array[file2_x_col+1])[i];
    if(application.verbose_state.verbose) {
      printf("Ratio of integrals of two histograms is %lf (should be very close to 1)\n", height1/height2);
    }
    if(height1/height2 > 1.001 || height2/height1 > 1.001) {
      printwarning(application.verbose_state.debug, "WARNING pstat: The two histograms appear to be normalised differently, so the derived numbers are unlikely to give useful results.");
    }
    chi2 = 0;
    dof = 0;
    nr_overlapping_bins = 0;
    for(i = -labs(offset_binnr)-10; i < number_values[0]+labs(offset_binnr)+10; i++) {
      int hist1_exist, hist2_exist;
      hist1_exist = hist2_exist = 0;
      if(i >= 0 && i < number_values[0]) {
 height1 = (input_array[1])[i];
 hist1_exist = 1;
      }else {
 height1 = 0;
      }
      i2 = i - offset_binnr;
      if(i2 >= 0 && i2 < number_values[file2_x_col]) {
 height2 = (input_array[file2_x_col+1])[i2];
 hist2_exist = 1;
 if(hist1_exist) {
   offset = (input_array[file2_x_col])[i2];
   offset -= (input_array[0])[i];
   offset /= binwidth;
   if(fabs(offset) > 0.001) {
     printerror(application.verbose_state.debug, "ERROR pstat: The bins of the two histograms do not appear to be aligned, but bins %ld and %ld have an offset of %lf.", i+1, i2+1, offset);
     return 0;
   }
   nr_overlapping_bins++;
 }
      }else {
 height2 = 0;
      }
      if(hist1_exist || hist2_exist) {
 if(height1 >= threshold1 && height2 >= threshold2 && (height1+height2) >= threshold3) {
   double delta_y, var;
   delta_y = height2 - height1;
   if(number_input_arrays == 4) {
     chi2 += delta_y*delta_y;
     dof++;
   }else {
     var = 0;
     if(hist1_exist)
       var += (input_array[2])[i] * (input_array[2])[i];
     if(hist2_exist)
       var += (input_array[file2_x_col+2])[i2] * (input_array[file2_x_col+2])[i2];
     chi2 += delta_y*delta_y/var;
     dof++;
   }
 }
      }
    }
    if(nr_overlapping_bins == 0) {
      printerror(application.verbose_state.debug, "ERROR pstat: There appears to be no overlap between the two histograms.");
      return 0;
    }
    if(application.verbose_state.verbose) {
      printf("Number of overlapping bins between two distributions:     %ld\n", nr_overlapping_bins);
      printf("Total number of bins considered with specified threshold: %ld\n", dof);
    }
    if(number_input_arrays != 6)
      fprintf(fout, "Total non-weighted chi square:   %f = %e\n", chi2, chi2);
    else if(application.verbose_state.verbose)
      fprintf(fout, "Total chi square:   %f = %e\n", chi2, chi2);
    if(number_input_arrays == 6) {
      fprintf(fout, "Reduced chi square: %f = %e\n", chi2/(double)dof, chi2/(double)dof);
    }
  }else if(typetest == CHI2TEST_CDF) {
    double x1, cdf1, cdf2, chi2;
    long start2, dof;
    if(number_input_arrays != 2) {
      printerror(application.verbose_state.debug, "ERROR pstat: The chi-square CDF test requires two columns of data to be specified.");
      return 0;
    }
    if(number_values[0] > number_values[1]) {
      i = number_values[0];
      number_values[0] = number_values[1];
      number_values[1] = i;
      input_array[2] = input_array[0];
      input_array[0] = input_array[1];
      input_array[1] = input_array[2];
    }
    gsl_sort(input_array[0], 1, number_values[0]);
    gsl_sort(input_array[1], 1, number_values[1]);
    start2 = 0;
    chi2 = 0;
    dof = 0;
    for(i = 0; i < number_values[0] - 1; i++) {
      cdf1 = (i+1)/(double)number_values[0];
      x1 = (input_array[0])[i] + 0.5*((input_array[0])[i+1]-(input_array[0])[i]);
      if(application.verbose_state.debug)
 printf("Going to find chi2 of observation distribution cdf point: (%e, %e)\n", x1, cdf1);
      if(x1 <= (input_array[1])[0]) {
 if(application.verbose_state.debug)
   printf("  model cdf starts after point of interest\n");
 cdf2 = 0;
 if(x1 == (input_array[1])[0]) {
   cdf2 = (j+0.5)/(double)number_values[1];
 }
      }else if(x1 >= (input_array[1])[number_values[1]-1]) {
 if(application.verbose_state.debug)
   printf("  model cdf ends before point of interest\n");
 cdf2 = 1;
 if(x1 == (input_array[1])[number_values[1]-1]) {
   cdf2 = (number_values[1]-1 +0.5)/(double)number_values[1];
 }
      }else {
 for(j = start2; j < number_values[1]; j++) {
   if((input_array[1])[j] >= x1) {
     if(application.verbose_state.debug)
       printf("  First model cdf point after point of interest: (%e %e)\n", (input_array[1])[j], (j+1)/(double)number_values[1]);
     long index1, index2;
     index1 = j-1;
     index2 = j;
     if(j == 0) {
       printerror(application.verbose_state.debug, "ERROR pstat: Bug!\n");
       return 0;
     }
     while((input_array[1])[index2] == (input_array[1])[index1]) {
       if(index2 < number_values[1] - 1) {
  index2++;
       }else if(index1 > 1) {
  index1--;
       }else {
  printerror(application.verbose_state.debug, "ERROR pstat: Something is wrong with the second input distribution, as all input values appear to be identical.\n");
  return 0;
       }
     }
     if(application.verbose_state.debug)
       printf("  Going to interpolate following model points: (%e %e) and (%e %e)\n", (input_array[1])[index1], (index1+1-0.5)/(double)number_values[1], (input_array[1])[index2], (index2+1-0.5)/(double)number_values[1]);
     if(index1 < 0 || index2 >= number_values[1]) {
       printerror(application.verbose_state.debug, "ERROR pstat: Interpolation failed - trying to access second array out of limits (index1=%ld, index2=%ld).\n", index1, index2);
       return 0;
     }
     cdf2 = ((index2+1-0.5)/(double)number_values[1] - (index1+1-0.5)/(double)number_values[1]) * (x1 - (input_array[1])[index1]) / ((input_array[1])[index2] - (input_array[1])[index1]) + (index1+1-0.5)/(double)number_values[1];
     if(!isfinite(cdf2)) {
       printerror(application.verbose_state.debug, "ERROR pstat: Interpolation failed: i=%ld, j=%ld (array2=[%e, %e, %e]), \n", i, j, (input_array[1])[j-1], (input_array[1])[j], (input_array[1])[j+1]);
       return 0;
     }
     start2 = j - 2;
     if(start2 < 0)
       start2 = 0;
     break;
   }
 }
      }
      chi2 += (cdf2-cdf1)*(cdf2-cdf1);
      dof++;
      if(application.verbose_state.debug) {
 printf("  cdf2=%e, cdf1=%e, diff=%e\n", cdf2, cdf1, (cdf2-cdf1));
 printf("  new chi2 = %e\n", chi2);
      }
    }
    if(application.verbose_state.verbose) {
      printf("\nTotal number of bins considered: %ld\n", dof);
    }
    fprintf(fout, "Non-weighted total chi square = %f = %e\n", chi2, chi2);
  }else {
    printerror(application.verbose_state.debug, "ERROR pstat: Bug!");
    return 0;
  }
  if(output_idx) {
    fclose(fout);
  }
  for(i = 0; i < number_input_arrays; i++) {
    free(input_array[i]);
  }
  terminateApplication(&application);
  return 0;
}

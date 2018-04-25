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
#include <psrsalsa.h>
int main(int argc, char **argv)
{
  long filesize, i;
  unsigned char byte1, byte2, result;
  FILE *fin1, *fin2, *fout;
  if(argc != 4) {
    printf("This program reads in two files which should be of the same size. A byte at the time is read in from both input files (assumed to be unsigned) and the average of the two is written to the output file. This can be used in combination with ImageMagick to obtain transparent figures. See the documentation of ppolFit.\n\n");
    printf("Usage: avrg_bin_files file1 file2 output_file\n\n");
    printCitationInfo();
    return 0;
  }
  fin1 = fopen(argv[1], "r");
  if(fin1 == NULL) {
    printf("Cannot open first file=%s\n", argv[1]);
    return 0;
  }
  fin2 = fopen(argv[2], "r");
  if(fin2 == NULL) {
    printf("Cannot open second file=%s\n", argv[2]);
    return 0;
  }
  fout = fopen(argv[3], "w");
  if(fout == NULL) {
    printf("Cannot open output file=%s\n", argv[3]);
    return 0;
  }
  fseek(fin1, 0, SEEK_END);
  filesize = ftell(fin1);
  fseek(fin2, 0, SEEK_END);
  if(filesize != ftell(fin2)) {
    printf("Input files have different sizes!\n");
    return 0;
  }
  rewind(fin1);
  rewind(fin2);
  for(i = 0; i < filesize; i++) {
    byte1 = fgetc(fin1);
    byte2 = fgetc(fin2);
    result = (byte1+byte2)/2;
    fputc(result, fout);
  }
  fclose(fout);
  fclose(fin1);
  fclose(fin2);
  return 0;
}

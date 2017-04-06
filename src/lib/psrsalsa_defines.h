/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/



#define DM_CONST 4.148808e3





#define maxNrRegions 200

#define maxNrVonMisesComponents 100

#define MaxPgplotDeviceLength 2000
#define PUMA_format 1
#define PSRCHIVE_ASCII_format 5
#define EPN_format 6
#define FITS_format 7
#define SIGPROC_format 8
#define PPOL_format 9
#define PPOL_SHORT_format 10
#define SIGPROC_ASCII_format 11
#define PSRSALSA_BINARY_format 20
#define MEMORY_format 99

#define PPGPLOT_GRAYSCALE 1
#define PPGPLOT_INVERTED_GRAYSCALE 2
#define PPGPLOT_RED 3
#define PPGPLOT_INVERTED_RED 4
#define PPGPLOT_GREEN 5
#define PPGPLOT_INVERTED_GREEN 6
#define PPGPLOT_BLUE 7
#define PPGPLOT_INVERTED_BLUE 8
#define PPGPLOT_CYAN 9
#define PPGPLOT_INVERTED_CYAN 10
#define PPGPLOT_HEAT 20
#define PPGPLOT_INVERTED_HEAT 21
#define PPGPLOT_COLD 22
#define PPGPLOT_INVERTED_COLD 23
#define PPGPLOT_HEAT2 30
#define PPGPLOT_INVERTED_HEAT2 31
#define PPGPLOT_HEAT3 32
#define PPGPLOT_INVERTED_HEAT3 33
#define PPGPLOT_HEAT4 34
#define PPGPLOT_INVERTED_HEAT4 35
#define PPGPLOT_DIVREDBLUE 36
#define PPGPLOT_INVERTED_DIVREDBLUE 37


#define GENTYPE_UNDEFINED 0
#define GENTYPE_PROFILE 1
#define GENTYPE_PULSESTACK 2
#define GENTYPE_SUBINTEGRATIONS 3

#define GENTYPE_SEARCHMODE 4
#define GENTYPE_BANDPASS 5
#define GENTYPE_DYNAMICSPECTRUM 6
#define GENTYPE_PENERGY 7
#define GENTYPE_POLNCAL 10
#define GENTYPE_LRFS 20
#define GENTYPE_2DFS 21
#define GENTYPE_S2DFSP3 22
#define GENTYPE_S2DFSP2 23
#define GENTYPE_P3FOLD 25
#define GENTYPE_HRFS_UNFOLDED 30
#define GENTYPE_HRFS 31
#define GENTYPE_LRCC 33
#define GENTYPE_LRAC 34
#define GENTYPE_RMMAP 50

#define GENTYPE_PADIST 101
#define GENTYPE_RECEIVERMODEL 200
#define GENTYPE_RECEIVERMODEL2 201

#define FOLDMODE_UNKNOWN -1
#define FOLDMODE_FIXEDPERIOD 1

#define TSAMPMODE_UNKNOWN -1
#define TSAMPMODE_FIXEDTSAMP 1
#define TSAMPMODE_LONGITUDELIST 2

#define TSUBMODE_UNKNOWN -1
#define TSUBMODE_FIXEDTSUB 1
#define TSUBMODE_TSUBLIST 2


#define POLTYPE_UNKNOWN -1
#define POLTYPE_STOKES 1
#define POLTYPE_COHERENCY 2
#define POLTYPE_ILVPAdPA 3
#define POLTYPE_PAdPA 4

#define FREQMODE_UNKNOWN -1
#define FREQMODE_UNIFORM 1


#define FEEDTYPE_UNKNOWN 0
#define FEEDTYPE_LINEAR 1
#define FEEDTYPE_CIRCULAR 2
#define FEEDTYPE_INV_LINEAR -1
#define FEEDTYPE_INV_CIRCULAR -2



#define MaxOutputNameLength 10000

#define MaxNrApplicationFilenames 1025

#define maxNrRotateStokes 10



#ifndef NAN
  #define NAN (0.0/0.0)
#endif

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif



#define printerror(debug_flag,...) \
  { fflush(stdout); fprintf_color(stderr, 2, __VA_ARGS__); \
  if(debug_flag) fprintf_color(stderr, 2, " (message generated in %s line %d)", __FILE__, __LINE__); \
  fprintf_color(stderr, 0, "\n"); }


#define printwarning(debug_flag,...) \
  fflush(stdout); fprintf_color(stderr, 7, __VA_ARGS__); \
  if(debug_flag) fprintf_color(stderr, 7, " (message generated in %s line %d)", __FILE__, __LINE__); \
  fprintf_color(stderr, 0, "\n")

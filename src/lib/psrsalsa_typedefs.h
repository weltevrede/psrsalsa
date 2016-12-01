/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


typedef struct {
  int nrRegions;

  int bins_defined[maxNrRegions];
  int left_bin[maxNrRegions], right_bin[maxNrRegions];

  int frac_defined[maxNrRegions];
  float left_frac[maxNrRegions], right_frac[maxNrRegions];
}regions_definition;

typedef struct {
  int verbose;
  int debug;
  int nocounters;
  int indent;
}verbose_definition;


typedef struct {

  double centre[maxNrVonMisesComponents], concentration[maxNrVonMisesComponents], height[maxNrVonMisesComponents];
  int nrcomponents;
}vonMises_def;


typedef struct {
  char plotDevice[MaxPgplotDeviceLength];
  int windowwidth, windowheight;
  float aspectratio;
  float dxplot;
  float xsize;
  float dyplot;
  float ysize;
  int noclear;
  int dontopen;
  int dontclose;
}pgplot_viewport_def;

typedef struct {
  int svp;
  float svp_x1, svp_x2, svp_y1, svp_y2;
  int swin;
  int swin_showtwice;
  float swin_x1, swin_x2, swin_y1, swin_y2;
  float TR[6];
}pgplot_frame_def_internal;

typedef struct {
  int drawbox;
  int box_lw;
  float box_labelsize;
  float box_xtick, box_ytick;
  int box_nxsub, box_nysub;
  char box_xopt[10], box_yopt[10];

  int drawtitle;
  float title_ch;
  int title_lw, title_f;
  char title[1000];

  int drawlabels;
  float label_ch;
  float dxlabel, dylabel;
  char xlabel[1000];
  char ylabel[1000];
  char wedgelabel[1000];
}pgplot_box_def;



typedef struct {
  char *timestamp;
  char *cmd;
  char *user;
  char *hostname;
  void *nextEntry;
}datafile_history_entry_definition;

typedef struct
{



  FILE *fptr, *fptr_hdr;
  fitsfile *fits_fptr;
  char *filename;
  int format;
  int version;
  int opened_flag, enable_write_flag;
  int dumpOnClose;






  char *psrname;
  char *observatory;
  char *instrument;
  char *scanID;

  char *institute;


  double telescope_X, telescope_Y, telescope_Z;

  int NrBits;
  char isDeDisp, isDeFarad, isDePar, isDebase;
  double dm, rm;
  double freq_ref;
  int feedtype;
  int poltype;
  long double mjd_start;
  char cableSwap;
  char cableSwapcor;


  long NrSubints, NrBins, NrPols, NrFreqChan;


  char isFolded;
  char foldMode;
  double fixedPeriod;
  char tsampMode;
  double fixedtsamp;
  double *tsamp_list;

  char tsubMode;
  double *tsub_list;


  double ra, dec;



  char freqMode;
  double uniform_freq_cent, uniform_bw;



  int gentype;
  char isTransposed;

  float xrange[2];
  float yrange[2];
  char xrangeset, yrangeset;
  datafile_history_entry_definition history;





  float *data;
  float *offpulse_rms;







  float *scales, *offsets, *weights;
  long long datastart;
}datafile_definition;
typedef struct {
  char progname[100], *genusage;
  int switch_verbose, switch_debug, switch_nocounters;
  verbose_definition verbose_state;
  int switch_formatlist;
  int switch_iformat, iformat;
  int switch_oformat, oformat;
  int switch_header;
  int switch_headerlist;
  int switch_onpulse, switch_onpulsef;
  int switch_polselect, polselectnr;
  int switch_itf, itf;
  int switch_rebin, dorebin; long rebin;
  int switch_nread; long nread;
  int switch_nskip; long nskip;
  int switch_conshift, doconshift;
  int switch_circshift, docircshift;
  int switch_rot, switch_rotdeg, doshiftphase; float shiftPhase_cmdline, shiftPhase;
  int switch_filelist, filelist;
  int switch_device; char pgplotdevice[MaxOutputNameLength];
  int switch_tscr; long dotscr;
  int switch_tscr_complete; int tscr_complete;
  int switch_TSCR, doTSCR;
  int switch_fscr; long dofscr;
  int switch_FSCR, doFSCR;
  int switch_dedisperse, do_dedisperse;
  int switch_deFaraday, do_deFaraday;
  int switch_changeRefFreq; double newRefFreq;
  int switch_stokes, dostokes;
  int switch_coherence, docoherence;
  int switch_noweights, noweights;
  int switch_scale, doscale; float scale_scale, scale_offset;
  int switch_debase, dodebase;
  int switch_onpulsegr, doonpulsegr;
  int switch_size, windowwidth, windowheight;
  int switch_macro; FILE *macro_ptr;
  int switch_noplotsubset, do_noplotsubset;
  int switch_cmap, cmap;
  int switch_cmaplist;
  int switch_insertparang, switch_deparang, do_parang_corr;
  int switch_history_cmd_only, history_cmd_only;
  int switch_norm, do_norm; float normvalue;
  int switch_normglobal, do_normglobal;
  int switch_clip, do_clip; float clipvalue;
  int switch_rotateQU, dorotateQU; float rotateQUangle;
  int switch_rotateUV, dorotateUV; float rotateUVangle;
  int switch_fchan, fchan_select;
  int switch_fixseed, fixseed;
  int switch_templatedata, template_data_index; datafile_definition template_file;
  int switch_template, template_specified;
  int switch_align, doalign;
  int switch_blocksize, blocksize;
  regions_definition onpulse;
  vonMises_def vonMises_components;
  int switch_ext; char *extension;
  int switch_output; char outputname[MaxOutputNameLength];
  int switch_shuffle, doshuffle;
  int doautot;
  int *fzapMask;
}psrsalsaApplication;

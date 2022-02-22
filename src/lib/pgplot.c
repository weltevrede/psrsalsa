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
#include <string.h>
#include <math.h>
#include "psrsalsa.h"
#include <cpgplot.h>
static double internal_pgplot_xmin = 0;
static double internal_pgplot_xmax = 1;
static double internal_pgplot_ymin = 0;
static double internal_pgplot_ymax = 1;
static int internal_pgplot_nrx = 1;
static int internal_pgplot_nry = 1;
void ppgplot_set_internal_mapping_coordinated(double xmin, double xmax, double ymin, double ymax, int nrx, int nry)
{
  internal_pgplot_xmin = xmin;
  internal_pgplot_xmax = xmax;
  internal_pgplot_ymin = ymin;
  internal_pgplot_ymax = ymax;
  internal_pgplot_nrx = nrx;
  internal_pgplot_nry = nry;
}
void print_pgplot_version_used(FILE *stream)
{
  char version[25];
  int length = 20;
  cpgqinf("VERSION", version, &length);
  fprintf(stream, "%s (library)", version);
}
int pgplotMapCoordinate_dbl(double x, double y, int *nx, int *ny)
{
  int clip;
  clip = 0;
  *nx = 0.5 + (x - internal_pgplot_xmin)*(internal_pgplot_nrx-1.0)/(internal_pgplot_xmax-internal_pgplot_xmin);
  *ny = 0.5 + (y - internal_pgplot_ymin)*(internal_pgplot_nry-1.0)/(internal_pgplot_ymax-internal_pgplot_ymin);
  if(*nx < 0) {
    *nx = 0;
    clip = 1;
  }
  if(*nx >= internal_pgplot_nrx) {
    *nx = internal_pgplot_nrx - 1;
    clip = 1;
  }
  if(*ny < 0) {
    *ny = 0;
    clip = 1;
  }
  if(*ny >= internal_pgplot_nry) {
    *ny = internal_pgplot_nry - 1;
    clip = 1;
  }
  return clip;
}
int pgplotMapCoordinate(float x, float y, int *nx, int *ny)
{
  double x2, y2;
  x2 = x;
  y2 = y;
  return pgplotMapCoordinate_dbl(x2, y2, nx, ny);
}
void pgplotMapCoordinateInverse_dbl(double *x, double *y, int nx, int ny)
{
  if(internal_pgplot_nrx != 1)
    *x = ((double)nx)*(internal_pgplot_xmax-internal_pgplot_xmin)/((double)internal_pgplot_nrx-1.0) + internal_pgplot_xmin;
  else
    *x = internal_pgplot_xmin;
  if(internal_pgplot_nry != 1)
    *y = ((double)ny)*(internal_pgplot_ymax-internal_pgplot_ymin)/((double)internal_pgplot_nry-1.0) + internal_pgplot_ymin;
  else
    *y = internal_pgplot_ymin;
}
void pgplotMapCoordinateInverse(float *x, float *y, int nx, int ny)
{
  double x2, y2;
  pgplotMapCoordinateInverse_dbl(&x2, &y2, nx, ny);
  *x = x2;
  *y = y2;
}
void pgplotMapCoordinateBinSize(float *dx, float *dy)
{
  *dx = (internal_pgplot_xmax-internal_pgplot_xmin)/((double)internal_pgplot_nrx-1.0);
  *dy = (internal_pgplot_ymax-internal_pgplot_ymin)/((double)internal_pgplot_nry-1.0);
}
void pgplot_setWindowsize(int windowwidth, int windowheight, float aspectratio)
{
  float x, y;
  if(windowwidth > 0 && windowheight > 0) {
    x = windowwidth*0.01175548589341692789994673739445429916373;
    y = (windowheight-1)/(float)windowwidth;
    ppgpap(x, y);
  }else if(aspectratio > 0) {
    ppgpap(0, aspectratio);
  }
}
void pgplot_makeframe(pgplot_frame_def_internal *frame)
{
  if(frame->svp)
    ppgsvp(frame->svp_x1, frame->svp_x2, frame->svp_y1, frame->svp_y2);
  if(frame->swin) {
    ppgswin(frame->swin_x1, frame->swin_x2, frame->swin_y1, frame->swin_y2);
    if(internal_pgplot_nrx != 1) {
      frame->TR[0] = internal_pgplot_xmin - (internal_pgplot_xmax-internal_pgplot_xmin)/(float)(internal_pgplot_nrx-1);
      frame->TR[1] = (internal_pgplot_xmax-internal_pgplot_xmin)/(float)(internal_pgplot_nrx-1);
    }else {
      if(internal_pgplot_xmin != internal_pgplot_xmax) {
 frame->TR[0] = internal_pgplot_xmin - (internal_pgplot_xmax-internal_pgplot_xmin);
 frame->TR[1] = internal_pgplot_xmax-internal_pgplot_xmin;
      }else {
 frame->TR[0] = internal_pgplot_xmin - (frame->swin_x2 - frame->swin_x1);
 frame->TR[1] = frame->swin_x2 - frame->swin_x1;
      }
    }
    frame->TR[2] = 0;
    frame->TR[3] = internal_pgplot_ymin - (internal_pgplot_ymax-internal_pgplot_ymin)/(float)(internal_pgplot_nry-1); frame->TR[4] = 0; frame->TR[5] = (internal_pgplot_ymax-internal_pgplot_ymin)/(float)(internal_pgplot_nry-1);
    if(internal_pgplot_nry == 1) {
      frame->TR[3] = internal_pgplot_ymin;
      frame->TR[5] = (internal_pgplot_ymax-internal_pgplot_ymin);
      if(frame->TR[5] == 0) {
 frame->TR[3] = -1;
 frame->TR[5] = 1;
      }
    }
  }
}
int pgplot_process_text_option(char *textoption, int outline_txt, int outline_lw, int outline_color, int argc, char **argv, int xunit_type, datafile_definition data, verbose_definition verbose)
{
  int i, j;
  float text_x, text_y, text_ch;
  int text_ci, text_lw, text_f;
  char *newtext;
  if(textoption != NULL && argv != NULL && argc != 0) {
    for(j = 0; j < argc; j++) {
      if(strcmp(argv[j], textoption) == 0) {
 i = sscanf(argv[j+1], "%f %f %f %d %d %d", &text_x, &text_y, &text_ch, &text_lw, &text_f, &text_ci);
 if(xunit_type == 1) {
   text_x /= 360.0;
 }
 if(i != 6) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR pgplot_process_text_option: Error parsing %s option. Expected 6 items, but could only parse %d items.", textoption, i);
   return 0;
 }
 ppgsch(text_ch);
 ppgscf(text_f);
 newtext = str_replace_header_params(data, argv[j+2], verbose);
 if(newtext == NULL) {
   printerror(verbose.debug, "ERROR pgplotPAplot: Cannot substiture header parameter in -text option");
   return 0;
 }
 if(outline_txt > 0) {
   ppgslw(text_lw+outline_lw);
   ppgsci(outline_color);
   ppgptxt(text_x, text_y, 0, 0, newtext);
 }
 ppgsci(text_ci);
 ppgslw(text_lw);
 ppgptxt(text_x, text_y, 0, 0, newtext);
 ppgscf(1);
 ppgsch(1);
 ppgslw(1);
 ppgsci(1);
 j += 2;
 free(newtext);
      }
    }
  }
  return 1;
}
int pgplot_process_errorbars_options(char *herrorbaroption, char *herrorbaroption2, char *verrorbaroption, char *verrorbaroption2, int argc, char **argv, int xunit_type, pgplot_frame_def_internal *frame, verbose_definition verbose)
{
  int i, j;
  if(herrorbaroption != NULL && argv != NULL && argc != 0) {
    for(j = 0; j < argc; j++) {
      float herr_x1, herr_x2, herr_x3, herr_y, herr_size, herr_lw;
      int herr_ci, ok;
      ok = 0;
      if(strcmp(argv[j], herrorbaroption) == 0) {
 ok = 1;
      }else if(herrorbaroption2 != NULL) {
 if(strcmp(argv[j], herrorbaroption2) == 0) {
   ok = 1;
 }
      }
      if(ok) {
 i = sscanf(argv[j+1], "%f %f %f %f %f %f %d", &herr_x1, &herr_x2, &herr_x3, &herr_y, &herr_size, &herr_lw, &herr_ci);
 if(xunit_type == 1) {
   herr_x1 /= 360.0;
   herr_x2 /= 360.0;
   herr_x3 /= 360.0;
 }
 if(i != 7) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR pgplot_process_errorbars_options: Error parsing %s option. Expected 7 arguments, but only parsed %d.", herrorbaroption, i);
   return 0;
 }
 ppgsch(herr_size);
 ppgsci(herr_ci);
 ppgslw(herr_lw);
 if(frame == NULL) {
   ppgpt1(herr_x2, herr_y, 2);
   ppgerr1(1, herr_x2, herr_y, herr_x3-herr_x2, herr_size);
   ppgerr1(3, herr_x2, herr_y, herr_x2-herr_x1, herr_size);
 }else {
   ppgpt1(herr_x2, herr_y*(frame->swin_y2-frame->swin_y1)+frame->swin_y1, 2);
   ppgerr1(1, herr_x2, herr_y*(frame->swin_y2-frame->swin_y1)+frame->swin_y1, herr_x3-herr_x2, herr_size);
   ppgerr1(3, herr_x2, herr_y*(frame->swin_y2-frame->swin_y1)+frame->swin_y1, herr_x2-herr_x1, herr_size);
 }
 ppgsch(1);
 ppgslw(1);
 ppgsci(1);
 j += 1;
      }
    }
  }
  if(verrorbaroption != NULL && argv != NULL && argc != 0) {
    for(j = 0; j < argc; j++) {
      float verr_y1, verr_y2, verr_y3, verr_x, verr_size, verr_lw;
      int verr_ci, ok;
      ok = 0;
      if(strcmp(argv[j], verrorbaroption) == 0) {
 ok = 1;
      }else if(verrorbaroption2 != NULL) {
 if(strcmp(argv[j], verrorbaroption2) == 0) {
   ok = 1;
 }
      }
      if(ok) {
 i = sscanf(argv[j+1], "%f %f %f %f %f %f %d", &verr_x, &verr_y1, &verr_y2, &verr_y3, &verr_size, &verr_lw, &verr_ci);
 if(xunit_type == 1) {
   verr_x /= 360.0;
 }
 if(i != 7) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR pgplot_process_errorbars_options: Error parsing %s option. Expected 7 arguments, but only parsed %d.", verrorbaroption, i);
   return 0;
 }
 ppgsch(verr_size);
 ppgsci(verr_ci);
 ppgslw(verr_lw);
 if(frame == NULL) {
   ppgpt1(verr_x, verr_y2, 2);
   ppgerr1(2, verr_x, verr_y2, verr_y3-verr_y2, verr_size);
   ppgerr1(4, verr_x, verr_y2, verr_y2-verr_y1, verr_size);
 }else {
   ppgpt1(verr_x, verr_y2*(frame->swin_y2-frame->swin_y1)+frame->swin_y1, 2);
   ppgerr1(2, verr_x, verr_y2*(frame->swin_y2-frame->swin_y1)+frame->swin_y1, verr_y3-verr_y2, verr_size);
   ppgerr1(4, verr_x, verr_y2*(frame->swin_y2-frame->swin_y1)+frame->swin_y1, verr_y2-verr_y1, verr_size);
 }
 ppgsch(1);
 ppgslw(1);
 ppgsci(1);
 j += 1;
      }
    }
  }
  return 1;
}
void pgplot_overlay_paswing(datafile_definition data, int xunit_type, float pa_bottom, float PAoffset, float overlayalpha, float overlaybeta, float overlaypa0, float overlayl0, int overlayPAfine, int nrJumps, float *jump_longitudes, float *jump_offsets, verbose_definition verbose)
{
  long i, j;
  float x, xold, I, Iold;
  for(j = 0; j < (data.NrBins); j++) {
    for(i = 0; i < overlayPAfine; i++) {
      x = get_pulse_longitude(data, 0, j, verbose);
      if(j == data.NrBins-1) {
 x += (get_pulse_longitude(data, 0, j, verbose)-get_pulse_longitude(data, 0, j-1, verbose))*i/(float)overlayPAfine;
      }else {
 x += (get_pulse_longitude(data, 0, j+1, verbose)-get_pulse_longitude(data, 0, j, verbose))*i/(float)overlayPAfine;
      }
      I = paswing(overlayalpha, overlaybeta, x, overlaypa0-PAoffset, overlayl0, nrJumps, jump_longitudes, jump_offsets, 0, 0, 0);
      if(xunit_type == 1) {
 x /= 360.0;
      }
      if(I > 180)
 I -= 180;
      if(I > 180)
 I -= 180;
      if(I < 0)
 I += 180;
      if(I < 0)
 I += 180;
      I += PAoffset;
      if(j != 0 || i != 0) {
 if(fabs(I-Iold-180) < fabs(I-Iold))
   Iold -= 180;
 if(fabs(I-Iold+180) < fabs(I-Iold))
   Iold += 180;
 if(fabs(I-Iold) < 30) {
   ppgmove(xold, Iold);
   ppgdraw(x, I);
   if(pa_bottom < -10) {
     ppgmove(xold, Iold-180);
     ppgdraw(x, I-180);
   }
 }
      }
      Iold = I;
      xold = x;
    }
  }
}
void draw_papoints(datafile_definition data, float sigma_limit, float loffset, float PAoffset, int xunit_type, verbose_definition verbose)
{
  long j;
  int ok;
  float I;
  for(j = 0; j < (data.NrBins); j++) {
    ok = 1;
    if(data.offpulse_rms != NULL) {
      if(data.data[j+data.NrBins] < sigma_limit*data.offpulse_rms[1])
 ok = 0;
    }
    if(data.data[j+4*data.NrBins] < 0)
      ok = 0;
    if(ok && !isnan(data.data[j+3*data.NrBins]) && !isnan(get_pulse_longitude(data, 0, j, verbose))) {
      I = derotate_180(data.data[j+3*data.NrBins]+PAoffset) + 180.0;
      float xpos;
      xpos = get_pulse_longitude(data, 0, j, verbose)+loffset;
      if(xunit_type == 1) {
 xpos /= 360.0;
      }
      ppgpt1(xpos, I, -2);
      ppgerr1(6, xpos, I, data.data[j+4*data.NrBins], 1.0);
      I = derotate_180(data.data[j+3*data.NrBins]+PAoffset);
      ppgpt1(xpos, I, -2);
      ppgerr1(6, xpos, I, data.data[j+4*data.NrBins], 1.0);
      I = derotate_180(data.data[j+3*data.NrBins]+PAoffset) - 180.0;
      ppgpt1(xpos, I, -2);
      ppgerr1(6, xpos, I, data.data[j+4*data.NrBins], 1.0);
      I = derotate_180(data.data[j+3*data.NrBins]+PAoffset) - 360.0;
      ppgpt1(xpos, I, -2);
      ppgerr1(6, xpos, I, data.data[j+4*data.NrBins], 1.0);
    }
  }
}
int pgplotPAplot(datafile_definition data, int showtotpol, int nopaswing, int showEll, pgplot_options_definition *pgplot, char *xlabel, char *ylabel, char *ylabel_pa, char *ylabel_ell, char *ylabel_spstatfrac, float longitude_left, float longitude_right, int xunit_type, float loffset, float Imin, float Imax, float pa_bottom, float pa_top, float PAoffset, float sigma_limit, float datalinewidth, float ysize2, int dashed, int noynumbers, char *textoption, char *textoption_pa, float ytick_pa, int nysub_pa, char *textoption_ell, float ytick_ell, int nysub_ell, char *textoption_padist, char *textoption_elldist, char *textoption_spstat, char *herrorbaroption, char *herrorbaroptionpa, char *herrorbaroptionpa2, char *herrorbaroptionell, char *herrorbaroptionpadist, char *herrorbaroptionelldist, char *verrorbaroption, char *verrorbaroptionpa, char *verrorbaroptionpa2, char *verrorbaroptionell, int argc, char **argv, int outline_txt, int outline_lw, int outline_color, int overlayPA, float overlayalpha, float overlaybeta, float overlaypa0, float overlayl0, int overlayPAfine, int nrJumps, float *jump_longitudes, float *jump_offsets, datafile_definition *padist, float padist_pamin, float padist_pamax, float padist_saturize, int padist_overlayavpa, int padist_paswing, datafile_definition *elldist, float elldist_saturize, int nowedge, datafile_definition *spstatfrac, verbose_definition verbose)
{
  int deviceID, ok, domove, showPAdist, showELLdist, showspstatfrac;
  float ymin, ymax, I, overallplotscaling;
  long i, j;
  char *newtext;
  pgplot_frame_def_internal frame;
  pgplot_options_definition *pgplot_backup;
  int showwedge = 1;
  if(nowedge)
    showwedge = 0;
  showPAdist = 0;
  if(padist != NULL) {
    if(padist->format != MEMORY_format) {
      printerror(verbose.debug, "ERROR pgplotPAplot: PA distribution does not appear to be loaded into memory\n");
      return 0;
    }
    if(padist->NrPols != 1) {
      printerror(verbose.debug, "ERROR pgplotPAplot: PA distribution is expected to have a single polarization channel\n");
      return 0;
    }
    if(padist->gentype != GENTYPE_PADIST) {
      printwarning(verbose.debug, "WARNING pgplotPAplot: The file opened as a PA distribution does not have a gentype set to be a PA distribution. It is currently set to %s.", returnGenType_str(padist->gentype));
      return 0;
    }
    if(padist->NrBins != data.NrBins) {
      printerror(verbose.debug, "ERROR pgplotPAplot: The number of pulse longitude bins is different in the PA distribution compared to the profile.");
      return 0;
    }
    if(padist_pamax - padist_pamin < 0) {
      float junk;
      junk = padist_pamax;
      padist_pamax = padist_pamin;
      padist_pamin = junk;
    }
    if(padist_pamax - padist_pamin > 360) {
      printerror(verbose.debug, "ERROR pgplotPAplot: The PA distribution cannot cover more than 360 degrees.");
      return 0;
    }
    showPAdist = 1;
  }
  showELLdist = 0;
  if(elldist != NULL) {
    if(elldist->format != MEMORY_format) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Ellipticity angle distribution does not appear to be loaded into memory\n");
      return 0;
    }
    if(elldist->NrPols != 1) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Ellipticity angle distribution is expected to have a single polarization channel\n");
      return 0;
    }
    if(elldist->gentype != GENTYPE_ELLDIST) {
      printwarning(verbose.debug, "WARNING pgplotPAplot: The file opened as an ellipticity angle distribution does not have a gentype set to be an ellipticity angle distribution. It is currently set to %s.", returnGenType_str(elldist->gentype));
      return 0;
    }
    if(elldist->NrBins != data.NrBins) {
      printerror(verbose.debug, "ERROR pgplotPAplot: The number of pulse longitude bins is different in the ellipticity angle distribution compared to the profile.");
      return 0;
    }
    showELLdist = 1;
  }
  showspstatfrac = 0;
  if(overlayPAfine)
    overlayPAfine = 100;
  else
    overlayPAfine = 1;
  overallplotscaling = 1;
  if(showEll && (showPAdist == 0 && showELLdist == 0)) {
    overallplotscaling = 0.78;
  }else if(showEll && (showPAdist || showELLdist)) {
    overallplotscaling = 0.53;
  }else if(showEll == 0 && (showPAdist || showELLdist)) {
    overallplotscaling = 0.63;
  }
  if(showPAdist && showELLdist) {
    overallplotscaling = 0.4;
  }
  if(showPAdist && showELLdist && showspstatfrac) {
    overallplotscaling = 0.35;
  }
  pgplot_backup = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
  if(pgplot_backup == NULL) {
    printerror(verbose.debug, "ERROR pgplotPAplot: Memory allocation error\n");
    return 0;
  }
  memcpy(pgplot_backup, pgplot, sizeof(pgplot_options_definition));
  pgplot->box.box_labelsize *= overallplotscaling;
  pgplot->box.label_ch *= overallplotscaling;
  float vp_left, vp_right, vp_profile_top, vp_profile_bottom, vp_spstatfrac_top, vp_spstatfrac_bottom, vp_paswing_top, vp_paswing_bottom, vp_ellswing_top, vp_ellswing_bottom, vp_padist_top, vp_padist_bottom, vp_elldist_top, vp_elldist_bottom;
  vp_left = 0.1+pgplot->viewport.dxplot;
  vp_right = 0.1+pgplot->viewport.dxplot+0.8*pgplot->viewport.xsize*overallplotscaling;
  vp_profile_bottom = 0.35+(1-overallplotscaling)*0.6+pgplot->viewport.dyplot;
  vp_profile_top = vp_profile_bottom +0.55*pgplot->viewport.ysize*overallplotscaling;
  if(showspstatfrac) {
    vp_spstatfrac_bottom = vp_profile_bottom -0.25*pgplot->viewport.ysize*ysize2*overallplotscaling;
    vp_spstatfrac_top = vp_profile_top - 0.55*pgplot->viewport.ysize*overallplotscaling;
  }else {
    vp_spstatfrac_bottom = vp_profile_bottom;
    vp_spstatfrac_top = vp_profile_bottom + overallplotscaling*0.25*pgplot->viewport.ysize*ysize2;
  }
  if(nopaswing == 0) {
    vp_paswing_bottom = vp_spstatfrac_bottom - overallplotscaling*0.25*pgplot->viewport.ysize*ysize2;
    vp_paswing_top = vp_spstatfrac_top - overallplotscaling*0.25*pgplot->viewport.ysize*ysize2;
  }else {
    vp_paswing_bottom = vp_spstatfrac_bottom;
    vp_paswing_top = vp_spstatfrac_bottom + overallplotscaling*0.25*pgplot->viewport.ysize*ysize2;
  }
  if(showEll) {
    vp_ellswing_bottom = vp_paswing_bottom - overallplotscaling*0.25*pgplot->viewport.ysize*ysize2;
    vp_ellswing_top = vp_paswing_top - overallplotscaling*0.25*pgplot->viewport.ysize*ysize2;
  }else {
    vp_ellswing_bottom = vp_paswing_bottom;
    vp_ellswing_top = vp_paswing_top;
  }
  if(showPAdist) {
    vp_padist_top = vp_ellswing_bottom;
    vp_padist_bottom = vp_padist_top - overallplotscaling*0.5*pgplot->viewport.ysize*ysize2;
  }else {
    vp_padist_top = vp_ellswing_top;
    vp_padist_bottom = vp_ellswing_bottom;
  }
  if(showELLdist) {
    vp_elldist_top = vp_padist_bottom;
    vp_elldist_bottom = vp_elldist_top - overallplotscaling*0.5*pgplot->viewport.ysize*ysize2;
  }else {
    vp_elldist_top = vp_padist_top;
    vp_elldist_bottom = vp_padist_bottom;
  }
  if(pgplot_opendevice(&(pgplot->viewport), &deviceID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Cannot open plot device");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  if(data.poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does only have PA points, but no profile");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  if(data.poltype != POLTYPE_ILVPAdPA && data.poltype != POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does not appear to contain PA data");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  if(data.poltype == POLTYPE_ILVPAdPA && data.NrPols != 5) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does not appear to contain PA data as NrPols != 5");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  if(data.poltype == POLTYPE_ILVPAdPATEldEl && data.NrPols != 8) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does not appear to contain PA data as NrPols != 8");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  if((showtotpol || showEll) && data.poltype != POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does not appear to contain the total polarization and/or ellipticity angle, but it was requested to be plotted");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  newtext = str_replace_header_params(data, pgplot->box.title, verbose);
  if(newtext == NULL) {
    printerror(verbose.debug, "ERROR pgplotPAplot: Cannot substuture header parameter in the tile");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  strcpy(pgplot->box.title, newtext);
  free(newtext);
  clear_pgplot_frame(&frame);
  ymin = ymax = NAN;
  for(j = 1; j < (data.NrBins); j++) {
    if(data.data[j] > ymax || isnan(ymax))
      ymax = data.data[j];
    if(data.data[j] < ymin || isnan(ymin))
      ymin = data.data[j];
    if(data.data[j+data.NrBins] > ymax)
      ymax = data.data[j+data.NrBins];
    if(data.data[j+data.NrBins] < ymin)
      ymin = data.data[j+data.NrBins];
    if(data.data[j+2*data.NrBins] > ymax)
      ymax = data.data[j+2*data.NrBins];
    if(data.data[j+2*data.NrBins] < ymin)
      ymin = data.data[j+2*data.NrBins];
    if(showtotpol) {
      if(data.data[j+5*data.NrBins] > ymax)
 ymax = data.data[j+5*data.NrBins];
      if(data.data[j+5*data.NrBins] < ymin)
 ymin = data.data[j+5*data.NrBins];
    }
  }
  ppgsvp(vp_left, vp_right, vp_profile_bottom, vp_profile_top);
  frame.swin_x1 = longitude_left;
  frame.swin_x2 = longitude_right;
  if(xunit_type == 1) {
    frame.swin_x1 /= 360.0;
    frame.swin_x2 /= 360.0;
  }
  if(Imin != Imax) {
    frame.swin_y1 = Imin;
    frame.swin_y2 = Imax;
    ppgswin(frame.swin_x1,frame.swin_x2,Imin,Imax);
  }else {
    frame.swin_y1 = ymin-(ymax-ymin)*0.05;
    frame.swin_y2 = ymax+(ymax-ymin)*0.05;
    ppgswin(frame.swin_x1,frame.swin_x2,ymin-(ymax-ymin)*0.05,ymax+(ymax-ymin)*0.05);
  }
  if(showspstatfrac == 0 && nopaswing == 1 && showEll == 0 && showPAdist == 0 && showELLdist == 0) {
    strcpy(pgplot->box.xlabel, xlabel);
    strcpy(pgplot->box.box_xopt, "bcnst");
  }else {
    pgplot->box.xlabel[0] = 0;
    strcpy(pgplot->box.box_xopt, "bcst");
  }
  if(noynumbers) {
    strcpy(pgplot->box.box_yopt, "bcts");
  }else {
    strcpy(pgplot->box.box_yopt, "bcnts");
  }
  if(!noynumbers) {
    pgplot->box.drawlabels = 1;
    strcpy(pgplot->box.ylabel, ylabel);
  }
  pgplot_drawbox(&(pgplot->box));
  pgplot->box.drawtitle = 0;
  ppgslw(datalinewidth);
  if(dashed) {
    ppgsls(4);
  }else
    ppgsci(3);
  domove = 1;
  for(j = 0; j < (data.NrBins); j++) {
    float xpos;
    xpos = get_pulse_longitude(data, 0, j, verbose)+loffset;
    if(xunit_type == 1) {
      xpos /= 360.0;
    }
    if(isnan(xpos) || isnan(data.data[j+2*data.NrBins])) {
      domove = 1;
    }else {
      if(domove) {
 ppgmove(xpos, data.data[j+2*data.NrBins]);
 domove = 0;
      }else {
 ppgdraw(xpos, data.data[j+2*data.NrBins]);
      }
    }
  }
  if(dashed)
    ppgsls(2);
  else
    ppgsci(2);
  domove = 1;
  for(j = 0; j < (data.NrBins); j++) {
    float xpos;
    xpos = get_pulse_longitude(data, 0, j, verbose)+loffset;
    if(xunit_type == 1) {
      xpos /= 360.0;
    }
    if(isnan(xpos) || isnan(data.data[j+1*data.NrBins])) {
      domove = 1;
    }else {
      if(domove) {
 ppgmove(xpos, data.data[j+1*data.NrBins]);
 domove = 0;
      }else {
 ppgdraw(xpos, data.data[j+1*data.NrBins]);
      }
    }
  }
  if(showtotpol) {
    if(dashed) {
      ppgsls(1);
    }
    if(!dashed) {
      ppgsci(4);
      ppgsls(2);
    }
    domove = 1;
    for(j = 0; j < (data.NrBins); j++) {
      float xpos;
      xpos = get_pulse_longitude(data, 0, j, verbose)+loffset;
      if(xunit_type == 1) {
 xpos /= 360.0;
      }
      if(isnan(xpos) || isnan(data.data[j+5*data.NrBins])) {
 domove = 1;
      }else {
 if(domove) {
   ppgmove(xpos, data.data[j+5*data.NrBins]);
   domove = 0;
 }else {
   ppgdraw(xpos, data.data[j+5*data.NrBins]);
 }
      }
    }
  }
  ppgsci(1);
  ppgsls(1);
  domove = 1;
  for(j = 0; j < (data.NrBins); j++) {
    float xpos;
    xpos = get_pulse_longitude(data, 0, j, verbose)+loffset;
    if(xunit_type == 1) {
      xpos /= 360.0;
    }
    if(isnan(xpos) || isnan(data.data[j])) {
      domove = 1;
    }else {
      if(domove) {
 ppgmove(xpos, data.data[j]);
 domove = 0;
      }else {
 ppgdraw(xpos, data.data[j]);
      }
    }
  }
  ppgsls(1);
  if(pgplot_process_errorbars_options(herrorbaroption, NULL, verrorbaroption, NULL, argc, argv, xunit_type, &frame, verbose) == 0) {
    printerror(verbose.debug, "ERROR pgplotPAplot: Error processing the %s and %s options", herrorbaroption, verrorbaroption);
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  if(pgplot_process_text_option(textoption, outline_txt, outline_lw, outline_color, argc, argv, xunit_type, data, verbose) == 0) {
    printerror(verbose.debug, "ERROR pgplotPAplot: Error processing %s option.", textoption);
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
  }
  if(nopaswing == 0) {
    ppgsls(1);
    ppgsci(1);
    ppgslw(1);
    ppgsvp(vp_left, vp_right, vp_paswing_bottom, vp_paswing_top);
    frame.swin_x1 = longitude_left;
    frame.swin_x2 = longitude_right;
    if(xunit_type == 1) {
      frame.swin_x1 /= 360.0;
      frame.swin_x2 /= 360.0;
    }
    frame.swin_y1 = pa_bottom;
    frame.swin_y2 = pa_top;
    ppgswin(frame.swin_x1,frame.swin_x2,pa_bottom,pa_top);
    ppgslw(pgplot->box.box_lw);
    pgplot->box.drawlabels = 1;
    if(showEll == 0 && showPAdist == 0 && showELLdist == 0) {
      strcpy(pgplot->box.xlabel, xlabel);
      strcpy(pgplot->box.box_xopt, "bcnst");
    }else {
      pgplot->box.xlabel[0] = 0;
      strcpy(pgplot->box.box_xopt, "bcst");
    }
    if(!noynumbers) {
      strcpy(pgplot->box.ylabel, ylabel_pa);
    }else {
      pgplot->box.ylabel[0] = 0;
    }
    float tmp_ytick = pgplot->box.box_ytick;
    int tmp_nysub = pgplot->box.box_nysub;
    pgplot->box.box_ytick = ytick_pa;
    pgplot->box.box_nysub = nysub_pa;
    pgplot_drawbox(&(pgplot->box));
    pgplot->box.box_ytick = tmp_ytick;
    pgplot->box.box_nysub = tmp_nysub;
    ppgslw(datalinewidth);
    if(overlayPA) {
      ppgsci(2);
      pgplot_overlay_paswing(data, xunit_type, pa_bottom, 0.0, overlayalpha, overlaybeta, overlaypa0, overlayl0, overlayPAfine, nrJumps, jump_longitudes, jump_offsets, verbose);
      ppgsci(1);
    }
    ppgsci(1);
    draw_papoints(data, sigma_limit, loffset, PAoffset, xunit_type, verbose);
    if(pgplot_process_errorbars_options(herrorbaroptionpa, herrorbaroptionpa2, verrorbaroptionpa, verrorbaroptionpa2, argc, argv, xunit_type, NULL, verbose) == 0) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Error processing the %s and %s options", herrorbaroptionpa, verrorbaroptionpa);
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    if(pgplot_process_text_option(textoption_pa, outline_txt, outline_lw, outline_color, argc, argv, xunit_type, data, verbose) == 0) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Error processing %s option.", textoption_pa);
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
    }
  }
  if(showEll) {
    ppgsls(1);
    ppgsci(1);
    ppgslw(1);
    ppgsvp(vp_left, vp_right, vp_ellswing_bottom, vp_ellswing_top);
    frame.swin_x1 = longitude_left;
    frame.swin_x2 = longitude_right;
    if(xunit_type == 1) {
      frame.swin_x1 /= 360.0;
      frame.swin_x2 /= 360.0;
    }
    frame.swin_y1 = -45;
    frame.swin_y2 = +45;
    ppgswin(frame.swin_x1,frame.swin_x2,-45,45);
    ppgslw(pgplot->box.box_lw);
    pgplot->box.drawlabels = 1;
    if(!noynumbers) {
      strcpy(pgplot->box.ylabel, ylabel_ell);
    }else {
      pgplot->box.ylabel[0] = 0;
    }
    if(showPAdist == 0 && showELLdist == 0) {
      strcpy(pgplot->box.xlabel, xlabel);
      strcpy(pgplot->box.box_xopt, "bcnst");
    }else {
      pgplot->box.xlabel[0] = 0;
      strcpy(pgplot->box.box_xopt, "bcst");
    }
    float tmp_ytick = pgplot->box.box_ytick;
    int tmp_nysub = pgplot->box.box_nysub;
    pgplot->box.box_ytick = ytick_ell;
    pgplot->box.box_nysub = nysub_ell;
    pgplot_drawbox(&(pgplot->box));
    pgplot->box.box_ytick = tmp_ytick;
    pgplot->box.box_nysub = tmp_nysub;
    ppgslw(datalinewidth);
    for(j = 0; j < (data.NrBins); j++) {
      ok = 1;
      if(data.offpulse_rms != NULL) {
 if(data.data[j+5*data.NrBins] < sigma_limit*data.offpulse_rms[5])
   ok = 0;
      }
      if(data.data[j+7*data.NrBins] < 0)
 ok = 0;
      if(ok && !isnan(data.data[j+6*data.NrBins]) && !isnan(get_pulse_longitude(data, 0, j, verbose))) {
 I = derotate_90(data.data[j+6*data.NrBins]) + 90.0;
 float xpos;
 xpos = get_pulse_longitude(data, 0, j, verbose)+loffset;
 if(xunit_type == 1) {
   xpos /= 360.0;
 }
 ppgpt1(xpos, I, -2);
 ppgerr1(6, xpos, I, data.data[j+7*data.NrBins], 1.0);
 I = derotate_90(data.data[j+6*data.NrBins]);
 ppgpt1(xpos, I, -2);
 ppgerr1(6, xpos, I, data.data[j+7*data.NrBins], 1.0);
 I = derotate_90(data.data[j+6*data.NrBins]) - 90.0;
 ppgpt1(xpos, I, -2);
 ppgerr1(6, xpos, I, data.data[j+7*data.NrBins], 1.0);
 I = derotate_90(data.data[j+6*data.NrBins]) - 180.0;
 ppgpt1(xpos, I, -2);
 ppgerr1(6, xpos, I, data.data[j+7*data.NrBins], 1.0);
      }
    }
    if(pgplot_process_errorbars_options(herrorbaroptionell, NULL, verrorbaroptionell, NULL, argc, argv, xunit_type, NULL, verbose) == 0) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Error processing the %s and %s options", herrorbaroptionell, verrorbaroptionell);
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    if(pgplot_process_text_option(textoption_ell, outline_txt, outline_lw, outline_color, argc, argv, xunit_type, data, verbose) == 0) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Error processing %s option.", textoption_ell);
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
    }
  }
  int cmaptype;
  cmaptype = pgplot_device_type(NULL, verbose);
  if(cmaptype <= 2)
    cmaptype = PPGPLOT_GRAYSCALE;
  else
    cmaptype = PPGPLOT_INVERTED_GRAYSCALE;
  if(showPAdist) {
    ppgsls(1);
    ppgsci(1);
    ppgslw(1);
    ppgsvp(vp_left, vp_right, vp_padist_bottom, vp_padist_top);
    ppgslw(pgplot->box.box_lw);
    pgplot->box.drawlabels = 1;
    if(!noynumbers) {
      strcpy(pgplot->box.ylabel, ylabel_pa);
    }else {
      pgplot->box.ylabel[0] = 0;
    }
    if(showELLdist == 0) {
      strcpy(pgplot->box.xlabel, xlabel);
      strcpy(pgplot->box.box_xopt, "bcnst");
    }else {
      pgplot->box.xlabel[0] = 0;
      strcpy(pgplot->box.box_xopt, "bcst");
    }
    ppgslw(datalinewidth);
    pgplot->viewport.dontopen = 1;
    pgplot->viewport.dontclose = 1;
    float xpos1, xpos2, xpos_start, xpos_end;
    xpos1 = get_pulse_longitude(data, 0, 0, verbose)+loffset;
    xpos2 = get_pulse_longitude(data, 0, padist->NrBins-1, verbose)+loffset;
    xpos_start = longitude_left;
    xpos_end = longitude_right;
    if(xunit_type) {
      xpos1 /= 360.0;
      xpos2 /= 360.0;
      xpos_start /= 360.0;
      xpos_end /= 360.0;
    }
    float ymin_data, ymax_data;
    ymin_data = -90;
    ymax_data = 90;
    long wraps;
    wraps = floor((padist_pamin-ymin_data)/180.0);
    ymin_data += 180.0*wraps;
    ymax_data += 180.0*wraps;
    wraps = ceil((padist_pamax-ymax_data)/180.0)+1;
    if(wraps < 1)
      wraps = 1;
    float *data_ptr;
    if(wraps == 1) {
      data_ptr = padist->data;
    }else if(wraps > 10) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Cannot plot pa-distribution more than 10 times above each other. Something appears to be wrong with the requested PA range for the PA distribution panel.\n");
      return 0;
    }else {
      long blocksize;
      blocksize = padist->NrBins*padist->NrSubints;
      data_ptr = (float *)malloc(wraps*blocksize*sizeof(float));
      if(data_ptr == NULL) {
 printerror(verbose.debug, "ERROR pgplotPAplot: Memory allocation error\n");
 return 0;
      }
      for(i = 0; i < wraps; i++) {
 memcpy(data_ptr+i*blocksize, padist->data, blocksize*sizeof(float));
      }
      ymax_data += 180.0*(wraps-1);
    }
    ymin_data += 0.5*(180.0/(float)padist->NrSubints);
    ymax_data -= 0.5*(180.0/(float)padist->NrSubints);
    if(pgplotMap(pgplot, data_ptr, padist->NrBins, wraps*padist->NrSubints, xpos1, xpos2, xpos_start, xpos_end, ymin_data, ymax_data, padist_pamin, padist_pamax, cmaptype, 0, 0, 0, NULL, 0, 0, padist_saturize, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, showwedge, 0, 1, 0, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplotPAplot: Plotting of the PA distribution failed.");
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    if(padist_overlayavpa) {
      ppgsci(2);
      draw_papoints(data, sigma_limit, loffset, PAoffset, xunit_type, verbose);
      if(padist_overlayavpa == 2) {
 ppgsci(3);
 draw_papoints(data, sigma_limit, loffset, PAoffset+90, xunit_type, verbose);
      }
      ppgsci(1);
    }
    ppgslw(datalinewidth);
    if(padist_paswing) {
      ppgsci(2);
      pgplot_overlay_paswing(data, xunit_type, pa_bottom, -90.0, overlayalpha, overlaybeta, overlaypa0, overlayl0, overlayPAfine, nrJumps, jump_longitudes, jump_offsets, verbose);
      ppgsci(1);
    }
    if(wraps > 1) {
      free(data_ptr);
    }
    pgplot_drawbox(&(pgplot->box));
    if(pgplot_process_errorbars_options(herrorbaroptionpadist, NULL, NULL, NULL, argc, argv, xunit_type, NULL, verbose) == 0) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Error processing the %s option", herrorbaroptionpadist);
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    if(pgplot_process_text_option(textoption_padist, outline_txt, outline_lw, outline_color, argc, argv, xunit_type, data, verbose) == 0) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Error processing %s option.", textoption_padist);
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
    }
  }
  if(showELLdist) {
    ppgsls(1);
    ppgsci(1);
    ppgslw(1);
    ppgsvp(vp_left, vp_right, vp_elldist_bottom, vp_elldist_top);
    ppgslw(pgplot->box.box_lw);
    pgplot->box.drawlabels = 1;
    strcpy(pgplot->box.xlabel, xlabel);
    if(!noynumbers) {
      strcpy(pgplot->box.ylabel, ylabel_ell);
    }else {
      pgplot->box.ylabel[0] = 0;
    }
    strcpy(pgplot->box.box_xopt, "bcnst");
    ppgslw(datalinewidth);
    pgplot->viewport.dontopen = 1;
    pgplot->viewport.dontclose = 1;
    float xpos1, xpos2, xpos_start, xpos_end;
    xpos1 = get_pulse_longitude(data, 0, 0, verbose)+loffset;
    xpos2 = get_pulse_longitude(data, 0, elldist->NrBins-1, verbose)+loffset;
    xpos_start = longitude_left;
    xpos_end = longitude_right;
    if(xunit_type) {
      xpos1 /= 360.0;
      xpos2 /= 360.0;
      xpos_start /= 360.0;
      xpos_end /= 360.0;
    }
    if(pgplotMap(pgplot, elldist->data, elldist->NrBins, elldist->NrSubints, xpos1, xpos2, xpos_start, xpos_end, -45+0.5*(90.0/(float)elldist->NrSubints), 45-0.5*(90.0/(float)elldist->NrSubints), -45, 45, cmaptype, 0, 0, 0, NULL, 0, 0, elldist_saturize, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, showwedge, 0, 1, 0, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplotPAplot: Plotting of the ellipticity angle distribution failed.");
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    pgplot_drawbox(&(pgplot->box));
    if(pgplot_process_errorbars_options(herrorbaroptionelldist, NULL, NULL, NULL, argc, argv, xunit_type, NULL, verbose) == 0) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Error processing the %s option", herrorbaroptionelldist);
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    if(pgplot_process_text_option(textoption_elldist, outline_txt, outline_lw, outline_color, argc, argv, xunit_type, data, verbose) == 0) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Error processing %s option.", textoption_elldist);
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
    }
  }
  if(pgplot->viewport.dontclose == 0)
    ppgclos();
  memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
  free(pgplot_backup);
  return 1;
}
int pgplotGraph1(pgplot_options_definition *pgplot, float *data, float *datax, float *sigma, int nrx, float xmin, float xmax, int dontsetranges, float xmin_show, float xmax_show, float ymin_show, float ymax_show, int forceMinZero, int hist, int noline, int linewidth, int pointtype, int color, int boxcolor, pulselongitude_regions_definition *regions, int onpulsecolor, verbose_definition verbose)
{
  float min, max;
  float x, xnext, y;
  int binnr, regnr, loopnr, deviceID, notset, color_prev, color_cur, curpencolor;
  pgplot_frame_def_internal frame;
  if(hist && datax != NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotGraph1: Cannot plot histogram when x-values are provided");
    return 0;
  }
  if(hist && nrx == 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotGraph1: Cannot plot histogram for one point");
    return 0;
  }
  if(hist == 0 && regions != NULL) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING pgplotGraph1: When showing selected regions, histogram mode is more precise");
  }
  if(regions != NULL) {
    for(regnr = 0; regnr < regions->nrRegions; regnr++) {
      if(regions->bins_defined[regnr] == 0) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR pgplotGraph1: onpulse region is not defined in bins");
 return 0;
      }
    }
  }
  xnext = 0;
  if(pgplot_opendevice(&(pgplot->viewport), &deviceID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotGraph1: Cannot open plot device");
    return 0;
  }
  notset = 1;
  min = max = 0;
  if(datax != NULL) {
    xmin = xmax = datax[0];
  }
  for(binnr = 0; binnr < nrx; binnr++) {
    if(datax == NULL) {
      if(nrx == 1) {
 x = xmin;
      }else {
 x = binnr*(xmax-xmin)/(float)(nrx-1) + xmin;
      }
    }else {
      x = datax[binnr];
    }
    if(xmin_show == xmax_show || (xmin_show < xmax_show && x >= xmin_show && x <= xmax_show) || (xmin_show > xmax_show && x <= xmin_show && x >= xmax_show)) {
      if(notset) {
 min = max = data[binnr];
 notset = 0;
      }
      if(data[binnr] > max)
 max = data[binnr];
      if(data[binnr] < min)
 min = data[binnr];
    }
    if(datax != NULL) {
      if(datax[binnr] > xmax)
 xmax = datax[binnr];
      if(datax[binnr] < xmin)
 xmin = datax[binnr];
    }
  }
  if(forceMinZero)
    min = 0;
  if(min == max) {
    min -= 1.0;
    max += 1.0;
  }
  if(ymin_show != ymax_show) {
    min = ymin_show;
    max = ymax_show;
  }else {
    if(forceMinZero == 0 || hist == 0)
      min -= (max-min)*0.01;
    max += (max-min)*0.01;
  }
  if(xmin_show == xmax_show) {
    xmin_show = xmin;
    xmax_show = xmax;
  }
  if(hist) {
    xmin_show -= 0.5*(xmax-xmin)/(float)(nrx-1);
    xmax_show += 0.5*(xmax-xmin)/(float)(nrx-1);
  }
  clear_pgplot_frame(&frame);
  frame.svp = 1;
  frame.svp_x1 = 0.1+pgplot->viewport.dxplot;
  frame.svp_x2 = 0.1+pgplot->viewport.dxplot+0.8*pgplot->viewport.xsize;
  frame.svp_y1 = 0.1+pgplot->viewport.dyplot;
  frame.svp_y2 = 0.1+pgplot->viewport.dyplot+0.8*pgplot->viewport.ysize;
  if(dontsetranges == 0)
    frame.swin = 1;
  else
    frame.swin = 0;
  frame.swin_x1 = xmin_show;
  frame.swin_x2 = xmax_show;
  frame.swin_y1 = min;
  frame.swin_y2 = max;
  pgplot_makeframe(&frame);
  ppgsci(boxcolor);
  curpencolor = boxcolor;
  pgplot_drawbox(&(pgplot->box));
  color_prev = color;
  ppgbbuf();
  ppgslw(linewidth);
  for(loopnr = 0; loopnr < 2; loopnr++) {
    for(binnr = 0; binnr < nrx; binnr++) {
      if(datax == NULL) {
 if(nrx == 1) {
   x = xmin;
   xnext = xmin;
 }else {
   x = binnr*(xmax-xmin)/(float)(nrx-1) + xmin;
   xnext = (binnr+1)*(xmax-xmin)/(float)(nrx-1) + xmin;
 }
 if(hist) {
   x -= 0.5*(xmax-xmin)/(float)(nrx-1);
   xnext -= 0.5*(xmax-xmin)/(float)(nrx-1);
 }
      }else {
 x = datax[binnr];
      }
      y = data[binnr];
      color_cur = color;
      if(regions != NULL) {
 for(regnr = 0; regnr < regions->nrRegions; regnr++) {
   if(hist == 0) {
     if(binnr > regions->left_bin[regnr] && binnr <= regions->right_bin[regnr]) {
       if(onpulsecolor < 0) {
  color_cur = color+1+regnr;
       }else {
  color_cur = onpulsecolor;
       }
     }
   }else {
     if(binnr >= regions->left_bin[regnr] && binnr <= regions->right_bin[regnr]) {
       if(onpulsecolor < 0) {
  color_cur = color+1+regnr;
       }else {
  color_cur = onpulsecolor;
       }
     }
   }
 }
      }
      if(loopnr == 0) {
 if(noline == 0) {
   if(hist == 0) {
     if(binnr == 0) {
       ppgmove(x, y);
     }else {
       if(color_cur != curpencolor) {
  ppgsci(color_cur);
  curpencolor = color_cur;
       }
       ppgdraw(x, y);
     }
   }else {
     if(binnr == 0) {
       if(forceMinZero) {
  ppgmove(x, 0);
  if(color_cur != curpencolor) {
    ppgsci(color_cur);
    curpencolor = color_cur;
  }
  ppgdraw(x, y);
       }else {
  ppgmove(x, y);
       }
     }else {
       if(color_cur == color || color_prev == color) {
  if(color != curpencolor) {
    ppgsci(color);
    curpencolor = color;
  }
       }else {
  if(color_cur == color_prev) {
    if(color_cur != curpencolor) {
      ppgsci(color_cur);
      curpencolor = color_cur;
    }
  }else {
    if(color_prev != curpencolor) {
      ppgsci(color_prev);
      curpencolor = color_prev;
    }
  }
       }
       ppgdraw(x, y);
     }
     ppgsci(color_cur);
     ppgdraw(xnext, y);
     if(binnr == nrx - 1) {
       if(forceMinZero) {
  ppgdraw(xnext, 0);
       }
     }
   }
 }
      }else {
 ppgsci(color);
 if(pointtype) {
   ppgpt1(x, y, pointtype);
 }
 if(sigma != NULL) {
   ppgerr1(6, x, y, sigma[binnr], 1.0);
 }
      }
      color_prev = color_cur;
    }
  }
  ppgebuf();
  ppgsci(1);
  if(!pgplot->viewport.dontclose)
    ppgclos();
  return 1;
}
int pgplot_opendevice(pgplot_viewport_definition *viewport, int *deviceID, verbose_definition verbose)
{
  if(viewport->dontopen == 0) {
    *deviceID = ppgopen(viewport->plotDevice);
    if(*deviceID <= 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplot_opendevice: Cannot open plot device called '%s'", viewport->plotDevice);
      return 0;
    }
    pgplot_setWindowsize(viewport->windowwidth, viewport->windowheight, viewport->aspectratio);
    ppgask(0);
  }
  if(viewport->noclear == 0)
    ppgpage();
  return 1;
}
void pgplot_clear_options(pgplot_options_definition *pgplot) {
  pgplot_clear_viewport_def(&(pgplot->viewport));
  clear_pgplot_box(&(pgplot->box));
}
void pgplot_clear_viewport_def(pgplot_viewport_definition *viewport)
{
  viewport->windowwidth = 0;
  viewport->windowheight = 0;
  viewport->aspectratio = -1;
  viewport->dxplot = 0;
  viewport->xsize = 1;
  viewport->dyplot = 0;
  viewport->ysize = 1;
  viewport->noclear = 0;
  viewport->dontopen = 0;
  viewport->dontclose = 0;
  sprintf(viewport->plotDevice, "?");
}
int pgplot_set_maptype(int maptype, verbose_definition verbose)
{
  float contrast, brightness;
  contrast = 1.0;
  brightness = 0.5;
  if(maptype == PPGPLOT_GRAYSCALE) {
    float datal[] = {0.0, 1.0};
    float datar[] = {1.0, 0.0};
    float datag[] = {1.0, 0.0};
    float datab[] = {1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_GRAYSCALE) {
    float datal[] = {0.0, 1.};
    float datar[] = {0.0, 1.};
    float datag[] = {0.0, 1.};
    float datab[] = {0.0, 1.};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_RED) {
    float datal[] = {0.0, 1.0};
    float datar[] = {1.0, 0.0};
    float datag[] = {0.0, 0.0};
    float datab[] = {0.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_RED) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 1.0};
    float datag[] = {0.0, 0.0};
    float datab[] = {0.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_GREEN) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {1.0, 0.0};
    float datab[] = {0.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_GREEN) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {0.0, 1.0};
    float datab[] = {0.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_BLUE) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {0.0, 0.0};
    float datab[] = {1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_BLUE) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {0.0, 0.0};
    float datab[] = {0.0, 1.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_CYAN) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {1.0, 0.0};
    float datab[] = {1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_CYAN) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {0.0, 1.0};
    float datab[] = {0.0, 1.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_HEAT) {
    float datal[] = {0.0, 0.4, 0.6, 0.8, 1.0};
    float datar[] = {1.0, 1.0, 1.0, 0.5, 0.0};
    float datag[] = {1.0, 1.0, 0.5, 0.0, 0.0};
    float datab[] = {1.0, 0.3, 0.0, 0.0, 0.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_HEAT) {
    float datal[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    float datar[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    float datag[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    float datab[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
  }else if(maptype == PPGPLOT_COLD) {
    float datal[] = {0.0, 0.4, 0.6, 0.8, 1.0};
    float datar[] = {1.0, 0.3, 0.0, 0.0, 0.0};
    float datag[] = {1.0, 1.0, 0.5, 0.0, 0.0};
    float datab[] = {1.0, 1.0, 1.0, 0.5, 0.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_COLD) {
    float datal[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    float datar[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    float datag[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    float datab[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
  }else if(maptype == PPGPLOT_HEAT2) {
    float datal[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
    float datar[] = { 1.0, 1.0, 1.0, 1.0, 0.6, 0.0, 0.0, 0.0, 0.0};
    float datag[] = { 1.0, 0.0, 0.6, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
    float datab[] = { 1.0, 0.0, 0.0, 0.0, 0.3, 1.0, 0.8, 0.3, 0.0};
    ppgctab(datal, datar, datag, datab, 9, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_HEAT2) {
    float datal[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
    float datar[] = { 0.0, 0.0, 0.0, 0.0, 0.6, 1.0, 1.0, 1.0, 1.0};
    float datag[] = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.6, 0.0, 1.0};
    float datab[] = { 0.0, 0.3, 0.8, 1.0, 0.3, 0.0, 0.0, 0.0, 1.0};
    ppgctab(datal, datar, datag, datab, 9, contrast, brightness);
  }else if(maptype == PPGPLOT_HEAT3) {
    float datal[] = {0.0, 0.2, 0.35, 0.5, 0.8, 1.0};
    float datar[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0};
    float datag[] = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    float datab[] = {1.0, 0.0, 0.0, 1.0, 1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 6, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_HEAT3) {
    float datal[] = {0.0, 0.2, 0.35, 0.5, 0.8, 1.0};
    float datar[] = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
    float datag[] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0};
    float datab[] = {0.0, 1.0, 1.0, 0.0, 0.0, 1.0};
    ppgctab(datal, datar, datag, datab, 6, contrast, brightness);
  }else if(maptype == PPGPLOT_HEAT4) {
    float datal[] = {0.0, 0.33, 0.75, 1.0};
    float datar[] = {1.0, 0.0, 1.0, 1.0};
    float datag[] = {1.0, 0.0, 0.0, 0.0};
    float datab[] = {1.0, 1.0, 1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 4, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_HEAT4) {
    float datal[] = {0.0, 0.33, 0.75, 1.0};
    float datar[] = {1.0, 1.0, 0.0, 1.0};
    float datag[] = {0.0, 0.0, 0.0, 1.0};
    float datab[] = {0.0, 1.0, 1.0, 1.0};
    ppgctab(datal, datar, datag, datab, 4, contrast, brightness);
  }else if(maptype == PPGPLOT_DIVREDBLUE) {
    float datal[] = {0.0, 0.5, 1.0};
    float datar[] = {1.0, 1.0, 0.0};
    float datag[] = {0.0, 1.0, 0.0};
    float datab[] = {0.0, 1.0, 1.0};
    ppgctab(datal, datar, datag, datab, 3, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_DIVREDBLUE) {
    float datal[] = {0.0, 0.5, 1.0};
    float datar[] = {0.0, 1.0, 1.0};
    float datag[] = {0.0, 1.0, 0.0};
    float datab[] = {1.0, 1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 3, contrast, brightness);
  }else if(maptype == PPGPLOT_INFERNO) {
    float datal[] = {0.000000, 0.003922, 0.007843, 0.011765, 0.015686, 0.019608, 0.023529, 0.027451,
       0.031373, 0.035294, 0.039216, 0.043137, 0.047059, 0.050980, 0.054902, 0.058824,
       0.062745, 0.066667, 0.070588, 0.074510, 0.078431, 0.082353, 0.086275, 0.090196,
       0.094118, 0.098039, 0.101961, 0.105882, 0.109804, 0.113725, 0.117647, 0.121569,
       0.125490, 0.129412, 0.133333, 0.137255, 0.141176, 0.145098, 0.149020, 0.152941,
       0.156863, 0.160784, 0.164706, 0.168627, 0.172549, 0.176471, 0.180392, 0.184314,
       0.188235, 0.192157, 0.196078, 0.200000, 0.203922, 0.207843, 0.211765, 0.215686,
       0.219608, 0.223529, 0.227451, 0.231373, 0.235294, 0.239216, 0.243137, 0.247059,
       0.250980, 0.254902, 0.258824, 0.262745, 0.266667, 0.270588, 0.274510, 0.278431,
       0.282353, 0.286275, 0.290196, 0.294118, 0.298039, 0.301961, 0.305882, 0.309804,
       0.313725, 0.317647, 0.321569, 0.325490, 0.329412, 0.333333, 0.337255, 0.341176,
       0.345098, 0.349020, 0.352941, 0.356863, 0.360784, 0.364706, 0.368627, 0.372549,
       0.376471, 0.380392, 0.384314, 0.388235, 0.392157, 0.396078, 0.400000, 0.403922,
       0.407843, 0.411765, 0.415686, 0.419608, 0.423529, 0.427451, 0.431373, 0.435294,
       0.439216, 0.443137, 0.447059, 0.450980, 0.454902, 0.458824, 0.462745, 0.466667,
       0.470588, 0.474510, 0.478431, 0.482353, 0.486275, 0.490196, 0.494118, 0.498039,
       0.501961, 0.505882, 0.509804, 0.513725, 0.517647, 0.521569, 0.525490, 0.529412,
       0.533333, 0.537255, 0.541176, 0.545098, 0.549020, 0.552941, 0.556863, 0.560784,
       0.564706, 0.568627, 0.572549, 0.576471, 0.580392, 0.584314, 0.588235, 0.592157,
       0.596078, 0.600000, 0.603922, 0.607843, 0.611765, 0.615686, 0.619608, 0.623529,
       0.627451, 0.631373, 0.635294, 0.639216, 0.643137, 0.647059, 0.650980, 0.654902,
       0.658824, 0.662745, 0.666667, 0.670588, 0.674510, 0.678431, 0.682353, 0.686275,
       0.690196, 0.694118, 0.698039, 0.701961, 0.705882, 0.709804, 0.713725, 0.717647,
       0.721569, 0.725490, 0.729412, 0.733333, 0.737255, 0.741176, 0.745098, 0.749020,
       0.752941, 0.756863, 0.760784, 0.764706, 0.768627, 0.772549, 0.776471, 0.780392,
       0.784314, 0.788235, 0.792157, 0.796078, 0.800000, 0.803922, 0.807843, 0.811765,
       0.815686, 0.819608, 0.823529, 0.827451, 0.831373, 0.835294, 0.839216, 0.843137,
       0.847059, 0.850980, 0.854902, 0.858824, 0.862745, 0.866667, 0.870588, 0.874510,
       0.878431, 0.882353, 0.886275, 0.890196, 0.894118, 0.898039, 0.901961, 0.905882,
       0.909804, 0.913725, 0.917647, 0.921569, 0.925490, 0.929412, 0.933333, 0.937255,
       0.941176, 0.945098, 0.949020, 0.952941, 0.956863, 0.960784, 0.964706, 0.968627,
       0.972549, 0.976471, 0.980392, 0.984314, 0.988235, 0.992157, 0.996078, 1.000000};
    float datar[] = {
      0.988362,0.982257,0.976511,0.971162,0.966249,0.961812,0.957896,0.954529,
      0.951740,0.949545,0.947937,0.946903,0.946403,0.946392,0.946809,0.947594,
      0.948683,0.950018,0.951546,0.953215,0.954997,0.956834,0.958720,0.960626,
      0.962517,0.964394,0.966243,0.968041,0.969783,0.971468,0.973088,0.974638,
      0.976108,0.977497,0.978806,0.980032,0.981173,0.982228,0.983196,0.984075,
      0.984865,0.985566,0.986175,0.986694,0.987124,0.987464,0.987714,0.987874,
      0.987945,0.987926,0.987819,0.987622,0.987337,0.986964,0.986502,0.985952,
      0.985315,0.984591,0.983779,0.982881,0.981895,0.980824,0.979666,0.978422,
      0.977092,0.975677,0.974176,0.972590,0.970919,0.969163,0.967322,0.965397,
      0.963387,0.961293,0.959114,0.956852,0.954506,0.952075,0.949562,0.946965,
      0.944285,0.941521,0.938675,0.935747,0.932737,0.929644,0.926470,0.923215,
      0.919879,0.916462,0.912966,0.909390,0.905735,0.902003,0.898192,0.894305,
      0.890341,0.886302,0.882188,0.878001,0.873741,0.869409,0.865006,0.860533,
      0.855992,0.851384,0.846709,0.841969,0.837165,0.832299,0.827372,0.822386,
      0.817341,0.812239,0.807082,0.801871,0.796607,0.791293,0.785929,0.780517,
      0.775059,0.769556,0.764010,0.758422,0.752794,0.747127,0.741423,0.735683,
      0.729909,0.724103,0.718264,0.712396,0.706500,0.700576,0.694627,0.688653,
      0.682656,0.676638,0.670599,0.664540,0.658463,0.652369,0.646260,0.640135,
      0.633998,0.627847,0.621685,0.615513,0.609330,0.603139,0.596940,0.590734,
      0.584521,0.578304,0.572081,0.565854,0.559624,0.553392,0.547157,0.540920,
      0.534683,0.528444,0.522206,0.515967,0.509730,0.503493,0.497257,0.491022,
      0.484789,0.478558,0.472328,0.466100,0.459875,0.453651,0.447428,0.441207,
      0.434987,0.428768,0.422549,0.416331,0.410113,0.403894,0.397674,0.391453,
      0.385228,0.379001,0.372768,0.366529,0.360284,0.354032,0.347771,0.341500,
      0.335217,0.328921,0.322610,0.316282,0.309935,0.303568,0.297178,0.290763,
      0.284321,0.277850,0.271347,0.264810,0.258234,0.251620,0.244967,0.238273,
      0.231538,0.224763,0.217949,0.211095,0.204209,0.197297,0.190367,0.183429,
      0.176493,0.169575,0.162689,0.155850,0.149073,0.142378,0.135778,0.129285,
      0.122908,0.116656,0.110536,0.104551,0.098702,0.092990,0.087411,0.081962,
      0.076637,0.071429,0.066331,0.061340,0.056449,0.051644,0.046915,0.042253,
      0.037668,0.033385,0.029432,0.025793,0.022447,0.019373,0.016561,0.013995,
      0.011663,0.009561,0.007676,0.006006,0.004547,0.003299,0.002267,0.001462};
    float datag[] = {
      0.998364,0.994109,0.989753,0.985282,0.980678,0.975924,0.971003,0.965896,
      0.960587,0.955063,0.949318,0.943348,0.937159,0.930761,0.924168,0.917399,
      0.910473,0.903409,0.896226,0.888942,0.881569,0.874129,0.866624,0.859069,
      0.851476,0.843848,0.836191,0.828515,0.820825,0.813122,0.805409,0.797692,
      0.789974,0.782258,0.774545,0.766837,0.759135,0.751442,0.743758,0.736087,
      0.728427,0.720782,0.713153,0.705540,0.697944,0.690366,0.682807,0.675267,
      0.667748,0.660250,0.652773,0.645320,0.637890,0.630485,0.623105,0.615750,
      0.608422,0.601122,0.593849,0.586606,0.579392,0.572209,0.565057,0.557937,
      0.550850,0.543798,0.536780,0.529798,0.522853,0.515946,0.509078,0.502249,
      0.495462,0.488716,0.482014,0.475356,0.468744,0.462178,0.455660,0.449191,
      0.442772,0.436405,0.430091,0.423831,0.417627,0.411479,0.405389,0.399359,
      0.393389,0.387481,0.381636,0.375856,0.370140,0.364492,0.358911,0.353399,
      0.347957,0.342586,0.337287,0.332060,0.326906,0.321827,0.316822,0.311892,
      0.307038,0.302260,0.297559,0.292933,0.288385,0.283913,0.279517,0.275197,
      0.270954,0.266786,0.262692,0.258674,0.254728,0.250856,0.247056,0.243327,
      0.239667,0.236077,0.232554,0.229097,0.225706,0.222378,0.219112,0.215906,
      0.212759,0.209670,0.206636,0.203656,0.200728,0.197851,0.195021,0.192239,
      0.189501,0.186807,0.184153,0.181539,0.178962,0.176421,0.173914,0.171438,
      0.168992,0.166575,0.164184,0.161817,0.159474,0.157151,0.154848,0.152563,
      0.150294,0.148039,0.145797,0.143567,0.141346,0.139134,0.136929,0.134729,
      0.132534,0.130341,0.128150,0.125960,0.123769,0.121575,0.119379,0.117179,
      0.114974,0.112764,0.110547,0.108322,0.106089,0.103848,0.101597,0.099338,
      0.097069,0.094790,0.092501,0.090203,0.087896,0.085580,0.083257,0.080927,
      0.078591,0.076253,0.073915,0.071579,0.069247,0.066925,0.064616,0.062325,
      0.060060,0.057827,0.055634,0.053490,0.051407,0.049396,0.047470,0.045644,
      0.043933,0.042353,0.040922,0.039647,0.038571,0.037705,0.037055,0.036621,
      0.036405,0.036405,0.036615,0.037030,0.037632,0.038400,0.039309,0.040329,
      0.041402,0.042489,0.043554,0.044559,0.045468,0.046242,0.046856,0.047293,
      0.047536,0.047574,0.047399,0.047008,0.046402,0.045583,0.044556,0.043328,
      0.041905,0.040294,0.038504,0.036590,0.034569,0.032474,0.030324,0.028139,
      0.025921,0.023702,0.021503,0.019331,0.017199,0.015133,0.013136,0.011225,
      0.009417,0.007713,0.006136,0.004692,0.003392,0.002249,0.001270,0.000466};
    float datab[] = {
      0.644924,0.631017,0.616760,0.602154,0.587206,0.571925,0.556275,0.540361,
      0.524203,0.507860,0.491426,0.474970,0.458592,0.442367,0.426373,0.410665,
      0.395289,0.380271,0.365627,0.351369,0.337475,0.323974,0.310820,0.298010,
      0.285546,0.273391,0.261534,0.249972,0.238686,0.227658,0.216877,0.206332,
      0.196018,0.185923,0.176037,0.166353,0.156863,0.147565,0.138453,0.129527,
      0.120785,0.112229,0.103863,0.095694,0.087731,0.079990,0.072489,0.065257,
      0.058329,0.051750,0.045581,0.039886,0.034916,0.030908,0.027814,0.025592,
      0.024202,0.023606,0.023770,0.024661,0.026250,0.028508,0.031409,0.034931,
      0.039050,0.043618,0.048392,0.053324,0.058367,0.063488,0.068659,0.073859,
      0.079073,0.084289,0.089499,0.094695,0.099874,0.105031,0.110164,0.115272,
      0.120354,0.125409,0.130438,0.135440,0.140417,0.145367,0.150292,0.155193,
      0.160070,0.164924,0.169755,0.174563,0.179350,0.184116,0.188860,0.193584,
      0.198286,0.202968,0.207628,0.212268,0.216886,0.221482,0.226055,0.230606,
      0.235133,0.239636,0.244113,0.248564,0.252988,0.257383,0.261750,0.266085,
      0.270390,0.274661,0.278898,0.283099,0.287264,0.291390,0.295477,0.299523,
      0.303526,0.307485,0.311399,0.315266,0.319085,0.322856,0.326576,0.330245,
      0.333861,0.337424,0.340931,0.344383,0.347777,0.351113,0.354388,0.357603,
      0.360757,0.363849,0.366879,0.369846,0.372748,0.375586,0.378359,0.381065,
      0.383704,0.386276,0.388781,0.391219,0.393589,0.395891,0.398125,0.400290,
      0.402385,0.404411,0.406369,0.408258,0.410078,0.411829,0.413511,0.415123,
      0.416667,0.418142,0.419549,0.420887,0.422156,0.423356,0.424488,0.425552,
      0.426548,0.427475,0.428334,0.429125,0.429846,0.430498,0.431080,0.431594,
      0.432039,0.432412,0.432714,0.432943,0.433098,0.433179,0.433183,0.433109,
      0.432955,0.432719,0.432400,0.431994,0.431497,0.430906,0.430217,0.429425,
      0.428524,0.427511,0.426377,0.425116,0.423721,0.422182,0.420491,0.418637,
      0.416608,0.414392,0.411976,0.409345,0.406485,0.403378,0.400007,0.396353,
      0.392400,0.388129,0.383522,0.378563,0.373238,0.367535,0.361447,0.354971,
      0.348111,0.340874,0.333277,0.325338,0.317085,0.308553,0.299776,0.290788,
      0.281624,0.272321,0.262912,0.253430,0.243904,0.234358,0.224813,0.215289,
      0.205799,0.196354,0.186962,0.177642,0.168414,0.159254,0.150164,0.141141,
      0.132232,0.123397,0.114621,0.105930,0.097327,0.088767,0.080282,0.071862,
      0.063460,0.055143,0.046836,0.038558,0.030909,0.024239,0.018570,0.013866};
    ppgctab(datal, datar, datag, datab, 256, contrast, brightness);
 }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplot_set_maptype: Maptype undefined");
    return 0;
  }
  return 1;
}
void clear_pgplot_frame(pgplot_frame_def_internal *frame)
{
  frame->svp = 0;
  frame->svp_x1 = 0;
  frame->svp_x2 = 0;
  frame->svp_y1 = 0;
  frame->svp_y2 = 0;
  frame->swin = 0;
  frame->swin_showtwice = 0;
  frame->swin_x1 = 0;
  frame->swin_x2 = 0;
  frame->swin_y1 = 0;
  frame->swin_y2 = 0;
  frame->TR[0] = 0;
  frame->TR[1] = 0;
  frame->TR[2] = 0;
  frame->TR[3] = 0;
  frame->TR[4] = 0;
  frame->TR[5] = 0;
}
void clear_pgplot_box(pgplot_box_definition *box)
{
  box->drawbox = 1;
  box->box_lw = 1;
  box->box_f = 1;
  box->box_labelsize = 1;
  box->box_xtick = 0;
  box->box_ytick = 0;
  box->box_nxsub= 0;
  box->box_nysub= 0;
  sprintf(box->box_xopt, "bcnsti");
  sprintf(box->box_yopt, "bcnsti");
  box->drawtitle = 1;
  box->title_ch = 1;
  box->title_lw = 1;
  box->title_f = 1;
  box->title[0] = 0;
  box->drawlabels = 1;
  box->label_ch = 1;
  box->label_lw = 1;
  box->label_f = 1;
  box->dxlabel = 0;
  box->dylabel = 0;
  box->xlabel[0] = 0;
  box->ylabel[0] = 0;
  box->wedgelabel[0] = 0;
}
void pgplot_drawbox(pgplot_box_definition *box)
{
  if(box->drawbox) {
    ppgslw(box->box_lw);
    ppgsch(box->box_labelsize);
    ppgscf(box->box_f);
    ppgbox(box->box_xopt,box->box_xtick,box->box_nxsub,box->box_yopt,box->box_ytick,box->box_nysub);
  }
  if(box->drawtitle) {
    ppgsch(box->title_ch);
    ppgslw(box->title_lw);
    ppgscf(box->title_f);
    ppgmtxt("T", (1.5*box->box_labelsize)/box->title_ch, 0.5, 0.5, box->title);
  }
  if(box->drawlabels) {
    ppgslw(box->label_lw);
    ppgsch(box->label_ch);
    ppgscf(box->label_f);
    ppgmtxt("L", ((2.5+box->dylabel)*box->box_labelsize)/box->label_ch, 0.5, 0.5, box->ylabel);
    ppgmtxt("B", ((1.5+box->dxlabel)*box->box_labelsize)/box->label_ch+1, 0.5, 0.5, box->xlabel);
  }
  ppgsch(1);
  ppgslw(1);
  ppgscf(1);
}
int pgplotMap(pgplot_options_definition *pgplot, float *cmap, int nrx, int nry, float xmin, float xmax, float xminshow, float xmaxshow, float ymin, float ymax, float yminshow, float ymaxshow, int maptype, int itf, int nogray, int nrcontours, float *contours, int contourlw, int forceMinZero, float saturize, int levelset, float levelmin, float levelmax, int levelInversion, int onlyData, int sideright, int forceMinZeroRight, int sidetop, int forceMinZeroTop, int sidelw, int showwedge, int showwedge_actualmax, int plotSubset, int showTwice, verbose_definition verbose)
{
  float datamin, datamax, xmin2, xmax2;
  float *collapse, x, y, remember2_x1, remember2_x2, remember2_y1, remember2_y2;
  int i, j, deviceID, firstc, lastc, firstpoint;
  float junk_f, lasty;
  pgplot_frame_def_internal pgplot_frame_internal;
  pgplot_options_definition *pgplot_backup;
  pgplot_backup = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
  if(pgplot_backup == NULL) {
    printerror(verbose.debug, "ERROR pgplotMap: Memory allocation error");
    return 0;
  }
  memcpy(pgplot_backup, pgplot, sizeof(pgplot_options_definition));
  clear_pgplot_frame(&pgplot_frame_internal);
  ppgplot_set_internal_mapping_coordinated(xmin, xmax, ymin, ymax, nrx, nry);
  lasty = 0;
  if(verbose.debug) {
    printf("pgplotMap -- dimensions: %dX%d points  xrange=%e %e yrange=%e %e\n", nrx, nry, xmin, xmax, ymin, ymax);
    printf("pgplotMap -- Horizontal dimensions %f - %f (shown %f - %f)\n", xmin, xmax, xminshow, xmaxshow);
    printf("pgplotMap -- Vertical dimensions %f - %f (shown %f - %f)\n", ymin, ymax, yminshow, ymaxshow);
  }
  if(pgplot_opendevice(&(pgplot->viewport), &deviceID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotMap: Cannot open plot device");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  if(onlyData == 0) {
    pgplot_setWindowsize(pgplot->viewport.windowwidth, pgplot->viewport.windowheight, pgplot->viewport.aspectratio);
  }
  ppgslw(1);
  ppgsci(1);
  ppgscf(1);
  ppgsch(1);
  if(onlyData == 0 || onlyData == 2) {
    pgplot_frame_internal.svp_x1 = pgplot->viewport.dxplot+0.15;
    pgplot_frame_internal.svp_y1 = pgplot->viewport.dyplot+0.15;
    if(sideright)
      pgplot_frame_internal.svp_x2 = (0.75-0.15)*pgplot->viewport.xsize+pgplot_frame_internal.svp_x1;
    else
      pgplot_frame_internal.svp_x2 = (0.90-0.15)*pgplot->viewport.xsize+pgplot_frame_internal.svp_x1;
    if(sidetop)
      pgplot_frame_internal.svp_y2 = (0.75-0.15)*pgplot->viewport.ysize+pgplot_frame_internal.svp_y1;
    else
      pgplot_frame_internal.svp_y2 = (0.90-0.15)*pgplot->viewport.ysize+pgplot_frame_internal.svp_y1;
    pgplot_frame_internal.svp = 1;
  }else {
    pgplot_frame_internal.svp = 0;
  }
  pgplot_frame_internal.swin = 1;
  pgplot_frame_internal.swin_showtwice = showTwice;
  pgplot_frame_internal.swin_x1 = xminshow;
  pgplot_frame_internal.swin_x2 = xmaxshow;
  pgplot_frame_internal.swin_y1 = yminshow;
  pgplot_frame_internal.swin_y2 = ymaxshow;
  if(showTwice)
    pgplot_frame_internal.swin_y2 += ymaxshow-yminshow;
  pgplot_makeframe(&pgplot_frame_internal);
  firstpoint = 1;
  datamin = datamax = 0;
  for(i = 0; i < nrx; i++) {
    for(j = 0; j <nry; j++) {
      pgplotMapCoordinateInverse(&x, &y, i, j);
      if(x >= xminshow && x <= xmaxshow && y >= yminshow && y <= ymaxshow) {
 if(firstpoint) {
   datamin = datamax = cmap[j*nrx+i];
   firstpoint = 0;
 }
 if(cmap[j*nrx+i] > datamax)
   datamax = cmap[j*nrx+i];
 if(cmap[j*nrx+i] < datamin) {
   datamin = cmap[j*nrx+i];
 }
      }
    }
  }
  if(verbose.debug) printf("pgplotMap -- Data range: %e to %e\n", datamin, datamax);
  if(forceMinZero)
    datamin = 0;
  double datamax_actual;
  datamax_actual = datamax;
  datamax = datamin + (datamax-datamin)/saturize;
  if(levelset) {
    datamax = levelmax;
    datamin = levelmin;
  }
  if(levelInversion) {
    junk_f = datamax;
    datamax = datamin;
    datamin = junk_f;
  }
  if(verbose.debug) printf("pgplotMap -- Color range: %e to %e\n", datamin, datamax);
  firstc = 17;
  lastc = 1000;
  ppgscir(firstc, lastc);
  cpgqcir(&firstc, &lastc);
  if(verbose.debug) printf("pgplotMap -- Color index range:  %d to %d\n", firstc, lastc);
  ppgsitf(itf);
  if(pgplot_set_maptype(maptype, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotMap: Cannot set colour map");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  int subset_nrx, subset_nry, subset_extrabefore, subset_extraafter;
  float *cmap_subset;
  if(plotSubset) {
    int subset_x0, subset_y0;
    subset_extrabefore = 1;
    subset_extraafter = 1;
    subset_x0 = (nrx-1)*(xminshow - xmin)/(xmax-xmin)-1;
    if(subset_x0 < 0)
      subset_x0 = 0;
    if(nrx != 1) {
      subset_nrx = (nrx-1)*(xmaxshow - xmin)/(xmax-xmin)+1 - subset_x0;
    }else {
      subset_nrx = 1 - subset_x0;
    }
    if(subset_nrx > nrx)
      subset_nrx = nrx;
    subset_y0 = (nry-1)*(yminshow - ymin)/(ymax-ymin)-1;
    if(subset_y0 < 0) {
      subset_y0 = 0;
      subset_extrabefore = 0;
    }
    if(nry != 1)
      subset_nry = (nry-1)*(ymaxshow - ymin)/(ymax-ymin)+1 - subset_y0;
    else
      subset_nry = +1 - subset_y0;
    if(subset_nry >= nry) {
      subset_nry = nry;
      subset_extraafter = 0;
    }
    if(subset_nrx == nrx && subset_nry == nry && subset_y0 == 0 && subset_x0 == 0)
      plotSubset = 0;
    if(verbose.debug && plotSubset) printf("pgplotMap -- changed to: %dX%d points  startx=%d starty=%d\n", subset_nrx, subset_nry, subset_x0, subset_y0);
    if(plotSubset) {
      cmap_subset = (float *)malloc(subset_nrx*subset_nry*sizeof(float));
      if(cmap_subset == NULL) {
 fflush(stdout);
 printwarning(verbose.debug, "WARNING pgplotMap: Cannot allocate memory to plot subset of the map");
 if(verbose.debug) {
   printwarning(verbose.debug, "Tried to allocate %ld x %ld floating points", subset_nrx, subset_nry);
 }
 plotSubset = 0;
      }
    }
    if(plotSubset) {
      for(i = 0; i < subset_nrx; i++) {
 for(j = 0; j < subset_nry; j++) {
   cmap_subset[j*subset_nrx+i] = cmap[(j+subset_y0)*nrx+i+subset_x0];
 }
      }
      pgplot_frame_internal.TR[0] += pgplot_frame_internal.TR[1]*subset_x0;
      pgplot_frame_internal.TR[3] += pgplot_frame_internal.TR[5]*subset_y0;
    }
  }
  if(plotSubset == 0) {
    subset_nrx = nrx;
    subset_nry = nry;
    subset_extrabefore = 0;
    subset_extraafter = 0;
    cmap_subset = cmap;
  }
  if(!nogray) {
    if(showTwice) {
      pgplot_frame_internal.TR[3] += pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
      ppgimag(cmap_subset, subset_nrx, subset_nry, 1, subset_nrx, 1, subset_nry, datamax, datamin, pgplot_frame_internal.TR);
      pgplot_frame_internal.TR[3] -= pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
    }
    ppgimag(cmap_subset, subset_nrx, subset_nry, 1, subset_nrx, 1, subset_nry, datamax, datamin, pgplot_frame_internal.TR);
  }
  if(showwedge && sideright == 0) {
    int free_txt;
    char *txt;
    free_txt = 0;
    if(showwedge_actualmax == 0) {
      txt = pgplot->box.wedgelabel;
    }else {
      txt = malloc(1000);
      if(txt == NULL) {
 printerror(verbose.debug, "ERROR pgplotMap: Cannot allocate memory");
 return 0;
      }
      free_txt = 1;
      if(fabs(datamax_actual) >= 1.0) {
 sprintf(txt, "Max. = %.2f", datamax_actual);
      }else if(fabs(datamax_actual) >= 0.1) {
 sprintf(txt, "Max. = %.3f", datamax_actual);
      }else if(fabs(datamax_actual) >= 0.01) {
 sprintf(txt, "Max. = %.4f", datamax_actual);
      }else {
 sprintf(txt, "Max. = %.2e", datamax_actual);
      }
    }
    ppgsch(pgplot->box.box_labelsize*0.5);
    ppgslw(ceil(0.5*pgplot->box.box_lw));
    ppgwedg("RI", 0, 6, datamax, datamin, txt);
    ppgslw(1);
    if(free_txt) {
      free(txt);
    }
  }
  if(nrcontours > 0) {
    int didallocate_contours = 0;
    if(contours == NULL) {
      contours = malloc(nrcontours*sizeof(float));
      if(contours == NULL) {
 printerror(verbose.debug, "ERROR pgplotMap: Cannot allocate memory");
 return 0;
      }
      didallocate_contours = 1;
      int contlevel;
      for(contlevel = 0; contlevel < nrcontours; contlevel++) {
 contours[contlevel] = datamin + (contlevel+1)*(datamax-datamin)/(float)(nrcontours+1);
      }
    }
    ppgslw(contourlw);
    if(showTwice) {
      pgplot_frame_internal.TR[3] += pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
      ppgcont(cmap_subset, subset_nrx, subset_nry, 1, subset_nrx, 1, subset_nry, contours, nrcontours, pgplot_frame_internal.TR);
      pgplot_frame_internal.TR[3] -= pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
    }
    ppgcont(cmap_subset, subset_nrx, subset_nry, 1, subset_nrx, 1, subset_nry, contours, nrcontours, pgplot_frame_internal.TR);
    if(didallocate_contours) {
      free(contours);
      contours = NULL;
    }
  }
  if(plotSubset)
    free(cmap_subset);
  if(onlyData != 0)
    pgplot->box.drawbox = 0;
  if(sidetop && onlyData != 1)
    pgplot->box.drawtitle = 0;
  if(onlyData != 0)
    pgplot->box.drawlabels = 0;
  pgplot_drawbox(&(pgplot->box));
  if(sideright) {
    collapse = (float *)calloc(nry, sizeof(float));
    if(collapse == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplotMap: Cannot allocate memory");
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    for(i = 0; i < nrx; i++) {
      x = i*(xmax-xmin)/(float)(nrx-1) + xmin;
      if(x >= xminshow && x <= xmaxshow) {
 for(j = 0; j <nry; j++) {
   collapse[j] += cmap[j*nrx+i];
 }
      }
    }
    xmin2 = xmax2 = collapse[0];
    i = 0;
    for(j = 0; j <nry; j++) {
      y = j*(ymax-ymin)/(float)(nry-1) + ymin;
      if(y >= yminshow && y <= ymaxshow) {
 if(i == 0)
   xmin2 = xmax2 = collapse[j];
 if(collapse[j] > xmax2)
   xmax2 = collapse[j];
 if(collapse[j] < xmin2)
   xmin2 = collapse[j];
 i++;
      }
    }
    if(xmin2 == xmax2) {
      float offset;
      offset = 0.05*xmax2;
      xmin2 -= offset;
      xmax2 += offset;
      if(xmin2 == xmax2) {
 xmin2 -= 0.5;
 xmax2 += 0.5;
      }
    }
    if(forceMinZeroRight)
      xmin2 = 0;
    if(sidetop)
      junk_f = (0.75-0.15)*pgplot->viewport.ysize+pgplot->viewport.dyplot+0.15;
    else
      junk_f = (0.90-0.15)*pgplot->viewport.ysize+pgplot->viewport.dyplot+0.15;
    if(onlyData != 0) {
      ppgqvp(0, &pgplot_frame_internal.svp_x1, &pgplot_frame_internal.svp_x2, &pgplot_frame_internal.svp_y1, &pgplot_frame_internal.svp_y2);
      ppgqwin(&remember2_x1, &remember2_x2, &remember2_y1, &remember2_y2);
    }
    ppgsvp(pgplot_frame_internal.svp_x2, pgplot_frame_internal.svp_x2+0.15, pgplot_frame_internal.svp_y1, junk_f);
    if(showTwice == 0)
      ppgswin(xmin2, xmax2, yminshow, ymaxshow);
    else
      ppgswin(xmin2, xmax2, yminshow, ymaxshow+ymaxshow-yminshow);
    if(showwedge) {
      ppgsch(pgplot->box.box_labelsize*0.5);
      ppgslw(ceil(0.5*pgplot->box.box_lw));
      ppgwedg("RI", 0, 6, datamax, datamin, pgplot->box.wedgelabel);
      ppgslw(1);
    }
    if(onlyData == 0 || onlyData == 2) {
      ppgslw(pgplot->box.box_lw);
      ppgsch(pgplot->box.box_labelsize*0.5);
      if(sideright == 4) {
 ppgbox("bcsti",0.0,0,"bc",0.0,0);
      }else if(sideright == 3) {
 ppgbox("bcnsti",0.0,0,"bc",0.0,0);
      }else if(sideright == 2) {
 ppgbox("bcsti",0.0,0,"bctsi",0.0,0);
      }else {
 ppgbox("bcnsti",0.0,0,"bctsi",0.0,0);
      }
      ppgsch(1);
      ppgslw(1);
    }
    ppgslw(sidelw);
    for(j = 0; j < nry; j++) {
      x = collapse[j];
      y = j*(ymax-ymin)/(float)(nry-1) + ymin;
      if(y >= yminshow && y <= ymaxshow) {
 if(j == 0)
   ppgmove(x, y);
 else {
   ppgdraw(x, y);
   lasty = y;
 }
      }else {
   ppgmove(x, y);
      }
    }
    if(showTwice) {
      for(j = 0; j < nry; j++) {
 x = collapse[j];
 if(plotSubset == 0)
   y = j*(ymax-ymin)/(float)(nry-1) + ymin + nry*(ymax-ymin)/(float)(nry-1);
 else
   y = j*(ymax-ymin)/(float)(nry-1) + ymin + pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
 if(y <= ymaxshow+ymaxshow-yminshow) {
   if(y > lasty)
     ppgdraw(x, y);
 }
 ppgsci(1);
      }
    }
    ppgslw(1);
    free(collapse);
    if(onlyData != 0) {
      ppgsvp(pgplot_frame_internal.svp_x1, pgplot_frame_internal.svp_x2, pgplot_frame_internal.svp_y1, pgplot_frame_internal.svp_y2);
      ppgswin(remember2_x1, remember2_x2, remember2_y1, remember2_y2);
    }
  }
  if(sidetop) {
    collapse = (float *)calloc(nrx, sizeof(float));
    if(collapse == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplotMap: Cannot allocate memory");
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    for(i = 0; i < nry; i++) {
      if(nry == 1)
 y = ymin;
      else
 y = i*(ymax-ymin)/(float)(nry-1) + ymin;
      if(y >= yminshow && y <= ymaxshow) {
 for(j = 0; j <nrx; j++) {
   collapse[j] += cmap[i*nrx+j];
 }
      }
    }
    i = 0;
    xmin2 = xmax2 = 0;
    for(j = 0; j <nrx; j++) {
      x = (j*(xmax-xmin)/(float)(nrx - 1)) + xmin;
      if(x >= xminshow && x <= xmaxshow) {
 if(i == 0) {
   xmin2 = xmax2 = collapse[j];
   i = 1;
 }
 if(collapse[j] > xmax2)
   xmax2 = collapse[j];
 if(collapse[j] < xmin2)
   xmin2 = collapse[j];
      }
    }
    if(xmin2 == xmax2) {
      float offset;
      offset = 0.05*xmax2;
      xmin2 -= offset;
      xmax2 += offset;
      if(xmin2 == xmax2) {
 xmin2 -= 0.5;
 xmax2 += 0.5;
      }
    }
    if(forceMinZeroRight)
      xmin2 = 0;
    if(sideright)
      junk_f = (0.75-0.15)*pgplot->viewport.xsize+pgplot->viewport.dxplot+0.15;
    else
      junk_f = (0.90-0.15)*pgplot->viewport.xsize+pgplot->viewport.dxplot+0.15;
    if(onlyData != 0) {
      ppgqvp(0, &pgplot_frame_internal.svp_x1, &pgplot_frame_internal.svp_x2, &pgplot_frame_internal.svp_y1, &pgplot_frame_internal.svp_y2);
      ppgqwin(&remember2_x1, &remember2_x2, &remember2_y1, &remember2_y2);
      junk_f = pgplot_frame_internal.svp_x2;
    }
    ppgsvp(pgplot_frame_internal.svp_x1, junk_f, pgplot_frame_internal.svp_y2, pgplot_frame_internal.svp_y2+0.15);
    ppgswin(xminshow, xmaxshow, xmin2, xmax2);
    if(onlyData == 0 || onlyData == 2) {
      ppgsch(pgplot->box.box_labelsize*0.5);
      ppgslw(pgplot->box.box_lw);
      if(sidetop == 2) {
 ppgbox("bcsti",0.0,0,"bctsi",0.0,0);
      }else {
 ppgbox("bcsti",0.0,0,"bcntsi",0.0,0);
      }
      ppgsch(1);
      ppgslw(1);
    }
    ppgslw(sidelw);
    for(j = 0; j < nrx; j++) {
      y = collapse[j];
      x = j*(xmax-xmin)/(float)(nrx-1) + xmin;
      if(j == 0)
 ppgmove(x, y);
      else
 ppgdraw(x, y);
    }
    ppgslw(1);
    if(onlyData == 0 || onlyData == 2) {
      ppgsch(pgplot->box.title_ch);
      ppgslw(pgplot->box.title_lw);
      ppgscf(pgplot->box.title_f);
      ppgptxt(0.5*(xmaxshow-xminshow)+xminshow, xmax2+0.3*(pgplot->box.label_ch/1.7)*(xmax2-xmin2), 0, 0.5, pgplot->box.title);
      ppgscf(1);
      ppgslw(1);
      ppgsch(pgplot->box.label_ch);
    }
    free(collapse);
    if(onlyData != 0) {
      ppgsvp(pgplot_frame_internal.svp_x1, pgplot_frame_internal.svp_x2, pgplot_frame_internal.svp_y1, pgplot_frame_internal.svp_y2);
      ppgswin(remember2_x1, remember2_x2, remember2_y1, remember2_y2);
    }
  }
  if(onlyData == 0 || onlyData == 2) {
    ppgsvp(pgplot_frame_internal.svp_x1, pgplot_frame_internal.svp_x2, pgplot_frame_internal.svp_y1, pgplot_frame_internal.svp_y2);
    if(showTwice == 0)
      ppgswin(xminshow, xmaxshow, yminshow, ymaxshow);
    else
      ppgswin(xminshow, xmaxshow, yminshow, ymaxshow+ymaxshow-yminshow);
  }
  memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
  free(pgplot_backup);
  if(!pgplot->viewport.dontclose)
    ppgclos();
  return 1;
}
int selectRegions(float *profileI, int nrBins, pgplot_options_definition *pgplot, int onlyOne, int powerTwo, int evenNumber, pulselongitude_regions_definition *regions, verbose_definition verbose)
{
  int i, j, k, bin1, bin2, replot;
  float x, y;
  char c;
  pgplot_options_definition *pgplot2;
  pgplot2 = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
  if(pgplot2 == NULL) {
    printerror(verbose.debug, "ERROR selectRegions: Memory allocation error");
    return 0;
  }
  x = y = 0;
  c = 0;
  bin1 = bin2 = 0;
  if(strcmp(pgplot->viewport.plotDevice, "?") == 0)
    printf("Select plotting device to show profile in order to select the onpulse region\n  ");
  memcpy(pgplot2, pgplot, sizeof(pgplot_options_definition));
  pgplot2->viewport.dontclose = 1;
  pgplotGraph1(pgplot2, profileI, NULL, NULL, nrBins, 0, nrBins-1, 0, 0, nrBins-1, 0, 0, 0, 1, 0, 1, 0, 1, 1, regions, -1, verbose);
  printf("Select onpulse region (click left and right of region, you can select multiple regions and pressing 'r' inside the PGPLOT window removes selections one by one).\nWhen finished press S inside the PGPLOT window.\n");
  i = 0;
  replot = 0;
  do {
    if(i == 0)
      cpgband(6, 0, 0.0, 0.0, &x, &y, &c);
    else
      cpgband(4, 0, bin1-0.5, 0.0, &x, &y, &c);
    if(c == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING selectRegions: PGPLOT device does not support cursor: no pulse phase selection can be made.");
      if(pgplot->viewport.dontclose == 0) {
 ppgclos();
 fflush(stdout);
      }
      free(pgplot2);
      return -1;
    }else if(c == 65) {
      if(i == 0)
 bin1 = (int)(x+0.5);
      else
 bin2 = (int)(x+0.5);
      i++;
    }else if(c == 115 || c == 83 ) {
      i = 3;
      printf("S is pressed\n");
    }else if(c == 'r' || c == 'R') {
      printf("R is pressed\n");
      i = 0;
      regions->nrRegions--;
      if(regions->nrRegions < 0)
 regions->nrRegions = 0;
      replot = 1;
    }
    if(i == 2) {
      if(bin2 < bin1) {
 i = bin1;
 bin1 = bin2;
 bin2 = i;
      }
      if(bin1 < 0)
 bin1 = 0;
      if(bin2 >= nrBins)
 bin2 = nrBins-1;
      if(verbose.verbose) printf("Selected region: %d - %d\n", bin1, bin2);
      if(powerTwo == 1) {
   y = log(bin2-bin1+1.0)/log(2.0);
 j = y;
 k = pow(2.0, j+1.0);
 bin1 -= 0.5*(k - (bin2-bin1+1.0));
 if(bin1 < 0)
   bin1 = 0;
 bin2 = bin1+k-1;
 if(verbose.verbose) printf("Selected region: %d - %d\n", bin1, bin2);
      }else if(evenNumber == 1) {
 j = bin2-bin1+1;
 j /= 2;
 j = 2*j-(bin2-bin1+1);
 if(j != 0) {
   bin2 += 1;
 }
 if(bin2 >= nrBins) {
   bin2 = nrBins-1;
   bin1 -= 1;
 }
 if(bin1 < 0)
   bin1 = 0;
 if(verbose.verbose) printf("Selected region: %d - %d\n", bin1, bin2);
      }
      if(onlyOne)
 regions->nrRegions = 0;
      regions->left_bin[regions->nrRegions] = bin1;
      regions->right_bin[regions->nrRegions] = bin2;
      regions->bins_defined[regions->nrRegions] = 1;
      regions->nrRegions += 1;
      if(regions->nrRegions == MAX_pulselongitude_regions) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR selectRegions: Too many regions selected.");
 ppgend();
 free(pgplot2);
 return regions->nrRegions;
      }
      replot = 1;
      i = 0;
    }
    if(replot) {
      memcpy(pgplot2, pgplot, sizeof(pgplot_options_definition));
      pgplot2->viewport.dontclose = 1;
      pgplot2->viewport.dontopen = 1;
      pgplotGraph1(pgplot2, profileI, NULL, NULL, nrBins, 0, nrBins-1, 0, 0, nrBins-1, 0, 0, 0, 1, 0, 1, 0, 1, 1, regions, -1, verbose);
      replot = 0;
    }
  }while(i < 3);
  if(pgplot->viewport.dontclose == 0)
    ppgclos();
  free(pgplot2);
  return regions->nrRegions;
}
int checkRegions(int bin, pulselongitude_regions_definition *regions, int whichregion, verbose_definition verbose)
{
  int i;
  if(regions == NULL)
    return 0;
  for(i = 0; i < regions->nrRegions; i++) {
    if(regions->bins_defined[i] == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING: checkRegions called without region defined in bins, call is ignored");
      return 0;
    }
    if(bin >= regions->left_bin[i] && bin <= regions->right_bin[i]) {
      if(i+1 == whichregion || whichregion == 0)
  return i+1;
    }
  }
  return 0;
}
int initPulselongitudeRegion(pulselongitude_regions_definition *region, verbose_definition verbose)
{
  region->bins_defined = malloc(MAX_pulselongitude_regions*sizeof(int));
  region->left_bin = malloc(MAX_pulselongitude_regions*sizeof(int));
  region->right_bin = malloc(MAX_pulselongitude_regions*sizeof(int));
  region->frac_defined = malloc(MAX_pulselongitude_regions*sizeof(int));
  region->left_frac = malloc(MAX_pulselongitude_regions*sizeof(float));
  region->right_frac = malloc(MAX_pulselongitude_regions*sizeof(float));
  if(region->bins_defined == NULL || region->left_bin == NULL || region->right_bin == NULL || region->frac_defined == NULL || region->left_frac == NULL || region->right_frac == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR initPulselongitudeRegion: Memory allocation error");
    return 0;
  }
  clearPulselongitudeRegion(region);
  return 1;
}
void freePulselongitudeRegion(pulselongitude_regions_definition *region)
{
  free(region->bins_defined);
  free(region->left_bin);
  free(region->right_bin);
  free(region->frac_defined);
  free(region->left_frac);
  free(region->right_frac);
}
void copyPulselongitudeRegion(pulselongitude_regions_definition source, pulselongitude_regions_definition *destination)
{
  int i;
  destination->nrRegions = source.nrRegions;
  if(source.nrRegions > 0) {
    for(i = 0; i < source.nrRegions; i++) {
      destination->bins_defined[i] = source.bins_defined[i];
      destination->left_bin[i] = source.left_bin[i];
      destination->right_bin[i] = source.right_bin[i];
      destination->frac_defined[i] = source.frac_defined[i];
      destination->left_frac[i] = source.left_frac[i];
      destination->right_frac[i] = source.right_frac[i];
    }
  }
}
void clearPulselongitudeRegion(pulselongitude_regions_definition *region)
{
  int i;
  region->nrRegions = 0;
  for(i = 0; i < MAX_pulselongitude_regions; i++) {
    region->frac_defined[i] = 0;
    region->bins_defined[i] = 0;
  }
}
void region_frac_to_int(pulselongitude_regions_definition *region, float scale, float offset)
{
  int i;
  for(i = 0; i < region->nrRegions; i++) {
    if(region->frac_defined[i]) {
      region->left_bin[i] = region->left_frac[i]*scale + offset;
      region->right_bin[i] = region->right_frac[i]*scale + offset;
      region->bins_defined[i] = 1;
    }
  }
}
void region_int_to_frac(pulselongitude_regions_definition *region, float scale, float offset)
{
  int i;
  for(i = 0; i < region->nrRegions; i++) {
    if(region->bins_defined[i]) {
      region->frac_defined[i] = 1;
      region->left_frac[i] = region->left_bin[i]*scale + offset;
      region->right_frac[i] = region->right_bin[i]*scale + offset;
    }
  }
}
int region_make_even(int regionnr, pulselongitude_regions_definition *region, int nrx)
{
  long nrx2, junk_int;
  nrx2 = region->right_bin[regionnr]-region->left_bin[regionnr] + 1;
  junk_int = nrx2 / 2;
  junk_int *= 2;
  junk_int -= nrx2;
  if(junk_int == 0) {
    return 1;
  }else {
    if(region->right_bin[regionnr] < nrx - 1) {
      region->right_bin[regionnr]++;
      return 2;
    }else if(region->left_bin[regionnr] > 0) {
      region->left_bin[regionnr]--;
      return 2;
    }
  }
  return 1;
}
void regionShowNextTimeUse(pulselongitude_regions_definition region, char *option, char *optionFrac, FILE *where)
{
  int i, ok;
  if(region.nrRegions > 0) {
    fprintf(where, "If you repeat this command, you could specify on the command line: ");
    for(i = 0; i < region.nrRegions; i++) {
      if(region.bins_defined[i])
 fprintf(where, "%s '%d %d' ", option, region.left_bin[i], region.right_bin[i]);
    }
    if(optionFrac != NULL) {
      ok = 0;
      for(i = 0; i < region.nrRegions; i++) {
 if(region.frac_defined[i])
   ok = 1;
      }
      if(ok) {
 fprintf(where, "\n  alternatively you can use the command line:    ");
 for(i = 0; i < region.nrRegions; i++) {
   if(region.frac_defined[i])
     fprintf(where, "%s '%f %f' ", optionFrac, region.left_frac[i], region.right_frac[i]);
   else if(region.bins_defined[i])
     fprintf(where, "%s '%d %d' ", option, region.left_bin[i], region.right_bin[i]);
 }
      }
    }
    fprintf(where, "\n");
  }
}
void printCMAPCommandlineOptions(FILE *printdevice)
{
  fprintf(printdevice, "Valid options for the -cmap option are:\n");
  fprintf(printdevice, "  GRAYSCALE\n");
  fprintf(printdevice, "  INVERTED_GRAYSCALE\n");
  fprintf(printdevice, "  RED\n");
  fprintf(printdevice, "  INVERTED_RED\n");
  fprintf(printdevice, "  GREEN\n");
  fprintf(printdevice, "  INVERTED_GREEN\n");
  fprintf(printdevice, "  BLUE\n");
  fprintf(printdevice, "  INVERTED_BLUE\n");
  fprintf(printdevice, "  CYAN\n");
  fprintf(printdevice, "  INVERTED_CYAN\n");
  fprintf(printdevice, "  HEAT\n");
  fprintf(printdevice, "  INVERTED_HEAT\n");
  fprintf(printdevice, "  HEAT2\n");
  fprintf(printdevice, "  INVERTED_HEAT2\n");
  fprintf(printdevice, "  HEAT3\n");
  fprintf(printdevice, "  INVERTED_HEAT3\n");
  fprintf(printdevice, "  HEAT4\n");
  fprintf(printdevice, "  INVERTED_HEAT4\n");
  fprintf(printdevice, "  COLD\n");
  fprintf(printdevice, "  INVERTED_COLD\n");
  fprintf(printdevice, "  DIVREDBLUE\n");
  fprintf(printdevice, "  INVERTED_DIVREDBLUE (a.k.a the Dutch flag)\n");
  fprintf(printdevice, "  INFERNO\n");
}
int cmap_parse_commandline(int argc, char **argv, int debug)
{
  int i, j, type;
  char identifier[100];
  type = 0;
  for(i = 1; i < argc - 1; i++) {
    if(strcmp(argv[i], "-cmap") == 0) {
      i++;
      j = sscanf(argv[i], "%s", identifier);
      if(j != 1) {
 fflush(stdout);
 printerror(debug, "ERROR cmap_parse_commandline:  Cannot parse -cmap option.");
 printCMAPCommandlineOptions(stderr);
 return 0;
      }else {
 if(strcasecmp(identifier,"GRAYSCALE") == 0) {
   type = PPGPLOT_GRAYSCALE;
 }else if(strcasecmp(identifier,"INVERTED_GRAYSCALE") == 0) {
   type = PPGPLOT_INVERTED_GRAYSCALE;
 }else if(strcasecmp(identifier,"RED") == 0) {
   type = PPGPLOT_RED;
 }else if(strcasecmp(identifier,"GREEN") == 0) {
   type = PPGPLOT_GREEN;
 }else if(strcasecmp(identifier,"BLUE") == 0) {
   type = PPGPLOT_BLUE;
 }else if(strcasecmp(identifier,"CYAN") == 0) {
   type = PPGPLOT_CYAN;
 }else if(strcasecmp(identifier,"INVERTED_RED") == 0) {
   type = PPGPLOT_INVERTED_RED;
 }else if(strcasecmp(identifier,"INVERTED_GREEN") == 0) {
   type = PPGPLOT_INVERTED_GREEN;
 }else if(strcasecmp(identifier,"INVERTED_BLUE") == 0) {
   type = PPGPLOT_INVERTED_BLUE;
 }else if(strcasecmp(identifier,"INVERTED_CYAN") == 0) {
   type = PPGPLOT_INVERTED_CYAN;
 }else if(strcasecmp(identifier,"HEAT") == 0) {
   type = PPGPLOT_HEAT;
 }else if(strcasecmp(identifier,"INVERTED_HEAT") == 0) {
   type = PPGPLOT_INVERTED_HEAT;
 }else if(strcasecmp(identifier,"COLD") == 0) {
   type = PPGPLOT_COLD;
 }else if(strcasecmp(identifier,"INVERTED_COLD") == 0) {
   type = PPGPLOT_INVERTED_COLD;
 }else if(strcasecmp(identifier,"DIVREDBLUE") == 0) {
   type = PPGPLOT_DIVREDBLUE;
 }else if(strcasecmp(identifier,"INVERTED_DIVREDBLUE") == 0 || strcasecmp(identifier,"NL") == 0) {
   type = PPGPLOT_INVERTED_DIVREDBLUE;
 }else if(strcasecmp(identifier,"HEAT2") == 0) {
   type = PPGPLOT_HEAT2;
 }else if(strcasecmp(identifier,"INVERTED_HEAT2") == 0) {
   type = PPGPLOT_INVERTED_HEAT2;
 }else if(strcasecmp(identifier,"HEAT3") == 0) {
   type = PPGPLOT_HEAT3;
 }else if(strcasecmp(identifier,"INVERTED_HEAT3") == 0) {
   type = PPGPLOT_INVERTED_HEAT3;
 }else if(strcasecmp(identifier,"HEAT4") == 0) {
   type = PPGPLOT_HEAT4;
 }else if(strcasecmp(identifier,"INVERTED_HEAT4") == 0) {
   type = PPGPLOT_INVERTED_HEAT4;
 }else if(strcasecmp(identifier,"INFERNO") == 0 || strcasecmp(identifier,"CRISTINA") == 0) {
   type = PPGPLOT_INFERNO;
 }else {
   fflush(stdout);
   printerror(debug, "ERROR cmap_parse_commandline:  '%s' not recognized as a color map.", identifier);
   printCMAPCommandlineOptions(stderr);
   return 0;
 }
      }
    }
  }
  return type;
}
void drawSphericalGrid(float dlat, float dlong, float rot_long, float rot_lat, int lw, int projection)
{
  float lat, lon, t, dt, x, y, sign, weight;
  int first, side, side2, oldside, oldside2, ok;
  ppgslw(lw);
  dlat *= M_PI/180.0;
  dlong *= M_PI/180.0;
  rot_long *= M_PI/180.0;
  rot_lat *= M_PI/180.0;
  dt = 0.01;
  if(projection != 3) {
    for(side = 0; side < 2; side++) {
      first = 1;
      for(t = 0; t <= 2.0*M_PI+2.0*dt; t += dt) {
 ok = 1;
 if(projection == 1) {
   x = 2*cos(t);
   y = sin(t);
   if(side == 1)
     ok = 0;
 }else if(projection == 2) {
   x = cos(t);
   y = sin(t);
   if(side == 0)
     x -= 1.125;
   else
     x += 1.125;
 }
 if(first && ok) {
   ppgmove(x, y);
   first = 0;
 }else if(ok) {
   ppgdraw(x, y);
 }
      }
    }
  }
  for(sign = -1; sign < 1.1; sign += 2) {
    for(lat = 0; lat <= 0.5*M_PI+2.0*dt; lat += dlat) {
      first = 1;
      oldside = -1;
      oldside2 = -1;
      for(t = -M_PI; t <= M_PI+2.0*dt; t += dt) {
 lon = t;
 if(lon > M_PI)
   lon = M_PI;
 if(projection == 1) {
   projectionHammerAitoff_xy(lon, sign*lat, rot_long, rot_lat, &x, &y);
   if(x < 0)
     side = 0;
   if(x > 0)
     side = 1;
   if(side != oldside)
     first = 1;
   oldside = side;
 }else if(projection == 2) {
   side = projection_sphere_xy(lon, sign*lat, rot_long, rot_lat, &x, &y, &weight);
   if(side != oldside)
     first = 1;
   oldside = side;
 }else if(projection == 3) {
   projection_longlat_xy(lon, sign*lat, rot_long, rot_lat, &x, &y);
   if(x < 0)
     side = 0;
   if(x > 0)
     side = 1;
   if(y < 0)
     side2 = 0;
   if(y > 0)
     side2 = 1;
   if(side != oldside || side2 != oldside2)
     first = 1;
   oldside = side;
   oldside2 = side2;
 }
 if(first) {
   ppgmove(x, y);
   first = 0;
 }else {
   ppgdraw(x, y);
 }
      }
    }
  }
  for(sign = -1; sign < 1.1; sign += 2) {
    for(lon = 0; lon <= M_PI; lon += dlong) {
      first = 1;
      oldside = -1;
      oldside2 = -1;
      for(t = -0.5*M_PI; t <= 0.5*M_PI+2.0*dt; t += dt) {
 lat = t;
 if(lat > 0.5*M_PI)
   lat = 0.5*M_PI;
 if(projection == 1) {
   projectionHammerAitoff_xy(sign*lon, lat, rot_long, rot_lat, &x, &y);
   if(x < 0)
     side = 0;
   if(x > 0)
     side = 1;
   if(side != oldside)
     first = 1;
   oldside = side;
 }else if(projection == 2) {
   side = projection_sphere_xy(sign*lon, lat, rot_long, rot_lat, &x, &y, &weight);
   if(side != oldside)
     first = 1;
   oldside = side;
 }else if(projection == 3) {
   projection_longlat_xy(sign*lon, lat, rot_long, rot_lat, &x, &y);
   if(x < 0)
     side = 0;
   if(x > 0)
     side = 1;
   if(y < 0)
     side2 = 0;
   if(y > 0)
     side2 = 1;
   if(side != oldside || side2 != oldside2)
     first = 1;
   oldside = side;
   oldside2 = side2;
 }
 if(first) {
   ppgmove(x, y);
   first = 0;
 }else {
   ppgdraw(x, y);
 }
      }
    }
  }
  ppgslw(1);
}
int pgplot_device_type(char *devicename, verbose_definition verbose)
{
  int i, n, found, type;
  char *curdev;
  if(devicename != NULL) {
    curdev = devicename;
  }else {
    curdev = calloc(1000, 1);
    if(curdev == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplot_device_type: Memory allocation error");
      return 0;
    }
    n = 998;
    curdev[0] = '/';
    ppgqinf("type", curdev+1, &n);
  }
  n = strlen(curdev);
  if(n <= 1)
    return 0;
  found = 0;
  for(i = n; i >= 0; i--) {
    if(curdev[i] == '/') {
      found = 1;
      break;
    }
  }
  type = 0;
  if(found) {
    if(strcasecmp(&curdev[i], "/xwindow") == 0) {
      type = 1;
    }else if(strcasecmp(&curdev[i], "/xw") == 0) {
      type = 1;
    }else if(strcasecmp(&curdev[i], "/xserve") == 0) {
      type = 2;
    }else if(strcasecmp(&curdev[i], "/xs") == 0) {
      type = 2;
    }else if(strcasecmp(&curdev[i], "/ps") == 0) {
      type = 3;
    }else if(strcasecmp(&curdev[i], "/vps") == 0) {
      type = 4;
    }else if(strcasecmp(&curdev[i], "/cps") == 0) {
      type = 5;
    }else if(strcasecmp(&curdev[i], "/vcps") == 0) {
      type = 6;
    }else if(strcasecmp(&curdev[i], "/latex") == 0) {
      type = 7;
    }else if(strcasecmp(&curdev[i], "/null") == 0) {
      type = 8;
    }else if(strcasecmp(&curdev[i], "/png") == 0) {
      type = 9;
    }else if(strcasecmp(&curdev[i], "/tpng") == 0) {
      type = 10;
    }else if(strcasecmp(&curdev[i], "/pw") == 0) {
      type = 100;
    }
  }
  if(devicename == NULL) {
    free(curdev);
  }
  return type;
}

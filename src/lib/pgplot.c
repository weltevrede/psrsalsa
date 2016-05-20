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
  *x = ((double)nx)*(internal_pgplot_xmax-internal_pgplot_xmin)/((double)internal_pgplot_nrx-1.0) + internal_pgplot_xmin;
  *y = ((double)ny)*(internal_pgplot_ymax-internal_pgplot_ymin)/((double)internal_pgplot_nry-1.0) + internal_pgplot_ymin;
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
    frame->TR[0] = internal_pgplot_xmin - (internal_pgplot_xmax-internal_pgplot_xmin)/(float)(internal_pgplot_nrx-1);
    frame->TR[1] = (internal_pgplot_xmax-internal_pgplot_xmin)/(float)(internal_pgplot_nrx-1); frame->TR[2] = 0;
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
int pgplotPAplot(datafile_definition data, float *Ppulse, pgplot_viewport_def viewport, pgplot_box_def box, char *xlabel, char *ylabel, char *ylabel_pa, float longitude_left, float longitude_right, float Imin, float Imax, float pa_bottom, float pa_top, float PAoffset, float sigma_limit, float datalinewidth, float ysize2, int dashed, int noynumbers, char *textoption, char *herrorbaroption, char *herrorbaroption2, int argc, char **argv, int outline_txt, int outline_lw, int outline_color, int overlayPA, float overlayalpha, float overlaybeta, float overlaypa0, float overlayl0, int overlayPAfine, int nrJumps, float *jump_longitudes, float *jump_offsets, verbose_definition verbose)
{
  int deviceID, text_ci, text_lw, text_f, ok, domove;
  float ymin, ymax, text_x, text_y, text_ch, I, Iold;
  long i, j;
  char *newtext;
  pgplot_frame_def_internal frame;
  if(pgplot_opendevice(viewport, &deviceID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Cannot open plot device");
    return 0;
  }
  if(data.poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does only have PA points, but no profile");
    return 0;
  }
  if(data.poltype != POLTYPE_ILVPAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does not appear to contain PA data");
    return 0;
  }
  if(data.NrPols != 5) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does not appear to contain PA data as NrPols != 5");
    return 0;
  }
  newtext = str_replace_header_params(data, box.title, verbose);
  if(newtext == NULL) {
    printerror(verbose.debug, "ERROR pgplotPAplot: Cannot substuture header parameter in the tile");
    return 0;
  }
  strcpy(box.title, newtext);
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
  }
  ppgsvp(0.1+viewport.dxplot, 0.1+viewport.dxplot+0.8*viewport.xsize, 0.35+viewport.dyplot, 0.35+viewport.dyplot+0.55*viewport.ysize);
  if(Imin != Imax) {
    frame.swin_x1 = longitude_left;
    frame.swin_x2 = longitude_right;
    frame.swin_y1 = Imin;
    frame.swin_y2 = Imax;
    ppgswin(longitude_left,longitude_right,Imin,Imax);
  }else {
    frame.swin_x1 = longitude_left;
    frame.swin_x2 = longitude_right;
    frame.swin_y1 = ymin-(ymax-ymin)*0.05;
    frame.swin_y2 = ymax+(ymax-ymin)*0.05;
    ppgswin(longitude_left,longitude_right,ymin-(ymax-ymin)*0.05,ymax+(ymax-ymin)*0.05);
  }
  if(noynumbers) {
    strcpy(box.box_xopt, "bcst");
    strcpy(box.box_yopt, "bcts");
  }else {
    strcpy(box.box_xopt, "bcst");
    strcpy(box.box_yopt, "bcnts");
  }
  if(!noynumbers) {
    box.drawlabels = 1;
    strcpy(box.ylabel, ylabel);
  }
  pgplot_drawbox(box);
  box.drawtitle = 0;
  ppgslw(datalinewidth);
  if(dashed) {
    ppgsls(4);
  }else
    ppgsci(3);
  domove = 1;
  for(j = 0; j < (data.NrBins); j++) {
    if(isnan(get_pulse_longitude(data, 0, j, verbose)) || isnan(data.data[j+2*data.NrBins])) {
      domove = 1;
    }else {
      if(domove) {
 ppgmove(get_pulse_longitude(data, 0, j, verbose), data.data[j+2*data.NrBins]);
 domove = 0;
      }else {
 ppgdraw(get_pulse_longitude(data, 0, j, verbose), data.data[j+2*data.NrBins]);
      }
    }
  }
  if(dashed)
    ppgsls(2);
  else
    ppgsci(2);
  domove = 1;
  for(j = 0; j < (data.NrBins); j++) {
    if(isnan(get_pulse_longitude(data, 0, j, verbose)) || isnan(data.data[j+1*data.NrBins])) {
      domove = 1;
    }else {
      if(domove) {
 ppgmove(get_pulse_longitude(data, 0, j, verbose), data.data[j+1*data.NrBins]);
 domove = 0;
      }else {
 ppgdraw(get_pulse_longitude(data, 0, j, verbose), data.data[j+1*data.NrBins]);
      }
    }
  }
  ppgsci(1);
  if(dashed) {
    ppgsls(1);
  }
  domove = 1;
  for(j = 0; j < (data.NrBins); j++) {
    if(isnan(get_pulse_longitude(data, 0, j, verbose)) || isnan(data.data[j])) {
      domove = 1;
    }else {
      if(domove) {
 ppgmove(get_pulse_longitude(data, 0, j, verbose), data.data[j]);
 domove = 0;
      }else {
 ppgdraw(get_pulse_longitude(data, 0, j, verbose), data.data[j]);
      }
    }
  }
  if(dashed) {
    ppgslw(4);
    ppgsls(1);
    ppgslw(datalinewidth);
  }
  if(!dashed) {
    ppgsci(4);
    ppgsls(2);
  }
  if(Ppulse != NULL) {
    ppgmove(get_pulse_longitude(data, 0, 0, verbose), Ppulse[0]);
    for(j = 1; j < (data.NrBins); j++) {
      ppgdraw(get_pulse_longitude(data, 0, j, verbose), Ppulse[j]);
    }
  }
  ppgsls(1);
  if(textoption != NULL && argv != NULL && argc != 0) {
    for(j = 0; j < argc; j++) {
      if(strcmp(argv[j], textoption) == 0) {
 i = sscanf(argv[j+1], "%f %f %f %d %d %d", &text_x, &text_y, &text_ch, &text_lw, &text_f, &text_ci);
 if(i != 6) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR pgplotPAplot: Error parsing %s option", textoption);
   return 0;
 }
 ppgsch(text_ch);
 ppgscf(text_f);
 newtext = str_replace_header_params(data, argv[j+2], verbose);
 if(newtext == NULL) {
   printerror(verbose.debug, "ERROR pgplotPAplot: Cannot substuture header parameter in -text option");
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
  if(herrorbaroption != NULL && argv != NULL && argc != 0) {
    for(j = 0; j < argc; j++) {
      float herr_x1, herr_x2, herr_x3, herr_y, herr_size, herr_lw;
      int herr_ci;
      if(strcmp(argv[j], herrorbaroption) == 0) {
 i = sscanf(argv[j+1], "%f %f %f %f %f %f %d", &herr_x1, &herr_x2, &herr_x3, &herr_y, &herr_size, &herr_lw, &herr_ci);
 if(i != 7) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR pgplotPAplot: Error parsing %s option", herrorbaroption);
   return 0;
 }
 ppgsch(herr_size);
 ppgsci(herr_ci);
 ppgslw(herr_lw);
 ppgpt1(herr_x2, herr_y*(frame.swin_y2-frame.swin_y1)+frame.swin_y1, 2);
 ppgerr1(1, herr_x2, herr_y*(frame.swin_y2-frame.swin_y1)+frame.swin_y1, herr_x3-herr_x2, herr_size);
 ppgerr1(3, herr_x2, herr_y*(frame.swin_y2-frame.swin_y1)+frame.swin_y1, herr_x2-herr_x1, herr_size);
 ppgsch(1);
 ppgslw(1);
 ppgsci(1);
 j += 1;
      }
    }
  }
  ppgsls(1);
  ppgsci(1);
  ppgslw(1);
  ppgsvp(0.1+viewport.dxplot,0.1+viewport.dxplot+0.8*viewport.xsize, 0.35+viewport.dyplot-0.25*viewport.ysize*ysize2, 0.35+viewport.dyplot);
  frame.swin_x1 = longitude_left;
  frame.swin_x2 = longitude_right;
  frame.swin_y1 = pa_bottom;
  frame.swin_y2 = pa_top;
  ppgswin(longitude_left,longitude_right,pa_bottom,pa_top);
  ppgslw(box.box_lw);
  box.drawlabels = 1;
  if(!noynumbers) {
    strcpy(box.xlabel, xlabel);
    strcpy(box.ylabel, ylabel_pa);
  }else {
    strcpy(box.xlabel, xlabel);
    box.ylabel[0] = 0;
  }
  strcpy(box.box_xopt, "bcnst");
  pgplot_drawbox(box);
  float x, xold;
  ppgslw(datalinewidth);
  if(overlayPAfine)
    overlayPAfine = 100;
  else
    overlayPAfine = 1;
  if(overlayPA) {
    ppgsci(2);
    for(j = 0; j < (data.NrBins); j++) {
      for(i = 0; i < overlayPAfine; i++) {
 x = get_pulse_longitude(data, 0, j, verbose);
 if(j == data.NrBins-1) {
   x += (get_pulse_longitude(data, 0, j, verbose)-get_pulse_longitude(data, 0, j-1, verbose))*i/(float)overlayPAfine;
 }else {
   x += (get_pulse_longitude(data, 0, j+1, verbose)-get_pulse_longitude(data, 0, j, verbose))*i/(float)overlayPAfine;
 }
 I = paswing(overlayalpha, overlaybeta, x, overlaypa0, overlayl0, nrJumps, jump_longitudes, jump_offsets, 0, 0);
 if(I > 180)
   I -= 180;
 if(I > 180)
   I -= 180;
 if(I < 0)
   I += 180;
 if(I < 0)
   I += 180;
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
    ppgsci(1);
  }
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
      ppgpt1(get_pulse_longitude(data, 0, j, verbose), I, -2);
      ymax = I+data.data[j+4*data.NrBins];
      ymin = I-data.data[j+4*data.NrBins];
      ppgerr1(6, get_pulse_longitude(data, 0, j, verbose), I, data.data[j+4*data.NrBins], 1.0);
      I = derotate_180(data.data[j+3*data.NrBins]+PAoffset);
      ppgpt1(get_pulse_longitude(data, 0, j, verbose), I, -2);
      ymax = I+data.data[j+4*data.NrBins];
      ymin = I-data.data[j+4*data.NrBins];
      ppgerr1(6, get_pulse_longitude(data, 0, j, verbose), I, data.data[j+4*data.NrBins], 1.0);
      I = derotate_180(data.data[j+3*data.NrBins]+PAoffset) - 180.0;
      ppgpt1(get_pulse_longitude(data, 0, j, verbose), I, -2);
      ymax = I+data.data[j+4*data.NrBins];
      ymin = I-data.data[j+4*data.NrBins];
      ppgerr1(6, get_pulse_longitude(data, 0, j, verbose), I, data.data[j+4*data.NrBins], 1.0);
      I = derotate_180(data.data[j+3*data.NrBins]+PAoffset) - 360.0;
      ppgpt1(get_pulse_longitude(data, 0, j, verbose), I, -2);
      ymax = I+data.data[j+4*data.NrBins];
      ymin = I-data.data[j+4*data.NrBins];
      ppgerr1(6, get_pulse_longitude(data, 0, j, verbose), I, data.data[j+4*data.NrBins], 1.0);
    }
  }
  if(herrorbaroption2 != NULL && argv != NULL && argc != 0) {
    for(j = 0; j < argc; j++) {
      float herr_x1, herr_x2, herr_x3, herr_y, herr_size, herr_lw;
      int herr_ci;
      if(strcmp(argv[j], herrorbaroption2) == 0) {
 i = sscanf(argv[j+1], "%f %f %f %f %f %f %d", &herr_x1, &herr_x2, &herr_x3, &herr_y, &herr_size, &herr_lw, &herr_ci);
 if(i != 7) {
   fflush(stdout);
   printerror(verbose.debug, "ERROR pgplotPAplot: Error parsing %s option", herrorbaroption2);
   return 0;
 }
 ppgsch(herr_size);
 ppgsci(herr_ci);
 ppgslw(herr_lw);
 ppgpt1(herr_x2, herr_y, 2);
 ppgerr1(1, herr_x2, herr_y, herr_x3-herr_x2, herr_size);
 ppgerr1(3, herr_x2, herr_y, herr_x2-herr_x1, herr_size);
 ppgsch(1);
 ppgslw(1);
 ppgsci(1);
 j += 1;
      }
    }
  }
  if(!viewport.dontclose)
    ppgclos();
  return 1;
}
int pgplotGraph1(pgplot_viewport_def viewport, float *data, float *datax, float *sigma, int nrx, float xmin, float xmax, int dontsetranges, float xmin_show, float xmax_show, float ymin_show, float ymax_show, int forceMinZero, pgplot_box_def pgplotbox, int hist, int noline, int pointtype, int color, int boxcolor, regions_definition *regions, verbose_definition verbose)
{
  float min, max;
  float x, xnext, y;
  int binnr, regnr, loopnr, deviceID, notset, color_prev, color_cur;
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
  if(pgplot_opendevice(viewport, &deviceID, verbose) == 0) {
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
  frame.svp_x1 = 0.1+viewport.dxplot;
  frame.svp_x2 = 0.1+viewport.dxplot+0.8*viewport.xsize;
  frame.svp_y1 = 0.1+viewport.dyplot;
  frame.svp_y2 = 0.1+viewport.dyplot+0.8*viewport.ysize;
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
  pgplot_drawbox(pgplotbox);
  ppgsci(1);
  color_prev = color;
  ppgbbuf();
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
     if(binnr > regions->left_bin[regnr] && binnr <= regions->right_bin[regnr])
       color_cur = color+1+regnr;
   }else {
     if(binnr >= regions->left_bin[regnr] && binnr <= regions->right_bin[regnr])
       color_cur = color+1+regnr;
   }
 }
      }
      if(loopnr == 0) {
 if(noline == 0) {
   if(hist == 0) {
     if(binnr == 0) {
       ppgmove(x, y);
     }else {
       ppgsci(color_cur);
       ppgdraw(x, y);
     }
   }else {
     if(binnr == 0) {
       if(forceMinZero) {
  ppgmove(x, 0);
  ppgdraw(x, y);
       }else {
  ppgmove(x, y);
       }
     }else {
       if(color_cur == color || color_prev == color) {
  ppgsci(color);
       }else {
  if(color_cur == color_prev)
    ppgsci(color_cur);
  else
    ppgsci(color_prev);
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
  if(!viewport.dontclose)
    ppgclos();
  return 1;
}
int pgplot_opendevice(pgplot_viewport_def viewport, int *deviceID, verbose_definition verbose)
{
  if(viewport.dontopen == 0) {
    *deviceID = ppgopen(viewport.plotDevice);
    if(*deviceID <= 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplot_opendevice: Cannot open plot device called '%s'", viewport.plotDevice);
      return 0;
    }
    pgplot_setWindowsize(viewport.windowwidth, viewport.windowheight, viewport.aspectratio);
    ppgask(0);
  }
  if(viewport.noclear == 0)
    ppgpage();
  return 1;
}
void pgplot_clear_viewport_def(pgplot_viewport_def *viewport)
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
void clear_pgplot_box(pgplot_box_def *box)
{
  box->drawbox = 1;
  box->box_lw = 1;
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
  box->dxlabel = 0;
  box->dylabel = 0;
  box->xlabel[0] = 0;
  box->ylabel[0] = 0;
  box->wedgelabel[0] = 0;
}
void pgplot_drawbox(pgplot_box_def box)
{
  if(box.drawbox) {
    ppgslw(box.box_lw);
    ppgsch(box.box_labelsize);
    ppgbox(box.box_xopt,box.box_xtick,box.box_nxsub,box.box_yopt,box.box_ytick,box.box_nysub);
  }
  if(box.drawtitle) {
    ppgsch(box.title_ch);
    ppgslw(box.title_lw);
    ppgscf(box.title_f);
    ppgmtxt("T", (1.5*box.box_labelsize)/box.title_ch, 0.5, 0.5, box.title);
  }
  if(box.drawlabels) {
    ppgslw(box.box_lw);
    ppgsch(box.label_ch);
    ppgmtxt("L", ((2.5+box.dylabel)*box.box_labelsize)/box.label_ch, 0.5, 0.5, box.ylabel);
    ppgmtxt("B", ((1.5+box.dxlabel)*box.box_labelsize)/box.label_ch+1, 0.5, 0.5, box.xlabel);
  }
  ppgsch(1);
  ppgslw(1);
}
int pgplotMap(pgplot_viewport_def viewport, float *cmap, int nrx, int nry, float xmin, float xmax, float xminshow, float xmaxshow, float ymin, float ymax, float yminshow, float ymaxshow, pgplot_box_def pgplotbox, int maptype, int itf, int nogray, int nrcontours, float *contours, int contourlw, int forceMinZero, float saturize, int levelset, float levelmin, float levelmax, int levelInversion, int onlyData, int sideright, int forceMinZeroRight, int sidetop, int forceMinZeroTop, int sidelw, int showwedge, int plotSubset, int showTwice, verbose_definition verbose)
{
  float datamin, datamax, xmin2, xmax2;
  float *collapse, x, y, remember2_x1, remember2_x2, remember2_y1, remember2_y2;
  int i, j, deviceID, firstc, lastc, firstpoint;
  float junk_f, lasty;
  pgplot_frame_def_internal pgplot_frame_internal;
  clear_pgplot_frame(&pgplot_frame_internal);
  internal_pgplot_xmin = xmin;
  internal_pgplot_xmax = xmax;
  internal_pgplot_ymin = ymin;
  internal_pgplot_ymax = ymax;
  internal_pgplot_nrx = nrx;
  internal_pgplot_nry = nry;
  lasty = 0;
  if(verbose.debug) {
    printf("pgplotMap -- dimensions: %dX%d points  xrange=%e %e yrange=%e %e\n", nrx, nry, xmin, xmax, ymin, ymax);
    printf("pgplotMap -- Horizontal dimensions %f - %f (shown %f - %f)\n", xmin, xmax, xminshow, xmaxshow);
    printf("pgplotMap -- Vertical dimensions %f - %f (shown %f - %f)\n", ymin, ymax, yminshow, ymaxshow);
  }
  if(pgplot_opendevice(viewport, &deviceID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotMap: Cannot open plot device");
    return 0;
  }
  if(onlyData == 0) {
    pgplot_setWindowsize(viewport.windowwidth, viewport.windowheight, viewport.aspectratio);
  }
  ppgslw(1);
  ppgsci(1);
  ppgscf(1);
  ppgsch(1);
  if(onlyData == 0 || onlyData == 2) {
    pgplot_frame_internal.svp_x1 = viewport.dxplot+0.15;
    pgplot_frame_internal.svp_y1 = viewport.dyplot+0.15;
    if(sideright)
      pgplot_frame_internal.svp_x2 = (0.75-0.15)*viewport.xsize+pgplot_frame_internal.svp_x1;
    else
      pgplot_frame_internal.svp_x2 = (0.90-0.15)*viewport.xsize+pgplot_frame_internal.svp_x1;
    if(sidetop)
      pgplot_frame_internal.svp_y2 = (0.75-0.15)*viewport.ysize+pgplot_frame_internal.svp_y1;
    else
      pgplot_frame_internal.svp_y2 = (0.90-0.15)*viewport.ysize+pgplot_frame_internal.svp_y1;
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
    return 0;
  }
  int subset_nrx, subset_nry, subset_extrabefore, subset_extraafter;
  float *cmap_subset;
  if(plotSubset) {
    int subset_x0, subset_y0;
    float subset_xmin, subset_xmax, subset_ymin, subset_ymax;
    subset_extrabefore = 1;
    subset_extraafter = 1;
    subset_x0 = (nrx-1)*(xminshow - xmin)/(xmax-xmin)-1;
    if(subset_x0 < 0)
      subset_x0 = 0;
    subset_nrx = (nrx-1)*(xmaxshow - xmin)/(xmax-xmin)+1 - subset_x0;
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
 plotSubset = 0;
      }
    }
    if(plotSubset) {
      for(i = 0; i < subset_nrx; i++) {
 for(j = 0; j < subset_nry; j++) {
   cmap_subset[j*subset_nrx+i] = cmap[(j+subset_y0)*nrx+i+subset_x0];
 }
      }
      subset_xmin = subset_x0*(xmax-xmin)/(float)(nrx-1) + xmin;
      subset_xmax = (subset_x0+subset_nrx)*(xmax-xmin)/(float)(nrx-1) + xmin;
      subset_ymin = subset_y0*(ymax-ymin)/(float)(nry-1) + ymin;
      subset_ymax = (subset_y0+subset_nry)*(ymax-ymin)/(float)(nry-1) + ymin;
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
    ppgsch(pgplotbox.box_labelsize*0.5);
    ppgwedg("RI", 0, 6, datamax, datamin, pgplotbox.wedgelabel);
  }
  if(nrcontours > 0) {
    ppgslw(contourlw);
    if(showTwice) {
      pgplot_frame_internal.TR[3] += pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
      ppgcont(cmap_subset, subset_nrx, subset_nry, 1, subset_nrx, 1, subset_nry, contours, nrcontours, pgplot_frame_internal.TR);
      pgplot_frame_internal.TR[3] -= pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
    }
    ppgcont(cmap_subset, subset_nrx, subset_nry, 1, subset_nrx, 1, subset_nry, contours, nrcontours, pgplot_frame_internal.TR);
  }
  if(plotSubset)
    free(cmap_subset);
  if(onlyData != 0)
    pgplotbox.drawbox = 0;
  if(sidetop && onlyData != 1)
    pgplotbox.drawtitle = 0;
  if(onlyData != 0)
    pgplotbox.drawlabels = 0;
  pgplot_drawbox(pgplotbox);
  if(sideright) {
    collapse = (float *)calloc(nry, sizeof(float));
    if(collapse == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplotMap: Cannot allocate memory");
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
      junk_f = (0.75-0.15)*viewport.ysize+viewport.dyplot+0.15;
    else
      junk_f = (0.90-0.15)*viewport.ysize+viewport.dyplot+0.15;
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
      ppgsch(pgplotbox.box_labelsize*0.5);
      ppgwedg("RI", 0, 6, datamax, datamin, pgplotbox.wedgelabel);
    }
    if(onlyData == 0 || onlyData == 2) {
      ppgslw(pgplotbox.box_lw);
      ppgsch(pgplotbox.box_labelsize*0.5);
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
      junk_f = (0.75-0.15)*viewport.xsize+viewport.dxplot+0.15;
    else
      junk_f = (0.90-0.15)*viewport.xsize+viewport.dxplot+0.15;
    if(onlyData != 0) {
      ppgqvp(0, &pgplot_frame_internal.svp_x1, &pgplot_frame_internal.svp_x2, &pgplot_frame_internal.svp_y1, &pgplot_frame_internal.svp_y2);
      ppgqwin(&remember2_x1, &remember2_x2, &remember2_y1, &remember2_y2);
      junk_f = pgplot_frame_internal.svp_x2;
    }
    ppgsvp(pgplot_frame_internal.svp_x1, junk_f, pgplot_frame_internal.svp_y2, pgplot_frame_internal.svp_y2+0.15);
    ppgswin(xminshow, xmaxshow, xmin2, xmax2);
    if(onlyData == 0 || onlyData == 2) {
      ppgsch(pgplotbox.box_labelsize*0.5);
      ppgslw(pgplotbox.box_lw);
      if(sidetop == 2)
 ppgbox("bcsti",0.0,0,"bctsi",0.0,0);
      else
 ppgbox("bcsti",0.0,0,"bcntsi",0.0,0);
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
      ppgsch(pgplotbox.title_ch);
      ppgslw(pgplotbox.title_lw);
      ppgscf(pgplotbox.title_f);
      ppgptxt(0.5*(xmaxshow-xminshow)+xminshow, xmax2+0.3*(pgplotbox.label_ch/1.7)*(xmax2-xmin2), 0, 0.5, pgplotbox.title);
      ppgscf(1);
      ppgslw(1);
      ppgsch(pgplotbox.label_ch);
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
  if(!viewport.dontclose)
    ppgclos();
  return 1;
}
int selectRegions(float *profileI, int nrBins, pgplot_viewport_def viewport, pgplot_box_def pgplotbox, int onlyOne, int powerTwo, int evenNumber, regions_definition *regions, verbose_definition verbose)
{
  int i, j, k, bin1, bin2, replot;
  float x, y;
  char c;
  pgplot_viewport_def viewport2;
  x = y = 0;
  c = 0;
  bin1 = bin2 = 0;
  if(strcmp(viewport.plotDevice, "?") == 0)
    printf("Select plotting device to show profile in order to select the onpulse region\n  ");
  memcpy(&viewport2, &viewport, sizeof(pgplot_viewport_def));
  viewport2.dontclose = 1;
  pgplotGraph1(viewport2, profileI, NULL, NULL, nrBins, 0, nrBins-1, 0, 0, nrBins-1, 0, 0, 0, pgplotbox, 1, 0, 0, 1, 1, regions, verbose);
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
      if(viewport.dontclose == 0) {
 ppgclos();
 fflush(stdout);
      }
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
      if(regions->nrRegions == maxNrRegions) {
 fflush(stdout);
 printerror(verbose.debug, "ERROR selectRegions: Too many regions selected.");
 ppgend();
 return regions->nrRegions;
      }
      replot = 1;
      i = 0;
    }
    if(replot) {
      memcpy(&viewport2, &viewport, sizeof(pgplot_viewport_def));
      viewport2.dontclose = 1;
      viewport2.dontopen = 1;
      pgplotGraph1(viewport2, profileI, NULL, NULL, nrBins, 0, nrBins-1, 0, 0, nrBins-1, 0, 0, 0, pgplotbox, 1, 0, 0, 1, 1, regions, verbose);
      replot = 0;
    }
  }while(i < 3);
  if(viewport.dontclose == 0)
    ppgclos();
  return regions->nrRegions;
}
int checkRegions(int bin, regions_definition *regions, int whichregion, verbose_definition verbose)
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
void clearRegion(regions_definition *region)
{
  int i;
  region->nrRegions = 0;
  for(i = 0; i < maxNrRegions; i++) {
    region->frac_defined[i] = 0;
    region->bins_defined[i] = 0;
  }
}
void region_frac_to_int(regions_definition *region, float scale, float offset)
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
void region_int_to_frac(regions_definition *region, float scale, float offset)
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
void regionShowNextTimeUse(regions_definition region, char *option, char *optionFrac, FILE *where)
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
    }
  }
  if(devicename == NULL) {
    free(curdev);
  }
  return type;
}

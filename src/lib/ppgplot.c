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
#include <stdlib.h>
#include "cpgplot.h"
static FILE *ppgdevice_internal = NULL;
static char ppgdevicename_internal[1000];
static int ppg_nextaux_internal;
int ppgopen(const char *device)
{
  fflush(stdout);
  if(ppgdevice_internal == NULL)
    return cpgopen(device);
  else {
    fprintf(ppgdevice_internal, "pgopen %s\n", device);
    return 1;
  }
}
void ppgclos(void)
{
  fflush(stdout);
  if(ppgdevice_internal == NULL)
    cpgclos();
  else {
    fprintf(ppgdevice_internal, "pgclos\n");
  }
}
void ppgend(void)
{
  fflush(stdout);
  if(ppgdevice_internal == NULL)
    cpgend();
  else {
    fprintf(ppgdevice_internal, "pgend\n");
  }
}
void ppgbbuf(void)
{
  if(ppgdevice_internal == NULL)
    cpgbbuf();
  else {
    fprintf(ppgdevice_internal, "pgbbuf\n");
  }
}
void ppgebuf(void)
{
  if(ppgdevice_internal == NULL)
    cpgebuf();
  else {
    fprintf(ppgdevice_internal, "pgebuf\n");
  }
}
void ppgpap(float width, float aspect)
{
  fflush(stdout);
  if(ppgdevice_internal == NULL)
    cpgpap(width, aspect);
  else {
    fprintf(ppgdevice_internal, "pgpap %f %f\n", width, aspect);
  }
}
void ppgsvp(float xleft, float xright, float ybot, float ytop)
{
  if(ppgdevice_internal == NULL)
    cpgsvp(xleft, xright, ybot, ytop);
  else {
    fprintf(ppgdevice_internal, "pgsvp %f %f %f %f\n", xleft, xright, ybot, ytop);
  }
}
int ppgqvp(int units, float *xleft, float *xright, float *ybot, float *ytop)
{
  if(ppgdevice_internal == NULL) {
    cpgqvp(units, xleft, xright, ybot, ytop);
    return 1;
  }else {
    fflush(stdout);
    fprintf(stderr, "ERROR ppgqvp not defined in pgplotinterpretor mode\n");
    return 0;
  }
}
void ppgqinf(const char *item, char *value, int *value_length)
{
  if(ppgdevice_internal == NULL) {
    cpgqinf(item, value, value_length);
    return;
  }else {
    fflush(stdout);
    fprintf(stderr, "ERROR ppgqinf not defined in pgplotinterpretor mode\n");
    return;
  }
}
int ppgqwin(float *xleft, float *xright, float *ybot, float *ytop)
{
  if(ppgdevice_internal == NULL) {
    cpgqwin(xleft, xright, ybot, ytop);
    return 1;
  }else {
    fflush(stdout);
    fprintf(stderr, "ERROR ppgqwin not defined in pgplotinterpretor mode\n");
    return 0;
  }
}
void ppgswin(float x1, float x2, float y1, float y2)
{
  if(ppgdevice_internal == NULL)
    cpgswin(x1, x2, y1, y2);
  else {
    fprintf(ppgdevice_internal, "pgswin %f %f %f %f\n", x1, x2, y1, y2);
  }
}
void ppglab(const char *xlbl, const char *ylbl, const char *toplbl)
{
  if(ppgdevice_internal == NULL)
    cpglab(xlbl, ylbl, toplbl);
  else {
    fprintf(ppgdevice_internal, "pglab '%s' '%s' '%s'\n", xlbl, ylbl, toplbl);
  }
}
void ppgbox(const char *xopt, float xtick, int nxsub, const char *yopt, float ytick, int nysub)
{
  if(ppgdevice_internal == NULL)
    cpgbox(xopt, xtick, nxsub, yopt, ytick, nysub);
  else {
    fprintf(ppgdevice_internal, "pgbox %s %f %d %s %f %d\n", xopt, xtick, nxsub, yopt, ytick, nysub);
  }
}
void ppgask(int flag)
{
  if(ppgdevice_internal == NULL)
    cpgask(flag);
  else {
    fprintf(ppgdevice_internal, "pgask %d\n", flag);
  }
}
void ppgsitf(int itf)
{
  if(ppgdevice_internal == NULL)
    cpgsitf(itf);
  else {
    fprintf(ppgdevice_internal, "pgsitf %d\n", itf);
  }
}
void ppgpage(void)
{
  if(ppgdevice_internal == NULL)
    cpgpage();
  else {
    fprintf(ppgdevice_internal, "pgpage\n");
  }
}
void ppgaxis(const char *opt, float x1, float y1, float x2, float y2, float v1, float v2, float step, int nsub, float dmajl, float dmajr, float fmin, float disp, float orient)
{
  if(ppgdevice_internal == NULL)
    cpgaxis(opt, x1, y1, x2, y2, v1, v2, step, nsub, dmajl, dmajr, fmin, disp, orient);
  else {
    fprintf(ppgdevice_internal, "pgaxis '%s' %f %f %f %f %f %f %f %d %f %f %f %f %f\n", opt, x1, y1, x2, y2, v1, v2, step, nsub, dmajl, dmajr, fmin, disp, orient);
  }
}
void ppgmtxt(const char *side, float disp, float coord, float fjust, const char *text)
{
  if(ppgdevice_internal == NULL)
    cpgmtxt(side, disp, coord, fjust, text);
  else {
    fprintf(ppgdevice_internal, "pgmtxt %s %f %f %f '%s'\n", side, disp, coord, fjust, text);
  }
}
void ppgtick(float x1, float y1, float x2, float y2, float v, float tikl, float tikr, float disp, float orient, const char *str)
{
  if(ppgdevice_internal == NULL)
    cpgtick(x1, y1, x2, y2, v, tikl, tikr, disp, orient, str);
  else {
    fprintf(ppgdevice_internal, "pgtick %f %f %f %f %f %f %f %f %f '%s'\n", x1, y1, x2, y2, v, tikl, tikr, disp, orient, str);
  }
}
void ppgslw(int lw)
{
  if(ppgdevice_internal == NULL)
    cpgslw(lw);
  else {
    fprintf(ppgdevice_internal, "pgslw %d\n", lw);
  }
}
void ppgsci(int ci)
{
  if(ppgdevice_internal == NULL)
    cpgsci(ci);
  else {
    fprintf(ppgdevice_internal, "pgsci %d\n", ci);
  }
}
void ppgscr(int ci, float cr, float cg, float cb)
{
  if(ppgdevice_internal == NULL)
    cpgscr(ci, cr, cg, cb);
  else {
    fprintf(ppgdevice_internal, "pgscr %d %f %f %f\n", ci, cr, cg, cb);
  }
}
void ppgsls(int ls)
{
  if(ppgdevice_internal == NULL)
    cpgsls(ls);
  else {
    fprintf(ppgdevice_internal, "pgsls %d\n", ls);
  }
}
void ppgsfs(int fs)
{
  if(ppgdevice_internal == NULL)
    cpgsfs(fs);
  else {
    fprintf(ppgdevice_internal, "pgsfs %d\n", fs);
  }
}
void ppgsch(float size)
{
  if(ppgdevice_internal == NULL)
    cpgsch(size);
  else {
    fprintf(ppgdevice_internal, "pgsch %f\n", size);
  }
}
void ppgqcs(int units, float *xch, float *ych)
{
  if(ppgdevice_internal == NULL) {
    cpgqcs(units, xch, ych);
    return;
  }else {
    fflush(stdout);
    fprintf(stderr, "ERROR ppgqcs not defined in pgplotinterpretor mode\n");
    return;
  }
}
void ppgscf(int font)
{
  if(ppgdevice_internal == NULL)
    cpgscf(font);
  else {
    fprintf(ppgdevice_internal, "pgscf %d\n", font);
  }
}
void ppgshs(float angle, float sepn, float phase)
{
  if(ppgdevice_internal == NULL)
    cpgshs(angle, sepn, phase);
  else {
    fprintf(ppgdevice_internal, "pgshs %f %f %f\n", angle, sepn, phase);
  }
}
void ppgarro(float x1, float y1, float x2, float y2)
{
  if(ppgdevice_internal == NULL)
    cpgarro(x1, y1, x2, y2);
  else {
    fprintf(ppgdevice_internal, "pgarro %f %f %f %f\n", x1, y1, x2, y2);
  }
}
void ppgmove(float x, float y)
{
  if(ppgdevice_internal == NULL)
    cpgmove(x, y);
  else {
    fprintf(ppgdevice_internal, "pgmove %f %f\n", x, y);
  }
}
void ppgdraw(float x, float y)
{
  if(ppgdevice_internal == NULL)
    cpgdraw(x, y);
  else {
    fprintf(ppgdevice_internal, "pgdraw %f %f\n", x, y);
  }
}
void ppgpt1(float xpt, float ypt, int symbol)
{
  if(ppgdevice_internal == NULL)
    cpgpt1(xpt, ypt, symbol);
  else {
    fprintf(ppgdevice_internal, "pgpt1 %f %f %d\n", xpt, ypt, symbol);
  }
}
void ppgerr1(int dir, float x, float y, float e, float t)
{
  if(ppgdevice_internal == NULL)
    cpgerr1(dir, x, y, e, t);
  else {
    fprintf(ppgdevice_internal, "pgerr1 %d %f %f %f %f\n", dir, x, y, e, t);
  }
}
void ppgcirc(float xcent, float ycent, float radius)
{
  if(ppgdevice_internal == NULL)
    cpgcirc(xcent, ycent, radius);
  else {
    fprintf(ppgdevice_internal, "pgcirc %f %f %f\n", xcent, ycent, radius);
  }
}
void ppgtext(float x, float y, const char *text)
{
  if(ppgdevice_internal == NULL)
    cpgtext(x, y, text);
  else {
    fprintf(ppgdevice_internal, "pgtext %f %f '%s'\n", x, y, text);
  }
}
void ppgptxt(float x, float y, float angle, float fjust, const char *text)
{
  if(ppgdevice_internal == NULL)
    cpgptxt(x, y, angle, fjust, text);
  else {
    fprintf(ppgdevice_internal, "pgptxt %f %f %f %f '%s'\n", x, y, angle, fjust, text);
  }
}
void ppgrect(float x1, float x2, float y1, float y2)
{
  if(ppgdevice_internal == NULL)
    cpgrect(x1, x2, y1, y2);
  else {
    fprintf(ppgdevice_internal, "pgrect %f %f %f %f\n", x1, x2, y1, y2);
  }
}
void ppgimag(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float a1, float a2, const float *tr)
{
  FILE *fout;
  char filename[2000];
  long ret;
  if(ppgdevice_internal == NULL)
    cpgimag(a, idim, jdim, i1, i2, j1, j2, a1, a2, tr);
  else {
    sprintf(filename, "%s.aux%03d", ppgdevicename_internal, ppg_nextaux_internal++);
    fout = fopen(filename, "wb");
    if(fout == NULL) {
      fflush(stdout);
      fprintf(stderr, "ERROR ppgimag: Cannot open %s\n", filename);
      exit(0);
    }
    ret = fwrite(a, sizeof(float), idim*jdim, fout);
    if(ret != idim*jdim) {
      fflush(stdout);
      fprintf(stderr, "ERROR ppgimag: Write error to %s (%ld != %ld)\n", filename, ret, (long)idim*jdim);
      exit(0);
    }
    fclose(fout);
    fprintf(ppgdevice_internal, "pgimag %s %d %d %d %d %d %d %f %f %f %f %f %f %f %f\n", filename, idim, jdim, i1, i2, j1, j2, a1, a2, tr[0], tr[1], tr[2], tr[3], tr[4], tr[5]);
  }
}
void ppgcont(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, const float *c, int nc, const float *tr)
{
  FILE *fout;
  char filename[2000];
  long ret;
  int i;
  if(ppgdevice_internal == NULL)
    cpgcont(a, idim, jdim, i1, i2, j1, j2, c, nc, tr);
  else {
    sprintf(filename, "%s.aux%03d", ppgdevicename_internal, ppg_nextaux_internal++);
    fout = fopen(filename, "wb");
    if(fout == NULL) {
      fflush(stdout);
      fprintf(stderr, "ERROR ppgcont: Cannot open %s\n", filename);
      exit(0);
    }
    ret = fwrite(a, sizeof(float), idim*jdim, fout);
    if(ret != idim*jdim) {
      fflush(stdout);
      fprintf(stderr, "ERROR ppgcont: Write error to %s (%ld != %ld)\n", filename, ret, (long)idim*jdim);
      exit(0);
    }
    fclose(fout);
    fprintf(ppgdevice_internal, "pgcont %s %d %d %d %d %d %d %d", filename, idim, jdim, i1, i2, j1, j2, nc);
    if(nc < 0)
      nc *= -1;
    for(i = 0; i < nc; i++)
      fprintf(ppgdevice_internal, " %f", c[i]);
    fprintf(ppgdevice_internal, " %f %f %f %f %f %f\n", tr[0], tr[1], tr[2], tr[3], tr[4], tr[5]);
  }
}
void ppgconl(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float c, const float *tr, const char *label, int intval, int minint)
{
  FILE *fout;
  char filename[2000];
  long ret;
  if(ppgdevice_internal == NULL)
    cpgconl(a, idim, jdim, i1, i2, j1, j2, c, tr, label, intval, minint);
  else {
    sprintf(filename, "%s.aux%03d", ppgdevicename_internal, ppg_nextaux_internal++);
    fout = fopen(filename, "wb");
    if(fout == NULL) {
      fflush(stdout);
      fprintf(stderr, "ERROR ppgconl: Cannot open %s\n", filename);
      exit(0);
    }
    ret = fwrite(a, sizeof(float), idim*jdim, fout);
    if(ret != idim*jdim) {
      fflush(stdout);
      fprintf(stderr, "ERROR ppgconl: Write error to %s (%ld != %ld)\n", filename, ret, (long)idim*jdim);
      exit(0);
    }
    fclose(fout);
    fprintf(ppgdevice_internal, "pgconl %s %d %d %d %d %d %d %f %f %f %f %f %f %f %d %d '%s'\n", filename, idim, jdim, i1, i2, j1, j2, c, tr[0], tr[1], tr[2], tr[3], tr[4], tr[5], intval, minint, label);
  }
}
void ppgpoly(int n, const float *xpts, const float *ypts)
{
  char filename[2000];
  long ret;
  FILE *fout;
  if(ppgdevice_internal == NULL)
    cpgpoly(n, xpts, ypts);
  else {
    sprintf(filename, "%s.aux%03d", ppgdevicename_internal, ppg_nextaux_internal++);
    fout = fopen(filename, "wb");
    if(fout == NULL) {
      fflush(stdout);
      fprintf(stderr, "ERROR ppgpoly: Cannot open %s\n", filename);
      exit(0);
    }
    ret = fwrite(xpts, sizeof(float), n, fout);
    if(ret != n) {
      fflush(stdout);
      fprintf(stderr, "ERROR ppgpoly: Write error to %s (%ld != %ld)\n", filename, ret, (long)n);
      exit(0);
    }
    ret = fwrite(ypts, sizeof(float), n, fout);
    if(ret != n) {
      fflush(stdout);
      fprintf(stderr, "ERROR ppgpoly: Write error to %s (%ld != %ld)\n", filename, ret, (long)n);
      exit(0);
    }
    fclose(fout);
    fprintf(ppgdevice_internal, "pgpoly %d %s\n", n, filename);
  }
}
int ppgcurs(float *x, float *y, char *ch)
{
  fflush(stdout);
  if(ppgdevice_internal == NULL)
    return cpgcurs(x, y, ch);
  else {
    fflush(stdout);
    fprintf(stderr, "ERROR ppgcurs not defined in pgplotinterpretor mode\n");
    return 0;
  }
}
int ppgband(int mode, int posn, float xref, float yref, float *x, float *y, char *ch_scalar)
{
  fflush(stdout);
  if(ppgdevice_internal == NULL) {
    return cpgband(mode, posn, xref, yref, x, y, ch_scalar);
  }else {
    fflush(stdout);
    fprintf(stderr, "ERROR ppgband not defined in pgplotinterpretor mode\n");
    return 0;
  }
}
void ppgscir(int icilo, int icihi)
{
  if(ppgdevice_internal == NULL)
    cpgscir(icilo, icihi);
  else {
    fprintf(ppgdevice_internal, "pgscir %d %d\n", icilo, icihi);
  }
}
void ppgctab(float *l, float *r, float *g, float *b, int nc, float contra, float bright)
{
  FILE *fout;
  char filename[2000];
  long ret;
  if(ppgdevice_internal == NULL)
    cpgctab(l, r, g, b, nc, contra, bright);
  else {
    sprintf(filename, "%s.aux%03d", ppgdevicename_internal, ppg_nextaux_internal++);
    fout = fopen(filename, "wb");
    if(fout == NULL) {
      fflush(stdout);
      fprintf(stderr, "ERROR ppgctab: Cannot open %s\n", filename);
      exit(0);
    }
    ret = fwrite(l, sizeof(float), nc, fout);
    ret += fwrite(r, sizeof(float), nc, fout);
    ret += fwrite(g, sizeof(float), nc, fout);
    ret += fwrite(b, sizeof(float), nc, fout);
    if(ret != 4*nc) {
      fflush(stdout);
      fprintf(stderr, "ERROR ppgctab: Write error to %s (%ld != %ld)\n", filename, ret, (long)4*nc);
      exit(0);
    }
    fclose(fout);
    fprintf(ppgdevice_internal, "pgctab %s %d %f %f\n", filename, nc, contra, bright);
  }
}
void ppgwedg(const char *side, float disp, float width, float fg, float bg, const char *label)
{
  if(ppgdevice_internal == NULL)
    cpgwedg(side, disp, width, fg, bg, label);
  else {
    fprintf(ppgdevice_internal, "pgwedg %s %f %f %f %f '%s'\n", side, disp, width, fg, bg, label);
  }
}
int ppgqid(int *id)
{
  if(ppgdevice_internal == NULL) {
    cpgqid(id);
    return 1;
  }else {
    fflush(stdout);
    fprintf(stderr, "ERROR ppgqid not defined in pgplotinterpretor mode\n");
    return 0;
  }
}
int ppgslct(int id)
{
  if(ppgdevice_internal == NULL) {
    cpgslct(id);
    return 1;
  }else {
    fflush(stdout);
    fprintf(stderr, "ERROR ppgslct not defined in pgplotinterpretor mode\n");
    return 0;
  }
}
void ppgsclp(int state)
{
  if(ppgdevice_internal == NULL) {
    cpgsclp(state);
    return;
  }else {
    fprintf(ppgdevice_internal, "pgsclp %d\n", state);
  }
}

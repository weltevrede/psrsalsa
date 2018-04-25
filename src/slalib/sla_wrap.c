/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <ctype.h>
extern void sla_altaz(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
void csla_altaz(double ha, double dec, double phi, double *az, double *azd, double *azdd, double *el, double *eld, double *eldd, double *pa, double *pad, double *padd)
{
  sla_altaz(&ha, &dec, &phi, az, azd, azdd, el, eld, eldd, pa, pad, padd);
}
extern void sla_dcc2s(double *v, double *ra, double *dec);
void csla_dcc2s(double v[3], double *a, double *b)
{
  sla_dcc2s(v, a, b);
}
extern void sla_dcs2c(double *ra, double *dec, double *v);
void csla_dcs2c(double a, double b, double v[3])
{
  sla_dcs2c(&a, &b, v);
}
extern void sla_dmxv(double *dm, double *va, double *vb);
void csla_dmxv(double dm[3][3], double va[3], double vb[3])
{
  sla_dmxv(&(dm[0][0]), va, vb);
}
extern double sla_dranrm(double *);
double csla_dranrm(double angle)
{
  return sla_dranrm(&angle);
}
extern double sla_dvdv(double *va, double *vb);
double csla_dvdv(double va[3], double vb[3])
{
  return sla_dvdv(va, vb);
}
extern double sla_dvn(double *v, double *uv, double *vm);
double csla_dvn(double v[3], double uv[3], double *vm)
{
  return sla_dvn(v, uv, vm);
}
extern double sla_epb(double *mjd);
double csla_epb(double date)
{
  return sla_epb(&date);
}
extern double sla_epj(double *mjd);
double csla_epj(double date)
{
  return sla_epj(&date);
}
extern void sla_evp(double *tdb, double *ep, double *dvb, double *dpb, double *dvh, double *dph);
void csla_evp(double date, double deqx, double dvb[3], double dpb[3], double dvh[3], double dph[3])
{
  sla_evp(&date, &deqx, dvb, dpb, dvh, dph);
}
extern double sla_gmst(double *mjd);
double csla_gmst(double ut1)
{
  return sla_gmst(&ut1);
}
extern void sla_nut(double *date, double *rmatn);
void csla_nut(double date, double rmatn[3][3])
{
  sla_nut(&date, &(rmatn[0][0]));
}
extern double sla_pa(double *HA, double *DEC, double *PHI);
double csla_pa(double ha, double dec, double phi)
{
  return sla_pa(&ha, &dec, &phi);
}
extern void sla_prebn(double *ep0, double *ep1, double *rmatp);
void csla_prebn(double bep0, double bep1, double rmatp[3][3])
{
  sla_prebn(&bep0, &bep1, &(rmatp[0][0]));
}
extern void sla_prec(double *ep0, double *ep1, double *rmatp);
void csla_prec(double ep0, double ep1, double rmatp[3][3])
{
  sla_prec(&ep0, &ep1, &(rmatp[0][0]));
}
extern void sla_prenut(double *epoch, double *date, double *rmatpn);
void csla_prenut(double ep0, double date, double rmatpn[3][3])
{
  sla_prenut(&ep0, &date, &(rmatpn[0][0]));
}

/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

void csla_altaz(double ha, double dec, double phi, double *az, double *azd, double *azdd, double *el, double *eld, double *eldd, double *pa, double *pad, double *padd);
void csla_dcc2s(double v[3], double *a, double *b);
void csla_dcs2c(double a, double b, double v[3]);
void csla_dmxv(double dm[3][3], double va[3], double vb[3] );
double csla_dranrm(double angle);
double csla_dvdv(double va[3], double vb[3]);
double csla_dvn(double v[3], double uv[3], double *vm);
double csla_epb(double date);
double csla_epj(double date );
void csla_evp(double date, double deqx, double dvb[3], double dpb[3], double dvh[3], double dph[3]);
double csla_gmst(double ut1);
void csla_nut(double date, double rmatn[3][3]);
double csla_pa(double ha, double dec, double phi);
void csla_prebn(double bep0, double bep1, double rmatp[3][3]);
void csla_prec(double ep0, double ep1, double rmatp[3][3]);
void csla_prenut(double ep0, double date, double rmatpn[3][3]);

/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <math.h>
float derotate_deg(float a)
{
  int i;
  i = fabs(a)/360.0;
  if(a > 0)
    a -= 360*i;
  else
    a += 360*i;
  if(a < 0)
    a += 360;
  return a;
}
float derotate_180(float a)
{
  int i;
  i = fabs(a)/180.0;
  if(a > 0)
    a -= 180*i;
  else
    a += 180*i;
  if(a < 0)
    a += 180;
  return a;
}
double derotate_180_double(double a)
{
  int i;
  i = fabs(a)/180.0;
  if(a > 0)
    a -= 180*i;
  else
    a += 180*i;
  if(a < 0)
    a += 180;
  return a;
}
double derotate_180_rad_double(double a)
{
  int i;
  i = fabs(a)/M_PI;
  if(a > 0)
    a -= M_PI*i;
  else
    a += M_PI*i;
  if(a < 0)
    a += M_PI;
  return a;
}
float derotate_90(float a)
{
  int i;
  i = fabs(a)/90.0;
  if(a > 0)
    a -= 90*i;
  else
    a += 90*i;
  if(a < 0)
    a += 90;
  return a;
}
double derotate_90_double(double a)
{
  int i;
  i = fabs(a)/90.0;
  if(a > 0)
    a -= 90*i;
  else
    a += 90*i;
  if(a < 0)
    a += 90;
  return a;
}
double derotate_180_small_double(double a)
{
  double y;
  y = derotate_180_double(a);
  if(y >= 90.0)
    y = y-180.0;
  return y;
}
float polar_angle_rad(float x, float y)
{
  float alpha;
  if(x == 0) {
    if(y <= 0)
      return 1.5*M_PI;
    else
      return 0.5*M_PI;
  }else {
    alpha = atan(y/x);
    if(x > 0 && y >= 0)
      return alpha;
    if(x > 0)
      return alpha + 2*M_PI;
    return alpha + M_PI;
  }
}
float paswing(float alpha, float beta, float l, float pa0, float l0, int nrJumps, float *jump_longitude, float *jump_offset, float add_height_longitude, float add_height_shift)
{
  float x1, y1, sa, dl, pa, dbcw, dha;
  int i;
  alpha *= M_PI/180.0;
  beta *= M_PI/180.0;
  if(l >= add_height_longitude) {
    dbcw = 2.0*add_height_shift*180.0/M_PI;
    dha = 10.0*add_height_shift*cos(alpha)*180.0/(3.0*M_PI);
  }else {
    dbcw = 0;
    dha = 0;
  }
  sa = sin(alpha);
  dl = (l-(l0+dbcw))*M_PI/180.0;
  y1 = sa*sin(dl);
  x1 = sin(alpha+beta)*cos(alpha);
  x1 -= cos(alpha+beta)*sa*cos(dl);
  pa = pa0+dha+atan2(y1,x1)*180.0/M_PI;
  if(nrJumps > 0) {
    for(i = 0; i < nrJumps; i++)
      if(l > jump_longitude[i])
 pa += jump_offset[i];
  }
  pa = derotate_180(pa);
  return pa;
}
double paswing_double(double alpha, double beta, double l, double pa0, double l0, int nrJumps, double *jump_longitude, double *jump_offset, double add_height_longitude, double add_height_shift)
{
  double x1, y1, sa, dl, pa, dbcw, dha;
  int i;
  alpha *= M_PI/180.0;
  beta *= M_PI/180.0;
  if(l >= add_height_longitude) {
    dbcw = 2.0*add_height_shift*180.0/M_PI;
    dha = 10.0*add_height_shift*cos(alpha)*180.0/(3.0*M_PI);
  }else {
    dbcw = 0;
    dha = 0;
  }
  sa = sin(alpha);
  dl = (l-(l0+dbcw))*M_PI/180.0;
  y1 = sa*sin(dl);
  x1 = sin(alpha+beta)*cos(alpha);
  x1 -= cos(alpha+beta)*sa*cos(dl);
  pa = pa0+dha+atan2(y1,x1)*180.0/M_PI;
  if(nrJumps > 0) {
    for(i = 0; i < nrJumps; i++)
      if(l > jump_longitude[i])
 pa += jump_offset[i];
  }
  pa = derotate_180_double(pa);
  return pa;
}

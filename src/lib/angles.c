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
double derotate_90_small_double(double a)
{
  double y;
  y = derotate_90_double(a);
  if(y >= 45.0)
    y = y-90.0;
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
void rotateZ_d(long nrpoints, double angle, double *x, double *y, double *z)
{
  long n;
  double x2, y2;
  angle *= M_PI/180.0;
  for(n = 0; n < nrpoints; n++) {
    x2 = x[n]*cos(angle);
    y2 = x[n]*sin(angle);
    x2 -= y[n]*sin(angle);
    y2 += y[n]*cos(angle);
    x[n] = x2;
    y[n] = y2;
  }
}
void rotateY(long nrpoints, float angle, float *x, float *y, float *z)
{
  long n;
  float x2, z2;
  angle *= M_PI/180.0;
  for(n = 0; n < nrpoints; n++) {
    x2 = x[n]*cos(angle);
    z2 = -x[n]*sin(angle);
    x2 += z[n]*sin(angle);
    z2 += z[n]*cos(angle);
    x[n] = x2;
    z[n] = z2;
  }
}
void rotateY_d(long nrpoints, double angle, double *x, double *y, double *z)
{
  long n;
  double x2, z2;
  angle *= M_PI/180.0;
  for(n = 0; n < nrpoints; n++) {
    x2 = x[n]*cos(angle);
    z2 = -x[n]*sin(angle);
    x2 += z[n]*sin(angle);
    z2 += z[n]*cos(angle);
    x[n] = x2;
    z[n] = z2;
  }
}
void rotateX(long nrpoints, float angle, float *x, float *y, float *z)
{
  long n;
  float y2, z2;
  angle *= M_PI/180.0;
  for(n = 0; n < nrpoints; n++) {
    y2 = y[n]*cos(angle);
    z2 = -y[n]*sin(angle);
    y2 += z[n]*sin(angle);
    z2 += z[n]*cos(angle);
    y[n] = y2;
    z[n] = z2;
  }
}
float paswing(float alpha, float beta, float l, float pa0, float l0, int nrJumps, float *jump_longitude, float *jump_offset, float add_height_longitude, float add_height_shift, int add_height_shift_bcw_only)
{
  float x1, y1, sa, dl, pa, dbcw, dha;
  int i;
  alpha *= M_PI/180.0;
  beta *= M_PI/180.0;
  if(l >= add_height_longitude) {
    dbcw = 2.0*add_height_shift*180.0/M_PI;
    if(add_height_shift_bcw_only == 0) {
      dha = 10.0*add_height_shift*cos(alpha)*180.0/(3.0*M_PI);
    }else {
      dha = 0;
    }
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
double paswing_double(double alpha, double beta, double l, double pa0, double l0, int nrJumps, double *jump_longitude, double *jump_offset, double add_height_longitude, double add_height_shift, int add_height_shift_bcw_only)
{
  double x1, y1, sa, dl, pa, dbcw, dha;
  int i;
  alpha *= M_PI/180.0;
  beta *= M_PI/180.0;
  if(l >= add_height_longitude) {
    dbcw = 2.0*add_height_shift*180.0/M_PI;
    if(add_height_shift_bcw_only == 0) {
      dha = 10.0*add_height_shift*cos(alpha)*180.0/(3.0*M_PI);
    }else {
      dha = 0;
    }
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
void spherical2Cartesian(double r, double theta, double phi, double *x, double *y, double *z)
{
  *x = r*sin(theta)*cos(phi);
  *y = r*sin(theta)*sin(phi);
  *z = r*cos(theta);
}
void cartesian2spherical(double *r, double *theta, double *phi, double x, double y, double z)
{
  float r2;
  *r = sqrt(x*x+y*y+z*z);
  r2 = sqrt(x*x+y*y);
  *theta = atan2(r2, z);
  *phi = atan2(y, x);
}
void projectionHammerAitoff_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y)
{
  float f;
  double x2, y2, z2, theta, phi, r;
  dlongitude *= 180.0/M_PI;
  dlatitude *= 180.0/M_PI;
  spherical2Cartesian(1, latitude+0.5*M_PI, longitude, &x2, &y2, &z2);
  rotateZ_d(1, dlongitude, &x2, &y2, &z2);
  rotateY_d(1, -dlatitude, &x2, &y2, &z2);
  cartesian2spherical(&r, &theta, &phi, x2, y2, z2);
  latitude = theta - 0.5*M_PI;
  longitude = phi;
  f = 1.0/sqrt(1.0+cos(latitude)*cos(0.5*longitude));
  *x = f*2.0*sqrt(2.0)*cos(latitude)*sin(0.5*longitude);
  *y = sqrt(2.0)*sin(latitude)*f;
  *x *= 0.5*sqrt(2);
  *y *= 0.5*sqrt(2);
}
int projection_sphere_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y, float *weight)
{
  float xE, yE, zE;
  int side;
  dlongitude *= 180.0/M_PI;
  dlatitude *= 180.0/M_PI;
  yE = sin(latitude);
  xE = cos(latitude)*sin(longitude);
  zE = -cos(latitude)*cos(longitude);
  rotateY(1, -dlongitude, &xE, &yE, &zE);
  rotateX(1, dlatitude, &xE, &yE, &zE);
  *weight = fabs(zE);
  if(zE > 0)
    side = 1;
  else
    side = 0;
  *x = xE;
  *y = yE;
  if(side == 0)
    *x -= 1.125;
  else
    *x = (-*x)+1.125;
  return side;
}
int projection_sphere_longlat(float x, float y, float dlongitude, float dlatitude, float *longitude, float *latitude)
{
  float z, r2;
  int side;
  if(x < 0) {
    x += 1.125;
    side = 0;
  }else {
    x = -(x - 1.125);
    side = 1;
  }
  if(x*x+y*y > 1)
    return 0;
  if(side == 0)
    z = -sqrt(1-x*x-y*y);
  else
    z = sqrt(1-x*x-y*y);
  rotateX(1, -dlatitude*180.0/M_PI, &x, &y, &z);
  rotateY(1, dlongitude*180.0/M_PI, &x, &y, &z);
  r2 = sqrt(z*z+x*x);
  *latitude = atan2(r2, y);
  *latitude = 0.5*M_PI-(*latitude);
  *longitude = atan2(x, -z);
  return 1;
}
void projection_longlat_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y)
{
  float xE, yE, zE, r2;
  dlongitude *= 180.0/M_PI;
  dlatitude *= 180.0/M_PI;
  yE = sin(latitude);
  xE = cos(latitude)*sin(longitude);
  zE = -cos(latitude)*cos(longitude);
  rotateX(1, dlatitude, &xE, &yE, &zE);
  rotateY(1, -dlongitude, &xE, &yE, &zE);
  r2 = sqrt(zE*zE+xE*xE);
  *y = atan2(r2, yE);
  *y = 0.5*M_PI-(*y);
  *y *= 180.0/M_PI;
  *x = atan2(xE, -zE)*180.0/M_PI;
}

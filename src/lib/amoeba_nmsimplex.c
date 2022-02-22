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
#include <math.h>
#include <stdio.h>
#define MAX_IT 5000
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0
float amoeba_nmsimplex(float (*objfunc)(float[]), float start[], float dx[], int n, float EPSILON, int *nfunc, float *reachedEpsilon, int verbose)
{
  int vs;
  int vh;
  int vg;
  int i,j,m,row;
  int k;
  float **v;
  float *f;
  float fr;
  float fe;
  float fc;
  float *vr;
  float *ve;
  float *vc;
  float *vm;
  float min;
  float s,cent;
  v = (float **) malloc ((n+1) * sizeof(float *));
  f = (float *) malloc ((n+1) * sizeof(float));
  vr = (float *) malloc (n * sizeof(float));
  ve = (float *) malloc (n * sizeof(float));
  vc = (float *) malloc (n * sizeof(float));
  vm = (float *) malloc (n * sizeof(float));
  for (i=0;i<=n;i++) {
    v[i] = (float *) malloc (n * sizeof(float));
  }
  for (i=0;i<n;i++) {
    v[0][i] = start[i];
  }
  for (i=1;i<=n;i++) {
    for (j=0;j<n;j++) {
      if (i-1 == j) {
 v[i][j] = start[j]+dx[j];
      }
      else {
 v[i][j] = start[j];
      }
    }
  }
  for (j=0;j<=n;j++) {
    f[j] = objfunc(v[j]);
  }
  k = n+1;
  if(verbose) {
    printf("Initial Values\n");
    for (j=0;j<=n;j++) {
      printf("Vertex %d:", j+1);
      for (i=0;i<n;i++) {
 printf(" %f",v[j][i]);
      }
      printf(": %f\n",f[j]);
    }
  }
  int nritterations;
  for (nritterations=1;nritterations<=MAX_IT;(nritterations)++) {
    vg=0;
    for (j=0;j<=n;j++) {
      if (f[j] > f[vg]) {
 vg = j;
      }
    }
    vs=0;
    for (j=0;j<=n;j++) {
      if (f[j] < f[vs]) {
 vs = j;
      }
    }
    vh=vs;
    for (j=0;j<=n;j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
 vh = j;
      }
    }
    for (j=0;j<=n-1;j++) {
      cent=0.0;
      for (m=0;m<=n;m++) {
 if (m!=vg) {
   cent += v[m][j];
 }
      }
      vm[j] = cent/n;
    }
    for (j=0;j<=n-1;j++) {
      vr[j] = vm[j]+ALPHA*(vm[j]-v[vg][j]);
    }
    fr = objfunc(vr);
    k++;
    if (fr < f[vh] && fr >= f[vs]) {
      for (j=0;j<=n-1;j++) {
 v[vg][j] = vr[j];
      }
      f[vg] = fr;
    }
    if ( fr < f[vs]) {
      for (j=0;j<=n-1;j++) {
 ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
      }
      fe = objfunc(ve);
      k++;
      if (fe < fr) {
 for (j=0;j<=n-1;j++) {
   v[vg][j] = ve[j];
 }
 f[vg] = fe;
      }
      else {
 for (j=0;j<=n-1;j++) {
   v[vg][j] = vr[j];
 }
 f[vg] = fr;
      }
    }
    if (fr >= f[vh]) {
      if (fr < f[vg] && fr >= f[vh]) {
 for (j=0;j<=n-1;j++) {
   vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
 }
 fc = objfunc(vc);
 k++;
      }
      else {
 for (j=0;j<=n-1;j++) {
   vc[j] = vm[j]-BETA*(vm[j]-v[vg][j]);
 }
 fc = objfunc(vc);
 k++;
      }
      if (fc < f[vg]) {
 for (j=0;j<=n-1;j++) {
   v[vg][j] = vc[j];
 }
 f[vg] = fc;
      }
      else {
 for (row=0;row<=n;row++) {
   if (row != vs) {
     for (j=0;j<=n-1;j++) {
       v[row][j] = v[vs][j]+(v[row][j]-v[vs][j])/2.0;
     }
   }
 }
 f[vg] = objfunc(v[vg]);
 k++;
 f[vh] = objfunc(v[vh]);
 k++;
      }
    }
    if(verbose) {
      printf("Iteration %d\n",nritterations);
      for (j=0;j<=n;j++) {
 printf("Vertex %d:", j+1);
 for (i=0;i<n;i++) {
   printf(" %f",v[j][i]);
 }
 printf(": %f\n",f[j]);
      }
    }
    float yhi, ylo;
    yhi = ylo = f[0];
    for (j=0;j<=n;j++) {
      if(f[j] > yhi)
 yhi = f[j];
      if(f[j] < ylo)
 ylo = f[j];
    }
    s = 2*fabs(yhi-ylo)/(fabs(yhi)+fabs(ylo));
    if (s < EPSILON) break;
  }
  *reachedEpsilon = s;
  vs=0;
  for (j=0;j<=n;j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  if(verbose) {
    printf("The minimum was found at epsilon=%e:\n  ", *reachedEpsilon);
    for (j=0;j<n;j++) {
      printf(" %e",v[vs][j]);
    }
    printf("\n");
  }
  for (j=0;j<n;j++) {
    start[j] = v[vs][j];
  }
  min=objfunc(v[vs]);
  k++;
  if(verbose) {
    printf("%d Function Evaluations\n",k);
    printf("%d Iterations through program\n",nritterations);
  }
  *nfunc = k;
  free(f);
  free(vr);
  free(ve);
  free(vc);
  free(vm);
  for (i=0;i<=n;i++) {
    free (v[i]);
  }
  free(v);
  return min;
}

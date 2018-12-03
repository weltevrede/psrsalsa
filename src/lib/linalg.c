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
#include <stdlib.h>
#include <string.h>
#include "psrsalsa.h"
#include <gsl/gsl_linalg.h>
int linalg_solve_matrix_eq_gauss_jordan(double *matrixa, double *matrixb, int n, int m, verbose_definition verbose)
{
  int *orig_pivot_row, *orig_pivot_column, *pivot_column_done;
  int pivotnr, pivot_col, pivot_row, rowi, coli, i;
  double *ptr1, *ptr2, tmpvalue;
  double largest_element_value;
  int currow;
  orig_pivot_row = malloc(n*sizeof(int));
  orig_pivot_column = malloc(n*sizeof(int));
  pivot_column_done = malloc(n*sizeof(int));
  if(orig_pivot_row == NULL || orig_pivot_column == NULL || pivot_column_done == NULL) {
    printerror(verbose.debug, "ERROR linalg_solve_matrix_eq_gauss_jordan: Memory allocation error.");
    return 1;
  }
  for(coli = 0; coli < n; coli++)
    pivot_column_done[coli] = 0;
  for(pivotnr = 0; pivotnr < n; pivotnr++) {
    largest_element_value = -1.0;
    for(rowi = 0; rowi < n; rowi++) {
      if(pivot_column_done[rowi] != 1) {
 for(coli = 0; coli < n; coli++) {
   if(pivot_column_done[coli] == 0) {
     tmpvalue = matrixa[rowi*n+coli];
     if(fabs(tmpvalue) >= largest_element_value) {
       largest_element_value = fabs(tmpvalue);
       pivot_row=rowi;
       pivot_col=coli;
     }
   }else if(pivot_column_done[coli] > 1) {
     printerror(verbose.debug, "ERROR linalg_solve_matrix_eq_gauss_jordan: The matrix equation is singular, no solution can be determined.");
     return 2;
   }
 }
      }
    }
    if(largest_element_value <= 0.0) {
      printerror(verbose.debug, "ERROR linalg_solve_matrix_eq_gauss_jordan: The matrix equation is singular, no solution can be determined.");
      return 2;
    }
    pivot_column_done[pivot_col] += 1;
    if(pivot_row != pivot_col) {
      ptr1 = matrixa + pivot_row*n;
      ptr2 = matrixa + pivot_col*n;
      for(i = 0; i < n; i++) {
 tmpvalue = *ptr1;
 *ptr1 = *ptr2;
 *ptr2 = tmpvalue;
      }
      ptr1 = matrixb + pivot_row*n;
      ptr2 = matrixb + pivot_col*n;
      for(i = 0; i < n; i++) {
 tmpvalue = *ptr1;
 *ptr1 = *ptr2;
 *ptr2 = tmpvalue;
      }
    }
    orig_pivot_row[pivotnr] = pivot_row;
    orig_pivot_column[pivotnr] = pivot_col;
    matrixa[pivot_col*n+pivot_col]=1.0;
    ptr1 = matrixa + pivot_col*n;
    for(i = 0; i < n; i++) {
      *ptr1 /= largest_element_value;
      ptr1 += 1;
    }
    ptr1 = matrixb + pivot_col*m;
    for(i = 0; i < m; i++) {
      *ptr1 /= largest_element_value;
      ptr1 += 1;
    }
    for(currow = 0; currow < n; currow++) {
      if(currow != pivot_col) {
 i = currow*n+pivot_col;
 tmpvalue = matrixa[i];
 matrixa[i] = 0.0;
 ptr1 = matrixa + n*currow;
 ptr2 = matrixa + n*pivot_col;
 for(i=1; i <= n; i++) {
   *ptr1 -= (*ptr2)*tmpvalue;
   ptr1++;
   ptr2++;
 }
 ptr1 = matrixb + currow*m;
 ptr2 = matrixb + pivot_col*m;
 for(i = 0; i < m; i++) {
   *ptr1 -= (*ptr2)*tmpvalue;
   ptr1++;
   ptr2++;
 }
      }
    }
  }
  for(coli = n-1; coli >= 0; coli--) {
    if(orig_pivot_row[coli] != orig_pivot_column[coli]) {
      for(rowi = 0; rowi <= n; rowi++) {
 tmpvalue = matrixa[rowi*n+orig_pivot_row[coli]];
 matrixa[rowi*n+orig_pivot_row[coli]] = matrixa[rowi*n+orig_pivot_column[coli]];
 matrixa[rowi*n+orig_pivot_column[coli]] = tmpvalue;
      }
    }
  }
  free(pivot_column_done);
  free(orig_pivot_row);
  free(orig_pivot_column);
  return 0;
}

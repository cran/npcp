/*#################################################################################
  ##
  ##   R package npcp by Ivan Kojadinovic Copyright (C) 2014
  ##
  ##   This file is part of the R package npcp.
  ##
  ##   The R package npcp is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 3 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package npcp is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with the R package npcp. If not, see <http://www.gnu.org/licenses/>.
  ##
  #################################################################################*/

/***********************************************************************

  Some utility functions

***********************************************************************/

#ifndef UTILITIES_H
#define UTILITIES_H

#define	MIN(x,y) ((x) > (y) ? (y) : (x))
#define	MAX(x,y) ((x) < (y) ? (y) : (x))

double parzen(double x);
double convrect(double x, int n);
double div_diff_xn_1_y_plus(int n, double *a, double y);
void pdf_sum_unif(int *n, double *y, int *ny, double *pdfy);
void gendepmult(int n, int M, int bw, int w, double *initseq, double *multipliers);
void makepseudoobs(double *x, int *index, int n, int d, int b, int e, double *V);
double zeroin(double ax,double bx,double (*f)(), double tol); // see zeroin.c

#endif

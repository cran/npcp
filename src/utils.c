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

#include <R.h>

/*****************************************************************************

  Computation of the divided difference [(x-y)^(n-1)_{+} : a_0,...,a_n]

*****************************************************************************/

double div_diff_xn_1_y_plus(int n, double *a, double y) {

  int i, j, k, r=0, s=0;
  double *dd, *b, *c, out;
  for (i=0; i<=n; i++)
      {
	  if (a[i] < y)
	      r++;
	  else
	      s++;
      }

  if (r==0 || s==0)
      return 0.0;

  b = Calloc(r, double);
  c = Calloc(s, double);

  j=0;
  k=0;
  for (i=0; i<=n; i++)
      if (a[i] < y)
	  b[j++]=a[i]-y;
      else
	  c[k++]=a[i]-y;

  dd = Calloc(s + 1, double);

  /* Initialize dd */
  dd[0] = 0.0;
  dd[1] = 1.0 / (c[0] - b[0]);

  /* Computation by induction */
  for (j=2; j<=s; j++)
      dd[j] = - b[0] * dd[j-1]/(c[j-1] - b[0]);

  for (i=2; i<=r; i++)
      for (j=1; j<=s; j++)
	  dd[j] = (c[j-1] * dd[j] - b[i-1] * dd[j-1])/(c[j-1] - b[i-1]);

  out = dd[s];

  Free(b); Free(c); Free(dd);

  return out;
}

/*****************************************************************************

  Pdf of sum of standard uniform random variables

*****************************************************************************/

void pdf_sum_unif(int *n, double *y, int *ny, double *pdfy)
{
    double *a = Calloc(*n+1,double);

    for(int i=0; i<=*n; i++)
	a[i] = *n - i;

    for (int i=0; i <*ny; i++)
	pdfy[i] = div_diff_xn_1_y_plus(*n,a,y[i]) * (double)(*n);

    Free(a);
}

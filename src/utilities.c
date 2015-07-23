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
#include <Rmath.h>

/***********************************************************************

  Make pseudo obs from portion b:(e-1) of the data in x; result in V
  b: line beginning
  e: line end (+1)

***********************************************************************/

void makepseudoobs(double *x, int *index, int n, int d, int b, int e, double *V)
{
    int i, j;
    if (e == b)
	return;
    for (j = 0; j < d; j++)
	{
	    for (i = 0; i < e - b; i++)
		index[i] = i;
	    rsort_with_index (&x[j * n + b], index, e - b);
	    for (i = 0; i < e - b; i++)
		V[j * n + b + index[i]] = (i + 1.0) / (e - b + 1.0);
	}
}

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


/***********************************************************************

  Parzen kernel

***********************************************************************/

double parzen(double x)
{
    double fabsx = fabs(x);
    if (fabsx <= 0.5)
	return 1.0 - 6.0 * R_pow_di(x,2) + 6.0 * R_pow_di(fabsx,3);
    else if (0.5 < fabsx && fabsx <= 1)
	return 2.0 * R_pow_di(1 - fabsx,3);
    else return 0.0;
}

/*****************************************************************************

  (Not exactly) a kernel obtained by convoluting standard uniform random variables

*****************************************************************************/

double convrect(double x, int n)
{
    double *a = Calloc(n+1,double), res;

    for(int i=0; i<=n; i++)
	a[i] = n - i;

    res = div_diff_xn_1_y_plus(n, a, x + n / 2.0)
	/ div_diff_xn_1_y_plus(n, a, n / 2.0);

    Free(a);

    return res;
}

/***********************************************************************

  Generate dependent multipliers

***********************************************************************/

void gendepmult(int n, int M, int bw, int w, double *initseq, double *multipliers)
{
    int i, j, m;
    double *weights = Calloc(2 * bw - 1, double);
    double norm;

    /* Bartlett weights */
    if (w == 1)
	{
	    norm = sqrt(3.0 * R_pow_di(bw, 3) / (2.0 * R_pow_di(bw, 2) + 1.0));
	    for (j = -(bw - 1); j <= bw - 1; j++)
		weights[j + bw - 1] = (1.0 - fabs(j) / bw) / bw * norm;
	}
    else /* Parzen weights */
	{
	    norm = 0.0;
	    for (j = -(bw - 1); j <= bw - 1; j++)
		{
		    weights[j + bw - 1] = parzen((double)j / (double)bw);
		    norm += R_pow_di(weights[j + bw - 1], 2);
		}
	    for (j = -(bw - 1); j <= bw - 1; j++)
		weights[j + bw - 1] /= sqrt(norm);
	}

    /* dependent multipliers */
    for (m = 0; m < M; m++)
	for (i = 0; i < n; i++)
	    {
		multipliers[i + m * n] = 0.0;
		for (j = 0; j < 2 * bw - 1; j++)
		    multipliers[i + m * n] += weights[j]
			* initseq[i + j + m * (n + 2*(bw-1))];
	    }

    Free(weights);
}

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
#include "utilities.h"

/////////////////////////////////////////////////////////////////////////
// U-STAT TEST
/////////////////////////////////////////////////////////////////////////


/***********************************************************************

  Change-point tests based on U-statistics
  bw: set bw to 1 for the iid case

***********************************************************************/

void cpTestU(double *h, int *n, double *influ, double *stat, int *M,
	     int *w, int *bw, int *method, double *stat0, double *avar,
	     double *initseq)
{
    int i, j, k, m, ln;
    double *multipliers = Calloc((*n) * (*M), double);
    double *u1 = Calloc(*n - 3, double);
    double *u2 = Calloc(*n - 3, double);
    double s, un, ukm, unkm, sqrtn = sqrt(*n);

    /* generate (dependent) multipliers */
    if (*method == 1 || *method == 2)
	gendepmult(*n, *M, *bw, *w, initseq, multipliers);

    /* global U-statistic */
    if (*method == 2 || *method == 3)
	{
	    un = 0.0;
	    for (j = 0; j < *n; j++)
		for (i = 0; i < j; i++)
		    un += h[i + j * (*n)];
	    un /= (*n) * (*n - 1) / 2.0;
	}

    /* for each possible breakpoint */
    /* test statistics */
    for (k = 2; k <= *n-2; k++)
	{
	    s = (double)k / (*n);

	    u1[k-2] = 0.0;
	    for (j = 0; j < k; j++)
		for (i = 0; i < j; i++)
		    u1[k-2] += h[i + j * (*n)];
	    u1[k-2] /= k * (k - 1) / 2.0;

	    u2[k-2] = 0.0;
	    for (j = k; j < *n; j++)
		for (i = k; i < j; i++)
		    u2[k-2] += h[i + j * (*n)];
	    u2[k-2] /= (*n - k) * (*n - k - 1) / 2.0;

	    stat[k-2] = sqrtn * s * (1.0 - s) * fabs(u1[k-2] - u2[k-2]);
	}

    if (*method == 1 || *method == 2)
	/* generate M approximate realizations */
	for (m = 0; m < *M; m++)
	    for (k = 2; k <= *n-2; k++)
		{
		    s = (double)k / (*n);

		    if (*method == 1)
			{
			    ukm = 0.0;
			    for (j = 0; j < k; j++)
				for (i = 0; i < k; i++)
				    if (i != j)
					ukm += multipliers[i + m * (*n)] * (h[i + j * (*n)] - u1[k-2]);
			    ukm /= k * (k - 1) / 2.0;

			    unkm = 0.0;
			    for (j = k; j < *n; j++)
				for (i = k; i < *n; i++)
				    if (i != j)
					unkm += multipliers[i + m * (*n)] * (h[i + j * (*n)] - u2[k-2]);
			    unkm /= (*n - k) * (*n - k - 1) / 2.0;

			    stat0[m + (k - 2) * (*M)] = sqrtn * s * (1.0 - s) * fabs(ukm - unkm);
			}
		    else
			{
			    ukm = 0.0;
			    for (i = 0; i < k; i++)
				ukm += multipliers[i + m * (*n)] * (influ[i] - un);

			    unkm = 0.0;
			    for (i = k; i < *n; i++)
				unkm += multipliers[i + m * (*n)] * (influ[i] - un);

			    stat0[m + (k - 2) * (*M)] = 2.0 / sqrtn * fabs( (1.0 - s) * ukm - s * unkm);
			}
		}

    /* asymptotic variance */
    if (*method == 3)
	{
	    ln = 2 * (*bw) - 1;
	    *avar = 0.0;
	    for (i = 0; i < *n; i++)
		for (j = MAX(0, i - ln + 1); j < MIN(*n, i + ln); j++)
		    //abs(i - j) < ln
		    if (*w == 1)
			*avar += parzen( (double)(i - j) / (double)ln)
			    * (influ[i] - un) * (influ[j] - un);
		    else
			*avar += convrect( 4.0 * (double)(i - j) / (double)ln, 8)
			    * (influ[i] - un) * (influ[j] - un);
	    //*avar /= *n;
	}

    Free(multipliers);
    Free(u1);
    Free(u2);
}

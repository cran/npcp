/*#################################################################################
  ##
  ##   R package npcp by Ivan Kojadinovic Copyright (C) 2017
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

/***********************************************************************

  CUSUM change-point tests based on the sample mean
  bw: set bw to 1 for the iid case

***********************************************************************/

void cpTestMean(double *x, int *n, double *stat, int *M,
		int *w, int *bw, int *method, double *stat0, double *avar,
		double *initseq) {

    int i, j, k, m, ln;
    double *multipliers = R_Calloc((*n) * (*M), double);
    double *sum = R_Calloc(*n - 1, double);
    double *mean1 = R_Calloc(*n - 1, double);
    double *mean2 = R_Calloc(*n - 1, double);
    double s, meann, prockm, procnkm, sqrtn = sqrt(*n);

    /* generate (dependent) multipliers */
    if (*method == 1 || *method == 2)
	gendepmult(*n, *M, *bw, *w, initseq, multipliers);

    /* partial sums and global mean */
    sum[0] = x[0];
    for (j = 1; j < *n-1; j++)
	sum[j] = sum[j-1] + x[j];
    meann = sum[*n-2] / (*n);

    /* for each possible breakpoint */
    /* test statistics */
    for (k = 1; k <= *n-1; k++) {

	s = (double)k / (*n);
	mean1[k-1] = sum[k-1] / k;
	mean2[k-1] = (sum[*n-2] - sum[k-1]) / (*n - k);
	stat[k-1] = sqrtn * s * (1.0 - s) * fabs(mean1[k-1] - mean2[k-1]);
    }

    if (*method == 1 || *method == 2)
	/* generate M approximate realizations */
	for (m = 0; m < *M; m++) {

	     /* partial sums for nonsequential method */
	    if (*method == 2) {
		sum[0] = multipliers[0 + m * (*n)] * (x[0] - meann);
		for (j = 1; j < *n-1; j++)
		    sum[j] = sum[j-1] + multipliers[j + m * (*n)] * (x[j] - meann);
	    }

	    /* for each possible breakpoint */
	    for (k = 1; k <= *n-1; k++) {

		s = (double)k / (*n);

		if (*method == 1) {
		    prockm = 0.0;
		    for (j = 0; j < k; j++)
			prockm += multipliers[j + m * (*n)] * (x[j] - mean1[k-1]);

		    procnkm = 0.0;
		    for (j = k; j < *n; j++)
			procnkm += multipliers[j + m * (*n)] * (x[j] - mean2[k-1]);
		}
		else {
		    prockm = sum[k-1];
		    procnkm = sum[*n-2] - sum[k-1];
		}
		stat0[m + (k - 1) * (*M)] = fabs( (1.0 - s) * prockm - s * procnkm )
		    / sqrtn;
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
			    * (x[i] - meann) * (x[j] - meann);
		    else
			*avar += convrect( 4.0 * (double)(i - j) / (double)ln, 8)
			    * (x[i] - meann) * (x[j] - meann);
	    //*avar /= *n;
	}

    R_Free(multipliers);
    R_Free(sum);
    R_Free(mean1);
    R_Free(mean2);
}

/*#################################################################################
  ##
  ##   R package npcp by Ivan Kojadinovic Copyright (C) 2020
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
// SEQUENTIAL CHANGE-POINT TEST BASED ON MEANS
/////////////////////////////////////////////////////////////////////////

/* Statistcs / detectors */
void seqCpMeanStat(double *X, int *m, int *n, double *r, double *s,
		   double *t, double *e, double *cs, int *wr, int *we)
{
    int i, j, k;
    double diff, term, rk, sk, tk, ek;
    int nm = (*n) - (*m);
    double *sum = R_Calloc(nm + 1, double);
    double sqrtm = sqrt(*m), msqrtm = (*m) * sqrtm,
	m2 = (*m) * (*m), m2sqrtm = msqrtm * (*m);

    /* compute \sum_{i=1}^j X_i for j = m, ..., n */
    sum[0] = 0.0;  /* first line */
    for (i = 0; i < *m; i++)
	sum[0] += X[i];
    for (i = *m; i < *n; i++) /* remaining lines */
	sum[i - *m + 1] = sum[i - *m] + X[i];

    /* detectors */
    for (k = *m + 1; k <= *n; k++)
	{
	    /* Rm, Sm, Tm, Em */
	    rk = 0.0;
	    sk = 0.0;
	    tk = 0.0;
	    ek = 0.0;
 	    for (j = *m; j <= k-1; j++)
		{
		    diff = fabs(sum[j - *m] / j - (sum[k - *m] - sum[j - *m]) / (k - j));
		    term = j * (k - j) * diff ;
		    if (term > rk)
			{
			    wr[k - *m - 1] = j;
			    rk = term;
			}
		    sk += term;
		    tk += term * term;
		    term = (k - j) * diff ;
		    if (term > ek)
			{
			    we[k - *m - 1] = j;
			    ek = term;
			}
		}
	    r[k - *m - 1] = rk / msqrtm;
	    s[k - *m - 1] = sk / m2sqrtm;
	    t[k - *m - 1] = sqrt(tk) / m2;
	    e[k - *m - 1] = ek / sqrtm;

	    /* Qm */
	    cs[k - *m - 1] = (k - *m) / sqrtm
		* fabs(sum[0] / (*m) - (sum[k - *m] - sum[0]) / (k - *m));
	}

    R_Free(sum);
}


/* Long-run variance estimation from the learning sample */
void LRVmean(double *x, int *m,  int *w, int *bw, double *avar) {

    int i, j, ln;
    double mean;

    /* global mean */
    mean = 0.0;
    for (i = 0; i < *m; i++)
	mean += x[i];
    mean /= *m;

    /* bandwidth */
    ln = 2 * (*bw) - 1;

    /* long-run variance estimation */
    *avar = 0.0;
    for (i = 0; i < *m; i++)
	for (j = imax2(0, i - ln + 1); j < imin2(*m, i + ln); j++) //abs(i - j) < ln
	    if (*w == 1)
		*avar += parzen( (double)(i - j) / (double)ln)
		    * (x[i] - mean) * (x[j] - mean);
	    else
		*avar += convrect( 4.0 * (double)(i - j) / (double)ln, 8)
		    * (x[i] - mean) * (x[j] - mean);
    *avar /= *m;
}




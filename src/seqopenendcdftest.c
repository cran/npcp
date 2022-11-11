/*#################################################################################
  ##
  ##   R package npcp by Ivan Kojadinovic Copyright (C) 2021
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
// SEQUENTIAL OPEN-END CHANGE-POINT TEST BASED ON CDFS
/////////////////////////////////////////////////////////////////////////

/* Statistcs / detectors */
void seqOpenEndCpDistStat(double *Y, int *m, int *n, int *p,
			   double *invsigma, double *r, int *wr)
{
    int i, j, k, l, q;
    double norm, rk;
    int nm = (*n) - (*m), nm1 = nm + 1;
    double *sum = Calloc(nm1 * (*p), double);
    double *diff = Calloc(*p, double);
    double *tmp = Calloc(*p, double);
    double msqrtm = (*m) * sqrt(*m);

    /* Precompute cumulative sums */
    for (l = 0; l < *p; l++)
	{
	    /* Compute sum_{i=1}^j y[i,l] for j = m, ..., n */
	    sum[0 + l * nm1] = 0.0;
	    /* First line */
	    for (i = 0; i < *m; i++)
		sum[0 + l * nm1] += Y[i + l * (*n)];
	    /* Remaining lines */
	    for (i = *m; i < *n; i++)
		sum[i - *m + 1 + l * nm1] = sum[i - *m + l * nm1] + Y[i + l * (*n)];
	}


    /* Detector */
    for (k = *m + 1; k <= *n; k++)
	{
	    rk = 0.0;
	    for (j = *m; j <= k-1; j++)
		{
		    for (l = 0; l < *p; l++)
			diff[l] = sum[j - *m + l * nm1] / j
			    - (sum[k - *m + l * nm1] - sum[j - *m + l * nm1]) / (k - j);
		    for (l = 0; l < *p; l++)
			{
			    tmp[l] = 0.0;
			    for (q = 0; q < *p; q++)
				tmp[l] += diff[q] * invsigma[q + l * (*p)];
			}
		    norm = 0.0;
		    for (l = 0; l < *p; l++)
			norm += tmp[l] * diff[l];
		    norm = j * (k - j) * sqrt(norm / *p); /* normalization */
		    if (norm > rk)
			{
			    wr[k - *m - 1] = j;
			    rk = norm;
			}
		}
	    r[k - *m - 1] = rk / msqrtm;
	}

    Free(sum);
    Free(diff);
    Free(tmp);
}






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
#include "utilities.h"

/////////////////////////////////////////////////////////////////////////
// CHANGE-POINT TEST BASED ON EMPIRICAL CDFS
/////////////////////////////////////////////////////////////////////////

/***********************************************************************

  Change-point tests based on the empirical cdfs
  bw: set bw to 1 for the iid case

***********************************************************************/

void cpTestF(double *X, int *n, int *d, double *cvm, double *ks,
	     int *M, int *w, int *bw, int *seq, double *cvm0, double *ks0,
	     double *initseq)
{

    int i, j, k, q;
    double procq, t, sqrtn = sqrt(*n), multk, multnk;
    int *ind = R_Calloc((*n) * (*n), int);
    double *ecdf = R_Calloc(*n, double);
    double *indk = R_Calloc(*n, double);
    double *indnk = R_Calloc(*n, double);
    double *sumk = R_Calloc(*n, double);
    double *sumnk = R_Calloc(*n, double);
    double *multipliers = R_Calloc((*n) * (*M), double);

    /* generate (dependent) multipliers */
    gendepmult(*n, *M, *bw, *w, initseq, multipliers);

    /* compute 1(X_i <= X_q) and ecdf(X_q) */
    for (q = 0; q < *n; q++)
	{
	    ecdf[q] = 0.0;
	    for (i = 0; i < *n; i++)
		{
		    ind[i + q * (*n)] = 1;
		    for (j = 0; j < *d; j++)
			ind[i + q * (*n)] *= (X[i + j * (*n)] <= X[q + j * (*n)]);
		    ecdf[q] += ind[i + q * (*n)];
		}
	ecdf[q] /= *n;
    }

    /* test statistics */
    for (q = 0; q < *n; q++)
	    sumk[q] = 0.0;
	    //sumnk[q] = ecdf[q] * (*n);

    for (k = 1; k <= *n - 1; k++)
	{
	    t = (double)k / (*n);
	    cvm[k-1] = 0.0;
	    ks[k-1] = 0.0;
	    /* eval at X_q */
	    for (q = 0; q < *n; q++)
		{
		    sumk[q] += ind[k - 1 + q * (*n)];
		    procq = ( sumk[q] - t * ecdf[q] * (*n) ) / sqrtn;
		    cvm[k-1] += procq * procq;
		    if (fabs(procq) > ks[k-1])
			ks[k-1] = fabs(procq);
		}
	}

    /* generate M approximate realizations */
    for (j = 0; j < *M; j++)
	{
	    /* realization number j */
	    multk = 0.0;
	    multnk = 0.0;
	    for (q = 0; q < *n; q++)
		{
		    sumk[q] = 0.0;
		    sumnk[q] = 0.0;
		    for (i = 0; i < *n; i++)
			sumnk[q] += multipliers[i + j * (*n)] * ind[i + q * (*n)];
		    if (*seq == 1)
			{
			    indk[q] = 0.0;
			    indnk[q] = ecdf[q] * (*n);
			}
		    multnk += multipliers[q + j * (*n)];
		}

	    for (k = 1; k <= *n - 1; k++)
		{
		    t = (double)k / *n;
		    cvm0[j + (k-1) * (*M)] = 0.0;
		    ks0[j + (k-1) * (*M)] = 0.0;
		    multk += multipliers[k - 1 + j * (*n)];
		    multnk -= multipliers[k - 1 + j * (*n)];
		    /* eval at X_q */
		    for (q = 0; q < *n; q++)
			{
			    sumk[q] += multipliers[k - 1 + j * (*n)] * ind[k - 1 + q * (*n)];
			    sumnk[q] -= multipliers[k - 1 + j * (*n)] * ind[k - 1 + q * (*n)];
			    if (*seq == 1)
				{
				    indk[q] += ind[k - 1 + q * (*n)];
				    indnk[q] -= ind[k - 1 + q * (*n)];
				    procq =
					((1 - t) * (sumk[q] - multk * indk[q] / k)
					 - t * (sumnk[q] -  multnk * indnk[q] / (*n - k)))
					/ sqrtn;
				}
			    else
				procq = ((1 - t) * (sumk[q] - multk * ecdf[q])
					 - t * (sumnk[q] -  multnk * ecdf[q])) / sqrtn;
			    cvm0[j + (k-1) * (*M)] += procq * procq; /* FIXME: missing division by n */
			    if (fabs(procq) > ks0[j + (k-1) * (*M)])
				ks0[j + (k-1) * (*M)] = fabs(procq);
			}
		}
	}

    R_Free(ind);
    R_Free(multipliers);
    R_Free(ecdf);
    R_Free(indk);
    R_Free(indnk);
    R_Free(sumk);
    R_Free(sumnk);
}

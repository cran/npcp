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
#include "set_utils.h"

/////////////////////////////////////////////////////////////////////////
// RHO TEST
/////////////////////////////////////////////////////////////////////////

#define POW -0.51

void statinflu_seq(int n, int d, int k, double *U, int a, double fa, 
		   double *proda, double *influa, double *rho, 
		   double *influ)
{
    int i, j, l, p, aj;
    double sumk, sumnk, prodj, hk = R_pow(n,POW), hnk = R_pow(n,POW), 
	Uip, Uim;

    /* compute proda[i] and influa[i] */
    for (i = 0; i < k; i++)
	{
	    proda[i] = 1.0;
	    influa[i] = 0.0;
	    for (j = 0; j < d; j++)
		if ((1<<j) & a)
		    {
			proda[i] *= (1.0 - U[i + n * j]);
			aj = a & ~(1 << j); // A \ {j}
			for (p = 0; p < k; p++)
			    {
				prodj = 1.0;
				for (l = 0; l < d; l++)
				    if ((1<<l) & aj)
					prodj *= (1.0 - U[p + n * l]);
				Uip = MIN(U[i + n * j] + hk, 1.0);
				Uim = MAX(U[i + n * j] - hk, 0.0);
				influa[i] += prodj *
				    ( MIN(Uip, U[p + n * j]) - MIN(Uim, U[p + n * j]) )
				    / (Uip - Uim);
			    }
		    }
	    influa[i] /= k;
	}
    for (i = k; i < n; i++)
	{
	    proda[i] = 1.0;
	    influa[i] = 0.0;
	    for (j = 0; j < d; j++)
		if ((1<<j) & a)
		    {
			proda[i] *= (1.0 - U[i + n * j]);
			aj = a & ~(1 << j); // A \ {j}
			for (p = k; p < n; p++)
			    {
				prodj = 1.0;
				for (l = 0; l < d; l++)
				    if ((1<<l) & aj)
					prodj *= (1.0 - U[p + n * l]);
				Uip = MIN(U[i + n * j] + hnk, 1.0);
				Uim = MAX(U[i + n * j] - hnk, 0.0);
				influa[i] += prodj *
				    ( MIN(Uip, U[p + n * j]) - MIN(Uim, U[p + n * j]) )
				    / (Uip - Uim);
			    }
		    }
	    influa[i] /= n - k;
	}

    /* compute statistic and influence matrix for subset a and update */
    sumk = 0.0;
    for (i = 0; i < k; i++)
	{
	    sumk += proda[i];
	    influ[i] += fa * (proda[i] - influa[i]);
	}
    sumnk = 0.0;
    for (i = k; i < n; i++)
	{
	    sumnk += proda[i];
	    influ[i] += fa * (proda[i] - influa[i]);
	}
    rho[k - 1] += fa * (sumk / k - sumnk / (n - k));

}

/***********************************************************************

  Given k (nonsequential setting), computes rho statistic for subset A

***********************************************************************/

void stat_nonseq(int n, int d, int k, double *U, int a, double fa,
		 double *proda, double *rho)
{
    int i, j;
    double sumk, sumnk;

    /* compute proda[i] */
    for (i = 0; i < n; i++)
	{
	    proda[i] = 1.0;
	    for (j = 0; j < d; j++)
		if ((1<<j) & a)
		    proda[i] *= (1.0 - U[i + n * j]);

	}

    /* compute statistic for subset a and update */
    sumk = 0.0;
    for (i = 0; i < k; i++)
	sumk += proda[i];

    sumnk = 0.0;
    for (i = k; i < n; i++)
	sumnk += proda[i];

    rho[k - 1] += fa * (sumk / k - sumnk / (n - k));

}

/***********************************************************************

  In the nonsequential setting, computes influence matrix for subset A

***********************************************************************/

void influ_nonseq(int n, int d, double *V, int a, double fa, 
		  double *proda, double *influa, double *influ)
{
    int i, j, l, p, aj;
    double meanproda,  meaninflua, prodj, Vip, Vim, hn = R_pow(n,POW);

     /* compute proda[i] and influa[i] */ 
    meanproda = 0.0;
    meaninflua = 0.0;
    for (i = 0; i < n; i++)
	{
	    proda[i] = 1.0;
	    influa[i] = 0.0;
	    for (j = 0; j < d; j++)
		if ((1<<j) & a)
		    {
			proda[i] *= (1.0 - V[i + n * j]);
			aj = a & ~(1 << j); // A \ {j}
			for (p = 0; p < n; p++)
			    {
				prodj = 1.0;
				for (l = 0; l < d; l++)
				    if ((1<<l) & aj)
					prodj *= (1.0 - V[p + n * l]);
				Vip = MIN(V[i + n * j] + hn, 1.0);
				Vim = MAX(V[i + n * j] - hn, 0.0);
				influa[i] += prodj *
				    ( MIN(Vip, V[p + n * j]) - MIN(Vim, V[p + n * j]) )
				    / (Vip - Vim);
			    }
		    }
	    influa[i] /= n; 
	    meanproda += proda[i];
	    meaninflua += influa[i];
	} 
    meanproda /= n;
    meaninflua /= n;

    /* compute influence matrix for subset a and update */
    for (i = 0; i < n; i++)
	influ[i] += fa * (proda[i] - meanproda - influa[i] + meaninflua);
}

/***********************************************************************

  Returns influence matrix for determing bw automatically

***********************************************************************/

void influRho(double *X, int *n, int *d, double *fbin, double *influ)
{
    int i;
    int *index = Calloc(*n, int);
    double *V = Calloc((*n) * (*d), double); // pseudo-obs not depending on k
    double *x = Calloc((*n) * (*d), double);
    double *influa = Calloc(*n, double);
    double *proda = Calloc(*n, double);
    
    /* compute pseudo-obs V not depending on k */
    for (i = 0; i < (*n) * (*d); i++)
	x[i] = X[i];
    makepseudoobs(x, index, *n, *d, 0, *n, V);

    for (i = 0; i < *n; i++)
	influ[i] = 0.0;
    /* update influence matrix for subset i */
    for (i = 0; i < (1<<*d); i++)
	if (fbin[i])
	    influ_nonseq(*n, *d, V, i, fbin[i],
			 proda, influa, influ);

    Free(index);
    Free(V);
    Free(x);
    Free(influa);
    Free(proda);
}



/***********************************************************************

  Change-point tests based on extensions of Spearman's rho
  bw: set bw to 1 for the iid case

***********************************************************************/

void cpTestRho(double *X, int *n, int *d, double *rho, double *fbin,
	       double *influ, int *influest, int *M, int *w, int *bw,
	       int *method, double *rho0, double *avar, double *initseq)
{
    int i, j, k, m, ln;
    double *U = Calloc((*n) * (*d), double); // pseudo-obs depending on k 
    int *index = Calloc(*n, int);
    double *V = Calloc((*n) * (*d), double); // pseudo-obs not depending on k
    double *x = Calloc((*n) * (*d), double);
    double *influa = Calloc(*n, double);
    double *proda = Calloc(*n, double);
    double *multipliers = Calloc((*n) * (*M), double);
    double s, sumk, meank, sumnk, meannk, sqrtn = sqrt(*n);

    /* generate (dependent) multipliers */
    if (*method == 1 || *method == 2)
	gendepmult(*n, *M, *bw, *w, initseq, multipliers);

    /* compute pseudo-obs V not depending on k */
    if (*method == 2  || *method == 3)
	{
	    for (i = 0; i < (*n) * (*d); i++)
		x[i] = X[i];
	    makepseudoobs(x, index, *n, *d, 0, *n, V);
	}

    /* for each possible breakpoint */
    for (k = 1; k <= *n-1; k++)
	{
	    s = (double)k / (*n);

	    /* compute pseudo-obs U depending on k */
	    for (i = 0; i < (*n) * (*d); i++)
		x[i] = X[i];
	    makepseudoobs(x, index, *n, *d, 0, k, U);
	    makepseudoobs(x, index, *n, *d, k, *n, U);

	    /* sequential multiplier method */
	    if (*method == 1)
		{
		    /* initialization */
		    rho[k - 1] = 0.0;
		    for (i = 0; i < *n; i++)
			influ[i] = 0.0;
		    /* update statistic and influence matrix for subset i */
		    for (i = 0; i < (1<<*d); i++)
			if (fbin[i])
			    statinflu_seq(*n, *d, k, U, i, fbin[i], proda, 
					  influa, rho, influ);
		    
		    /* statistic for breakpoint k */
		    rho[k-1] = sqrtn * s * (1.0 - s) * fabs(rho[k-1]);

		    /* generate M approximate realizations */
		    for (m = 0; m < *M; m++)
			{
			    meank = 0.0;
			    for (i = 0; i < k; i++)
				meank += multipliers[i + m * (*n)];
			    meank /= k;

			    sumk = 0.0;
			    for (i = 0; i < k; i++)
				sumk += (multipliers[i + m * (*n)] - meank)
				    * influ[i];

			    meannk = 0.0;
			    for (i = k; i < *n; i++)
				meannk += multipliers[i + m * (*n)];
			    meannk /= *n - k;

			    sumnk = 0.0;
			    for (i = k; i < *n; i++)
				sumnk += (multipliers[i + m * (*n)] - meannk)
					* influ[i];

			    rho0[m + (k - 1) * (*M)] =
				fabs( (1.0 - s) * sumk - s * sumnk ) / sqrtn;
			}
		}
	    else /* nonsequential multiplier method or asymptotic variance */
 		{
		    rho[k - 1] = 0.0;
		    /* update statistic and influence matrix for subset i */
		    for (i = 0; i < (1<<*d); i++)
			if (fbin[i])
			    stat_nonseq(*n, *d, k, U, i, fbin[i], proda, rho);
		    /* statistic for breakpoint k */
		    rho[k-1] = sqrtn * s * (1.0 - s) * fabs(rho[k-1]);
		}
	}

    /* nonsequential multiplier method or asymptotic variance */
    if ((*method == 2 || *method == 3) && *influest == 0)
	{
	    /* update influence matrix for subset i */
	    for (i = 0; i < *n; i++)
		influ[i] = 0.0;
	    for (i = 0; i < (1<<*d); i++)
		if (fbin[i])
		    influ_nonseq(*n, *d, V, i, fbin[i],
				 proda, influa, influ);
	}

    /* nonsequential multiplier method */
    if (*method == 2)
	{
	     /* generate M approximate realizations */
	    for (m = 0; m < *M; m++)
		{
		    sumk = 0.0;
		    sumnk = 0.0;
		    for (i = 0; i < *n; i++)
			sumnk += multipliers[i + m * (*n)] * influ[i];

		    for (k = 1; k <= *n-1; k++)
			{
			    s = (double)k / (*n);
			    rho0[m + (k - 1) * (*M)] = 0.0;
			    sumk += multipliers[k - 1 + m * (*n)] * influ[k - 1];
			    rho0[m + (k - 1) * (*M)] =
				fabs(sumk - s * sumnk) / sqrtn;
			}
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
			    * influ[i] * influ[j];
		    else
			*avar += convrect( 4.0 * (double)(i - j) / (double)ln, 8)
			    * influ[i] * influ[j];
	    //*avar /= *n;
	}



    Free(U);
    Free(index);
    Free(V);
    Free(x);
    Free(influa);
    Free(proda);
    Free(multipliers);
}

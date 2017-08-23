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
// EMPIRICAL COPULA TESTS
/////////////////////////////////////////////////////////////////////////

/***********************************************************************

 Empirical copula at u
 U: pseudo-obs n * d
 b: line beginning
 e: line end (+1)

***********************************************************************/

double ec(double *U, int n, int d, int b, int e, double *u)
{
    int i,j;
    double ind, res = 0.0;

    for (i = b; i < e; i++)
	{
	    ind = 1.0;
	    for (j = 0; j < d; j++)
		ind *= (U[i + n * j] <= u[j]);
	    res += ind;
	}
    return res/(e-b);
}

/***********************************************************************

 Derivative of the empirical copula

***********************************************************************/

double derec(double *U, int n, int d, double *u, double *v, double denom, int b, int e)
{
    return (ec(U, n, d, b, e, u) - ec(U, n, d, b, e, v)) / denom;
}


/***********************************************************************

 Univariate empirical c.d.f. from portion b:(e-1)
 x: univariate data
 b: line beginning
 e: line end (+1)

***********************************************************************/

double uecdf(double *x, int b, int e, double u)
{
    int i;
    double res = 0.0;
    for (i = b; i < e; i++)
	    res += (x[i] <= u);
    return res/(e-b);
}

/***********************************************************************

  Compute influence matrix at every V pseudo-obs from portion b:(e-1)
  of pseudo-obs U

***********************************************************************/

void makeinflumat(int n, int d, int b, int e, double *U, double *V,
		  double *u, double *v, double *w, double *der,
		  double *influ)
{

    int i, j, l;
    double h = 1.0 / sqrt(e - b), ind, deriv, empc, denom;

    for (l = 0; l < n; l++)
	{
	    /* derivatives */
	    for (j = 0; j < d; j++)
		{
		    w[j] = V[j * n + l]; // pseudo-obs V line l
		    u[j] = w[j];
		    v[j] = w[j];
		}
	    for (j = 0; j < d; j++)
		{
		    u[j] += h; v[j] -= h;
		    denom = MIN(u[j], 1.0) - MAX(v[j], 0.0);
		    //der[j] = derec(U, n, d, u, v, 2.0 * h, b, e);
		    der[j] = derec(U, n, d, u, v, denom, b, e);
		    u[j] = w[j]; v[j] = w[j];
		}

	    empc = ec(U, n, d, b, e, w);

	    /* for each pseudo-obs U[i] at pseudo-obs V[l] */
	    for (i = b; i < e; i++)
		{
		    ind = 1.0;
		    deriv = 0.0;
		    for (j = 0; j < d; j++)
			{
			    ind *= (U[i + j * n] <= w[j]);
			    deriv += der[j] * ((U[i + j * n] <= w[j])
			    		       - uecdf(&U[j * n], b, e, w[j]));
			    // uecdf replaces V[l + j * n] * (e - b + 1.0) / (e - b));
			}
		    influ[i + l * n] = (ind - empc - deriv) / sqrt(n);
		}
	}
}




/***********************************************************************

  Change-point tests based on the empirical copula
  bw: set bw to 1 for the iid case

***********************************************************************/

void cpTestC(double *X, int *n, int *d, double *cvm, int *M,
	     int *we, int *bw, int *seq, double *cvm0, double *initseq)
{
    int i, j, k, l, m;
    int *index = Calloc(*n, int);
    double *U = Calloc((*n) * (*d), double); // pseudo-obs depending on k
    double *V = Calloc((*n) * (*d), double); // pseudo-obs not depending on k
    double *x = Calloc((*n) * (*d), double);
    double *u = Calloc(*d, double);
    double *v = Calloc(*d, double);
    double *w = Calloc(*d, double);
    double *sumk = Calloc(*n, double);
    double *sumnk = Calloc(*n, double);
    double *der = Calloc(*d, double);
    double *influ = Calloc((*n) * (*n), double);
    double *multipliers = Calloc((*n) * (*M), double);
    double s, diff;

    /* generate (dependent) multipliers */
    gendepmult(*n, *M, *bw, *we, initseq, multipliers);

    /* pseudo-obs in V */
    for (i = 0; i < (*n) * (*d); i++)
	x[i] = X[i];
    makepseudoobs(x, index, *n, *d, 0, *n, V);

    /* compute influence matrix nonseq at every V pseudo-obs */
    /* does not depend on k */
    if (*seq == 0)
	makeinflumat(*n, *d, 0, *n, V, V, u, v, w, der, influ);

    /* for each possible breakpoint */
    for (k = 1; k <= *n-1; k++)
	{
	    s = (double)k / (*n); // lambda_n(0,s)

	    /* compute pseudo-obs U depending on k */
	    for (i = 0; i < (*n) * (*d); i++)
		x[i] = X[i];
	    makepseudoobs(x, index, *n, *d, 0, k, U);
	    makepseudoobs(x, index, *n, *d, k, *n, U);

	    /* compute influence matrix seq at every V pseudo-obs */
	    /* does depend on k */
	    if (*seq == 1)
		{
		    makeinflumat(*n, *d, 0, k, U, V, u, v, w, der, influ);
		    makeinflumat(*n, *d, k, *n, U, V, u, v, w, der, influ);
		}

	    /* compute statistics */
	    cvm[k - 1] = 0.0;
	    for (l = 0; l < *n; l++) // at every pseudo-obs V
		{
		    for (j = 0; j < *d; j++)
			u[j] = V[j * (*n) + l]; // at pseudo-obs V line l
		    diff = ec(U, *n, *d, 0, k, u) - ec(U, *n, *d, k, *n, u);
		    cvm[k - 1] += diff * diff;
		}
	    cvm[k - 1] *= (*n) * s * s * (1.0 - s) * (1.0 - s);

	    /* generate M approximate realizations */
	    /* if sequential multiplier approach */
	    if (*seq == 1)
		for (m = 0; m < *M; m++)
		    {
			cvm0[m + (k - 1) * (*M)] = 0.0;
			for (l = 0; l < *n; l++)  /* for every pseudo-obs V */
			    {
				sumk[l] = 0.0;
				for (i = 0; i < k; i++)
				    sumk[l] += multipliers[i + m * (*n)]
					* influ[i + l * (*n)];

				sumnk[l] = 0.0;
				for (i = k; i < *n; i++)
				    sumnk[l] += multipliers[i + m * (*n)]
					* influ[i + l * (*n)];

				diff = (1 - s) * sumk[l] - s * sumnk[l];
				cvm0[m + (k - 1) * (*M)] += diff * diff;
			    }
		    }
	}

    /* generate M approximate realizations */
    /* if nonsequential multiplier approach */
    if (*seq == 0)
	for (m = 0; m < *M; m++)
	    {
		for (l = 0; l < *n; l++)
		    {
			sumk[l] = 0.0;
			sumnk[l] = 0.0;
			for (i = 0; i < *n; i++)
			    sumnk[l] += multipliers[i + m * (*n)]
				* influ[i + l * (*n)];

		    }
		for (k = 1; k <= *n-1; k++)
		    {
			s = (double)k / (*n);
			cvm0[m + (k - 1) * (*M)] = 0.0;
			for (l = 0; l < *n; l++)  /* for every pseudo-obs V */
			    {
				sumk[l] += multipliers[k - 1 + m * (*n)]
				    * influ[k - 1 + l * (*n)];

				diff = sumk[l] - s * sumnk[l];
				cvm0[m + (k - 1) * (*M)] += diff * diff;
			    }
		    }
	    }

    Free(index);
    Free(U);
    Free(V);
    Free(x);
    Free(u);
    Free(v);
    Free(w);
    Free(sumk);
    Free(sumnk);
    Free(der);
    Free(influ);
    Free(multipliers);
}

/***********************************************************************

 Form h-lagged data
 in: (n + h - 1) x d
 out: n x h x d
 pairwise == 1: only the first and last set of d columns

***********************************************************************/

void lagged(int n, int d, int h, double *in, double *out, int b, int e,
	    int pairwise)
{
    int i, j, l;
    if (pairwise) {
	for (i = b; i < e; i++)
	    for (j = 0; j < d; j++) {
		/* first set of d columns */
		out[i + j * n] = in[i + j * (n + h - 1)];
		/* second set of columns */
		out[i + (j + d) * n] = in[i + h - 1 + j * (n + h - 1)];
	    }
    }
    else {
	for (i = b; i < e; i++)
	    for (l = 0; l < h; l++)
		for (j = 0; j < d; j++)
		    out[i + (j + l * d) * n] = in[i + l + j * (n + h - 1)];
    }
}

/***********************************************************************

  Change-point tests based on the empirical autocopula
  bw: set bw to 1 for the iid case

***********************************************************************/

void cpTestAutocop(double *X, int *n, int *d, int *h, double *cvm, int *M,
		   int *we, int *bw, int *seq, double *cvm0, double *initseq,
		   int *pairwise)
{
    int i, j, k, l, m;
    int nh = *n - *h + 1; /* first dim of lagged data */
    int dh; /* second dim of lagged data */
    if (*pairwise == 1)
	dh = (*d) * 2; /* only first and last set of d columns */
    else
	dh = (*d) * (*h); /* all sets */
    int *index = Calloc(*n, int);
    double *U = Calloc((*n) * (*d), double); // pseudo-obs depending on k
    double *Uh = Calloc(nh * dh, double); // lagged pseudo-obs depending on k
    double *V = Calloc((*n) * (*d), double); // pseudo-obs not depending on k
    double *Vh = Calloc(nh * dh, double); // lagged pseudo-obs not depending on k
    double *x = Calloc((*n) * (*d), double);
    double *u = Calloc(dh, double);
    double *v = Calloc(dh, double);
    double *w = Calloc(dh, double);
    double *sumk = Calloc(nh, double);
    double *sumnk = Calloc(nh, double);
    double *der = Calloc(dh, double);
    double *influ = Calloc(nh * nh, double);
    double *multipliers = Calloc(nh * (*M), double);
    double s, diff;

    /* generate (dependent) multipliers */
    gendepmult(nh, *M, *bw, *we, initseq, multipliers);

    /* pseudo-obs in V */
    for (i = 0; i < (*n) * (*d); i++)
	x[i] = X[i];
    makepseudoobs(x, index, *n, *d, 0, *n, V);

    /* form h-lagged pseudo-obs not depending on k */
    lagged(nh, *d, *h, V, Vh, 0, nh, *pairwise);

    /* compute influence matrix nonseq at every V pseudo-obs */
    /* does not depend on k */
    if (*seq == 0)
	makeinflumat(nh, dh, 0, nh, Vh, Vh, u, v, w, der, influ);

    /* for each possible breakpoint */
    for (k = 1; k <= nh-1; k++)
	{
	    s = (double)k / nh; // lambda_n(0,s)

	    /* compute pseudo-obs U depending on k */
	    for (i = 0; i < (*n) * (*d); i++)
		x[i] = X[i];
	    makepseudoobs(x, index, *n, *d, 0, k + *h - 1, U);
	    lagged(nh, *d, *h, U, Uh, 0, k, *pairwise);
	    for (i = 0; i < (*n) * (*d); i++)
		x[i] = X[i];
	    makepseudoobs(x, index, *n, *d, k, *n, U);
	    lagged(nh, *d, *h, U, Uh, k, nh, *pairwise);

	    /* compute influence matrix seq at every V pseudo-obs */
	    /* does depend on k */
	    if (*seq == 1)
		{
		    makeinflumat(nh, dh, 0, k, Uh, Vh, u, v, w, der, influ);
		    makeinflumat(nh, dh, k, nh, Uh, Vh, u, v, w, der, influ);
		}

	    /* compute statistics */
	    cvm[k - 1] = 0.0;
	    for (l = 0; l < nh; l++) // at every pseudo-obs Vh
		{
		    for (j = 0; j < dh; j++)
			u[j] = Vh[j * nh + l]; // at pseudo-obs Vh line l
		    diff = ec(Uh, nh, dh, 0, k, u) - ec(Uh, nh, dh, k, nh, u);
		    cvm[k - 1] += diff * diff;
		}
	    cvm[k - 1] *= nh * s * s * (1.0 - s) * (1.0 - s);

	    /* generate M approximate realizations */
	    /* if sequential multiplier approach */
	    if (*seq == 1)
		for (m = 0; m < *M; m++)
		    {
			cvm0[m + (k - 1) * (*M)] = 0.0;
			for (l = 0; l < nh; l++)  /* for every pseudo-obs Vh */
			    {
				sumk[l] = 0.0;
				for (i = 0; i < k; i++)
				    sumk[l] += multipliers[i + m * nh]
					* influ[i + l * nh];

				sumnk[l] = 0.0;
				for (i = k; i < nh; i++)
				    sumnk[l] += multipliers[i + m * nh]
					* influ[i + l * nh];

				diff = (1 - s) * sumk[l] - s * sumnk[l];
				cvm0[m + (k - 1) * (*M)] += diff * diff;
			    }
		    }
	}

    /* generate M approximate realizations */
    /* if nonsequential multiplier approach */
    if (*seq == 0)
	for (m = 0; m < *M; m++)
	    {
		for (l = 0; l < nh; l++)
		    {
			sumk[l] = 0.0;
			sumnk[l] = 0.0;
			for (i = 0; i < nh; i++)
			    sumnk[l] += multipliers[i + m * nh]
				* influ[i + l * nh];

		    }
		for (k = 1; k <= nh-1; k++)
		    {
			s = (double)k / nh;
			cvm0[m + (k - 1) * (*M)] = 0.0;
			for (l = 0; l < nh; l++)  /* for every pseudo-obs Vh */
			    {
				sumk[l] += multipliers[k - 1 + m * nh]
				    * influ[k - 1 + l * nh];

				diff = sumk[l] - s * sumnk[l];
				cvm0[m + (k - 1) * (*M)] += diff * diff;
			    }
		    }
	    }

    Free(index);
    Free(U);
    Free(Uh);
    Free(V);
    Free(Vh);
    Free(x);
    Free(u);
    Free(v);
    Free(w);
    Free(sumk);
    Free(sumnk);
    Free(der);
    Free(influ);
    Free(multipliers);
}

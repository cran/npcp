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
//#include <Rmath.h>

#define	MIN(x,y) ((x) > (y) ? (y) : (x))
#define	MAX(x,y) ((x) < (y) ? (y) : (x))

/***********************************************************************

  Generate dependent multipliers

***********************************************************************/

void gendepmult(int n, int M, int bw, double *weights,
		double *initseq, double *multipliers)
{
    int i, j, m;

    GetRNGstate();
    for (m = 0; m < M; m++)
	for (i = 0; i < n; i++)
	    {
		multipliers[i + m * n] = 0.0;
		for (j = 0; j < 2 * bw - 1; j++)
		    multipliers[i + m * n] += weights[j]
			* initseq[i + j + m * (n + 2*(bw-1))];
	    }
    PutRNGstate();
}

/***********************************************************************

  Change-point tests based on the empirical cdfs
  bw: set bw to 1 for the iid case

***********************************************************************/

void cptestF(double *X, int *n, int *d, double *cvm, double *ks, int *M,
	     double *weights, int *bw, int *seq, double *cvm0,
	     double *ks0, double *initseq)
{

    int i, j, k, q;
    double procq, t, sqrtn = sqrt(*n), multk, multnk;
    int *ind = Calloc((*n) * (*n), int);
    double *ecdf = Calloc(*n, double);
    double *indk = Calloc(*n, double);
    double *indnk = Calloc(*n, double);
    double *sumk = Calloc(*n, double);
    double *sumnk = Calloc(*n, double);
    double *multipliers = Calloc((*n) * (*M), double);

    /* generate (dependent) multipliers */
    gendepmult(*n, *M, *bw, weights, initseq, multipliers);

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

    /* generate N approximate realizations */
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
			    cvm0[j + (k-1) * (*M)] += procq * procq;
			    if (fabs(procq) > ks0[j + (k-1) * (*M)])
				ks0[j + (k-1) * (*M)] = fabs(procq);
			}
		}
	}

    Free(ind);
    Free(multipliers);
    Free(ecdf);
    Free(indk);
    Free(indnk);
    Free(sumk);
    Free(sumnk);
}


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
	     double *weights, int *bw, int *seq, double *cvm0,
	     double *initseq)
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
    gendepmult(*n, *M, *bw, weights, initseq, multipliers);

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

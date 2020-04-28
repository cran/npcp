/*#################################################################################
  ##
  ##   R package npcp by Ivan Kojadinovic Copyright (C) 2019
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
// SEQUENTIAL CHANGE-POINT TEST BASED ON EMPIRICAL CDFS
/////////////////////////////////////////////////////////////////////////

/***********************************************************************

  Sequential change-point tests based on the empirical cdfs
  bw: set bw to 1 for the iid case

***********************************************************************/

/* Statistcs / detectors */
void seqCpDistStat(double *Y, int *m, int *n, int *d, double *mac,
		   double *mmc, double *mmk, double *mc, double *mk,
		   double *gamma, double *delta, int *wmc, int *wmk)
{

    int i, j, k, l;
    double d1k, d2k, d3k, proc1, proc2, term;
    int ind, nm = (*n) - (*m);
    double *sum = Calloc((nm + 1) * (*n), double);
    double *q = Calloc(*n + 1, double);
    double f = R_pow((double)(*m), 1.5);

    /* compute \sum_{i=1}^j 1(X_i <= X_l) for i = m, ..., n*/
    for (l = 0; l < *n; l++)
	{
	    /* first line (line 0) */
	    sum[0 + l * (nm + 1)] = 0.0;
	    for (i = 0; i < *m; i++) {
		ind = 1;
		for (j = 0; j < *d; j++)
		    ind *= (Y[i + j * (*n)] <= Y[l + j * (*n)]);
		sum[0 + l * (nm + 1)] += ind;
	    }

	    /* remaining lines */
	    for (i = *m; i < *n; i++)
		{
		    ind = 1;
		    for (j = 0; j < *d; j++)
			ind *= (Y[i + j * (*n)] <= Y[l + j * (*n)]);
		    sum[i - (*m) + 1 + l * (nm+1)] = sum[i - (*m) + l * (nm+1)] + ind;
		}
	}

    /* precompute weights */
    for (i = 0; i <= *n; i++)
	q[i] = R_pow((double)i/(*m), *gamma);


    /* detectors */
    for (k = *m + 1; k <= *n; k++)
	{
	    /* Rm,q, Sm,q, Tm,q */
	    d1k = 0.0;
	    d2k = 0.0;
	    d3k = 0.0;
	    for (j = *m; j <= k-1; j++)
		{
		    proc1 = 0.0;
		    proc2 = 0.0;
		    for (l = 0; l < k; l++)
			{
			    term = j * (k - j) / fmax2(q[j] *  q[k-j], *delta)
				* (sum[j - *m + l * (nm+1)] / j
				   - (sum[k - *m + l * (nm+1)]
				      - sum[j - *m + l * (nm+1)]) / (k - j));
			    proc1 = fmax2(proc1, fabs(term));
			    proc2 += term * term;
			}
		    proc2 /= k;
		    if (proc1 > d1k)
			{
			    wmc[k - *m - 1] = j;
			    d1k = proc1;
			}
		    if (proc2 > d2k)
			{
			    wmk[k - *m - 1] = j;
			    d2k = proc2;
			}
		    d3k += proc2;
		}
	    mmk[k - *m - 1] = d1k / f;
	    mmc[k - *m - 1] = d2k / (f * f);
	    mac[k - *m - 1] = d3k / (f * f * (*m));

	    /* Pm, Qm */
	    proc1 = 0.0;
	    proc2 = 0.0;
	    for (l = 0; l < k; l++)
		{
		    term = (k - *m) * sum[0 + l * (nm+1)]
			- (*m) * (sum[k - *m + l * (nm+1)] - sum[0 + l * (nm+1)]);
		    proc1 = fmax2(proc1, fabs(term));
		    proc2 += term * term;
		}
	    mk[k - *m - 1] =  proc1 / f;
	    mc[k - *m - 1] =  proc2 / (f * f * k);
	}

    Free(sum);
    Free(q);
}

/////////////////////////////////////////////////////////////////////////
/* Multiplier replicates sequential version 1 */
/////////////////////////////////////////////////////////////////////////

void seqCpDistMultSeq1(double *X, int *m, int *n, int *d, int *B, int *w,
		       int *bw, double *mac0, double *mmc0, double *mmk0,
		       double *mc0, double *mk0, double *gamma,
		       double *delta, double *initseq)
{
    int i, j, k, l, b;
    int *I = Calloc((*m) * (*m), int); /* indicators m x m */
    int mm = (int)((*m) * (*m) / (double)(*n)), mmm1 = *m - mm + 1;
    double *sumb = Calloc(mmm1 * (*m), double);
    double *F = Calloc(mmm1 * (*m), double); /* sequential edfs */
    double *multipliers = Calloc((*m) * (*B), double);
    double *summult = Calloc(mmm1, double);
    double d1k, d2k, d3k, proc1, proc2, term;
    double f = sqrt(mm) * mm;
    double *q = Calloc(*m + 1, double);

    /* generate (dependent) multipliers */
    gendepmult(*m, *B, *bw, *w, initseq, multipliers);

    /* compute 1(X_i <= X_l) and ecdf(X_l) */
    for (l = 0; l < *m; l++)
    	{
    	    F[0 + l * mmm1] = 0.0;
    	    for (i = 0; i < mm; i++)
    		{
    		    I[i + l * (*m)] = 1;
    		    for (j = 0; j < *d; j++)
    			I[i + l * (*m)] *= (X[i + j * (*m)] <= X[l + j * (*m)]);
		    F[0 + l * mmm1] += I[i + l * (*m)];
    		}
	    for (i = mm; i < *m; i++)
    		{
    		    I[i + l * (*m)] = 1;
    		    for (j = 0; j < *d; j++)
    			I[i + l * (*m)] *= (X[i + j * (*m)] <= X[l + j * (*m)]);
		    F[i - mm + 1 + l * mmm1] = F[i - mm + l * mmm1] + I[i + l * (*m)];
		    F[i - mm + l * mmm1] /= i;
    		}
	    F[mmm1 - 1 + l * mmm1] /= *m;
	}

    /* precompute weights */
    for (i = 0; i <= *m; i++)
	q[i] = R_pow((double)i / mm, *gamma);

    /* generate B approximate realizations of maxima of detectors */
    for (b = 0; b < *B; b++)
    	{
    	    /* first line (line 0) of matrix sumb */
    	    for (l = 0; l < *m; l++)
    		{
		    sumb[0 + l * mmm1] = 0.0;
    		    for (i = 0; i < mm; i++)
    			sumb[0 + l * mmm1]
			    += multipliers[i + b * (*m)] * I[i + l * (*m)] ;
    		}

	    /* summult[0] */
	    summult[0] = 0.0;
	    for (i = 0; i < mm; i++)
		summult[0] += multipliers[i + b * (*m)];

    	    /* multiplier replicates of maxima of detectors */
    	    for (k = mm + 1; k <= *m; k++)
    	    	{
		    /* line k - mm of matrix sumb */
		    for (l = 0; l < *m; l++)
			sumb[k - mm + l * mmm1] = sumb[k - mm - 1 + l * mmm1]
			    + multipliers[k - 1 + b * (*m)] * I[k - 1 + l * (*m)];

		    /* summult[k - mm] */
		    summult[k - mm] = summult[k - mm - 1]
			+ multipliers[k - 1 + b * (*m)];

		    /* Rm,q, Sm,q, Tm,q */
		    d1k = 0.0;
		    d2k = 0.0;
		    d3k = 0.0;
		    for (j = mm; j <= k-1; j++)
			{
			    proc1 = 0.0;
			    proc2 = 0.0;
			    for (l = 0; l < k; l++)
				{
				    term = (k * (sumb[j - mm + l * mmm1]
						 - summult[j - mm]
						 * F[j - mm + l * mmm1])
				    	    - j * (sumb[k - mm + l * mmm1]
						   - summult[k - mm]
						   * F[k - mm + l * mmm1]))
					/ fmax2(q[j] * q[k - j], *delta);
				    proc1 = fmax2(proc1, fabs(term));
				    proc2 += term * term;
				}
			    d1k = fmax2(d1k, proc1);
			    d2k = fmax2(d2k, proc2 / k);
			    d3k += proc2 / k;
			}

		    mmk0[k - mm - 1 +  b * ((*m) - mm)] = d1k / f;
		    mmc0[k - mm - 1 +  b * ((*m) - mm)] = d2k / (f * f);
		    mac0[k - mm - 1 +  b * ((*m) - mm)] = d3k / (f * f * mm);


		    /* Pm, Qm */
		    proc1 = 0.0;
		    proc2 = 0.0;
		    for (l = 0; l < k; l++)
			{
			    term = k * (sumb[0 + l * mmm1] - summult[0] * F[0 + l * mmm1])
			    	- mm * (sumb[k - mm + l * mmm1] - summult[k - mm] * F[k - mm + l * mmm1]);
			    proc1 = fmax2(proc1, fabs(term));
			    proc2 += term * term;
			}
		    mk0[k - mm - 1 +  b * ((*m) - mm)] = proc1 / f;
		    mc0[k - mm - 1 +  b * ((*m) - mm)] = proc2 / (f * f * k);
    	    	}
    	}

    Free(I);
    Free(sumb);
    Free(F);
    Free(summult);
    Free(multipliers);
    Free(q);
}

/////////////////////////////////////////////////////////////////////////
/* Multiplier replicates sequential version 2 */
/////////////////////////////////////////////////////////////////////////
void seqCpDistMultSeq2(double *X, int *m, int *n, int *d, int *B, int *w,
		       int *bw, double *mac0, double *mmc0, double *mmk0,
		       double *mc0, double *mk0, double *gamma,
		       double *delta, double *initseq)
{
    int i, j, k, l, b;
    int *I = Calloc((*m) * (*m), int); /* indicators m x m */
    int mm = (int)((*m) * (*m) / (double)(*n)), mmm1 = *m - mm + 1;
    double *sumb = Calloc(mmm1 * (*m), double);
    double *sum = Calloc(mmm1 * (*m), double); /* sequential edfs */
    double *multipliers = Calloc((*m) * (*B), double);
    double *summult = Calloc(mmm1, double);
    double d1k, d2k, d3k, proc1, proc2, term;
    double f = sqrt(mm) * mm;
    double *q = Calloc(*m + 1, double);

    /* generate (dependent) multipliers */
    gendepmult(*m, *B, *bw, *w, initseq, multipliers);

    /* compute 1(X_i <= X_l) and ecdf(X_l) */
    for (l = 0; l < *m; l++)
    	{
    	    sum[0 + l * mmm1] = 0.0;
    	    for (i = 0; i < mm; i++)
    		{
    		    I[i + l * (*m)] = 1;
    		    for (j = 0; j < *d; j++)
    			I[i + l * (*m)] *= (X[i + j * (*m)] <= X[l + j * (*m)]);
		    sum[0 + l * mmm1] += I[i + l * (*m)];
    		}
	    for (i = mm; i < *m; i++)
    		{
    		    I[i + l * (*m)] = 1;
    		    for (j = 0; j < *d; j++)
    			I[i + l * (*m)] *= (X[i + j * (*m)] <= X[l + j * (*m)]);
		    sum[i - mm + 1 + l * mmm1] = sum[i - mm + l * mmm1] + I[i + l * (*m)];
    		}
	}

    /* precompute weights */
    for (i = 0; i <= *m; i++)
	q[i] = R_pow((double)i / mm, *gamma);

    /* generate B approximate realizations of maxima of detectors */
    for (b = 0; b < *B; b++)
    	{
    	    /* first line (line 0) of matrix sumb */
    	    for (l = 0; l < *m; l++)
    		{
		    sumb[0 + l * mmm1] = 0.0;
    		    for (i = 0; i < mm; i++)
    			sumb[0 + l * mmm1]
			    += multipliers[i + b * (*m)] * I[i + l * (*m)] ;
    		}

	    /* summult[0] */
	    summult[0] = 0.0;
	    for (i = 0; i < mm; i++)
		summult[0] += multipliers[i + b * (*m)];

    	    /* multiplier replicates of maxima of detectors */
    	    for (k = mm + 1; k <= *m; k++)
    	    	{
		    /* line k - mm of matrix sumb */
		    for (l = 0; l < *m; l++)
			sumb[k - mm + l * mmm1] = sumb[k - mm - 1 + l * mmm1]
			    + multipliers[k - 1 + b * (*m)] * I[k - 1 + l * (*m)];

		    /* summult[k - mm] */
		    summult[k - mm] = summult[k - mm - 1]
			+ multipliers[k - 1 + b * (*m)];

		    /* Rm,q, Sm,q, Tm,q */
		    d1k = 0.0;
		    d2k = 0.0;
		    d3k = 0.0;
		    for (j = mm; j <= k-1; j++)
			{
			    proc1 = 0.0;
			    proc2 = 0.0;
			    for (l = 0; l < k; l++)
				{
				    term = ((k - j)
					    * (sumb[j - mm + l * mmm1]
					       - summult[j - mm]
					       / j * sum[j - mm + l * mmm1])
				    	    - j * (sumb[k - mm + l * mmm1]
						   - sumb[j - mm + l * mmm1]
						   - (summult[k - mm] - summult[j - mm])
						   / (k - j)
						   * (sum[k - mm + l * mmm1]
						      - sum[j - mm + l * mmm1])))
					/ fmax2(q[j] * q[k - j], *delta);
				    proc1 = fmax2(proc1, fabs(term));
				    proc2 += term * term;
				}
			    d1k = fmax2(d1k, proc1);
			    d2k = fmax2(d2k, proc2 / k);
			    d3k += proc2 / k;
			}

		    mmk0[k - mm - 1 +  b * ((*m) - mm)] = d1k / f;
		    mmc0[k - mm - 1 +  b * ((*m) - mm)] = d2k / (f * f);
		    mac0[k - mm - 1 +  b * ((*m) - mm)] = d3k / (f * f * mm);

		    /* Pm, Qm */
		    proc1 = 0.0;
		    proc2 = 0.0;
		    for (l = 0; l < k; l++)
			{
			    term = k * (sumb[0 + l * mmm1]
					- summult[0] / mm * sum[0 + l * mmm1])
			    	- mm * (sumb[k - mm + l * mmm1]
					- summult[k - mm] / k * sum[k - mm + l * mmm1]);
			    proc1 = fmax2(proc1, fabs(term));
			    proc2 += term * term;
			}
		    mk0[k - mm - 1 +  b * ((*m) - mm)] = proc1 / f;
		    mc0[k - mm - 1 +  b * ((*m) - mm)] = proc2 / (f * f * k);
    	    	}
    	}

    Free(I);
    Free(sumb);
    Free(sum);
    Free(summult);
    Free(multipliers);
    Free(q);
}

/////////////////////////////////////////////////////////////////////////
/* Multiplier replicates non sequential version */
/////////////////////////////////////////////////////////////////////////
void seqCpDistMultNonSeq(double *X, int *m, int *n, int *d, int *B,
			 int *w, int *bw, double *mac0, double *mmc0,
			 double *mmk0, double *mc0, double *mk0,
			 double *gamma, double *delta, double *initseq)
{
    int i, j, k, l, b;
    int *I = Calloc((*m) * (*m), int); /* indicators m x m */
    int mm = (int)((*m) * (*m) / (double)(*n)), mmm1 = *m - mm + 1;
    double *sumb = Calloc(mmm1 * (*m), double);
    double *F = Calloc(*m, double); /* ecdfs */
    double *multipliers = Calloc((*m) * (*B), double);
    double d1k, d2k, d3k, proc1, proc2, term;
    double f = sqrt(mm) * mm;
    double *q = Calloc(*m + 1, double);

    /* generate (dependent) multipliers */
    gendepmult(*m, *B, *bw, *w, initseq, multipliers);

    /* compute 1(X_i <= X_l) and ecdf(X_l) */
    for (l = 0; l < *m; l++)
    	{
    	    F[l] = 0.0;
    	    for (i = 0; i < *m; i++)
    		{
    		    I[i + l * (*m)] = 1;
    		    for (j = 0; j < *d; j++)
    			I[i + l * (*m)] *= (X[i + j * (*m)] <= X[l + j * (*m)]);
		    F[l] += I[i + l * (*m)];
    		}
	    F[l] /= *m;
	}

    /* precompute weights */
    for (i = 0; i <= *m; i++)
	q[i] = R_pow((double)i / mm, *gamma);

    /* generate B approximate realizations of maxima of detectors */
    for (b = 0; b < *B; b++)
    	{
    	    /* first line (line 0) of matrix sumb */
    	    for (l = 0; l < *m; l++)
    		{
    		    sumb[0 + l * mmm1] = 0.0;
    		    for (i = 0; i < mm; i++)
    			sumb[0 + l * mmm1] += multipliers[i + b * (*m)]
    			    * (I[i + l * (*m)] - F[l]);
    		}

    	    /* multiplier replicates of maxima of detectors */
    	    for (k = mm + 1; k <= *m; k++)
    	    	{
		    /* line k - mm of matrix sumb */
		    for (l = 0; l < *m; l++)
			sumb[k - mm + l * mmm1] = sumb[k - mm - 1 + l * mmm1]
			    + multipliers[k - 1 + b * (*m)]
			    * (I[k - 1 + l * (*m)] - F[l]);

		    /* Rm,q, Sm,q, Tm,q */
		    d1k = 0.0;
		    d2k = 0.0;
		    d3k = 0.0;
		    for (j = mm; j <= k-1; j++)
			{
			    proc1 = 0.0;
			    proc2 = 0.0;
			    for (l = 0; l < k; l++)
				{
				    term = (k * sumb[j - mm + l * mmm1]
					    - j * sumb[k - mm + l * mmm1])
					/ fmax2(q[j] * q[k - j], *delta);
				    proc1 = fmax2(proc1, fabs(term));
				    proc2 += term * term;
				}
			    d1k = fmax2(d1k, proc1);
			    d2k = fmax2(d2k, proc2 / k);
			    d3k += proc2 / k;
			}

		    mmk0[k - mm - 1 +  b * ((*m) - mm)] = d1k / f;
		    mmc0[k - mm - 1 +  b * ((*m) - mm)] = d2k / (f * f);
		    mac0[k - mm - 1 +  b * ((*m) - mm)] = d3k / (f * f * mm);

		    /* Pm, Qm */
		    proc1 = 0.0;
		    proc2 = 0.0;
		    for (l = 0; l < k; l++)
			{
			    term = (k * sumb[0 + l * mmm1]
				    - mm * sumb[k - mm + l * mmm1]);
			    proc1 = fmax2(proc1, fabs(term));
			    proc2 += term * term;
			}
		    mk0[k - mm - 1 +  b * ((*m) - mm)] = proc1 / f;
		    mc0[k - mm - 1 +  b * ((*m) - mm)] = proc2 / (f * f * k);
    	    	}
    	}

    Free(I);
    Free(sumb);
    Free(F);
    Free(multipliers);
    Free(q);
}


/////////////////////////////////////////////////////////////////////////
/* Generates iid sample from the beta copula */
/////////////////////////////////////////////////////////////////////////

void rBetaCopula(int *r, int *m, int *d, int *n, double *x)
{
    int i, j, I, rr;

    GetRNGstate();

    for (i = 0; i < *n; i++)
	{
	    I = (int)runif(0.0, *m);
	    for (j = 0; j < *d; j++)
		{
		    rr = r[I + j * (*m)];
		    x[i + j * (*n)] = rbeta(rr, *m + 1 - rr);
		}
	}

    PutRNGstate();
}


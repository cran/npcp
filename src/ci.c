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
// U-STAT confidence intervals
/////////////////////////////////////////////////////////////////////////


/***********************************************************************

  U-STAT confidence intervals
  bw: set bw to 1 for the iid case

***********************************************************************/

/* void ciU(double *h, int *n, double *influ, double *stat, int *M, */
/* 	 int *w, int *bw, double *stat0, double *avar, double *initseq) */
/* { */
/*     int i, j, m, ln; */
/*     double *multipliers = Calloc((*n) * (*M), double); */
/*     double un, sqrtn = sqrt(*n); */

/*     /\* generate (dependent) multipliers *\/ */
/*     //if (*method == 1) */
/*     gendepmult(*n, *M, *bw, *w, initseq, multipliers); */

/*     /\* global U-statistic *\/ */
/*     un = 0.0; */
/*     for (j = 0; j < *n; j++) */
/* 	for (i = 0; i < j; i++) */
/* 	    un += h[i + j * (*n)]; */
/*     un /= (*n) * (*n - 1) / 2.0; */
/*     *stat = un; */


/*     //if (*method == 1) /\* generate M approximate realizations *\/ */
/*     for (m = 0; m < *M; m++) */
/* 	{ */
/* 	    un = 0.0; */
/* 	    for (j = 0; j < *n; j++) */
/* 		for (i = 0; i < *n; i++) */
/* 		    if (i != j) */
/* 			un += multipliers[i + m * (*n)] * (h[i + j * (*n)] - *stat); */
/* 	    un /= (*n) * (*n - 1) / 2.0; */
/* 		stat0[m] = sqrtn * un; */
/* 	} */
/*     //else /\* asymptotic variance *\/ */
/*     //	{ */
/*     ln = 2 * (*bw) - 1; */
/*     *avar = 0.0; */
/*     for (i = 0; i < *n; i++) */
/* 	for (j = MAX(0, i - ln + 1); j < MIN(*n, i + ln); j++) */
/* 	    //abs(i - j) < ln */
/* 	    if (*w == 1) */
/* 		*avar += parzen( (double)(i - j) / (double)ln) */
/* 		    * (influ[i] - *stat) * (influ[j] - *stat); */
/* 	    else */
/* 		*avar += convrect( 4.0 * (double)(i - j) / (double)ln, 8) */
/* 		    * (influ[i] - *stat) * (influ[j] - *stat); */
/*     *avar /= *n; */
/*     // } */

/*     Free(multipliers); */
/* } */

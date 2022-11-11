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

#ifndef NPCP_DEFS_H
#define NPCP_DEFS_H

#include <R.h>

// ./seqopenendcdftest.c ///////////////////////////////////////////////////////////////////

void seqOpenEndCpDistStat(double *Y, int *m, int *n, int *p,
			  double *invsigma, double *r, int *wr);

// ./seqmeantest.c ///////////////////////////////////////////////////////////////////

void seqCpMeanStat(double *X, int *m, int *n, double *r, double *s,
		   double *t, double *e, double *cs, int *wr, int *we);

void LRVmean(double *x, int *m,  int *w, int *bw, double *avar);

// ./seqcdftest.c ///////////////////////////////////////////////////////////////////

void seqCpDistStat(double *Y, int *m, int *n, int *d, double *mac,
		   double *mmc, double *mmk, double *mc, double *mk,
		   double *gamma, double *delta, int *wmc, int *wmk, int *scale);

void seqCpDistMultSeq1(double *X, int *m, int *n, int *d, int *B, int *w,
		       int *bw, double *mac0, double *mmc0, double *mmk0,
		       double *mc0, double *mk0, double *gamma,
		       double *delta, double *initseq, int *scale);

void seqCpDistMultSeq2(double *X, int *m, int *n, int *d, int *B, int *w,
		       int *bw, double *mac0, double *mmc0, double *mmk0,
		       double *mc0, double *mk0, double *gamma,
		       double *delta, double *initseq, int *scale);

void seqCpDistMultNonSeq(double *X, int *m, int *n, int *d, int *B,
			 int *w, int *bw, double *mac0, double *mmc0,
			 double *mmk0, double *mc0, double *mk0,
			 double *gamma, double *delta, double *initseq, int *scale);

void rBetaCopula(int *r, int *m, int *d, int *n, double *x);

// ./cdftest.c ///////////////////////////////////////////////////////////////////

void cpTestF(double *X, int *n, int *d, double *cvm, double *ks, int *M,
	     int *w, int *bw, int *seq, double *cvm0, double *ks0,
	     double *initseq);

// ./ectest.c ///////////////////////////////////////////////////////////////////

void cpTestC(double *X, int *n, int *d, double *cvm, int *M,
	     int *we, int *bw, int *seq, double *cvm0, double *initseq);

void cpTestAutocop(double *X, int *n, int *d, int *h, double *cvm, int *M,
		   int *we, int *bw, int *seq, double *cvm0, double *initseq,
		   int *pairwise);

// ./set_utils.c ///////////////////////////////////////////////////////////////

void k_power_set(int *n, int *k, int *power_set);

void natural2binary(int *n, double *sf, int *power_set, double *sf_out);

// ./rhotest.c ///////////////////////////////////////////////////////////////////

void cpTestRho(double *X, int *n, int *d, double *rho, double *fbin,
	       double *influ, int *influest, int *M, int *w, int *bw,
	       int *method, double *rho0, double *avar, double *initseq);

void influRho(double *X, int *n, int *d, double *fbin, double *influ);

// ./meantest.c ///////////////////////////////////////////////////////////////////

void cpTestMean(double *x, int *n, double *stat, int *M,
		int *w, int *bw, int *method, double *stat0, double *avar,
		double *initseq);

// ./Utest.c ///////////////////////////////////////////////////////////////////

void cpTestU(double *h, int *n, double *influ, double *stat, int *M,
	     int *w, int *bw, int *method, double *stat0, double *avar,
	     double *initseq);

// ./bmtest.c ///////////////////////////////////////////////////////////////////

void cpTestBM(double *X, int *n, int *r, double *stat, double *gamman,
	      double *deltan, double *gamma, double *delta,
	      int *meth, int *landwehr, int *noties, int* center,
	      double *param, double *avar);

void fitGEV(double *X, int *n, double *gamma, double *delta, int *meth,
	    int *landwehr, int *noties, double *param, double *avar);

// ./utilities.c ///////////////////////////////////////////////////////////////////

void pdf_sum_unif(int *n, double *y, int *ny, double *pdfy);



#endif

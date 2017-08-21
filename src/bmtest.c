/*#################################################################################
  ##
  ##   R package npcp by Ivan Kojadinovic Copyright (C) 2015
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
#include <Rinternals.h>
#include "utilities.h"

/////////////////////////////////////////////////////////////////////////
// Block maxima PWM tests
/////////////////////////////////////////////////////////////////////////

/***********************************************************************

  Compute scaled ranks (ecdfs) from portion b:(e-1)
  of the univariate sequence in x.
  b: line beginning
  e: line end (+1)
  KEEPS INDEX

***********************************************************************/

void ecdfs(double *x, int *index, int n, int b, int e, double *ecdf,
	   double gamma, double delta, int noties)
{
    int i,j;
    if (e == b)
	return;
    if (noties) // if no ties, use sorting
	{
	    for (i = 0; i < e - b; i++)
		index[b + i] = i;
	    rsort_with_index (&x[b], &index[b], e - b);
	    for (i = 0; i < e - b; i++)
		ecdf[b + index[b + i]] = (i + 1.0 + gamma) / (e - b + delta);
	}
    else // else maximal ranks
	for (i = 0; i < e - b; i++)
	    {
		ecdf[b + i] = 0.0;
		for (j = 0; j < e - b; j++)
		    if (x[b + j] <= x[b + i])
			ecdf[b + i] += 1.0;
		ecdf[b + i] = (ecdf[b + i] + gamma) / (e - b + delta);
	    }
}

/***********************************************************************

  C.d.f. of the GEV(mu, sigma, xi)

***********************************************************************/

double pgev(double x, double mu, double sigma, double xi)
{
    if (sigma <= 0.0)
	{
	    Rprintf("Error: invalid sigma in pgev\n");
	    return(NAN);
	}
    x = (x - mu) / sigma;
    if (xi == 0.0)
        return(exp(-exp(-x)));
    else
	return(exp(-R_pow(MAX(1.0 + xi * x, 0.0), -1.0 / xi)));

}

/***********************************************************************

  Function and derivative involved in probability weighted moments

***********************************************************************/

double omega(double x, double a, double b, int survival)
{
    if (survival)
	x = 1.0 - x;

    if (a == 0.0 && b == 0.0)
	return(1.0);
    else if (a >= 0.0 && b == 0.0)
	return(R_pow(x, a));
    else if (a >= 0.0 && b >= 0.0)
	return(R_pow(x, a) * R_pow(-log(x), b));
    else
	{
	    Rprintf("Wrong combination of powers a and b in omega\n");
	    return(NAN);
	}
}

double domega(double x, double a, double b, int survival)
{
    double res;

    if (survival)
	x = 1.0 - x;

    if (a == 0.0 && b == 0.0)
	res = 0.0;
    else if (a >= 1.0 && b == 0.0)
	res = a * R_pow(x, a - 1.0);
    else if (a >= 1.0 && b >= 1.0)
	res = a * R_pow(x, a - 1.0) * R_pow(-log(x), b) - R_pow(x, a) * b / x * R_pow(-log(x), b - 1.0);
    else
	{
	    Rprintf("Wrong combination of powers a and b in domega\n");
	    return(NAN);
	}

    if (survival)
	return(-res);
    return (res);
}

/***********************************************************************

  Estimate PWMs for each omega function

***********************************************************************/

void estimate_pwm(int p, int s, int e, double *beta, double *X,
		  double *cdf, double *a, double *b, int survival)
{
    int i,l;
    for (l = 0; l < p; l++)
	{
	    beta[l] = 0.0;
	    for (i = s; i < e; i++)
		beta[l] += X[i] * omega(cdf[i], a[l], b[l], survival);
	    beta[l] /= e - s;
	}
}

/***********************************************************************

  Landwehr's PWMs estimates

***********************************************************************/

void estimate_landwehr(int p, int s, int e, double *beta, double *x)
{
    int i;
    R_rsort (&x[s], e - s);
    beta[0] = 0.0;
    beta[1] = 0.0;
    beta[2] = 0.0;
    for (i = 0; i < e - s; i++)
	{
	    beta[0] += x[i + s];
	    beta[1] += i / (e - s - 1.0) * x[i + s]; // zero based
	    beta[2] += i / (e - s - 1.0) * (i - 1.0) / (e - s - 2.0) * x[i + s];
	}
    beta[0] /= e - s;
    beta[1] /= e - s;
    beta[2] /= e - s;
}



/***********************************************************************

  GEV PWM

***********************************************************************/

/* parameters of the GEV (PWM) */
int gev_pwm(double *x, double *param)//, int approx)
{
    double c = (2.0 * x[1] - x[0]) / (3.0 * x[2] - x[0]) - log(2) / log(3);
    param[2] = -7.8590 * c - 2.9554 * c * c; //  xi
    /* double b = (3.0 * x[2] - x[0]) / (2.0 * x[1] - x[0]); */
    /* double f(double x) */
    /* { */
    /* 	return ( (R_pow(3.0,x) - 1.0) / (R_pow(2.0,x) - 1.0) - b ); */
    /* } */
    /* param[2]=zeroin(-5.0,5.0,f,0.0); */
    if (ISNAN(param[2]) || param[2] >= 1.0) // feasibility condition
	return(1);
    double gamma1xi = gammafn(1.0 - param[2]);
    double pow2xi = R_pow(2.0, param[2]);
    param[1] = (2.0 * x[1] - x[0]) * param[2] / (gamma1xi * (pow2xi - 1.0)); // sigma
    if (ISNAN(param[1]) || param[1] <= 0.0) // feasibility condition
	return(1);
    param[0] =  x[0] + (2.0 * x[1] - x[0]) * ( gamma1xi - 1.0 ) / (gamma1xi * (1.0 - pow2xi)); //mu
    if (ISNAN(param[0]))
	return(1);
    else
	return(0);
}

/***********************************************************************

  c <- "( (2 * x1 - x0) / (3 * x2 - x0) - log(2)/log(3) )"
  xi <- parse( text = paste( "-7.8590 * ", c, " - 2.9554 * ", c, "^2" ) )
  x <- paste("x", 0:2, sep="")
  deriv.xi <- deriv(xi, x)

***********************************************************************/

/* gradient */
void grad_xi_gev_pwm(double *x, double *grad)
{
    double expr3 = 2 * x[1] - x[0];
    double expr5 = 3 * x[2] - x[0];
    double expr10 = expr3/expr5 - log(2)/log(3);
    double expr16 = expr5 * expr5;
    double expr18 = 1/expr5 - expr3/expr16;
    double expr24 = 2/expr5;
    double expr32 = expr3 * 3/expr16;

    grad[0] = 7.859 * expr18 + 2.9554 * (2 * (expr18 * expr10));

    grad[1] = -(7.859 * expr24 + 2.9554 * (2 * (expr24 * expr10)));

    grad[2] = 7.859 * expr32 + 2.9554 * (2 * (expr32 * expr10));
}

/***********************************************************************

  c <- "( (2 * x1 - x0) / (3 * x2 - x0) - log(2)/log(3) )"
  xi <- paste( "( -7.8590 * ", c, " - 2.9554 * ", c, "^2 )" )
  denom <- paste("( gamma(1 -", xi, ") * (1 - 2^", xi, ") )")
  sigma <- parse( text = paste(" - (2 * x1 - x0) * ", xi, "/", denom) )
  x <- paste("x", 0:2, sep="")
  deriv.sigma <- deriv(sigma, x)

***********************************************************************/

/* gradient */
void grad_sigma_gev_pwm(double *x, double *grad)
{
    double expr2 = 2 * x[1] - x[0];
    double expr6 = 3 * x[2] - x[0];
    double expr8 = log(2);
    double expr11 = expr2/expr6 - expr8/log(3);
    double expr15 = -7.859 * expr11 - 2.9554 * expr11 * expr11;
    double expr16 = -expr2 * expr15;
    double expr17 = 1 - expr15;
    double expr18 = gammafn(expr17);
    double expr19 = R_pow(2.0,expr15);
    double expr20 = 1 - expr19;
    double expr21 = expr18 * expr20;
    double expr24 = expr6 * expr6;
    double expr26 = 1/expr6 - expr2/expr24;
    double expr31 = 7.859 * expr26 + 2.9554 * (2 * (expr26 * expr11));
    double expr39 = expr18 * digamma(expr17);
    double expr44 = expr21 * expr21;
    double expr47 = 2/expr6;
    double expr52 = 7.859 * expr47 + 2.9554 * (2 * (expr47 * expr11));
    double expr67 = expr2 * 3/expr24;
    double expr72 = 7.859 * expr67 + 2.9554 * (2 * (expr67 * expr11));

    grad[0] = (expr15 - expr2 * expr31)/expr21 + expr16
	* (expr18 * (expr19 * (expr8 * expr31)) + expr31 * expr39 * expr20)/expr44;

    grad[1] = (expr2 * expr52 - 2 * expr15)/expr21 - expr16
	* (expr52 * expr39 * expr20 + expr18 * (expr19 * (expr8 * expr52)))/expr44;

    grad[2] = -(expr2 * expr72/expr21 - expr16
		* (expr18 * (expr19 * (expr8 * expr72)) + expr72 * expr39 * expr20)/expr44);

}

/***********************************************************************

  c <- "( (2 * x1 - x0) / (3 * x2 - x0) - log(2)/log(3) )"
  xi <- paste( "( -7.8590 * ", c, " - 2.9554 * ", c, "^2 )" )
  denom <- paste("( gamma(1 -", xi, ") * (1 - 2^", xi, ") )")
  mu <- parse( text = paste(" x0 + (2 * x1 - x0) * ( gamma(1 -", xi, ") - 1) /", denom) )
  x <- paste("x", 0:2, sep="")
  deriv.mu <- deriv(mu, x)

***********************************************************************/

/* gradient */
void grad_mu_gev_pwm(double *x, double *grad)
{
    double expr2 = 2 * x[1] - x[0];
    double expr5 = 3 * x[2] - x[0];
    double expr7 = log(2);
    double expr10 = expr2/expr5 - expr7/log(3);
    double expr14 = -7.859 * expr10 - 2.9554 * expr10 * expr10;
    double expr15 = 1 - expr14;
    double expr16 = gammafn(expr15);
    double expr17 = expr16 - 1;
    double expr18 = expr2 * expr17;
    double expr19 = R_pow(2.0, expr14);
    double expr20 = 1 - expr19;
    double expr21 = expr16 * expr20;
    double expr25 = expr5 * expr5;
    double expr27 = 1/expr5 - expr2/expr25;
    double expr32 = 7.859 * expr27 + 2.9554 * (2 * (expr27 * expr10));
    double expr34 = expr16 * digamma(expr15);
    double expr35 = expr32 * expr34;
    double expr45 = expr21 * expr21;
    double expr50 = 2/expr5;
    double expr55 = 7.859 * expr50 + 2.9554 * (2 * (expr50 * expr10));
    double expr56 = expr55 * expr34;
    double expr69 = expr2 * 3/expr25;
    double expr74 = 7.859 * expr69 + 2.9554 * (2 * (expr69 * expr10));
    double expr75 = expr74 * expr34;

    grad[0] = 1.0 - ((expr2 * expr35 + expr17)/expr21 - expr18
		     * (expr16 * (expr19 * (expr7 * expr32)) + expr35 * expr20)/expr45);

    grad[1] = (2 * expr17 + expr2 * expr56)/expr21 - expr18
	* (expr56 * expr20 + expr16 * (expr19 * (expr7 * expr55)))/expr45;

    grad[2] = -(expr2 * expr75/expr21 - expr18
		* (expr16 * (expr19 * (expr7 * expr74)) + expr75 * expr20)/expr45);
}


/***********************************************************************

  GEV GPWM

***********************************************************************/

/* parameters of the GEV (GPWM) */
int gev_gpwm(double *x, double *param)
{
    double y = 2.0 * (x[0] - x[1])  / (x[0] - 9.0 * x[2] / 4.0);
    param[2] = (1.442853 - R_pow(-y, 0.4054651) ) / 0.1183375; // xi
    /* double f(double x) */
    /* { */
    /* 	return ( x / (1.0 - R_pow(1.5,x)) - y); */
    /* } */
    /* param[2]=zeroin(-5.0,5.0,f,0.0); */
    if (ISNAN(param[2]) || param[2] >= 2.0) // feasibility condition
	return(1);
    double gamma2xi = gammafn(2.0 - param[2]);
    param[1] = R_pow(2.0, 3.0 - param[2]) * (x[0] - x[1]) / gamma2xi; // sigma
    if (ISNAN(param[1]) || param[1] <= 0.0) // feasibility condition
	return(1);
    param[0] = param[1] / param[2]  * (1.0 - R_pow(2.0, param[2]) * gamma2xi) + 4.0 * x[0]; // mu
    if (ISNAN(param[0]))
	return(1);
    else
	return(0);
}

/***********************************************************************

  y <- "( 2.0 * (x0 - x1)  / (x0 - 9.0 * x2 / 4.0) )"
  xi <- parse( text = paste( "(1.442853 - (-",y,")^0.4054651) / 0.1183375 " ) )
  x <- paste("x", 0:2, sep="")
  deriv.xi <- deriv(xi, x)

***********************************************************************/

/* gradient */
void grad_xi_gev_gpwm(double *x, double *grad)
{
    double expr2 = 2.0 * (x[0] - x[1]);
    double expr5 = x[0] - 9.0 * x[2] / 4.0;
    double expr7 = -(expr2/expr5);
    double expr11 = 2.0/expr5;
    double expr12 = expr5 * expr5;
    double expr15 = R_pow(expr7,-0.5945349);

    grad[0] = 0.4054651 * ((expr11 - expr2/expr12) * expr15)/0.1183375;

    grad[1] = -(0.4054651 * (expr11 * expr15)/0.1183375);

    grad[2] = 0.4054651 * (expr2 * (9.0/4.0)/expr12 * expr15)/0.1183375;
}

/***********************************************************************

  y <- "( 2.0 * (x0 - x1)  / (x0 - 9.0 * x2 / 4.0) )"
  xi <- paste( "( (1.442853 - (-",y,")^0.4054651)/0.1183375 )" )
  sigma <- parse( text = paste("2^(3 -", xi, ") * (x0 - x1) / gamma(2 -", xi, ")") )
  x <- paste("x", 0:2, sep="")
  deriv.sigma <- deriv(sigma, x)

***********************************************************************/

/* gradient */
void grad_sigma_gev_gpwm(double *x, double *grad)
{
    double expr1 = x[0] - x[1];
    double expr2 = 2.0 * expr1;
    double expr5 = x[0] - 9.0 * x[2]/4.0;
    double expr7 = -(expr2/expr5);
    double expr10 = (1.442853 - R_pow(expr7,0.4054651))/0.1183375;
    double expr12 = R_pow(2.0,3 - expr10);
    double expr13 = expr12 * expr1;
    double expr14 = 2 - expr10;
    double expr15 = gammafn(expr14);
    double expr17 = log(2);
    double expr18 = 2/expr5;
    double expr19 = expr5 * expr5;
    double expr22 = R_pow(expr7,-0.5945349);
    double expr25 = 0.4054651 * ((expr18 - expr2/expr19) * expr22)/0.1183375;
    double expr32 = expr15 * digamma(expr14);
    double expr35 = expr15 * expr15;
    double expr40 = 0.4054651 * (expr18 * expr22)/0.1183375;
    double expr55 = 0.4054651 * (expr2 * (9.0/4.0)/expr19 * expr22)/0.1183375;

    grad[0] =  (expr12 - expr12 * (expr17 * expr25) * expr1)/expr15
	+ expr13 * (expr25 * expr32)/expr35;

    grad[1] =  (expr12 * (expr17 * expr40) * expr1 - expr12)/expr15
	- expr13 * (expr40 * expr32)/expr35;

    grad[2] = -(expr12 * (expr17 * expr55) * expr1/expr15
		- expr13 * (expr55 * expr32)/expr35);
}


/***********************************************************************

  y <- "( 2.0 * (x0 - x1)  / (x0 - 9.0 * x2 / 4.0) )"
  xi <- paste( "( (1.442853 - (-",y,")^0.4054651)/0.1183375 )" )
  sigma <- paste("( 2^(3 -", xi, ") * (x0 - x1) / gamma(2 -", xi, ") )")
  mu <- parse( text = paste( sigma, "/", xi, "* (1 - 2^", xi, "* gamma(2 -",xi,")) + 4 * x0"))
  x <- paste("x", 0:2, sep="")
  deriv.mu <- deriv(mu, x)

***********************************************************************/

/* gradient */
void grad_mu_gev_gpwm(double *x, double *grad)
{
    double expr1 = x[0] - x[1];
    double expr2 = 2.0 * expr1;
    double expr5 = x[0] - 9.0 * x[2]/4.0;
    double expr7 = -(expr2/expr5);
    double expr10 = (1.442853 - R_pow(expr7,0.4054651))/0.1183375;
    double expr12 = R_pow(2.0,3.0 - expr10);
    double expr13 = expr12 * expr1;
    double expr14 = 2 - expr10;
    double expr15 = gammafn(expr14);
    double expr16 = expr13/expr15;
    double expr17 = expr16/expr10;
    double expr18 = R_pow(2,expr10);
    double expr20 = 1 - expr18 * expr15;
    double expr24 = log(2);
    double expr25 = 2/expr5;
    double expr26 = expr5 * expr5;
    double expr29 = R_pow(expr7,-0.5945349);
    double expr32 = 0.4054651 * ((expr25 - expr2/expr26) * expr29)/0.1183375;
    double expr33 = expr24 * expr32;
    double expr39 = expr15 * digamma(expr14);
    double expr40 = expr32 * expr39;
    double expr42 = expr15 * expr15;
    double expr47 = expr10 * expr10;
    double expr60 = 0.4054651 * (expr25 * expr29)/0.1183375;
    double expr61 = expr24 * expr60;
    double expr66 = expr60 * expr39;
    double expr86 = 0.4054651 * (expr2 * (9.0/4.0)/expr26 * expr29)/0.1183375;
    double expr87 = expr24 * expr86;
    double expr90 = expr86 * expr39;

    grad[0] = (((expr12 - expr12 * expr33 * expr1)/expr15 + expr13 * expr40/expr42)/expr10
	       - expr16 * expr32/expr47) * expr20 - expr17
	* (expr18 * expr33 * expr15 - expr18 * expr40) + 4.0;

    grad[1] = (((expr12 * expr61 * expr1 - expr12)/expr15 - expr13 * expr66/expr42)/expr10
	       + expr16 * expr60/expr47) * expr20 - expr17
	* (expr18 * expr66 - expr18 * expr61 * expr15);

    grad[2] = -(expr17 * (expr18 * expr87 * expr15 - expr18 * expr90)
		+ ((expr12 * expr87 * expr1/expr15 - expr13 * expr90/expr42)/expr10
		   + expr16 * expr86/expr47) * expr20);
}


/***********************************************************************
  Wrappers
***********************************************************************/

int g(double *x, double *param, int method)
{
    switch (method)
	    {
	    case 1 :
		return(gev_pwm(x, param));

	    case 2 :
		return(gev_gpwm(x, param));

	    default:
		Rprintf("Wrong method in g\n");
		return(1);
	    }
}

/* wrapper for gradient*/
void grad_g(double *x, int method, int statistic, double *grad)
{
    if (method == 1) //PWM
	switch (statistic)
	    {
	    case 0 :
		grad_mu_gev_pwm(x, grad);
		return;

	    case 1 :
		grad_sigma_gev_pwm(x, grad);
		return;

	    case 2 :
		grad_xi_gev_pwm(x, grad);
		return;

	    default:
		Rprintf("Error in grad_g\n");
		return;
	}
    else if (method == 2) //GPWM
	switch (statistic)
	    {
	    case 0 :
		grad_mu_gev_gpwm(x, grad);
		return;

	    case 1 :
		grad_sigma_gev_gpwm(x, grad);
		return;

	    case 2 :
		grad_xi_gev_gpwm(x, grad);
		return;

	    default:
		Rprintf("Error in grad_g\n");
		return;
	    }
    else
	{
	    Rprintf("Error in grad_g\n");
	    return;
	}
}

/***********************************************************************

  Change-point test based on PWM for block maxima with the GEV in mind
  gamman, deltan: for cdf from entire sample
  gamma, delta: for cdfs from subsamples

***********************************************************************/

void cpTestBM(double *X, int *n, int *r, double *stat, double *gamman,
	      double *deltan, double *gamma, double *delta,
	      int *meth, int *landwehr, int *noties, int* center,
	      double *param, double *avar)
{
    double *cdfn = Calloc(*n, double); // estimated cdfs from entire sample
    double *cdf = Calloc(*n, double); // estimated cdfs from subsamples
    double *x = Calloc(*n, double); // copy of data for sorting
    int *index = Calloc(*n, int); // for sorting
    int p = 3, survival = 0; // for omega, domega functions
    double *Y = Calloc((*n) * p, double); // influence functions
    double *beta = Calloc(p, double); // pwm estimates
    double *betak = Calloc(p, double); // pwm estimates
    double *betank = Calloc(p, double); // pwm estimates
    double *grad = Calloc(p * p, double); // for gradient of g
    double *cova = Calloc(p * p, double); // covariance matrix of influence terms
    double *a = Calloc(p, double); // powers of pwms
    double *b = Calloc(p, double); // powers of pwms
    double *paramk = Calloc(p, double); // parameter estimates from subsamples
    double *paramnk = Calloc(p, double); // parameter estimates from subsamples
    int i, j, k, l, m, failed, failedk, failednk;
    double meanj, meanl, s, sqrtn = sqrt(*n);

    /* define a and b according to st */
    switch (*meth)
	{

	case 1 : // PWM
	    a[0] = 0; a[1] = 1; a[2] = 2;
	    b[0] = 0; b[1] = 0; b[2] = 0;
	    break;

	case 2 : // GPWM
	    a[0] = 1; a[1] = 1; a[2] = 2;
	    b[0] = 1; b[1] = 2; b[2] = 1;
	    break;

	default:
	    Rprintf("Wrong method in cptestBM\n");
	    return;
	}

    /* copy data for sorting */
    for (i = 0; i < *n; i++)
	x[i] = X[i];
    /* estimate empirical cdfs from subsamples (does not depend on k) */
    ecdfs(x, index, *n, 0, *n, cdfn, *gamman, *deltan, *noties);

    if (*meth == 1 && *landwehr) // PWM with Landwehr
	{
	    /* copy data for sorting */
	    for (i = 0; i < *n; i++)
		x[i] = X[i];
	    estimate_landwehr(p, 0, *n, beta, x);
	}

    else /* estimate PWMs for each omega function */
	estimate_pwm(p, 0, *n, beta, X, cdfn, a, b, survival);

    /* estimate corresponding GEV parameters */
    failed = g(beta, param, *meth);
    /* feasibility criterion */
    if (failed)
	{
	    param[0] = 0.0; param[1] = 0.0; param[2] = 0.0;
	    Rprintf("Warning: invalid estimates from entire sample 1,...,%d\n",*n);
	}
    /* center with respect to estimated location parameter */
    if (*center)
	for (i = 0; i < *n; i++)
	    X[i] -= param[0];

    /* split each sample into k and n-k */
    for (k = *r; k <= *n - *r; k++)
	{
	    s = (double)k / *n;

	    /* copy data for sorting */
	    for (i = 0; i < *n; i++)
		x[i] = X[i];

	    if (*meth == 1 && *landwehr) // PWM with Landwehr
		{
		    estimate_landwehr(p, 0, k, betak, x);
		    estimate_landwehr(p, k, *n, betank, x);
		}
	    else
		{
		    /* estimate PWMs from k and n-k, for each omega function */
		    ecdfs(x, index, *n, 0, k, cdf, *gamma, *delta, *noties);
		    ecdfs(x, index, *n, k, *n, cdf, *gamma, *delta, *noties);
		    estimate_pwm(p, 0, k, betak, X, cdf, a, b, survival);
		    estimate_pwm(p, k, *n, betank, X, cdf, a, b, survival);
		}

	    /* level k statistics */
	    failedk = g(betak, paramk, *meth);
	    failednk = g(betank, paramnk, *meth);
	    /* feasibility criterion */
	    if (failedk)
		Rprintf("Warning: invalid estimates from subsample 1,...,%d\n",k);
	    /* feasibility criterion */
	    if (failednk)
		Rprintf("Warning: invalid estimates from subsample %d,...,%d\n",k+1,*n);
	    if (failedk || failednk)
		{
		    paramk[0] = 0.0; paramk[1] = 0.0; paramk[2] = 0.0;
		    paramnk[0] = 0.0; paramnk[1] = 0.0; paramnk[2] = 0.0;
		}
	    for (l = 0; l < p; l++)
		stat[k - *r + l * (*n - (2 * (*r) - 1))] = sqrtn * s * (1.0 - s)
		    * fabs(paramk[l] - paramnk[l]);
	}

    /* estimate influence functions from n, for each omega function */
    for (l = 0; l < p; l++)
	{
	    for (i = 0; i < *n; i++)
		{
		    Y[i + l * (*n)] = 0.0;
		    for (j = 0; j < *n; j++)
			if (X[i] <= X[j])
			    Y[i + l * (*n)] += X[j] * domega(cdfn[j],
							     a[l], b[l], survival);
		    Y[i + l * (*n)] /= *n;
		    Y[i + l * (*n)] += X[i] * omega(cdfn[i], a[l], b[l], survival);
		}
	}

    /* covariances */
    for (l = 0; l < p; l++)
	for (j = 0; j <= l; j++)
	    {
		meanl = 0.0;
		meanj = 0.0;
		for (i = 0; i < *n; i++)
		    {
			meanl += Y[i + l * (*n)];
			meanj += Y[i + j * (*n)];
		    }
		meanl /= *n;
		meanj /= *n;

		cova[l + j * p] = 0.0;
		for (i = 0; i < *n; i++)
		    cova[l + j * p] += (Y[i + j * (*n)] - meanj) * (Y[i + l * (*n)] - meanl);
		//cova[l + j * p] /= *n; // not divided because of pkolmogorov1x
		cova[j + l * p] = cova[l + j * p];
	    }

    /* compute gradients from entire sample */
    for (l = 0; l < p; l++) /* for 1, 2, 3 PWMs */
	grad_g(beta, *meth, l, &grad[l * p]);

	    /* asymptotic variances */
    for (l = 0; l < p; l++)
	{
	    avar[l] = 0.0;
	    for (j = 0; j < p; j++)
		for (m = 0; m < p; m++)
		    avar[l] += grad[j + l * p] * grad[m + l * p] * cova[j + m * p];
	}


    Free(cdf);
    Free(cdfn);
    Free(x);
    Free(index);
    Free(Y);
    Free(beta);
    Free(betak);
    Free(betank);
    Free(grad);
    Free(cova);
    Free(a);
    Free(b);
    Free(paramk);
    Free(paramnk);
}

/***********************************************************************

  Fit GEV by PWMs or GPWMs

***********************************************************************/

void fitGEV(double *X, int *n, double *gamma, double *delta, int *meth,
	    int *landwehr, int *noties, double *param, double *avar)
{
    double *cdfn = Calloc(*n, double); // estimated cdfs from entire sample
    double *x = Calloc(*n, double); // copy of data for sorting
    int *index = Calloc(*n, int); // for sorting
    int p = 3, survival = 0; // for omega, domega functions
    double *Y = Calloc((*n) * p, double); // influence functions
    double *beta = Calloc(p, double); // pwm estimates
    double *grad = Calloc(p * p, double); // for gradient of g
    double *cova = Calloc(p * p, double); // covariance matrix of influence terms
    double *a = Calloc(p, double); // powers of pwms
    double *b = Calloc(p, double); // powers of pwms
    int i, j, l, m, failed;
    double meanj, meanl;

    /* define a and b according to st */
    switch (*meth)
	{

	case 1 : // PWM
	    a[0] = 0; a[1] = 1; a[2] = 2;
	    b[0] = 0; b[1] = 0; b[2] = 0;
	    break;

	case 2 : // GPWM
	    a[0] = 1; a[1] = 1; a[2] = 2;
	    b[0] = 1; b[1] = 2; b[2] = 1;
	    break;

	default:
	    Rprintf("Wrong statistics in cpTestBM\n");
	    return;
	}

    /* copy data for sorting */
    for (i = 0; i < *n; i++)
	x[i] = X[i];
    /* estimate empirical cdf from sample */
    ecdfs(x, index, *n, 0, *n, cdfn, *gamma, *delta, *noties);

    if (*meth == 1 && *landwehr) // PWM with Landwehr
	{
	    /* copy data for sorting */
	    for (i = 0; i < *n; i++)
		x[i] = X[i];
	    estimate_landwehr(p, 0, *n, beta, x);
	}
    else // estimate PWMs for each omega function
	estimate_pwm(p, 0, *n, beta, X, cdfn, a, b, survival);

    /* estimate corresponding GEV parameters */
    failed = g(beta, param, *meth);
    /* feasibility criterion */
    if (failed)
	{
	    param[0] = 0.0; param[1] = 0.0; param[2] = 0.0;
	    Rprintf("Warning: invalid estimates\n");
	}

    /* estimate influence functions from n, for each omega function */
    for (l = 0; l < p; l++)
	{
	    for (i = 0; i < *n; i++)
		{
		    Y[i + l * (*n)] = 0.0;
		    for (j = 0; j < *n; j++)
			if (X[i] <= X[j])
			    Y[i + l * (*n)] += X[j] * domega(cdfn[j],
							     a[l], b[l], survival);
		    Y[i + l * (*n)] /= *n;
		    Y[i + l * (*n)] += X[i] * omega(cdfn[i], a[l], b[l], survival);
		}
	}

    /* covariances */
    for (l = 0; l < p; l++)
	for (j = 0; j <= l; j++)
	    {
		meanl = 0.0;
		meanj = 0.0;
		for (i = 0; i < *n; i++)
		    {
			meanl += Y[i + l * (*n)];
			meanj += Y[i + j * (*n)];
		    }
		meanl /= *n;
		meanj /= *n;

		cova[l + j * p] = 0.0;
		for (i = 0; i < *n; i++)
		    cova[l + j * p] += (Y[i + j * (*n)] - meanj) * (Y[i + l * (*n)] - meanl);
		cova[l + j * p] /= *n;
		cova[j + l * p] = cova[l + j * p];
	    }

    /* compute gradients from entire sample */
    for (l = 0; l < p; l++) /* for 1, 2, 3 PWMs */
	grad_g(beta, *meth, l, &grad[l * p]);

    /* asymptotic variances */
    for (l = 0; l < p; l++)
	{
	    avar[l] = 0.0;
	    for (j = 0; j < p; j++)
		for (m = 0; m < p; m++)
		    avar[l] += grad[j + l * p] * grad[m + l * p] * cova[j + m * p];
	}


    Free(cdfn);
    Free(x);
    Free(index);
    Free(Y);
    Free(beta);
    Free(grad);
    Free(cova);
    Free(a);
    Free(b);
}


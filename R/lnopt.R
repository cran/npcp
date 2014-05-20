#################################################################################
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
#################################################################################

#################################################################################
## Utility functions
#################################################################################

bartlett <- function(x)
{
    pmax(1 - abs(x), 0)
}

parzen <- function(x)
{
    ifelse(abs(x) <= 1/2, 1 - 6 * x^2 + 6 * abs(x)^3,
           ifelse(1/2 <= abs(x) & abs(x) <= 1, 2 * (1 - abs(x))^3, 0))
}

pdfsumunif <- function(x,n)
{
    nx <- length(x)

    .C("pdf_sum_unif",
       as.integer(n),
       as.double(x),
       as.integer(nx),
       pdf = double(nx),
       PACKAGE = "npcp")$pdf
}

convrect <- function(x, n)
{
    pdfsumunif(x + n/2, n) / pdfsumunif(n / 2, n)
}

flattop <- function(x, a=0.5)
{
    pmin( pmax((1-abs(x))/(1-a), 0), 1)
}

#################################################################################
## Functions for determining L
#################################################################################

## Adapted from Matlab code by A. Patton and the R translation
## by C. Parmeter and J. Racine
mval <- function(rho, lagmax, kn, rho.crit)
{
    ## Compute the number of insignificant runs following each rho(k),
    ## k=1,...,lagmax.
    num.ins <- sapply(1:(lagmax-kn+1),
                      function(j) sum((abs(rho) < rho.crit)[j:(j+kn-1)]))

    ## If there are any values of rho(k) for which the kn proceeding
    ## values of rho(k+j), j=1,...,kn are all insignificant, take the
    ## smallest rho(k) such that this holds (see footnote c of
    ## Politis and White for further details).
    if(any(num.ins == kn))
        return( which(num.ins == kn)[1] )
    else
    {
      ## If no runs of length kn are insignificant, take the smallest
      ## value of rho(k) that is significant.
      if(any(abs(rho) > rho.crit))
      {
          lag.sig <- which(abs(rho) > rho.crit)
          k.sig <- length(lag.sig)

          if(k.sig == 1)
              ## When only one lag is significant, mhat is the sole
              ## significant rho(k).
              return( lag.sig )
          else
              ## If there are more than one significant lags but no runs
              ## of length kn, take the largest value of rho(k) that is
              ## significant.
              return( max(lag.sig) )
      }
      else
          ## When there are no significant lags, mhat must be the
          ## smallest positive integer (footnote c), hence mhat is set
          ## to one.
          return( 1 )
  }
}

Lval <- function(x, method=mean)
{
    n <- nrow(x)
    d <- ncol(x)

    ## parameters for adapting the approach of Politis and White (2004)
    kn <- max(5, ceiling(log10(n)))
    lagmax <- ceiling(sqrt(n)) + kn
    rho.crit <- 1.96 * sqrt(log10(n)/n)

    m <- numeric(d)
    for (i in 1:d)
    {
        rho <- acf(x[,i], lag.max = lagmax, type = "correlation",
                   plot = FALSE)$acf[-1]
        m[i] <- mval(rho, lagmax, kn, rho.crit)
    }
    return( 2 * method(m) )
}

#################################################################################
## Function for estimating b.opt = round( (ln.opt + 1)/ 2)
#################################################################################

b.opt <- function(x, m=5, weights = c("parzen", "bartlett"),
                  combine.method=c("max","median","mean","min"))
{
    weights <- match.arg(weights)
    combine.method <- match.arg(combine.method)
    kernel <- switch(weights,
                     bartlett = parzen,
                     parzen = function(x) convrect(x*4,8))
    method <- switch(combine.method,
                     min = min,
                     median = median,
                     mean = mean,
                     max = max)
    stopifnot(m > 0)
    stopifnot(is.matrix(x))
    n <- nrow(x)
    d <- ncol(x)
    U <- apply(x,2,rank)/(n+1)

    ## parameters for adapting the approach of Politis and White (2004)
    kn <- max(5, ceiling(log10(n)))
    lagmax <- ceiling(sqrt(n)) + kn
    rho.crit <- 1.96 * sqrt(log10(n)/n)

    ## make grid -- suboptimal
    z <- seq(1/(m+1), 1 - 1/(m+1), len = m)
    v <- vector("list",d)
    for (i in 1:d)
        v[[i]] <- z
    g <- as.matrix(expand.grid(v)) ## grid on [0,1]^d
    ng <- nrow(g)


    ## compute gamma.n
    gamma.n <- array(NA, c(ng, ng, 2*lagmax+1))
    for (i in 1:ng)
        for (j in 1:i)
        {
            gamma.n[i,j,] <- as.numeric(ccf(apply(ifelse(U <= g[i,],1,0),1,prod),
                                            apply(ifelse(U <= g[j,],1,0),1,prod),
                                            lag.max = lagmax, type = "covariance",
                                            plot = FALSE)$acf)
            ## use "symmetry" of ccf
            gamma.n[j,i,] <- gamma.n[i,j,(2*lagmax+1):1]
        }

    ## determine L
    L <- Lval(x, method=method)

    ## compute K.n and sigma.n for all u, v on the grid
    K.n <- sigma.n <- matrix(0,ng,ng)
    for (i in 1:ng)
        for (j in 1:ng)
        {
            ft <- flattop(-lagmax:lagmax/L)
            sigma.n[i,j] <- sum(ft * gamma.n[i,j,])
            K.n[i,j] <- sum(ft * (-lagmax:lagmax)^2 * gamma.n[i,j,])
        }

    ## compute Gamma.n and Delta.n
    Gamma.n.2 <- as.numeric(hessian(kernel,0))^2 / 4 * mean(K.n^2)

    kernel2 <- function(x) kernel(x)^2
    Delta.n <- integrate(kernel2,-1,1)$value *
        (mean(diag(sigma.n))^2 + mean(sigma.n^2))
    ln.opt <- (4 * Gamma.n.2 / Delta.n * n)^(1/5)
    return( round((ln.opt + 1) / 2) )
}


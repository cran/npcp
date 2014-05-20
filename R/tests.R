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
## Functions necessary for generating dependent multipliers
#################################################################################

bartlettweights <- function(b)
{
    (1 - abs(-(b-1):(b-1))/b)/b * sqrt(3*b^3/(2*b^2+1))
}

parzenweights <- function(b)
{
    w <- parzen(-(b-1):(b-1)/b)
    w / sqrt(sum(w^2)) ## variance correction
}

#################################################################################
## Some change-point tests based on the empirical cdfs
#################################################################################

cpTestFn <- function(x, statistic = c("cvmmax", "cvmmean", "ksmax", "ksmean"),
                     mult.method = c("nonseq", "seq"), b = 1,
                     weights = c("parzen", "bartlett"),
                     m = 5, combine.method=c("max","median","mean","min"),
                     N = 1000, init.seq = NULL)
{
    statistic <- match.arg(statistic)
    mult.method <- match.arg(mult.method)
    weights <- match.arg(weights)
    combine.method <- match.arg(combine.method)

    stopifnot(is.matrix(x))
    d <- ncol(x)
    n <- nrow(x)
    npb <- n - 1 # number of possible breakpoints

    if (is.null(b))
        b <- b.opt(x, m=m, weights=weights,
                   combine.method=combine.method)
    stopifnot(b >= 1)
    w <- switch(weights,
                bartlett = bartlettweights(b),
                parzen = parzenweights(b))

    ## initial standard normal sequence for generating dependent multipliers
    if (is.null(init.seq))
        init.seq <- rnorm(N * (n + 2 * (b - 1)))
    else
        stopifnot(length(init.seq) == N * (n + 2 * (b - 1)))

    ## test
    out <- .C("cptestF",
              as.double(x),
              as.integer(n),
              as.integer(d),
              cvm = double(npb),
              ks = double(npb),
              as.integer(N),
              as.double(w),
              as.integer(b),
              as.integer(mult.method == "seq"),
              cvm0 = double(N * npb),
              ks0 = double(N * npb),
              as.double(init.seq),
              PACKAGE = "npcp")

    ## for CVM statistics
    cvm <- out$cvm
    cvm0 <- matrix(out$cvm0,N,npb)

    ## for KS statistics
    ks <- out$ks
    ks0 <- matrix(out$ks0,N,npb)

    pval <- function(phi,s0,s)
        ( sum( apply(s0,1,phi) >= phi(s) ) + 0.5 ) / (N + 1)

    statistics <- c(cvmmax = max(cvm), cvmmean = sum(cvm)/n,
                    ksmax = max(ks), ksmean = sum(ks)/n)
    p.values <- c(cvmmax = pval(max,cvm0,cvm),
                  cvmmean = pval(sum,cvm0,cvm),
                  ksmax = pval(max,ks0,ks),
                  ksmean = pval(sum,ks0,ks))

    structure(class = "htest",
              list(method = sprintf("Test for change-point detection based on the multivariate empirical c.d.f. with 'mult.method'=\"%s\"", mult.method),
                   statistic = statistics[statistic],
                   p.value =  p.values[statistic],
                   cvm = c(cvm=cvm), ks = c(ks=ks),
                   all.statistics = statistics,
                   all.p.values = p.values, b = c(b=b),
                   data.name = deparse(substitute(x))))
}


#################################################################################
## Change-point tests based on the empirical copula
#################################################################################

cpTestCn <- function(x, mult.method = c("seq", "nonseq"), b = 1,
                     weights = c("parzen", "bartlett"), m = 5,
                     combine.method=c("max","median","mean","min"),
                     N = 1000, init.seq = NULL)
{
    mult.method <- match.arg(mult.method)
    weights <- match.arg(weights)
    combine.method <- match.arg(combine.method)

    stopifnot(is.matrix(x))
    d <- ncol(x)
    n <- nrow(x)
    npb <- n - 1 # number of possible breakpoints

    if (is.null(b))
        b <- b.opt(x, m=m, weights=weights,
                   combine.method=combine.method)
    stopifnot(b >= 1)
    w <- switch(weights,
                bartlett = bartlettweights(b),
                parzen = parzenweights(b))

    ## initial standard normal sequence for generating dependent multipliers
    if (is.null(init.seq))
        init.seq <- rnorm(N * (n + 2 * (b - 1)))
    else
        stopifnot(length(init.seq) == N * (n + 2 * (b - 1)))

    ## test
    out <- .C("cpTestC",
              as.double(x),
              as.integer(n),
              as.integer(d),
              cvm = double(npb),
              as.integer(N),
              as.double(w),
              as.integer(b),
              as.integer(mult.method == "seq"),
              cvm0 = double(N * npb),
              as.double(init.seq),
              PACKAGE = "npcp")

    cvm <- out$cvm
    cvm0 <- matrix(out$cvm0,N,npb)
    statistic <- max(cvm)

    structure(class = "htest",
              list(method = sprintf("Test for change-point detection based on the empirical copula with 'mult.method'=\"%s\"", mult.method),
                   statistic = statistic,
                   p.value =  ( sum( apply(cvm0,1,max) >= statistic ) + 0.5 ) / (N + 1),
                   cvm = c(cvm=cvm), b = c(b=b),
                   data.name = deparse(substitute(x))))
}

#################################################################################
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
#################################################################################

#################################################################################
## Change-point tests based on the empirical dfs
#################################################################################

cpDist <- function(x, statistic = c("cvmmax", "cvmmean", "ksmax", "ksmean"),
                   method = c("nonseq", "seq"), b = NULL,
                   gamma = 0, delta = 1e-4, weights = c("parzen", "bartlett"),
                   m = 5, L.method = c("max","median","mean","min"),
                   N = 1000, init.seq = NULL, include.replicates = FALSE) {

    statistic <- match.arg(statistic)
    method <- match.arg(method)
    weights <- match.arg(weights)
    L.method <- match.arg(L.method)

    ## Power and constant in weight function
    stopifnot(gamma >= 0 && gamma <= 0.5)
    stopifnot(delta >= 0 && delta <= 1)

    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    d <- ncol(x)
    n <- nrow(x)
    npb <- n - 1 # number of possible breakpoints

    ## Bandwith parameter
    if (is.null(b))
        b <- bOptEmpProc(x, m = m, weights = weights,
                         L.method = L.method)
    stopifnot(b >= 1L)

    ## initial standard normal sequence for generating dependent multipliers
    if (is.null(init.seq))
        init.seq <- rnorm(N * (n + 2 * (b - 1)))
    else
        stopifnot(length(init.seq) >= N * (n + 2 * (b - 1)))

    ## test
    out <- .C("cpTestF",
              as.double(x),
              as.integer(n),
              as.integer(d),
              cvm = double(npb),
              ks = double(npb),
              as.integer(N),
              as.integer(weights == "bartlett"),
              as.integer(b),
              as.integer(method == "seq"),
              cvm0 = double(N * npb),
              ks0 = double(N * npb),
              as.double(init.seq),
              PACKAGE = "npcp")

    ## Weights
    s <- seq_len(npb) / n
    w <- 1 / pmax((s * (1-s))^gamma, delta)
    w.mat <- matrix(w, N, npb, byrow = TRUE)

    ## Weighted (replicates of) CVM statistics
    cvm <- w * out$cvm
    cvm0 <- w.mat * matrix(out$cvm0, N, npb)

    ## Weighted (replicates of KS) statistics
    ks <- w * out$ks
    ks0 <- w.mat * matrix(out$ks0, N, npb)

    pval <- function(s0, s) ( sum( s0 >= s ) + 0.5 ) / (N + 1)


    statistics <- c(cvmmax = max(cvm), cvmmean = sum(cvm)/n,
                    ksmax = max(ks), ksmean = sum(ks)/n)

    replicates <- cbind(cvmmax = apply(cvm0,1,max),
                        cvmmean = apply(cvm0,1,sum)/n,
                        ksmax = apply(ks0,1,max),
                        ksmean = apply(ks0,1,sum)/n)

    p.values <- c(cvmmax = pval(replicates[,"cvmmax"], statistics["cvmmax"]),
                  cvmmean = pval(replicates[,"cvmmean"], statistics["cvmmean"]),
                  ksmax = pval(replicates[,"ksmax"], statistics["ksmax"]),
                  ksmean = pval(replicates[,"ksmean"], statistics["ksmean"]))

    structure(class = "htest",
              list(method = sprintf("Test for change-point detection sensitive to changes in the distribution function with 'method'=\"%s\"", method),
                   statistic = statistics[statistic],
                   p.value =  p.values[statistic],
                   cvm = c(cvm=cvm), ks = c(ks=ks),
                   all.statistics = statistics,
                   all.p.values = p.values, b = c(b=b),
                   gamma = gamma, delta = delta,
                   all.replicates = if (include.replicates) replicates else NULL,
                   replicates = if (include.replicates) replicates[,statistic] else NULL,
                   data.name = deparse(substitute(x))))
}


#################################################################################
## Change-point tests based on the empirical copula
#################################################################################

cpCopula <- function(x, method = c("seq", "nonseq"), b = NULL,
                     weights = c("parzen", "bartlett"), m = 5,
                     L.method=c("max","median","mean","min"),
                     N = 1000, init.seq = NULL, include.replicates = FALSE) {

    method <- match.arg(method)
    weights <- match.arg(weights)
    L.method <- match.arg(L.method)

    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    d <- ncol(x)
    stopifnot(d > 1L)
    n <- nrow(x)
    npb <- n - 1 # number of possible breakpoints

    if (is.null(b))
        b <- bOptEmpProc(x, m=m, weights=weights,
                         L.method=L.method)
    stopifnot(b >= 1L)

    ## initial standard normal sequence for generating dependent multipliers
    if (is.null(init.seq))
        init.seq <- rnorm(N * (n + 2 * (b - 1)))
    else
        stopifnot(length(init.seq) >= N * (n + 2 * (b - 1)))

    ## test
    out <- .C("cpTestC",
              as.double(x),
              as.integer(n),
              as.integer(d),
              cvm = double(npb),
              as.integer(N),
              as.integer(weights == "bartlett"),
              as.integer(b),
              as.integer(method == "seq"),
              cvm0 = double(N * npb),
              as.double(init.seq),
              PACKAGE = "npcp")

    cvm <- out$cvm
    cvm0 <- matrix(out$cvm0,N,npb)
    statistic <- c(cvmmax=max(cvm))
    replicates <- apply(cvm0,1,max)

    structure(class = "htest",
              list(method = sprintf("Test for change-point detection sensitive to changes in the copula with 'method'=\"%s\"", method),
                   statistic = statistic,
                   p.value =  ( sum( replicates >= statistic ) + 0.5 ) / (N + 1),
                   cvm = c(cvm=cvm), b = c(b=b),
                   replicates = if (include.replicates) replicates else NULL,
                   data.name = deparse(substitute(x))))
}


#################################################################################
## Related to the df of the KS statistic
## From the source of ks.test -- credit to Rcore
#################################################################################

pkolmogorov1x <- function(x, n) {

    if (x <= 0)
        return(0)
    if (x >= 1)
        return(1)
    j <- seq.int(from = 0, to = floor(n * (1 - x)))
    1 - x * sum(exp(lchoose(n, j) + (n - j) * log(1 - x - j/n) + (j - 1) *
                        log(x + j/n)))
}

#################################################################################
## Change-point tests based on multivariate extensions of Spearman's rho
#################################################################################

cpRho <- function(x, method = c("mult", "asym.var"),
                      ## method = c("seq", "nonseq", "asym.var"),
                      statistic = c("pairwise", "global"),
                      ## L.method = c("pseudo","max","median","mean","min"),
                      b = NULL, weights = c("parzen", "bartlett"),
                      N = 1000, init.seq = NULL, include.replicates = FALSE) {

    method <- match.arg(method)
    statistic <- match.arg(statistic)
    weights <- match.arg(weights)
    ## L.method <- match.arg(L.method)
    L.method <- "pseudo"

    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    d <- ncol(x)
    stopifnot(d > 1L)
    n <- nrow(x)
    npb <- n - 1 # number of possible breakpoints

    if (is.null(b)) {

        res <- bOptRho(x, statistic=statistic, weights=weights, L.method=L.method)
        b <- res$b
        influ <- res$influnonseq
        fbin <- res$fbin
        influest <- TRUE
    }
    else {
        ## f in natural order
        f <- switch(statistic,
                    global = c(rep(0,2^d - 1),1),
                    pairwise = c(rep(0, d + 1), rep(1, choose(d,2)),
                    rep(0, 2^d - choose(d,2) - d - 1)))

        ## convert f into binary order
        powerset <-  .C("k_power_set",
                        as.integer(d),
                        as.integer(d),
                        powerset = integer(2^d),
                        PACKAGE="npcp")$powerset

        fbin <- .C("natural2binary",
                   as.integer(d),
                   as.double(f),
                   as.integer(powerset),
                   fbin = double(2^d),
                   PACKAGE="npcp")$fbin

        influ <- double(n)
        influest <- FALSE
    }
    stopifnot(b >= 1L)

    m <- switch(method,
                "mult" = 1, ##"seq" = 1,
                ## "nonseq" = 2,
                "asym.var" = 3)

    ## if (method %in% c("seq", "nonseq"))
    if (method == "mult") {
        ## initial standard normal sequence for generating dependent multipliers
        if (is.null(init.seq))
            init.seq <- rnorm(N * (n + 2 * (b - 1)))
        else
            stopifnot(length(init.seq) >= N * (n + 2 * (b - 1)))
    }

    ## test
    out <- .C("cpTestRho",
              as.double(x),
              as.integer(n),
              as.integer(d),
              rho = double(npb),
              as.double(fbin),
              influ = as.double(influ),
              as.integer(influest),
              as.integer(N),
              as.integer(weights == "bartlett"),
              as.integer(b),
              as.integer(m),
              rho0 = double(N * npb),
              avar = double(1),
              as.double(init.seq),
              PACKAGE = "npcp")

    rho <- out$rho
    stat <- max(rho)
    replicates <- NULL

    ## if (method %in% c("seq", "nonseq"))
    if (method == "mult") {

        rho0 <- matrix(out$rho0,N,npb)
        replicates <- apply(rho0,1,max)
        p.value <- ( sum( replicates >= stat ) + 0.5 ) / (N + 1)
    }
    else {

        if (!(out$avar > .Machine$double.eps)) { ## OK?

            cat("b =",b,"\n")
            cat("influ",out$influ,"\n")
            stop("The asymptotic variance is numerically equal to zero.")
        }
        else {
            #if (n <= 100)
            p.value <- 2 * (1 - pkolmogorov1x(stat / sqrt(out$avar), n))
            #else
            #p.value <- 2 * exp(-2 * n * stat^2 / out$avar)
            p.value <- min(1, max(0, p.value)) ## guard as in ks.test
        }
    }

    structure(class = "htest",
              list(method = sprintf("Test for change-point detection sensitive to changes in Spearman's rho with 'statistic'=\"%s\" and 'method'=\"%s\"", statistic, method),
                   statistic = c(rhomax=stat),
                   p.value = p.value,
                   rho = c(rho=rho), b = c(b=b),
                   replicates = if (include.replicates) replicates else NULL,
                   data.name = deparse(substitute(x))))
}

#################################################################################
## Change-point tests based on the sample mean
#################################################################################

cpMean <- function(x, method = c("nonseq", "seq", "asym.var"),
                   b = NULL, weights = c("parzen", "bartlett"),
                   N = 1000, init.seq = NULL, include.replicates = FALSE) {

    method <- match.arg(method)
    weights <- match.arg(weights)

    if(!is.vector(x, "numeric")) {
        warning("coercing 'x' to a numeric.")
        stopifnot(is.double(x <- as.double(x)))
    }

    n <- length(x)
    npb <- n - 1 # number of possible breakpoints

    if (is.null(b))
        b <- bOpt(x, weights=weights)
    else
        stopifnot((b <- as.integer(b)) >= 1L)

    m <- switch(method,
                "seq" = 1,
                "nonseq" = 2,
                "asym.var" = 3)

    if (method %in% c("seq", "nonseq")) {
        ## initial standard normal sequence for generating dependent multipliers
        if (is.null(init.seq))
            init.seq <- rnorm(N * (n + 2 * (b - 1)))
        else
            stopifnot(length(init.seq) >= N * (n + 2 * (b - 1)))
    }

    ## test
    out <- .C("cpTestMean",
              as.double(x),
              as.integer(n),
              stat = double(npb),
              as.integer(N),
              as.integer(weights == "bartlett"),
              as.integer(b),
              as.integer(m),
              stat0 = double(N * npb),
              avar = double(1),
              as.double(init.seq),
              PACKAGE = "npcp")

    stat <- out$stat
    statistic <- max(stat)
    replicates <- NULL

    if (method %in% c("seq", "nonseq")) {
        stat0 <- matrix(out$stat0,N,npb)
        replicates <- apply(stat0,1,max)
        p.value <- ( sum( replicates >= statistic ) + 0.5 ) / (N + 1)
    }
    else {
        if (!(out$avar > .Machine$double.eps)) { ## OK?
            cat("b =", b, "\n")
            cat("x", x, "\n")
            stop("The asymptotic variance is numerically equal to zero")
        }
        else {
            #if (n <= 100)
            p.value <- 2 * (1 - pkolmogorov1x(statistic / sqrt(out$avar), n))
            #else
            #p.value <- 2 * exp(-2 * n * stat^2 / out$avar)
            p.value <- min(1, max(0, p.value)) ## guard as in ks.test
        }
    }

    structure(class = "htest",
              list(method = sprintf("Test for change-point detection sensitive to changes in the expectation with 'method'=\"%s\"", method),
                   statistic = c(statistic=statistic),
                   p.value = p.value,
                   statistics = c(stat=stat), b = c(b=b),
                   replicates = if (include.replicates) replicates else NULL,
                   data.name = deparse(substitute(x))))

}

#################################################################################
## Change-point tests based on U-statistics
#################################################################################

## internal function
.cpU <- function(x, h.func, name, method = c("nonseq", "seq", "asym.var"),
                     b = NULL, weights = c("parzen", "bartlett"),
                     N = 1000, init.seq = NULL, include.replicates = FALSE) {

    method <- match.arg(method)
    weights <- match.arg(weights)

    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    n <- nrow(x)
    d <- ncol(x)

    npb <- n - 3 # number of possible breakpoints

    h <- matrix(0,n,n)
    for (i in seq_len(n))
        for (j in  seq_len(i))
            if (i != j) {
                h[i,j] <- h.func(x[i,],x[j,])
                h[j,i] <- h[i,j]
            }

    influ <- colSums(h) / (n-1) ## h1.n without centering term

    if (is.null(b))
        b <- bOpt(influ, weights=weights)
    else
        stopifnot((b <- as.integer(b)) >= 1L)

    m <- switch(method,
                "seq" = 1,
                "nonseq" = 2,
                "asym.var" = 3)

    if (method %in% c("seq", "nonseq")) {
        ## initial standard normal sequence for generating dependent multipliers
        if (is.null(init.seq))
            init.seq <- rnorm(N * (n + 2 * (b - 1)))
        else
            stopifnot(length(init.seq) >= N * (n + 2 * (b - 1)))
    }

    ## test
    out <- .C("cpTestU",
              as.double(h),
              as.integer(n),
              as.double(influ),
              u = double(npb),
              as.integer(N),
              as.integer(weights == "bartlett"),
              as.integer(b),
              as.integer(m),
              u0 = double(N * npb),
              avar = double(1),
              as.double(init.seq),
              PACKAGE = "npcp")

    u <- out$u
    names(u) <- paste0("u",seq.int(2, n-2))
    stat <- max(u)
    replicates <- NULL

    if (method %in% c("seq", "nonseq")) {
        u0 <- matrix(out$u0,N,npb)
        replicates <- apply(u0,1,max)
        p.value <- ( sum( replicates >= stat ) + 0.5 ) / (N + 1)
    }
    else {
        if (!(out$avar > .Machine$double.eps)) { ## OK?
            cat("b =",b,"\n")
            cat("h1.n",influ,"\n")
            stop("The asymptotic variance is numerically equal to zero.")
        }
        else {
            #if (n <= 100)
            p.value <- 2 * (1 - pkolmogorov1x(stat / (2 * sqrt(out$avar)), n))
            #else
            #p.value <- 2 * exp(-2 * n * stat^2 / out$avar)
            p.value <- min(1, max(0, p.value)) ## guard as in ks.test
        }
    }

    structure(class = "htest",
              list(method = sprintf("Test for change-point detection sensitive to changes in %s with 'method'=\"%s\"", name, method),
                   statistic = c(umax=stat),
                   p.value = p.value,
                   u = u, b = c(b=b),
                   replicates = if (include.replicates) replicates else NULL,
                   data.name = deparse(substitute(x))))

}

#################################################################################
## Wrappers
#################################################################################

## Kendall's tau
cpTau <- function(x, method = c("seq", "nonseq", "asym.var"),
                      b = NULL, weights = c("parzen", "bartlett"),
                      N = 1000, init.seq = NULL, include.replicates = FALSE) {

    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot(ncol(x) > 1L)
    method <- match.arg(method)
    weights <- match.arg(weights)

    .cpU(x, h.func = function(x, y) prod(x < y) + prod(y < x),
         name = "Kendall's tau", method = method,
         b = b, weights = weights, N = N, init.seq = init.seq,
         include.replicates = include.replicates)
}

## Covariance
cpCov <- function(x, method = c("nonseq", "seq", "asym.var"),
                      b = NULL, weights = c("parzen", "bartlett"),
                      N = 1000, init.seq = NULL, include.replicates = FALSE) {

    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot(ncol(x) == 2L)

    .cpU(x, h.func = function(x, y) (x[1] - y[1]) * (x[2] - y[2]) / 2,
         name = "the covariance", method = method,
         b = b, weights = weights, N = N, init.seq = init.seq,
         include.replicates = include.replicates)
}

## Variance
cpVar <- function(x, method = c("nonseq", "seq", "asym.var"),
                      b = NULL, weights = c("parzen", "bartlett"),
                      N = 1000, init.seq = NULL, include.replicates = FALSE) {

    if(!is.vector(x, "numeric")) {
        warning("coercing 'x' to a numeric.")
        stopifnot(is.double(x <- as.double(x)))
    }
    x <- matrix(x)
    .cpU(x, h.func = function(x, y) (x - y)^2 / 2,
         name = "the variance", method = method,
         b = b, weights = weights, N = N, init.seq = init.seq,
         include.replicates = include.replicates)
}

## Autocovariance
cpAutocov <- function(x, lag = 1, method = c("nonseq", "seq", "asym.var"), b = NULL,
                          weights = c("parzen", "bartlett"), N = 1000, init.seq = NULL,
                          include.replicates = FALSE) {

    if(!is.vector(x, "numeric")) {
        warning("coercing 'x' to a numeric.")
        stopifnot(is.double(x <- as.double(x)))
    }
    stopifnot((lag <- as.integer(lag)) >= 1L)

    ## lagged data
    x <- cbind(x[-(1:lag)], x[1:(length(x) - lag)])

    .cpU(x, h.func = function(x, y) (x[1] - y[1]) * (x[2] - y[2]) / 2,
         name = paste("the autocovariance at lag", lag),
         method = method, b = b, weights = weights, N = N, init.seq = init.seq,
         include.replicates = include.replicates)
}

## Gini's mean difference
cpGini <- function(x, method = c("nonseq", "seq", "asym.var"),
                   b = NULL, weights = c("parzen", "bartlett"),
                   N = 1000, init.seq = NULL, include.replicates = FALSE) {

    if(!is.vector(x, "numeric")) {
        warning("coercing 'x' to a numeric.")
        stopifnot(is.double(x <- as.double(x)))
    }
    x <- matrix(x)
    .cpU(x, h.func = function(x, y) abs(x - y),
         name = "Gini's mean difference", method = method,
         b = b, weights = weights, N = N, init.seq = init.seq,
         include.replicates = include.replicates)
}

#################################################################################
## Change-point tests based on probability weighted moments
## for block maxima with the GEV in mind
#################################################################################

cpBlockMax <- function(x, method = c("pwm", "gpwm"), r=10) {

    stopifnot(is.vector(x, "numeric"))
    n <- length(x)
    r <- as.integer(r)
    stopifnot(r > 0L && n - (2 * r - 1) >= 1L)
    method <- match.arg(method)

    ## internal parameters (see paper)
    r <- 10 # omit breakpoints at the beginning and the end
    gamma.var <- -0.35 # for cdf from entire sample
    delta.var <- 0 # for cdf from entire sample
    gamma.stat <- 0 # for cdfs from subsamples
    delta.stat <- 0 # for cdfs from subsamples
    landwehr <- TRUE # use Landwehr's PWM estimates if method=="pwm"
    noties <- TRUE # continous distribution assumed
    center <- TRUE # center wrt to location parameter estimate

    npb <- n - (2 * r - 1) # number of possible breakpoints


    meth <- switch(method, "pwm" = 1, "gpwm" = 2)
    p <- 3 # number of statistics

    out <- .C("cpTestBM",
              as.double(x),
              as.integer(n),
              as.integer(r), # omit breakpoints at the beginning and the end
              stat = double(npb * p),
              as.double(gamma.var),
              as.double(delta.var),
              as.double(gamma.stat),
              as.double(delta.stat),
              as.integer(meth),
              as.integer(landwehr),
              as.integer(noties),
              as.integer(center),
              param = double(p),
              avar = double(p),
              PACKAGE = "npcp")

    ## statistics
    stat <- matrix(out$stat, npb, p)
    stat.loc <- stat[,1]
    names(stat.loc) <- paste0("loc",r:(n-r))
    stat.scale <- stat[,2]
    names(stat.scale) <- paste0("scale",r:(n-r))
    stat.shape <- stat[,3]
    names(stat.shape) <- paste0("shape",r:(n-r))
    maxstat.loc <- max(stat.loc)
    maxstat.scale <- max(stat.scale)
    maxstat.shape <- max(stat.shape)

    ## corresponding p-values
    if (any(out$avar <= .Machine$double.eps)) ## OK?
        stop("Some of the asymptotic variances are numerically equal to zero.")
    else {
        var.offset <- if (method=="pwm") c(0,10,20) else rep(0,3) ## variance offset
        out$avar <- out$avar * (var.offset + n)/n
        p.value.loc <- 2 * (1 - pkolmogorov1x(maxstat.loc / sqrt(out$avar[1]), n))
        p.value.loc <- min(1, max(0, p.value.loc)) ## guard as in ks.test
        p.value.scale <- 2 * (1 - pkolmogorov1x(maxstat.scale / sqrt(out$avar[2]), n))
        p.value.scale <- min(1, max(0, p.value.scale)) ## guard as in ks.test
        p.value.shape <- 2 * (1 - pkolmogorov1x(maxstat.shape / sqrt(out$avar[3]), n))
        p.value.shape <- min(1, max(0, p.value.shape)) ## guard as in ks.test
    }


    structure(class = "htest",
              list(method = sprintf("Test for change-point detection in the distribution of block maxima with 'method'=\"%s\"", method),
                   statistic = c(loc.stat = maxstat.loc, scale.stat = maxstat.scale, shape.stat = maxstat.shape),
                   #parameter = c(loc.param = out$param[1], scale.param = out$param[2], shape.param = out$param[3]),
                   parameter = c(p.value.loc = p.value.loc, p.value.scale = p.value.scale, p.value.shape = p.value.shape), ## improve
                   pvalues = c(p.value.loc = p.value.loc, p.value.scale = p.value.scale, p.value.shape = p.value.shape),
                   stats.loc = stat.loc,
                   stats.scale = stat.scale,
                   stats.shape = stat.shape,
                   data.name = deparse(substitute(x))))
}


#################################################################################
## NOT EXPORTED
#################################################################################

fitGEV <- function(x, method = c("pwm", "gpwm"), gamma=-0.35, delta=0,
                   landwehr = TRUE, noties = TRUE) {
    method <- match.arg(method)
    stopifnot(is.vector(x, "numeric"))
    stopifnot(is.vector(gamma, "numeric"))
    stopifnot(is.vector(delta, "numeric"))
    n <- length(x)
    meth <- switch(method, "pwm" = 1, "gpwm" = 2)
    p <- 3 # number of statistics

    out <- .C("fitGEV",
              as.double(x),
              as.integer(n),
              as.double(gamma),
              as.double(delta),
              as.integer(meth),
              as.integer(landwehr),
              as.integer(noties),
              param = double(p),
              avar = double(p),
              PACKAGE = "npcp")

    sderr <- sqrt(out$avar/n)

    list(parameters=c(loc = out$param[1], scale = out$param[2], shape = out$param[3]),
         sderrs = c(loc = sderr[1], scale = sderr[2], shape = sderr[3]))
}

#################################################################################
## Change-point tests in the mean
#################################################################################

cpMeanIID <- function(x, method = c("cusum", "permut"), N = 1000) {

    ## checks
    stopifnot(is.vector(x, "numeric") & all(!is.na(x)))
    method <- match.arg(method)

    n <- length(x)
    m <- mean(x)

    ## returns vector of n - 1 intermediate k-statistics
    statistics <- function(x) {
        sapply(1:(n-1), function(k) {
            abs(sum(x[1:k] - m)) * if (method == "permut") sqrt(n / (k * (n - k))) else 1 / sqrt(n)
        })
    }
    stats <- statistics(x)

    ## test statistic
    stat <- max(stats)

    ## simulate based on N permutations
    if (method == "permut") {
        stat0 <- replicate(N, max(statistics(sample(x))))
        p.value <- ( sum( stat0 >= stat ) + 0.5 ) / (N + 1)
    } else { ## use "estimated" asymptotic distribution
        p.value <- 2 * (1 - pkolmogorov1x(stat / sqrt(sum((x - m)^2)), n))
        p.value <- min(1, max(0, p.value)) ## guard as in ks.test
    }


    structure(class = "htest",
              list(method = "Test for change-point detection based on the sample mean",
                   statistic = c(statistic = stat),
                   statistics = c(statistic = stats),
                   p.value =  p.value,
                   data.name = deparse(substitute(x))))
}

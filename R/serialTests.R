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
## Change-point tests based on the serial empirical copula
#################################################################################

cpAutocop <- function(x, lag = 1, b = NULL, bivariate = FALSE,
                      weights = c("parzen", "bartlett"), m = 5,
                      N = 1000, init.seq = NULL, include.replicates = FALSE) {

    ## settings
    method <- "nonseq" # c("nonseq", "seq")
    L.method <- "max" #c("max","median","mean","min")

    ## checks
    weights <- match.arg(weights)
    stopifnot((lag <- as.integer(lag)) >= 1L)
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot((d <- ncol(x)) == 1L) # one-column matrix
    n <- nrow(x)

    ## bandwidth parameter
    if (is.null(b))
        b <- bOptEmpProc(x, m=m, weights=weights,
                         L.method=L.method)
    b <- max(b,(h <- lag + 1))
    stopifnot(b >= 1L)

    ## initial standard normal sequence for generating dependent multipliers
    if (is.null(init.seq))
        init.seq <- rnorm(N * (n - h + 1 + 2 * (b - 1)))
    else
        stopifnot(length(init.seq) >= N * (n - h + 1 + 2 * (b - 1)))

    npb <- n - h # number of possible breakpoints

    ## test
    out <- .C("cpTestAutocop",
              as.double(x),
              as.integer(n),
              as.integer(d),
              as.integer(h),
              cvm = double(npb),
              as.integer(N),
              as.integer(weights == "bartlett"),
              as.integer(b),
              as.integer(method == "seq"),
              cvm0 = double(N * npb),
              as.double(init.seq),
              as.integer(bivariate),
              PACKAGE = "npcp")

    cvm <- out$cvm
    cvm0 <- matrix(out$cvm0,N,npb)
    statistic <- c(cvmmax=max(cvm))
    replicates <- apply(cvm0,1,max)

    structure(class = "htest",
              list(method = if (bivariate)
                                sprintf("Test of change-point detection sensitive to changes in the bivariate serial copula at lag %s", lag)
                            else
                                sprintf("Test of change-point detection sensitive to changes in the %s-dimensional autocopula", h),
                   statistic = statistic,
                   p.value =  ( sum( replicates >= statistic ) + 0.5 ) / (N + 1),
                   cvm = c(cvm=cvm), b = c(b=b), h = h,
                   replicates = if (include.replicates) replicates else NULL,
                   data.name = deparse(substitute(x))))
}


#################################################################################
## Stationarity tests based on df and serial copula CUSUM tests
#################################################################################

## For univariate time series
stDistAutocop <- function(x, lag = 1, b = NULL, pairwise = FALSE,
                          weights = c("parzen", "bartlett"), m = 5, N = 1000) {
    ## checks
    weights <- match.arg(weights)
    stopifnot((lag <- as.integer(lag)) >= 1L)
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot(ncol(x) == 1L) # one-column matrix

    ## Initial random sequence; the constant 20 is non-optimal
    init.seq <- rnorm(N * 20 * length(x))

    ## Bandwith parameter
    if (is.null(b))
        b <- bOptEmpProc(x, m = m, weights = weights)
    stopifnot(b >= 1L)

    ## Df test
    res <- cpDist(x, b = b, weights = weights,
                  init.seq = init.seq, include.replicates = TRUE)

    p.values <- res$p.value
    replicates <- res$replicates

    ## Serial empirical copula tests
    if (pairwise) {
        for (i in 1:lag) {
            res <- cpAutocop(x, lag = i,
                             b = b, bivariate = TRUE,
                             weights = weights, init.seq = init.seq,
                             include.replicates = TRUE)
            p.values <- c(p.values, res$p.value)
            replicates <- cbind(replicates, res$replicates)
        }
        names(p.values) <- c("Dist", paste0("sCop", 1:lag))
    }
    else {
        res <- cpAutocop(x, lag=lag,
                         b = b, bivariate = FALSE,
                         weights = weights, init.seq = init.seq,
                         include.replicates = TRUE)
        p.values <- c(p.values, res$p.value)
        replicates <- cbind(replicates, res$replicates)
        names(p.values) <- c("Dist", paste0("Autocop", lag+1))
    }

    ## Fisher combined test
    combined <- combinePValues(p.values = p.values, replicates = replicates,
                               w = if (pairwise) c(1, rep(1/lag, lag)) else c(1, 1))

    structure(class = "htest",
              list(method = sprintf("Combined test of stationarity sensitive to changes in the distribution function and the %s-dimensional autocopula with 'pairwise'=\"%s\"", lag+1, pairwise),
                   statistic = c(statistic = combined[1]),
                   p.value =  c(p.value = combined[2]),
                   component.p.values = p.values, b = c(b=b), h = lag+1,
                   data.name = deparse(substitute(x))))
}



#################################################################################
## NOT EXPORTED BELOW
#################################################################################

#################################################################################
## Combination of dependent p-values with corresponding multiplier replicates
#################################################################################

combinePValues <- function(p.values, replicates, w = rep(1, length(p.values))) {

    ## Checks
    stopifnot(is.matrix(replicates))
    stopifnot(is.numeric(p.values))
    stopifnot(all(p.values < 1L) && all(p.values > 0L))
    m <- nrow(replicates)
    p <- ncol(replicates)
    stopifnot(p == length(p.values))
    if (any(apply(replicates, 2, anyDuplicated)))
        warning("the component samples of 'replicates' contain ties")

    ## Combination of dependent p-values
    pv0 <- ((m + 1 - apply(replicates, 2, rank)) + 0.5) / (m + 1)
    ## Weighted Fisher
    f <- apply(rbind(p.values, pv0), 1, function(x) - sum(w * log(x)))
    ## Return statistic and p-value
    c(f[1], ( sum(f[-1] >= f[1]) + 0.5 ) / (m + 1), use.names = FALSE)
    ## Weighted Stouffer
    ## s <- apply(rbind(p.values, pv0), 1, function(x) sum(w * qnorm(1 - x)))
    ## stouffer = ( sum(s[-1] >= s[1]) + 0.5 ) / (m + 1))
}

## For multivariate time series
multStTestFnCn <- function(x, #method = c("nonseq", "seq"),
                           lag = 1, b = NULL,
                           weights = c("parzen", "bartlett"),
                           m = 5, L.method = c("max","median","mean","min"),
                           N=1000) {

    ## Checks
    #method <- match.arg(method)
    weights <- match.arg(weights)
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot((d <- ncol(x)) > 1L)

    ## Bandwith parameter
    if (is.null(b))
        b <- bOptEmpProc(x, m = m, weights = weights,
                         L.method = L.method)
    stopifnot(b >= 1L)

    ## Initial random sequence
    init.seq <- rnorm(N * 20 * (n <- nrow(x)))

    ## Df test
    res <- cpDist(x, #method = method,
                  b = b, weights = weights,
                  init.seq = init.seq, include.replicates = TRUE)

    p.values <- res$p.value
    replicates <- res$replicates

    ## Copula test
    res <- cpCopula(x, #method = method,
                    b = b, weights = weights,
                    init.seq = init.seq, include.replicates = TRUE)
    p.values <- c(p.values, res$p.value)
    replicates <- cbind(replicates,res$replicates)

    names(p.values) <- c("Fn", "Cn")

    ## Pairwise serial empirical copula tests
    for (j in 1:d)
        for (k in 1:d) {
            if (j == k) {
                ## "Auto-correlation" terms
                for (i in 1:lag) {
                    res <- cpAutocop(x[, j, drop = FALSE], #method = method,
                                     lag = i,
                                     b = b, weights = weights, init.seq = init.seq,
                                     include.replicates = TRUE)
                    p.values <- c(p.values, res$p.value)
                    replicates <- cbind(replicates,res$replicates)
                }
            }
            else { ## "Cross-correlation" terms
                for (i in 1:lag) {
                    res <- cpCopula(cbind(x[-(1:lag), j], x[1:(n - lag), k]),
                                        #method = method,
                                    b = b, weights = weights,
                                    init.seq = init.seq, include.replicates = TRUE)
                    p.values <- c(p.values, res$p.value)
                    replicates <- cbind(replicates,res$replicates)
                }
            }
            names(p.values) <- c(names(p.values), paste0("sCn",j,".",k,"lag",1:lag))
        }

    ## Fisher combined test
    combpval <- combinePValues(p.values = p.values, replicates = replicates)
    c(p.values, combpval)
}


#################################################################################
## Stationarity test based on mean, variance, autocovariances CUSUM tests
#################################################################################

stAutocov <- function(x, method = c("nonseq", "seq"), lag = 1, b = NULL,
                          weights = c("parzen", "bartlett"),
                          include.mean = TRUE, N = 1000) {

    ## Checks
    method <- match.arg(method)
    weights <- match.arg(weights)
    if(!is.vector(x, "numeric")) {
        warning("coercing 'x' to a numeric.")
        stopifnot(is.double(x <- as.double(x)))
    }

    ## Initial random sequence
    init.seq <- rnorm(N * 20 * length(x))

    ## Estimate bandwidth
    if (is.null(b))
        b <- bOpt(x, weights = weights)
    stopifnot(b >= 1L)

    ## Variance test
    res <- cpVar(x, method = method, b = b, weights = weights,
                     init.seq = init.seq, include.replicates = TRUE)
    p.values <- res$p.value
    replicates <- res$replicates

    ## Autocovariance tests
    for (i in 1:lag) {
        res <- cpAutocov(x, method = method, lag = i, b = b, weights = weights,
                             init.seq = init.seq, include.replicates = TRUE)
        p.values <- c(p.values, res$p.value)
        replicates <- cbind(replicates,res$replicates)
    }

    ## Combined test: variance + autocovariances
    combpval.autocov <- combinePValues(p.values = p.values, replicates = replicates,
                                       w = c(1, rep(1/lag, lag)))

    ## Mean test
    if (include.mean) {
        res <- cpMean(x, method = method, b = b, weights = weights,
                          init.seq = init.seq, include.replicates = TRUE)
        p.values <- c(res$p.value, p.values)
        replicates <- cbind(res$replicates, replicates)
        ## Combined test; mean + variance + autocovariances
        combpval.all <- combinePValues(p.values = p.values, replicates = replicates,
                                       w = c(1, 1, rep(1/lag, lag)))
        names(p.values) <- c("mean", "var", paste0("acv",1:lag))
        c(p.values, acv = combpval.autocov, all = combpval.all)
    }
    else {
        names(p.values) <- c("var", paste0("acv",1:lag))
        c(p.values, fisher.acv = combpval.autocov)
    }
}

#################################################################################
## Form h-lagged data
## in: (n + h - 1) x d
## out: n x h x d
#################################################################################

lagged <- function(n, d, h, x, pairwise = FALSE) {

    res <- matrix(NA, n, d * h)
    for (i in 1:n)
        for (l in if (pairwise) c(0,(h-1)) else 0:(h-1))
            for (j in 0:(d-1))
                res[i + (j + l * d) * n] <- x[i + l + j * (n + h - 1)]
    res
}


#################################################################################
##  Stationarity tests based on Kendall's tau
#################################################################################

cpSerialTau <- function(x, method = c("seq", "nonseq", "asym.var"),
                            lag = 1, b = NULL, weights = c("parzen", "bartlett"),
                            N = 1000, init.seq = NULL) {

    method <- match.arg(method)
    weights <- match.arg(weights)
    stopifnot((lag <- as.integer(lag)) >= 1L)
    h <- lag + 1

    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    n <- nrow(x)
    d <- ncol(x)
    ## form lagged data
    x <- lagged(n-h+1, d, h, x)
    res <- cpTau(x, method = method,
                     b = b, weights = weights, N = N, init.seq = init.seq)

    structure(class = "htest",
              list(method = sprintf("Test of change-point detection based on the serial version of Kendall's tau with 'method'=\"%s\" and 'h'=\"%s\"", method, h),
                   statistic = res$statistic,
                   p.value = res$p.value,
                   u = res$u, b = res$b, h = h,
                   data.name = deparse(substitute(x))))
}

#################################################################################
##  Stationarity tests based on Spearman's rho
#################################################################################

cpSerialRho <- function(x, method = c("mult", "asym.var"),
                      statistic = c("pairwise", "global"),
                      h = 2, b = NULL, weights = c("parzen", "bartlett"),
                      N = 1000, init.seq = NULL) {

    method <- match.arg(method)
    statistic <- match.arg(statistic)
    weights <- match.arg(weights)
    stopifnot((lag <- as.integer(lag)) >= 1L)
    h <- lag + 1

    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    n <- nrow(x)
    d <- ncol(x)
    ## form lagged data
    x <- lagged(n-h+1, d, h, x)
    res <- cpRho(x, method = method, statistic = statistic,
                     b = b, weights = weights, N = N, init.seq = init.seq)

    structure(class = "htest",
              list(method = sprintf("Test of change-point detection based on the serial version of Spearman's rho with 'method'=\"%s\", 'statistic'=\"%s\" and 'h'=\"%s\"", statistic, method, h),
                   statistic = res$statistic,
                   p.value = res$p.value,
                   rho = res$rho, b = res$b, h = h,
                   data.name = deparse(substitute(x))))
}

#################################################################################
##  Stationarity test based on cpDist on lagged time-series
#################################################################################

stDist <- function(x, method = c("nonseq", "seq"),
                   lag = 1, b = NULL, weights = c("parzen", "bartlett"),
                   N = 1000, init.seq = NULL) {

    method <- match.arg(method)
    weights <- match.arg(weights)
    stopifnot((lag <- as.integer(lag)) >= 1L)
    h <- lag + 1

    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    n <- nrow(x)
    d <- ncol(x)
    ## form lagged data
    y <- lagged(n-h+1, d, h, x)
    res <- cpDist(y, method = method,
                  b = b, weights = weights, N = N, init.seq = init.seq)

    structure(class = "htest",
              list(method = sprintf("Test of change-point detection based on the serial version of empirical distribution function with 'method'=\"%s\" and 'h'=\"%s\"", method, h),
                   statistic = res$statistic,
                   p.value = res$p.value,
                   b = res$b, h = h,
                   data.name = deparse(substitute(y))))
}

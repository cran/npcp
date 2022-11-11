#################################################################################
##
##   R package npcp by Ivan Kojadinovic, Alex Verhoijsen Copyright (C) 2022
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


##' @title Data-adaptive point selection procedure from an initial grid of size r^d in [0,1]^d
##' @param x n x d matrix of data in [0,1]^d
##' @param r is the number of initial points per dimension
##' @param plot logical; only used if d = 2
##' @param kappa constant for box probabilities
##' @return a matrix of selected points (rows)

selectPoints <- function(x, r, kappa = 1.5, plot = FALSE) {

    if(!is.matrix(x))
        stopifnot(is.matrix(x <- as.matrix(x)))

    ## Transposed pseudo-observations (to ease comparison below) corresonding to x
    m <- nrow(x)
    d <- ncol(x)
    U <- t(apply(x, 2, rank) / (m + 1))

    ## We will consider boxes of volume s^d
    s <- 1 / (r + 1)

    ## Set cutoff probability q
    q <- s^d / kappa # 1 / ((r + 1)^d * kappa)

    ## Consider an initial uniform grid covering of [0,1]^d
    u <- matrix(rep(1:r / (r + 1), d), nrow = r, ncol = d)
    grid <- as.matrix(do.call(expand.grid, split(u, col(u))))

    ## Count the proportion of pseudo-observations in each box of volume s^d
    ## and compare to cutoff probability q
    keep.pts <- rep(NA, r^d)
    for (k in 1:r^d)
        keep.pts[k] <- mean(apply(U > grid[k,] - s & U <= grid[k,], 2, prod)) >= q
    grid.pts <- grid[keep.pts, ] # selected grid points

    ## Check grid.pts
    if (nrow(grid.pts) == 0)
        stop("No points selected. Consider changing the parameters 'r' or 'kappa'.")

    ## Transform to original margins
    pts <- matrix(NA, nrow(grid.pts), d)
    for (j in 1:d)
        pts[,j] = quantile(x = x[,j], probs = grid.pts[,j])

    ## Plot the selected points in dimension two
    if (plot && d == 2) {

        par(mfrow = c(1, 2))

        plot(t(U), main = paste("Min. # of pseudo-obs. per box:", ceiling(m * q)), xlab = "", ylab = "")
        points(grid, col = "blue", pch = 16)
        for (i in 1:r) {
            abline(v = i / (r + 1), lty = 2, col = "blue")
            abline(h = i / (r + 1), lty = 2, col = "blue")
        }
        points(grid.pts, col = "red", pch = 16)

        plot(x, main = "Original scale", xlab = "", ylab = "")
        points(pts, col = "red", pch = 16)

        par(mfrow = c(1, 1))
    }

    pts
}



##' @title Sequential open-end change-point detection tests based on the empirical
##'        dfs evaluted at p number of fixed points
##' @param x.learn learning sample (length m)
##' @param x sample (length nm)
##' @param pts evaluation points (rows). If not provided by user,
##'            chosen automatically from the learning sample using parameter r.
##' @param r number of evaluation points per dimension (if 'pts = NULL')
##'        to be chosen from the learning sample; r >= 2.
##' @param sigma covariance matrix to be used; if NULL, estimated by lrvar
##' @param kappa constant involved in the point selection procedure; default is 2
##' @return detector and additional arguments used in the computation of the detector

detOpenEndCpDist <- function(x.learn, x, pts = NULL, r = NULL, sigma = NULL, kappa = 1.5, ...) {

    ## Technical constant
    eta <- 0.001

    ## Checks
    if(!is.matrix(x.learn)) {
        warning("coercing 'x.learn' to a matrix.")
        stopifnot(is.matrix(x.learn <- as.matrix(x.learn)))
    }
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot(kappa > 1)

    ## Form all data
    x.all <- rbind(x.learn, x)

    ## Size of the learning sample
    m <- nrow(x.learn)

    ## Size of the sequence to be monitored
    nm <- nrow(x) # n - m
    ## Index of last observation
    n <- m + nm

    ## Dimension of the data: x.learn and x need to have the same number of columns
    stopifnot((d <- ncol(x.learn)) == ncol(x))

    ## Check r
    if (!is.null(r)) stopifnot(r >= 2L && r^d <= 50L)

    ## Check 'pts' if provided
    if(!is.null(pts) && !is.matrix(pts)) {
        warning("coercing 'pts' to a matrix.")
        stopifnot(is.matrix(pts <- as.matrix(pts)))
        stopifnot(d == ncol(pts))
    }

    ## If 'pts' and number of points per dimension are not provided, some default values for r
    if(is.null(pts) && is.null(r)) {
        if (d == 1) {
            warning("Setting 'r' to 5 in the univariate case.")
            r <- 5
        }
        else if (d == 2) {
            warning("Starting form a 4^2 grid in the bivariate case")
            r <- 4
        }
        else if (d == 3) {
            warning("Starting form a 3^3 grid in the trivariate case")
            r <- 3
        }
        else
            stop("The case 'd > 3' has not been dealt with yet")
    }

    ## If 'pts' is not provided by the user, get the points from the learning sample
    if(is.null(pts))
        pts <- if (d == 1)  ## if d = 1, use uniformly spaced quantiles from the learning sample
                   matrix(quantile(x = x.learn, probs = seq(0, 1, length.out = r+2)[-c(1, r+2)]))
               else ## if d > 1, select points surrounded by sufficiently many observations first
                   selectPoints(x.learn, r, kappa = kappa, plot = FALSE)

    ## Number of evaluation points
    p <- nrow(pts)

    ## Prepare data: apply indicator function, for each of the points
    ## in 'pts', and for each of the dimensions
    y <- matrix(1, n, p)
    for (j in 1:p)
        for (k in 1:d)
            y[,j] <- y[,j] * as.integer(x.all[,k] <= pts[j,k])

    ## The long-run covariance matrix estimated from the learning data set
    if (is.null(sigma))
        sigma <- try(sandwich::lrvar(y[1:m,], ...) * m)

    ## Inverse of the estimate
    sigma.inv <- try(solve(sigma))

    ## Warning if computation or inversion of the estimate of the long-run covariance matrix fails
    if (inherits(sigma, "try-error") || inherits(sigma.inv, "try-error")) {
        warning("Problem with the computation or the inversion of the estimate of the long-run covariance matrix. Consider passing alternative parameters to 'sandwich::lrvar()' through '...' or changing the parameters 'r' or 'kappa' passed to 'selectPoints()'. Alternatively, provide your own evalution points.")

        return(structure(class = "det.OpenEndCpDist",
                         list(det = NA, w = NA, m = m, n = n, p = p, pts = pts,
                              sigma = NA, eta = eta)))
    }

    ## Time grid
    tg <- seq.int(m+1, n) / m

    ## Compute detector
    det <- .C("seqOpenEndCpDistStat",
              as.double(y),
              as.integer(m),
              as.integer(n),
              as.integer(p),
              as.double(solve(sigma)), # inverse of the estimated long-run covariance matrix
              det = double(nm),
              w = integer(nm),
              PACKAGE = "npcp")

    structure(class = "det.OpenEndCpDist",
              list(det = det$det / tg^(1.5 + eta),
                   w = det$w + 1, # offset necessary: see formula
                   m = m, n = n, p = p, pts = pts, sigma = sigma, eta = eta))
}

##' @title Monitoring of sequential open-end change-point tests based on empirical dfs evluated
##'        at p fixed points
##' @param det output from detOpenEndCpDist()
##' @param alpha significance level (0.1, 0.05, 0.01)
##' @param plot if TRUE, plot detector and threshold
##' @return alarm, time of alarm, and additional arguments.

monOpenEndCpDist <- function(det, alpha = 0.05, plot = TRUE) {

    if(!inherits(det, "det.OpenEndCpDist"))
        warning("class of 'det' should be 'det.OpenEndCpDist'.")

    ## Check whether the detector contains NAs
    if (any(is.na(det$det))) {
        warning("Montinoring cannot be carried out because the detector contains NAs")
        return(list(alarm = NA, time.alarm = NA,
                    times.max = NA, time.change = NA,
                    alpha = alpha, sigma = det$sigma, eta = det$eta, detector = det$det, threshold = NA))
    }

    stopifnot(alpha %in% c(0.1, 0.05, 0.01))

    p <- det$p
    m <- det$m
    n <- det$n

    ## Load estimated quantile if p = 2, 5, 10, 20 or estimate quantile
    ## using the asymptotic regression model
    if (p %in% quantiles$d[[paste0((1 - alpha) * 100, "%")]]$p)
        thresh <- quantiles$d[[paste0((1 - alpha) * 100, "%")]]$quantiles[paste(p)]
    else {
        x <- log(p)
        c <- quantiles$d[[paste0((1 - alpha) * 100, "%")]]$c
        d <- quantiles$d[[paste0((1 - alpha) * 100, "%")]]$d
        e <- quantiles$d[[paste0((1 - alpha) * 100, "%")]]$e
        thresh <- 2 - (c + (d - c) * (1 - exp(-x / e)))
    }

    ## Look at exceedance and time of alarm if mean
    conds <- (det$det <= thresh)
    alarm <- !all(conds)
    ta <- if (alarm) which.max(1 - as.double(conds)) else NA # time of alarm
    tm <- det$w

    ## Plot the detector and the threshold
    if(plot) {
        plot(seq.int(m + 1, n), det$det, type = "l",
             xlab = "Monitoring period", ylab = "",
             ylim = c(0, 1.2 * max(det$det, thresh)))
        abline(thresh, 0, lty = 2)
        legend("topleft", c("threshold", "detector"), lty = c(2, 1))
    }

    list(alarm = alarm, time.alarm = if (alarm) ta + m else NA,
         times.max = tm, time.change = if (alarm) tm[ta] else NA,
         alpha = alpha, sigma = det$sigma, eta = det$eta, detector = det$det, threshold = thresh)
}


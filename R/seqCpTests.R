#################################################################################
##
##   R package npcp by Ivan Kojadinovic Copyright (C) 2020
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
## Sequential change-point tests based on the empirical dfs
#################################################################################

rBetaCopula <- function(x, n) {

    ## Checks
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot(n >= 1L)

    m <- nrow(x)
    d <- ncol(x)

    matrix(.C("rBetaCopula",
              as.integer(apply(x, 2, rank)),
              as.integer(m),
              as.integer(d),
              as.integer(n),
              x = double(n * d),
              PACKAGE = "npcp")$x, n, d)
}

rBetaCopulaRanks <- function(r, n) {

    m <- nrow(r)
    d <- ncol(r)

    .C("rBetaCopula",
       as.integer(r),
       as.integer(m),
       as.integer(d),
       as.integer(n),
       x = double(n * d),
       PACKAGE = "npcp")$x
}


#################################################################################
## Simulation / bootstrap
#################################################################################

simCpDist <- function(x.learn = NULL, m = NULL, n, gamma = 0.25, delta = 1e-4,
                      B = 1000, method = c("sim", "mult"), b = NULL,
                      weights = c("parzen", "bartlett"),
                      g = 5, L.method = c("max","median","mean","min")) {

    ## Checks
    method <- match.arg(method) # "beta", "mult.seq1", "mult.seq2", "mult.nonseq"
    weights <- match.arg(weights)
    L.method <- match.arg(L.method)

    ## Checks related to x.learn and m
    if (is.null(x.learn)) {
        if (is.null(m))
            stop("either 'x.learn' or 'm' needs to be specified")
        else {
            d <- 1
            if (method != "sim")
                stop("if 'x.learn' is not specified, only 'method = \"sim\"' can be used")
        }
    }
    else {
        if(!is.matrix(x.learn)) {
            warning("coercing 'x.learn' to a matrix.")
            stopifnot(is.matrix(x.learn <- as.matrix(x.learn)))
        }
        d <- ncol(x.learn)
        if (is.null(m))
            m <- nrow(x.learn)
        else
            stopifnot(nrow(x.learn) == m)
    }

    ## Other checks
    stopifnot(m > 1L && n > m)
    stopifnot(gamma >= 0 && gamma <= 0.5)
    stopifnot(delta >= 0 && delta <= 1)

    nm <- n - m

    ########## sim ###################################

    if (method == "sim") {

        if (d > 1L)
            stop("Setting 'method = \"sim\"' is possible only for univariate (independent) observations")

        do1 <- function() {
            stat <- .C("seqCpDistStat",
                       as.double(runif(n)),
                       as.integer(m),
                       as.integer(n),
                       as.integer(d),
                       mac = double(nm),
                       mmc = double(nm),
                       mmk = double(nm),
                       mc = double(nm),
                       mk = double(nm),
                       as.double(gamma),
                       as.double(delta),
                       wmc = integer(nm),
                       wmk = integer(nm),
                       PACKAGE = "npcp")
            c(stat$mac, stat$mmc, stat$mmk, stat$mc, stat$mk)
        }

        ## Replicates
        rep <- t(replicate(B, do1()))
        mac0 <- rep[,1:nm]
        mmc0 <- rep[,(nm+1):(2*nm)]
        mmk0 <- rep[,(2*nm+1):(3*nm)]
        mc0 <- rep[,(3*nm+1):(4*nm)]
        mk0 <- rep[,(4*nm+1):(5*nm)]
        pmax <- nm
    }

    ########## beta ###################################

    else if (method == "beta") {

        if (d <= 1L)
            stop("Setting 'method = \"beta\"' is possible only for multivariate (independent) observations")

        r <- apply(x.learn, 2, rank)

        do1 <- function() {
            stat <- .C("seqCpDistStat",
                       as.double(rBetaCopulaRanks(r, n)),
                       as.integer(m),
                       as.integer(n),
                       as.integer(d),
                       mac = double(nm),
                       mmc = double(nm),
                       mmk = double(nm),
                       mc = double(nm),
                       mk = double(nm),
                       as.double(gamma),
                       as.double(delta),
                       wmc = integer(nm),
                       wmk = integer(nm),
                       PACKAGE = "npcp")
            c(stat$mac, stat$mmc, stat$mmk, stat$mc, stat$mk)
        }

        ## Replicates
        rep <- t(replicate(B, do1()))
        mac0 <- rep[,1:nm]
        mmc0 <- rep[,(nm+1):(2*nm)]
        mmk0 <- rep[,(2*nm+1):(3*nm)]
        mc0 <- rep[,(3*nm+1):(4*nm)]
        mk0 <- rep[,(4*nm+1):(5*nm)]
        pmax <- nm
    }

    ########## mult ###################################

    else {

        mm <- floor(m * m / n)
        mmm <- m - mm

        ## Bandwith parameter
        if (is.null(b))
            b <- bOptEmpProc(x.learn, m = g, weights = weights,
                             L.method = L.method)
        stopifnot(b >= 1L)
        init.seq <- rnorm(B * (m + 2 * (b - 1)))

        ## CFuncName <- switch(method,
        ##                     "mult.seq1" = "seqCpDistMultSeq1",
        ##                     "mult.seq2" = "seqCpDistMultSeq2",
        ##                     "mult.nonseq" = "seqCpDistMultNonSeq")

        rep <- .C("seqCpDistMultNonSeq",
                  as.double(x.learn),
                  as.integer(m),
                  as.integer(n),
                  as.integer(d),
                  as.integer(B),
                  as.integer(1),
                  as.integer(b),
                  mac0 = double(B * mmm),
                  mmc0 = double(B * mmm),
                  mmk0 = double(B * mmm),
                  mc0 = double(B * mmm),
                  mk0 = double(B * mmm),
                  as.double(gamma),
                  as.double(delta),
                  as.double(init.seq),
                  PACKAGE = "npcp")

        ## Replicates
        mac0 <- matrix(rep$mac0, B, mmm, byrow = TRUE)
        mmc0 <- matrix(rep$mmc0, B, mmm, byrow = TRUE)
        mmk0 <- matrix(rep$mmk0, B, mmm, byrow = TRUE)
        mc0 <- matrix(rep$mc0, B, mmm, byrow = TRUE)
        mk0 <- matrix(rep$mk0, B, mmm, byrow = TRUE)
        pmax <- mmm

    }

    structure(class = "sims.cpDist",
              list(mac = mac0, mmc = mmc0, mmk = mmk0, mc = mc0, mk = mk0,
                   d = d, m = m, n = n, gamma = gamma, delta = delta,
                   B = B, method = method, pmax = pmax))
}

#################################################################################
## Threshold functions
#################################################################################

threshCpDist <- function(sims, p = 1, alpha = 0.05, type = 7) {

    ## Checks on 'sims'
    if (!inherits(sims, "sims.cpDist"))
        stop("'sims' should be obtained by 'simCpDist()'")

    mac0 <- sims$mac
    mmc0 <- sims$mmc
    mmk0 <- sims$mmk
    mc0 <- sims$mc
    mk0 <- sims$mk
    d <- sims$d
    m <- sims$m
    n <- sims$n
    nm <- n - m
    B <- sims$B
    method <- sims$method
    pmax <- sims$pmax

    stopifnot(p >= 1L)
    if (p > pmax) stop("The maximum possible value for 'p' is ", pmax)

    stopifnot(alpha > 0 && alpha <= 0.5)

    ## Threshold functions

    ## Blocks for computing block maxima
    bs <- pmax %/% p # ideal size of blocks; blocks of size bs+1 may exist
    s <- rep(bs, p) # sizes of blocks initialised at bs
    r <- pmax - p * bs # remaining number of lines
    if (r > 0) s[seq_len(r)] <- bs + 1 # size of first r blocks incremented if r > 0
    bl <- c(0, cumsum(s)) # 0 + ending line of each block

    ## Blocks for computing thresholds
    if (method %in% c("sim", "beta"))
        st <- s # sizes of blocks for computing thresholds
    else  { # multiplier replicates
        bs <- nm %/% p # ideal size of blocks; blocks of size bs+1 may exist
        st <- rep(bs, p) # sizes of blocks initialised at bs
        r <- nm - p * bs # remaining number of lines
        if (r > 0) st[seq_len(r)] <- bs + 1 # size of first r blocks incr. if r > 0
    }

    ## Quantile order
    q.prob <- (1 - alpha)^(1/p)

    ## Threshold functions
    computeThreshFunc <- function(rep) {

        ## Compute block maxima
        rep.max <- matrix(0,B,0)
        for (i in 1:p)
            rep.max <- cbind(rep.max,
                             apply(rep[,(bl[i] + 1):bl[i+1],drop=FALSE], 1, max))

        ## Compute thresholds
        threshold <- numeric(p)
        threshold[1] <- quantile(rep.max[,1], probs = q.prob, type = type)
        if (p > 1)
            for (i in 2:p) {
                rep.max <- rep.max[rep.max[,i-1] <= threshold[i-1],]
                threshold[i] <- quantile(rep.max[,i], probs = q.prob, type = type)
            }

        ## Compute threshold function so that it is of length nm
        rep(threshold, times = st) # threshold function
    }

    ## Compute threshold functions for the 5 detectors
    thresh.mac <- computeThreshFunc(mac0)
    thresh.mmc <- computeThreshFunc(mmc0)
    thresh.mmk <- computeThreshFunc(mmk0)
    thresh.mc <- computeThreshFunc(mc0)
    thresh.mk <- computeThreshFunc(mk0)

    structure(class = "thresh.cpDist",
              list(mac = computeThreshFunc(mac0),
                   mmc = computeThreshFunc(mmc0),
                   mmk = computeThreshFunc(mmk0),
                   mc = computeThreshFunc(mc0),
                   mk = computeThreshFunc(mk0),
                   d = d, m = m, n = n, gamma = sims$gamma, delta = sims$delta,
                   B = B, method = method, p = p, alpha = alpha, type = type))
}

#################################################################################
## Detector functions
#################################################################################

detCpDist <- function(x.learn, x, gamma = 0.25, delta = 1e-4) {

    ## Checks
    if(!is.matrix(x.learn)) {
        warning("coercing 'x.learn' to a matrix.")
        stopifnot(is.matrix(x.learn <- as.matrix(x.learn)))
    }
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }

    stopifnot(gamma >= 0 && gamma <= 0.5)
    stopifnot(delta >= 0 && delta <= 1)
    stopifnot(ncol(x) == (d <- ncol(x.learn)))

    m <- nrow(x.learn)
    nm <- nrow(x) # n - m
    n <- m + nm

    ## Detector function #############################

    det <- .C("seqCpDistStat",
              as.double(rbind(x.learn, x)),
              as.integer(m),
              as.integer(n),
              as.integer(d),
              mac = double(nm),
              mmc = double(nm),
              mmk = double(nm),
              mc = double(nm),
              mk = double(nm),
              as.double(gamma),
              as.double(delta),
              wmc = integer(nm),
              wmk = integer(nm),
              PACKAGE = "npcp")

    structure(class = "det.cpDist",
              list(mac = det$mac, mmc = det$mmc, mmk = det$mmk,
                   mc = det$mc, mk = det$mk,
                   wmc = det$wmc + 1, wmk = det$wmk + 1, # offset necessary: see formula
                   d = d, m = m, gamma = gamma, delta = delta))
}

#################################################################################
## Monitoring step
#################################################################################

monCpDist <- function(det, thresh,
                      statistic = c("mac", "mmc", "mmk", "mk", "mc"),
                      plot = TRUE) {

    ## Checks on 'det' and 'thresh'
    if (!inherits(det, "det.cpDist"))
        stop("'det' should be obtained by 'detCpDist()'")

    if (!inherits(thresh, "thresh.cpDist"))
        stop("'thresh' should be obtained by 'threshCpDist()'")

    if (det$d != thresh$d || det$m != thresh$m)
        stop("'det' and 'thresh' have not been computed from the same learning sample")

    if (det$gamma != thresh$gamma)
        stop("'det' and 'thresh' have not been computed with the same value of 'gamma'")

    if (det$delta != thresh$delta)
        stop("'det' and 'thresh' have not been computed with the same value of 'delta'")

    ## Other checks
    statistic <- match.arg(statistic)
    ds <- det[[statistic]] # detector values
    ts <- thresh[[statistic]] # threshold values
    if ((l <- length(ds)) > length(ts))
        stop("the number of detector values is greater than the number of threshold values")
    ## Check for too large gamma value
    if (statistic != "mac" && thresh$method == "mult" && thresh$gamma > 0.25)
        warning("the test might be too conservative with these settings; consider decreasing gamma")

    ## Alarm ?
    conds <- (ds <= ts[seq_len(l)])
    alarm <- !all(conds)
    ta <- if (alarm) which.max(1 - as.double(conds)) else NA # time of alarm
    tm <- if (statistic %in% c("mac", "mmc")) det$wmc
          else if (statistic == "mmk") det$wmk else NA

    if (plot) {
        mon.period <- (thresh$m+1):thresh$n # monitoring period
        plot(mon.period, ts, type = "l", lty = 1,
             xlab = "Monitoring period", ylab = "",
             ylim = c(0, 1.1 * max(ts, ds)))
        points(mon.period[1:l], ds, type = "b", lty = 3)
        legend("topleft", c("threshhold function", "detector function"), lty = c(1, 3))
    }
    list(alarm = alarm, time.alarm = if (alarm) ta + thresh$m else NA,
         times.max = tm, time.change = if (alarm) tm[ta] else NA)
}



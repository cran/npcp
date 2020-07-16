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
## Sequential change-point tests based on the differences of means
#################################################################################

#################################################################################
## Compute detectors
#################################################################################

detCpMean <- function(x.learn, x, sigma = NULL, b = NULL,
                      weights = c("parzen", "bartlett")) {
    ## Constants
    epsilon <- 1e-10
    gamma <- list(r = c(0, 0.1, 0.25),
                  s = c(0, 0.25, 0.45, 0.65, 0.85),
                  t = c(0, 0.25, 0.45),
                  e = c(0, 0.25, 0.45),
                  cs = c(0, 0.25, 0.45))

    eta <- c(0.1, 0.05, 0.01, 0.005, 0.001)#, 0.0005, 0.0001, 0)
    ne <- length(eta)

    ## Checks
    if(!is.double(x.learn)) {
        warning("coercing 'x.learn' to a double.")
        stopifnot(is.double(x.learn <- as.double(x.learn)))
    }
    if(!is.double(x)) {
        warning("coercing 'x' to a double.")
        stopifnot(is.double(x <- as.double(x)))
    }
    stopifnot(is.double(epsilon) && epsilon > 0)

    m <- length(x.learn)
    nm <- length(x) # n - m
    n <- m + nm

    ## Bandwith for long run variance estimation
    if (is.null(sigma)) {
        weights <- match.arg(weights)
        if (is.null(b))
            b <- bOpt(x.learn, weights=weights)
        else
            stopifnot((b <- as.integer(b)) >= 1L)

        ## Estimate long run variance from learning sample
        sigma <- sqrt(.C("LRVmean",
                         as.double(x.learn),
                         as.integer(m),
                         as.integer(weights == "bartlett"),
                         as.integer(b),
                         avar = double(1),
                         PACKAGE = "npcp")$avar)
    }
    else
        stopifnot(is.double(sigma) && sigma > 0)

    ## Scale information
    tg <- seq.int(m+1, n) / m # time grid

    ## Scale parameter
    sr <- 1.5
    ss <- 2.5
    st <- 2
    se <- 1
    scs <- 1

    ## Function for scaling trajectories of r, s, t
    scale.traj.rst <- function(x, sx, gamma) {
        out <- matrix(NA, length(x), 0)
        cnames <- c()
        for (i in 1:ne)
            for (j in 1:length(gamma)) {
                cnames <- c(cnames, paste0("eta.", eta[i], ".gamma.", gamma[j]))
                out <- cbind(out, x / (tg^(sx + eta[i]) * pmax(((tg - 1) / tg)^gamma[j], epsilon)))
            }
        colnames(out) <- cnames
        out
    }

    ## Function for scaling trajectories of e, cs
    scale.traj.ecs <- function(x, sx, gamma) {
        out <- matrix(NA, length(x), 0)
        cnames <- c()
        for (j in 1:length(gamma)) {
            cnames <- c(cnames, paste0("gamma.", gamma[j]))
            out <- cbind(out, x / (tg^sx * pmax(((tg - 1) / tg)^gamma[j], epsilon)))
        }
        colnames(out) <- cnames
        out
    }

    ## Compute detector
    det <- .C("seqCpMeanStat",
              as.double(c(x.learn, x)),
              as.integer(m),
              as.integer(n),
              r = double(nm),
              s = double(nm),
              t = double(nm),
              e = double(nm),
              cs = double(nm),
              wr = integer(nm),
              we = integer(nm),
              PACKAGE = "npcp")

    structure(class = "det.cpMean",
              list(r = scale.traj.rst(det$r, sr, gamma$r),
                   s = scale.traj.rst(det$s, ss, gamma$s),
                   t = scale.traj.rst(det$t, st, gamma$t),
                   e = scale.traj.ecs(det$e, se, gamma$e),
                   cs = scale.traj.ecs(det$cs, scs, gamma$cs),
                   wr = det$wr + 1, # offset necessary: see formula
                   we = det$we + 1, # offset necessary: see formula
                   m = m, n = n, b = b, weights = weights,
                   sigma = sigma, gamma = gamma, eta = eta))
}

#################################################################################
## Monitor
#################################################################################

monCpMean <- function(det, statistic = c("t", "s", "r", "e", "cs"), eta = 0.001,
                      gamma = 0.45, alpha = 0.05, sigma = NULL, plot = TRUE) {

    ## Checks on 'det'
    if (!inherits(det, "det.cpMean"))
        stop("'det' should be obtained by 'detCpMean()'")

    statistic <- match.arg(statistic)
    stopifnot(is.double(gamma) && gamma >= 0)
    stopifnot(alpha %in% c(0.1, 0.05, 0.01))

    m <- det$m
    n <- det$n

    ## Value of sigma
    if (is.null(sigma))
        sigma <- det$sigma
    else
        stopifnot(is.double(sigma) && sigma > 0)

    ## Is the chosen value of gamma possible for the chosen statisic ?
    if (min(abs(det$gamma[[statistic]] - gamma)) > 0) {
        gamma <- det$gamma[[statistic]][which.min(abs(det$gamma[[statistic]] - gamma))]
        warning("The chosen value of gamma cannot be used with statistic '", statistic,
                "'. Using gamma = ", gamma, " instead.")
    }

    ## Eta correctly chosen ? Set threshold and detector
    if (statistic %in% c("t", "s", "r")) {
        if (min(abs(det$eta - eta)) > 0) {
            eta <- det$eta[which.min(abs(det$eta - eta))]
            warning("The chosen value of eta cannot be used. Using eta = ", eta, " instead.")
        }
        ## Threshold for the chosen statistic
        ts <- quantiles[[statistic]][paste0((1-alpha)*100,"%"), paste0("gamma.", gamma), paste0("eta.", eta)]
        ## Detector
        ds <- det[[statistic]][,paste0("eta.", eta, ".gamma.", gamma)] / sigma

    } else {
        ## Threshold for the chosen statistic
        ts <- quantiles[[statistic]][paste0((1-alpha)*100,"%"), paste0("gamma.", gamma)] # quantiles for the chosen statistic
        ## Detector
        ds <- det[[statistic]][,paste0("gamma.", gamma)] / sigma
    }

    ## Alarm ?
    conds <- (ds <= ts)
    alarm <- !all(conds)
    ta <- if (alarm) which.max(1 - as.double(conds)) else NA # time of alarm
    tm <- if (statistic %in% c("r", "s", "t")) det$wr
          else if (statistic == "e") det$we else NA

    if (plot) {
        mon.period <- seq.int(m+1, n) # monitoring period
        plot(mon.period, rep(ts, times = n - m),  type = "l", lty = 1,
             xlab = "Monitoring period", ylab = "",
             ylim = c(0, 1.2 * max(ts, ds)))
        points(mon.period, ds, type = "b", lty = 3)
        legend("topleft", c("threshhold", "detector"), lty = c(1, 3))
    }

    list(alarm = alarm, time.alarm = if (alarm) ta + m else NA,
         times.max = tm, time.change = if (alarm) tm[ta] else NA,
         statistic = statistic, eta = eta, gamma = gamma, alpha = alpha,
         sigma = sigma, detector = ds, threshold = ts)
}



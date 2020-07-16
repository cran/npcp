\name{seqCpDist}
\alias{seqCpDist}
\alias{simCpDist}
\alias{threshCpDist}
\alias{detCpDist}
\alias{monCpDist}
\title{Closed-end Sequential Test for Change-Point Detection in Possibly
  Multivariate Time Series Sensitive to Changes in the Contemporary
  Distribution Function}
\description{
  Closed-end nonparametric sequential test for change-point detection
  based on the (multivariate) empirical distribution function. The
  observations can be continuous univariate or multivariate, and
  serially independent or dependent (strongly mixing). To carry out the
  test, four steps are required. The first step consists of simulating
  under the null many trajectories of the detector function. The second
  step consists of estimating a piecewise constant threshold function
  from these trajectories. The third step consists of computing the
  detector function from the data to be monitored. The fourth and last
  step consists of comparing the detector function with the estimated
  threshold function. Each of these steps corresponds to one of the
  functions in the usage section below. The current implementation is
  preliminary and not optimized for real-time monitoring (but could
  still be used for that). If the observations to be monitored are
  univariate and can be assumed serially independent, the simulation of
  the trajectories of the detector functions can be carried out using
  Monte Carlo simulation. In all other cases, the test relies on a
  \emph{dependent multiplier bootstrap}. Details can be found in the
  second reference.
}
\usage{
simCpDist(x.learn = NULL, m = NULL, n, gamma = 0.25, delta = 1e-4,
          B = 1000, method = c("sim", "mult"), b = NULL,
          weights = c("parzen", "bartlett"), g = 5,
          L.method = c("max","median","mean","min"))

threshCpDist(sims, p = 1, alpha = 0.05, type = 7)

detCpDist(x.learn, x, gamma = 0.25, delta = 1e-4)

monCpDist(det, thresh, statistic = c("mac", "mmc", "mmk", "mk", "mc"),
          plot = TRUE)
}
\arguments{
  \item{x.learn}{a data matrix whose rows are continuous observations,
    representing the learning sample.}
  \item{m}{a strictly positive integer specifying the size of the
    learning sample if \code{x.learn} is not specified; the latter
    implies that the observations are univariate and assumed to be
    independent; if \code{m} is not specified, it is taken equal to
    \code{nrow(x.learn)}.}
  \item{n}{a strictly positive integer specifying the monitoring horizon;
    the monitoring period is \code{m+1}, ..., \code{n}.}
  \item{gamma}{a real parameter between 0 and 0.5 appearing in the definition
    of the weight function used in the detector function.}
  \item{delta}{a real parameter between 0 and 1 appearing in the definition
    of the weight function used in the detector function.}
  \item{B}{the number of trajectories of the detector function to simulate
    under the null.}
  \item{method}{a string specifying the trajectory simulation method;
    can be either \code{"sim"} (Monte Carlo simulation -- only in
    the univariate case under the assumption of serial independence)
    or \code{"mult"} (the dependent multiplier bootstrap).}
  \item{b}{strictly positive integer specifying the value of the
    bandwidth parameter determining the serial dependence when
    generating dependent multiplier sequences using the 'moving average
    approach'; see Section 5 of the first reference. The value 1 will
    create i.i.d. multiplier sequences suitable for serially independent
    observations. If set to \code{NULL}, \code{b} will be estimated from
    \code{x.learn} using the function \code{\link{bOptEmpProc}()}; see
    the procedure described in Section 5 of the first reference.}
  \item{weights}{a string specifying the kernel for creating the
   weights used in the generation of dependent multiplier sequences within the
   'moving average approach'; see Section 5 of the first reference.}
 \item{g}{a strictly positive integer specifying the number of points of
   the uniform grid on \eqn{(0,1)^d}{(0,1)^d} (where \eqn{d}{d} is
   \code{ncol(x)}) involved in the estimation of the bandwidth
   parameter; see Section 5 of the first reference. The number of points
   of the grid is given by \code{g^ncol(x)} so that \code{g} needs to be
   decreased as \eqn{d}{d} increases.}
 \item{L.method}{a string specifying how the parameter \eqn{L}{L} involved
   in the estimation of the bandwidth parameter is computed; see
   Section 5 of the first reference.}
 \item{sims}{an object of class \code{sims.cpDist} containing simulated
   trajectories of the detector function under the null.}
 \item{p}{a strictly positive integer specifying the number of steps of
   the piece constant threshold function; \code{p} should not be
   taken too large (say, smaller than 4) if \code{method = "mult"}.}
 \item{alpha}{the value of the desired significance level for the
   sequential test.}
 \item{type}{an integer between 1 and 9 selecting one of the nine quantile
   algorithms detailed in the help of the function \code{\link{quantile}()}.}
 \item{x}{a data matrix whose rows are continuous observations
   corresponding to the new observations to be monitored for a change in
   contemporary distribution.}
 \item{det}{an object of class \code{det.cpDist} representing a detector
   function computed using \code{detCpDist()}.}
 \item{thresh}{an object of class \code{thresh.cpDist} representing a threshold
   function estimated using \code{threshCpDist()}.}
 \item{statistic}{a string specifying the statistic/detector to be used
   for the monitoring; can be either \code{"mac"}, \code{"mmc"},
   \code{"mmk"}, \code{"mc"} or \code{"mk"}; the last letter
   specifies whether it is a Cram\enc{é}{e}r-von Mises-like
   statistic (letter \code{"c"}) or a Kolmogorov-Smirnov-like
   statistic (letter \code{"k"}); the letters before specify the
   type of aggregation steps used to compute the detectors
   (\code{"m"} for maximum, \code{"a"} for average); \code{"mac"} corresponds to the
   detector \eqn{T_{m,q}}{T[m][','][q]} in the second reference,
   \code{"mmc"} to the detector \eqn{S_{m,q}}{S[m][','][q]},
   \code{"mmk"} to the detector \eqn{R_{m,q}}{R[m][','][q]},
   \code{"mc"} to the detector \eqn{Q_m}{Q[m]} and \code{"mk"} to
   the detector \eqn{P_m}{P[m]}.}
 \item{plot}{logical indicating whether the monitoring should be
   plotted.}
}
\details{
  The testing procedure is described in detail in the second reference.
}
\value{
  All functions return lists whose components have explicit names. The
  function \code{monCpDist()} in particular returns a list whose
  components are
  \item{alarm}{a logical indicating whether the detector function has
    exceeded the threshold function.}
  \item{time.alarm}{an integer corresponding to the time at
    which the detector function has exceeded the threshold function or
    \code{NA}.}
  \item{times.max}{a vector of times at which the successive detectors
    \code{"mmc"} (if \code{statistic = "mac"} or \code{statistic =
    "mmc"}) or \code{"mmk"} (if \code{statistic = "mmk"}) have reached
    their maximum; a vector of \code{NA}'s if \code{statistic = "mc"} or
    \code{statistic = "mk"}; this sequence of times can be used to
    estimate the time of change from the time of alarm.}
  \item{time.change}{an integer giving the estimated time of change if
    \code{alarm} is \code{TRUE}; the latter is simply the value in
    \code{times.max} which corresponds to \code{time.alarm}.}
}
\references{
  A. B\enc{ü}{u}cher and I. Kojadinovic (2016), A dependent multiplier
  bootstrap for the sequential empirical copula process under strong
  mixing, \emph{Bernoulli} \bold{22:2}, pages 927-968,
  \url{http://arxiv.org/abs/1306.3930}.

  I. Kojadinovic and G. Verdier (2020), Nonparametric sequential
  change-point detection for multivariate time series based on empirical
  distribution functions, 52 pages,
  \url{http://arxiv.org/abs/2004.12322}.
}
\note{
  This is a test for continuous (multivariate) time series.
}
\seealso{
  see \code{\link{cpDist}()} for the corresponding a posteriori (offline) test.
}
\examples{
\dontrun{
## Example of montoring for the period m+1, ..., n
m <- 100 # size of the learning sample
n <- 150 # monitoring horizon

## The learning sample
set.seed(123)
x.learn <- matrix(rnorm(m))

## New observations with a large change in mean
## to simulate monitoring for the period m+1, ..., n
k <- 125 ## the true change-point
x <- matrix(c(rnorm(k-m), rnorm(n-k, mean = 2)))

## Step 1: Simulation of B trajectories of the detector functions under the null
B <- 1e4

## Under the assumption of serial independence
## (no need to specify the learning sample)
traj.sim <- simCpDist(m = m, n = n, B = B, method = "sim")

## Without the assumption of serial independence
## (the learning sample is compulsory; the larger it is, the better;
## the monitoring horizon n should not be too large)
traj.mult <- simCpDist(x.learn = x.learn, n = n, B = B, method = "mult")

## Step 2: Compute threshold functions with p steps
p <- 2
tf.sim <- threshCpDist(traj.sim, p = p) # p can be taken large
                                        # if B is very large
tf.mult <- threshCpDist(traj.mult, p = p) # p should not be taken too
                                          # large unless both m and B
                                          # are very large

## Step 3: Compute the detectors for the monitoring period m+1, ... , n
det <- detCpDist(x.learn = x.learn, x = x)

## Step 4: Monitoring

## Simulate the monitoring with the first threshold function
monCpDist(det, tf.sim)

## Simulate the monitoring with the second threshold function
monCpDist(det, tf.mult)

## Simulate the monitoring with the first threshold function
## and another detector function
monCpDist(det, tf.sim, statistic = "mmk")

## Alternative steps 3 and 4:

## Compute the detectors for the monitoring period m+1, ... , m+20 only
det <- detCpDist(x.learn = x.learn, x = x[1:20,,drop = FALSE])

## Simulate the monitoring with the first threshold function
monCpDist(det, tf.sim)

## Simulate the monitoring with the second threshold function
monCpDist(det, tf.mult)
}
}
\keyword{htest}
\keyword{nonparametric}
\keyword{ts}
\keyword{multivariate}

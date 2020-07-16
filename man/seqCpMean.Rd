\name{seqCpMean}
\alias{seqCpMean}
\alias{detCpMean}
\alias{monCpMean}
\title{Open-end Nonparametric Sequential Change-Point Detection Test
  for Univariate Time Series Sensitive to Changes in the Mean}
\description{
  Open-end nonparametric sequential test for change-point detection
  based on the retrospective CUSUM statistic. The observations need to
  be univariate but can be serially dependent. To carry out the test,
  two steps are required. The first step consists of computing a
  detector function. The second step consists of comparing the detector
  function to a suitable constant threshold function. Each of these
  steps corresponds to one of the functions in the usage section
  below. The current implementation is preliminary and not optimized for
  real-time monitoring (but could still be used for that). Details can
  be found in the third reference.
}
\usage{
detCpMean(x.learn, x, sigma = NULL, b = NULL,
          weights = c("parzen", "bartlett"))

monCpMean(det, statistic = c("t", "s", "r", "e", "cs"), eta = 0.001,
          gamma = 0.45, alpha = 0.05, sigma = NULL, plot = TRUE)
}
\arguments{
 \item{x.learn}{a numeric vector representing the learning sample.}
 \item{x}{a numeric vector representing the observations collected after the
   beginning of the monitoring for a change in mean.}
 \item{sigma}{an estimate of the long-run variance of the time series of
   which \code{x.learn} is a stretch. If set to \code{NULL},
   \code{sigma} will be estimated using an approach similar to
   those described in the fourth reference.}
 \item{b}{strictly positive integer specifying the value of the
   bandwidth for the estimation of the long-run variance if
   \code{sigma} is not provided. If set to \code{NULL}, \code{b} will
   be estimated from \code{x.learn} using the function
   \code{\link{bOpt}()}.}
 \item{weights}{a string specifying the kernel for creating the weights
   used for the estimation of the long-run variance if \code{sigma} is
   not provided; see Section 5 of the first reference.}
 \item{det}{an object of class \code{det.cpMean} representing a detector
   function computed using \code{detCpMean()}.}
 \item{statistic}{a string specifying the statistic/detector to be used
   for the monitoring; can be either \code{"t"}, \code{"s"},
   \code{"r"}, \code{"e"} or \code{"cs"}; \code{"t"} corresponds to the
   detector \eqn{T_{m}}{T[m]} in the third reference, \code{"s"} to the
   detector \eqn{S_{m}}{S[m]}, \code{"r"} to the detector \eqn{R_{m}}{R[m]},
   \code{"e"} to the detector \eqn{E_m}{E[m]} and \code{"cs"} to
   so-called ordinary CUSUM detector denoted by \eqn{Q_m}{Q[m]} in the
   third reference. Note that the detector \eqn{E_m}{E[m]}
   was proposed in the second reference.}
 \item{eta}{a real parameter whose role is described in detail in
   the third reference.}
 \item{gamma}{a real parameter that can improve the power of the
   sequential test at the beginning of the monitoring; possible values
   are 0, 0.1, 0.25, 0.45, 0.65 and 0.85, but not for all statistics; see
   the third reference.}
\item{alpha}{the value of the desired significance level for the
   sequential test.}
\item{plot}{logical indicating whether the monitoring should be
   plotted.}
}
\details{
  The testing procedure is described in detail in the third
  reference. An alternative way of estimating the long-run variance is
  to use the function \code{\link[sandwich]{lrvar}()} of the package
  \pkg{sandwich} and to pass it through the argument \code{sigma}.
}
\value{
  Both functions return lists whose components have explicit names. The
  function \code{monCpMean()} in particular returns a list whose
  components are
  \item{alarm}{a logical indicating whether the detector function has
    exceeded the threshold function.}
  \item{time.alarm}{an integer corresponding to the time at
    which the detector function has exceeded the threshold function or
    \code{NA}.}
  \item{times.max}{a vector of times at which the successive detectors
    \code{"r"} (if \code{statistic = "r"}, \code{statistic = "s"} or
    \code{statistic = "t"}) or \code{"e"} (if \code{statistic = "e"})
    have reached their maximum; a vector of \code{NA}'s if
    \code{statistic = "cs"}; this sequence of times can be used to
    estimate the time of change from the time of alarm.}
  \item{time.change}{an integer giving the estimated time of change if
    \code{alarm} is \code{TRUE}; the latter is simply the value in
    \code{times.max} which corresponds to \code{time.alarm}.}
  \item{statistic}{the value of \code{statistic} in the call of the function.}
  \item{eta}{the value of \code{eta} in the call of the function.}
  \item{gamma}{the value of \code{gamma} in the call of the function.}
  \item{alpha}{the value of \code{alpha} in the call of the function.}
  \item{sigma}{the value of \code{sigma} in the call of the function.}
  \item{detector}{the successive values of the chosen detector.}
  \item{threshold}{the value of the constant threshold for the chosen detector.}
}
\references{
  A. B\enc{ü}{u}cher and I. Kojadinovic (2016), A dependent multiplier
  bootstrap for the sequential empirical copula process under strong
  mixing, \emph{Bernoulli} \bold{22:2}, pages 927-968,
  \url{http://arxiv.org/abs/1306.3930}.

  J. G\enc{ö}{o}smann, T. Kley and H. Dette (2020), A new approach for
  open-end sequential change point monitoring, 49 pages,
  \url{http://arxiv.org/abs/1906.03225}.

  M. Holmes and I. Kojadinovic (2020), Open-end nonparametric sequential
  change-point detection based on the retrospective CUSUM statistic, 41
  pages, available on the arXiv.

  D.N. Politis and H. White (2004), Automatic block-length selection for the
  dependent bootstrap, \emph{Econometric Reviews} \bold{23(1)}, pages 53-70.
}
%\note{}
\seealso{
  see \code{\link{cpMean}()} for the corresponding a posteriori (offline) test.
}
\examples{
\dontrun{
## Example of open-end montoring
m <- 100 # size of the learning sample

## The learning sample
set.seed(123)
x.learn <- rnorm(m)

## New observations with a change in mean
## to simulate monitoring for the period m+1, ..., n
n <- 5000
k <- 2500 ## the true change-point
x <- c(rnorm(k-m), rnorm(n-k, mean = 0.2))

## Step 1: Compute the detector
det <- detCpMean(x.learn = x.learn, x = x)

## Step 2: Monitoring with the default detector
m1 <- monCpMean(det)
str(m1)

## Monitoring with another detector
m2 <- monCpMean(det, statistic = "s", gamma = 0.85)
str(m2)
}
}
\keyword{htest}
\keyword{nonparametric}
\keyword{ts}
\keyword{multivariate}


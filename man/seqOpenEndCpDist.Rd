\name{seqOpenEndCpDist}
\alias{seqOpenEndCpDist}
\alias{detOpenEndCpDist}
\alias{monOpenEndCpDist}

\title{
Open-end Nonparametric Sequential Change-Point Detection Test
  for (Possibly) Multivariate Time Series Sensitive to
  Changes in the Distribution Function
}
\description{
  Open-end nonparametric sequential test for change-point detection based
  on a retrospective CUSUM statistic constructed from differences of
  empirical distribution functions. The observations can be univariate
  or multivariate (low-dimensional), and serially dependent.  To carry
  out the test, two steps are required. The first step consists of
  computing a detector function. The second step consists of comparing
  the detector function to a suitable constant threshold function. Each
  of these steps corresponds to one of the functions in the usage
  section below. The current implementation is preliminary and not
  optimized for real-time monitoring (but could still be used for
  that). Details can be found in the first reference.
}
\usage{
detOpenEndCpDist(x.learn, x, pts = NULL, r = NULL, sigma = NULL, kappa = 1.5, ...)

monOpenEndCpDist(det, alpha = 0.05, plot = TRUE)
}

\arguments{
  \item{x.learn}{a numeric matrix representing the learning sample.}
  \item{x}{a numeric matrix representing the observations collected
    after the beginning of the monitoring.}  
  \item{pts}{a numeric matrix whose rows represent the evaluation
             points; if not provided by user, chosen automatically from
             the learning sample using parameter \code{r}.}	   
  \item{r}{integer greater or equal than 2 representing the number of
    evaluation points per dimension to be chosen from the learning
    sample; used only if \code{pts = NULL}.}
  \item{sigma}{a numeric matrix representing the covariance matrix to be
    used; if \code{NULL}, estimated by \code{sandwich::lrvar()}.}  
  \item{kappa}{constant involved in the point selection procedure; used
  only if the multivariate case; should be larger than 1.}
\item{...}{optional arguments passed to \code{sandwich::lrvar()}.}
  \item{det}{an object of class \code{det.OpenEndCpDist} representing a
   detector function computed using \code{detOpenEndCpDist()}.}
 \item{alpha}{the value of the desired significance level for the
   sequential test.}
  \item{plot}{logical indicating whether the monitoring should be
   plotted.}
}
\details{
  The testing procedure is described in detail in the first reference.
}
\value{
  Both functions return lists whose components have explicit names. The
  function \code{monOpenEndCpDist()} in particular returns a list whose
  components are
  \item{alarm}{a logical indicating whether the detector function has
    exceeded the threshold function.}
  \item{time.alarm}{an integer corresponding to the time at
    which the detector function has exceeded the threshold function or
    \code{NA}.}
  \item{times.max}{a vector of times at which the successive detectors
    have reached their maximum; this sequence of times can be used to
    estimate the time of change from the time of alarm.}
  \item{time.change}{an integer giving the estimated time of change if
    \code{alarm} is \code{TRUE}; the latter is simply the value in
    \code{times.max} which corresponds to \code{time.alarm}.}
  \item{statistic}{the value of \code{statistic} in the call of the function.}
  \item{eta}{the value of \code{eta} in the call of the function.}
  \item{p}{number of evaluations points of the empirical distribution functions.}
  \item{pts}{evaluation points of the empirical distribution functions.}
  \item{alpha}{the value of \code{alpha} in the call of the function.}
  \item{sigma}{the value of \code{sigma} in the call of the function.}
  \item{detector}{the successive values of the detector.}
  \item{threshold}{the value of the constant threshold for the detector.}
}
\references{
  M. Holmes, I. Kojadinovic and A. Verhoijsen (2022), Multi-purpose
  open-end monitoring procedures for multivariate observations based on
  the empirical distribution function, 45 pages,
  \url{https://arxiv.org/abs/2201.10311}.

  M. Holmes and I. Kojadinovic (2021), Open-end nonparametric sequential
  change-point detection based on the retrospective CUSUM statistic,
  \emph{Electronic Journal of Statistics} \bold{15:1}, pages 2288-2335,
  \doi{10.1214/21-EJS1840}.
}

\seealso{
  See \code{\link{detOpenEndCpMean}()} for the corresponding test
  sensitive to changes in the mean, \code{\link{selectPoints}()} for the
  underlying point selection procedure used in the multivariate case
  and \code{\link[sandwich]{lrvar}()} for information on the estimation
  of the underlying long-run covariance matrix.
}
\examples{
\dontrun{
## Example of open-end monitoring
m <- 800 # size of the learning sample
nm <- 5000 # number of collected observations after the start
n <- nm + m # total number of observations

set.seed(456)

## Univariate, no change in distribution
r <- 5 # number of evaluation points
x <- rnorm(n)
## Step 1: Compute the detector
det <- detOpenEndCpDist(x.learn = matrix(x[1:m]),
                        x = matrix(x[(m + 1):n]), r = r)
## Step 2: Monitoring
mon <- monOpenEndCpDist(det = det, alpha = 0.05, plot = TRUE)

## Univariate, change in distribution
k <- 2000 # m + k + 1 is the time of change
x[(m + k + 1):n] <- rt(nm - k, df = 3)
det <- detOpenEndCpDist(x.learn = matrix(x[1:m]),
                        x = matrix(x[(m + 1):n]), r = r)
mon <- monOpenEndCpDist(det = det, alpha = 0.05, plot = TRUE)

## Bivariate, no change
d <- 2
r <- 4 # number of evaluation points per dimension
x <- matrix(rnorm(n * d), nrow = n, ncol = d)
det <- detOpenEndCpDist(x.learn = x[1:m, ], x = x[(m + 1):n, ], r = r)
mon <- monOpenEndCpDist(det = det, alpha = 0.05, plot = TRUE)

## Bivariate, change in the mean of the first margin
x[(m + k + 1):n, 1] <- x[(m + k + 1):n, 1] + 0.3
det <- detOpenEndCpDist(x.learn = x[1:m, ], x = x[(m + 1):n, ], r = r)
mon <- monOpenEndCpDist(det = det, alpha = 0.05, plot = TRUE)

## Bivariate, change in the dependence structure
x1 <- rnorm(n)
x2 <- c(rnorm(m + k, 0.2 * x1[1:(m + k)], sqrt((1 - 0.2^2))),
        rnorm(nm - k, 0.7 * x1[(m + k + 1):n], sqrt((1 - 0.7^2))))
x <- cbind(x1, x2)
det <- detOpenEndCpDist(x.learn = x[1:m, ], x = x[(m + 1):n, ], r = r)
mon <- monOpenEndCpDist(det = det, alpha = 0.05, plot = TRUE)
}
}
\keyword{htest}
\keyword{nonparametric}
\keyword{ts}
\keyword{univariate}
\keyword{multivariate}
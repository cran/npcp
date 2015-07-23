\name{cpTestBM}
\alias{cpTestBM}
\title{Nonparametric tests for change-point detection in the
  distribution of independent block maxima}
\description{
  Nonparametric tests for change-point detection in the distribution of
  independent block maxima based either on the probability weighted moment
  method (see the second reference) or on the generalized probability
  weighted moment method (see the first reference) for
  estimating the parameters of the generalized extreme value (GEV)
  distribution. It is assumed that the block maxima are independent and
  that their unknown c.d.f.s are continuous, but not necessarily that they are GEV
  distributed. Three statistics are computed. Under the assumption that
  the block maxima are GEV distributed, these are statistics
  particularly sensitive to changes in the location, scale and shape
  parameters of the GEV. Details can be found in third reference.
}
\usage{
cpTestBM(x, method = c("pwm", "gpwm"), r=10)
}
\arguments{
  \item{x}{a numeric vector representing independent block maxima whose
    unknown c.d.f.s are assumed continuous.}
  \item{method}{a string specifying how statistics will be defined; can
    be either \code{"pwm"} (the probability weighted moment method) or
    \code{"gpwm"} (the generalized probability weighted moment method).
    The method \code{"pwm"} is suggested for climate block maxima
    that are typically not too heavy tailed, more precisely, whose
    distributions are in
    the maximum domains of attraction of GEV distributions with shape parameters
    smaller than a half. The method \code{"gpwm"} should be prefered
    otherwise.}
  \item{r}{strictly positive integer specifying the set of breakpoints
    that will be tested; more precisely, starting from the initial
    sample of block maxima, the tests compare subsamples formed by the
    \code{k} first maxima and \code{n-k} last maxima
    for \code{k} in the set \code{{r,...,n-r}}, where \code{n}
    is the sample size.}
}
\details{
  Approximate p-values are computed
  from the estimated asymptotic null distributions, which involve the
  Kolmogorov distribution. The latter is dealt with reusing code from
  the \code{\link{ks.test}()} function; credit to RCore.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{value of the three test statistics.}
  \item{pvalues}{corresponding approximate p-values.}
  \item{stats.loc}{the values of the \code{n - (2 * r - 1)} intermediate
    change-point statistics sensitive to changes in the location;
    the first test statistic is defined as the maximum of those.}
  \item{stats.scale}{the values of the \code{n - (2 * r - 1)} intermediate
    change-point statistics sensitive to changes in the scale;
    the second test statistic is defined as the maximum of those.}
  \item{stats.shape}{the values of the \code{n - (2 * r - 1)} intermediate
    change-point statistics sensitive to changes in the shape;
    the third test statistic is defined as the maximum of those.}
}
\references{
  J. Diebolt, A. Guillou, P. Naveau and P. Ribereau (2008), Improving
  probability-weighted moment methods for the generalized extreme-value
  distribution, \emph{REVSTAT} 6, pages 33-50.

  J.R.M. Hosking, J.R. Wallis and E.F. Wood (1985), Estimation of the
  generalized extreme-value distribution by the method of
  probability-weighted moments, \emph{Technometrics} 27, pages 251-261.

  I. Kojadinovic and P. Naveau (2015), Nonparametric tests for
  change-point detection in the distribution of block maxima based on
  probability weighted moments, \url{http://arxiv.org/abs/1507.06121}.
}
\note{
  The tests were derived under the assumption of block maxima whose
  unknown c.d.f.s are continuous, which implies that ties occur with
  probability zero. A way to deal with ties based on randomization
  is proposed in the third reference.
}
\seealso{
  \code{\link{cpTestFn}()} for a related test based on the empirical c.d.f.
}
\examples{
\dontrun{
require(evd)
n <- 100
k <- 50 ## the true change-point

## change in the shape parameter of a GEV
x <- rgev(k,loc=0,scale=1,shape=-0.8)
y <- rgev(k,loc=0,scale=1,shape=0.4)
cp <- cpTestBM(c(x,y))
cp
## estimated change-point
which(cp$stats.shape == max(cp$stats.shape))

## change in the scale parameter of a GEV
x <- rgev(k,loc=0,scale=0.5,shape=0)
y <- rgev(k,loc=0,scale=1,shape=0)
cp <- cpTestBM(c(x,y))
cp
## estimated change-point
which(cp$stats.scale == max(cp$stats.scale))

## change in the location parameter of a GEV
x <- rgev(k,loc=0,scale=1,shape=0)
y <- rgev(k,loc=0.5,scale=1,shape=0)
cp <- cpTestBM(c(x,y))
cp
## estimated change-point
which(cp$stats.loc == max(cp$stats.loc))
}
}
\keyword{htest}
\keyword{nonparametric}
\keyword{ts}

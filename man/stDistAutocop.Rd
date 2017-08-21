\name{stDistAutocop}
\alias{stDistAutocop}
\title{Combined Test of Stationarity for Univariate Continuous Time Series Sensitive
  to Changes in the Distribution Function and the Autocopula}
\description{
  A nonparametric test of stationarity for univariate continuous time
  series resulting from a combination \enc{à}{a} la Fisher of the
  change-point test sensitive to changes in the distribution function
  implemented in \code{\link{cpDist}()} and the change-point test
  sensitive to changes in the autcopula implemented in
  \code{\link{cpAutocop}()}. Approximate p-values are obtained by
  combining two \emph{multiplier} resampling schemes. Details can be
  found in the first reference.
}
\usage{
stDistAutocop(x, lag = 1, b = NULL, pairwise = FALSE,
              weights = c("parzen", "bartlett"), m = 5, N = 1000)
}
\arguments{
  \item{x}{a one-column matrix containing continuous observations.}
  \item{lag}{an integer specifying at which lag to consider the
    autocopula; the autcopula is a (\code{lag+1})-dimensional copula.}
  \item{b}{strictly positive integer specifying the value of the
    bandwidth parameter determining the serial dependence when
    generating dependent multiplier sequences using the 'moving average
    approach'; see Section 5 of the second reference. If set to
    \code{NULL}, \code{b} will be estimated using the
    function \code{\link{bOptEmpProc}()}; see the first reference.}
  \item{pairwise}{a logical specifying whether the test should focus
    only on the bivariate margins of the (\code{lag+1})-dimensional
    autocopula.}
  \item{weights}{a string specifying the kernel for creating the
   weights used in the generation of dependent multiplier sequences within the
   'moving average approach'; see Section 5 of the second reference.}
  \item{m}{a strictly positive integer specifying the number of points of the
   uniform grid on \eqn{(0,1)}{(0,1)} involved in the estimation of the bandwidth
   parameter; see Section 5 of the third reference. }
 \item{N}{number of multiplier replications.}
}
\details{
  The testing procedure is described in detail in the second section of
  the first reference.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{value of the test statistic.}
  \item{p.value}{corresponding approximate p-value \enc{à}{a} Fisher.}
  \item{component.p.values}{p-values of the component tests arising the combination.}
  \item{b}{the value of parameter \code{b}.}
}
\references{
  A. B\enc{ü}{u}cher, J.-D. Fermanian and I. Kojadinovic (2017),
  Combining cumulative sum change-point detection tests for assessing
  the stationarity of continuous univariate time series,
  \url{http://arxiv.org/}.

  A. B\enc{ü}{u}cher and I. Kojadinovic (2016), A dependent multiplier
  bootstrap for the sequential empirical copula process under strong
  mixing, \emph{Bernoulli} \bold{22:2}, pages 927-968,
  \url{http://arxiv.org/abs/1306.3930}.
}
\note{
  This is a tests for a continuous multivariate time series.
}
\seealso{
  see \code{\link{cpDist}()} and \code{\link{cpAutocop}()} for the
  component tests.
}
\examples{
## AR1 example
n <- 200
k <- n/2 ## the true change-point
x <- matrix(c(arima.sim(list(ar = -0.1), n = k),
              arima.sim(list(ar = 0.5), n = n - k)))
stDistAutocop(x)

## AR2 example
n <- 200
k <- n/2 ## the true change-point
x <- matrix(c(arima.sim(list(ar = c(0,-0.1)), n = k),
              arima.sim(list(ar = c(0,0.5)), n = n - k)))
\dontrun{
stDistAutocop(x)
stDistAutocop(x, lag = 2)}
stDistAutocop(x, lag = 2, pairwise = TRUE)
}

\keyword{htest}
\keyword{nonparametric}
\keyword{ts}

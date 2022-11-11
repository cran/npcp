\name{cpAutocop}
\alias{cpAutocop}
\title{Test for Change-Point Detection in Univariate Observations
  Sensitive to Changes in the Autocopula}
\description{
  Nonparametric test for change-point detection particularly sensitive
  to changes in the autocopula of univariate continuous observations.
  Approximate p-values for the test statistic are obtained by means of a
  \emph{multiplier} approach. Details can be found in the first reference.
}
\usage{
cpAutocop(x, lag = 1, b = NULL, bivariate = FALSE,
          weights = c("parzen", "bartlett"), m = 5,
          N = 1000, init.seq = NULL, include.replicates = FALSE)
}
\arguments{
  \item{x}{a one-column matrix containing continuous observations.}
  \item{lag}{an integer specifying at which lag to consider the
    autocopula; the autocopula is a (\code{lag+1})-dimensional copula.}
  \item{b}{strictly positive integer specifying the value of the
    bandwidth parameter determining the serial dependence when
    generating dependent multiplier sequences using the 'moving average
    approach'; see Section 5 of the second reference. If set to
    \code{NULL}, \code{b} will be estimated using the
    function \code{\link{bOptEmpProc}()}; see the first reference.}
  \item{bivariate}{a logical specifying whether the test should focus
    only on the bivariate margin of the (\code{lag+1})-dimensional autocopula
    obtained from the first and the last dimension.}
  \item{weights}{a string specifying the kernel for creating the
   weights used in the generation of dependent multiplier sequences within the
   'moving average approach'; see Section 5 of the second reference.}
  \item{m}{a strictly positive integer specifying the number of points of the
   uniform grid on \eqn{(0,1)}{(0,1)} involved in the estimation of the bandwidth
   parameter; see Section 5 of the second reference. }
 \item{N}{number of multiplier replications.}
 \item{init.seq}{a sequence of independent standard normal variates of
   length \code{N * (nrow(x) - lag  + 2 * (b - 1))} used to generate dependent
   multiplier sequences.}
 \item{include.replicates}{a logical specifying whether the
   object of \code{\link{class}} \code{htest} returned by the function
   (see below) will include the multiplier replicates.}
}
\details{
  The approximate p-value is computed as
  \deqn{(0.5 +\sum_{i=1}^N\mathbf{1}_{\{S_i\ge S\}})/(N+1),}{%
        (0.5 + sum(S[i] >= S, i=1, .., N)) / (N+1),}
  where \eqn{S} and \eqn{S_i}{S[i]} denote the test statistic and
  a multiplier replication, respectively. This ensures that the
  approximate p-value is a number strictly between 0 and 1, which is
  sometimes necessary for further treatments.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{value of the test statistic.}
  \item{p.value}{corresponding approximate p-value.}
  \item{cvm}{the values of the \code{length(x)-lag-1} intermediate
    Cram\enc{é}{e}r-von Mises change-point statistics; the test
    statistic is defined as the maximum of those.}
  \item{b}{the value of parameter \code{b}.}
}
\references{
   A. B\enc{ü}{u}cher, J.-D. Fermanian and I. Kojadinovic (2019),
  Combining cumulative sum change-point detection tests for assessing
  the stationarity of univariate time series, \emph{Journal of Time
  Series Analysis} \bold{40}, pages 124-150,
  \url{https://arxiv.org/abs/1709.02673}.

  A. B\enc{ü}{u}cher and I. Kojadinovic (2016), A dependent multiplier
  bootstrap for the sequential empirical copula process under strong
  mixing, \emph{Bernoulli} \bold{22:2}, pages 927-968,
  \url{https://arxiv.org/abs/1306.3930}.
}
\note{
  This is a tests for a continuous univariate time series.
}
\seealso{
  \code{\link{cpAutocov}()} for a related test based on
  the autocovariance.
}
\examples{
## AR1 example
n <- 200
k <- n/2 ## the true change-point
x <- matrix(c(arima.sim(list(ar = -0.5), n = k),
              arima.sim(list(ar = 0.5), n = n - k)))
cp <- cpAutocop(x)
cp
## Estimated change-point
which(cp$cvm == max(cp$cvm))

## AR2 example
n <- 200
k <- n/2 ## the true change-point
x <- matrix(c(arima.sim(list(ar = c(0,-0.5)), n = k),
              arima.sim(list(ar = c(0,0.5)), n = n - k)))
cpAutocop(x)
cpAutocop(x, lag = 2)
cpAutocop(x, lag = 2, bivariate = TRUE)
}
\keyword{htest}
\keyword{univariate}
\keyword{nonparametric}
\keyword{ts}

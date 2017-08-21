\name{cpU}
\alias{cpMean}
\alias{cpVar}
\alias{cpGini}
\alias{cpAutocov}
\alias{cpCov}
\alias{cpTau}
\title{Some CUSUM Tests for Change-Point Detection Based on U-statistics}
\description{
  Nonparametric CUSUM tests for change-point detection particularly sensitive
  to changes in certain quantities that can be estimated using
  one-sample U-statistics of order one or two. So far, the quantities under
  consideration are the expectation (thus corresponding to the standard
  CUSUM test based on the sample mean), the variance, Gini's mean difference, the
  autocovariance at a specified lag, the covariance for bivariate data
  and Kendall's tau for multivariate data.
  The observations can be serially independent or dependent (strongly mixing).
  Approximate p-values for the test statistic are obtained by means of a
  \emph{multiplier} approach or by estimating the asymptotic null
  distribution. Details can be found in the first reference.
}
\usage{
cpMean(x, method = c("nonseq", "seq", "asym.var"),
       b = NULL, weights = c("parzen", "bartlett"),
       N = 1000, init.seq = NULL, include.replicates = FALSE)

cpVar(x, method = c("nonseq", "seq", "asym.var"),
      b = NULL, weights = c("parzen", "bartlett"),
      N = 1000, init.seq = NULL, include.replicates = FALSE)

cpGini(x, method = c("nonseq", "seq", "asym.var"),
       b = NULL, weights = c("parzen", "bartlett"),
       N = 1000, init.seq = NULL, include.replicates = FALSE)

cpAutocov(x, lag = 1, method = c("nonseq", "seq", "asym.var"),
          b = NULL, weights = c("parzen", "bartlett"),
          N = 1000, init.seq = NULL, include.replicates = FALSE)

cpCov(x, method = c("nonseq", "seq", "asym.var"),
      b = NULL, weights = c("parzen", "bartlett"),
      N = 1000, init.seq = NULL, include.replicates = FALSE)

cpTau(x, method = c("seq", "nonseq", "asym.var"),
      b = NULL, weights = c("parzen", "bartlett"),
      N = 1000, init.seq = NULL, include.replicates = FALSE)

}
\arguments{
  \item{x}{a numeric vector or a data matrix containing continuous observations.}
  \item{lag}{an integer specifying at which lag to consider the autocovariance.}
  \item{method}{a string specifying the method for computing the
    approximate p-value for the test statistic; can be either
    \code{"seq"} (the 'check' approach in the first reference),
    \code{"nonseq"} (the 'hat' approach in the first reference), or
    \code{"asym.var"} (the approach based on the estimation of the
    asymptotic null distribution of the test statistic described
    in the first reference). The 'seq' approach appears overall to lead to
    better behaved tests for \code{\link{cpTau}()}. More experiments
    are necessary for the other functions.}
  \item{b}{strictly positive integer specifying the value of the
    bandwidth parameter determining the serial dependence when
    generating dependent multiplier sequences using the 'moving average
    approach'; see Section 5 of the second reference. The value 1
    will create i.i.d. multiplier
    sequences suitable for serially independent observations. If set to
    \code{NULL}, \code{b} will be estimated from \code{x} using the
    procedure described in the first reference.}
 \item{weights}{a string specifying the kernel for creating the
   weights used in the generation of dependent multiplier sequences within the
   'moving average approach'; see Section 5 of the second reference.}
 \item{N}{number of multiplier replications.}
 \item{init.seq}{a sequence of independent standard normal variates of
   length \code{N * (nrow(x) + 2 * (b - 1))} used to generate dependent
   multiplier sequences.}
 \item{include.replicates}{a logical specifying whether the
   object of \code{\link{class}} \code{htest} returned by the function
   (see below) will include the multiplier replicates, if generated.}

}
\details{
  When \code{method} is either \code{"seq"} or \code{"nonseq"},
  the approximate p-value is computed as
  \deqn{(0.5 +\sum_{i=1}^N\mathbf{1}_{\{S_i\ge S\}})/(N+1),}{%
        (0.5 + sum(S[i] >= S, i=1, .., N)) / (N+1),}
  where \eqn{S} and \eqn{S_i}{S[i]} denote the test statistic and
  a multiplier replication, respectively. This ensures that the
  approximate p-value is a number strictly between 0 and 1, which is
  sometimes necessary for further treatments.

  When \code{method = "asym.var"}, the approximate p-value is computed
  from the estimated asymptotic null distribution, which involves the
  Kolmogorov distribution. The latter is dealt with reusing code from
  the \code{\link{ks.test}()} function; credit to RCore.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{value of the test statistic.}
  \item{p.value}{corresponding approximate p-value.}
  \item{u}{the values of the \code{nrow(x)-3} intermediate
    change-point statistics; the test statistic is defined as
    the maximum of those.}
  \item{b}{the value of parameter \code{b}.}
}
\references{
  A. B\enc{ü}{u}cher and I. Kojadinovic (2016), Dependent multiplier
  bootstraps for non-degenerate U-statistics under mixing conditions
  with applications, \emph{Journal of Statistical Planning and
    Inference} \bold{170}, pages 83-105, \url{http://arxiv.org/abs/1412.5875}.

  A. B\enc{ü}{u}cher and I. Kojadinovic (2016), A dependent multiplier
  bootstrap for the sequential empirical copula process under strong
  mixing, \emph{Bernoulli} \bold{22:2}, pages 927-968,
  \url{http://arxiv.org/abs/1306.3930}.
}
%\note{}
\seealso{
  \code{\link{cpDist}()} for a related test based on the multivariate
  empirical d.f., \code{\link{cpCopula}()} for a related test based on
  the empirical copula, \code{\link{cpAutocop}()} for a related test based on
  the empirical autocopula, \code{\link{cpRho}()} for a related test based on
  Spearman's rho, \code{\link{bOpt}()} for the function used to
  estimate \code{b} from \code{x} if \code{b = NULL}.
}
\examples{
## The standard CUSUM test based on the sample mean
cp <- cpMean(c(rnorm(50), rnorm(50, mean=1)), b=1)
cp
## Estimated change-point
which(cp$statistics == cp$statistic)

## Testing for changes in the autocovariance
n <- 200
k <- n/2 ## the true change-point
x <- c(arima.sim(list(ar = -0.5), n = k),
       arima.sim(list(ar = 0.5), n = n - k))
cp <- cpAutocov(x)
cp
## Estimated change-point
which(cp$u == cp$statistic)
## Another example
x <- c(arima.sim(list(ar = c(0,-0.5)), n = k),
       arima.sim(list(ar = c(0,0.5)), n = n - k))
cpAutocov(x)
cp <- cpAutocov(x, lag = 2)
cp
## Estimated change-point
which(cp$u == cp$statistic)

\dontrun{
## Testing for changes in Kendall's tau
require(copula)
n <- 100
k <- 50 ## the true change-point
u <- rCopula(k,gumbelCopula(1.5))
v <- rCopula(n-k,gumbelCopula(3))
x <- rbind(u,v)
cp <- cpTau(x)
cp
## Estimated change-point
which(cp$u == cp$statistic)

## Testing for changes in the covariance
cp <- cpCov(x)
cp
## Estimated change-point
which(cp$u == cp$statistic)}
}
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{ts}

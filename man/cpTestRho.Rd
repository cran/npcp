\name{cpTestRho}
\alias{cpTestRho}
\title{Test for change-point detection based on Spearman's rho}
\description{
  Nonparametric test for change-point detection particularly sensitive
  to changes in Spearman's rho in multivariate time series. The
  observations can be serially independent or dependent (strongly mixing).
  Approximate p-values for the test statistic are obtained by means of a
  \emph{multiplier} approach or by estimating the asymptotic null
  distribution. Details can be found in first reference.
}
\usage{
cpTestRho(x, method = c("mult", "asym.var"),
          statistic = c("pairwise", "global"),
          b = 1, weights = c("parzen", "bartlett"),
          N = 1000, init.seq = NULL)
}
\arguments{
  \item{x}{a data matrix whose rows are multivariate continuous
    observations.}
  \item{method}{a string specifying the method for computing the
    approximate p-value for the test statistic; can be either
    \code{"mult"} (the multiplier approach 'tilde' in the first
    reference) or \code{"asym.var"} (the approach based on the
    estimation of the asymptotic null distribution of the test statistic
    described in the first reference). The 'mult' approach appears to lead to
    better behaved tests.}
  \item{statistic}{a string specifying the test statistic; can be either
    \code{"pairwise"} (the statistic \eqn{S_{n,3}}{S[n,3]} in the first
    reference) or \code{"global"} (the statistic \eqn{S_{n,1}}{S[n,1]}
    in the first reference).}
  \item{b}{strictly positive integer specifying the value of the
    bandwidth parameter determining the serial dependence when
    generating dependent multiplier sequences using the 'moving average
    approach'; see Section 5 of the second reference. The default
    value is 1, which will create i.i.d. multiplier
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
}
\details{
  When \code{method == "mult"}, the approximate p-value is computed as
  \deqn{(0.5 +\sum_{i=1}^N\mathbf{1}_{\{S_i\ge S\}})/(N+1),}{%
        (0.5 + sum(S[i] >= S, i=1, .., N)) / (N+1),}
  where \eqn{S} and \eqn{S_i}{S[i]} denote the test statistic and
  a multiplier replication, respectively. This ensures that the
  approximate p-value is a number strictly between 0 and 1, which is
  sometimes necessary for further treatments.

  When \code{method == "asym.var"}, the approximate p-value is computed
  from the estimated asymptotic null distribution, which involves the
  Kolmogorov distribution. The latter is dealt with reusing code from
  the \code{\link{ks.test}()} function; credit to RCore.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{value of the test statistic.}
  \item{p.value}{corresponding approximate p-value.}
  \item{rho}{the values of the \code{nrow(x)-1} intermediate
    change-point statistics; the test statistic is defined as
    the maximum of those.}
  \item{b}{the value of parameter \code{b}.}
}
\references{
  I. Kojadinovic, J-F. Quessy and T. Rohmer (2014),
  Testing the constancy of Spearman's rho in multivariate time series,
  \emph{Annals of the Institute of Statistical Mathematics}, in press,
  \url{http://arxiv.org/abs/1407.1624}.

  A. B\enc{ü}{u}cher and I. Kojadinovic (2014), A dependent multiplier
  bootstrap for the sequential empirical copula process under strong
  mixing, \emph{Bernoulli}, in press, \url{http://arxiv.org/abs/1306.3930}.
}
\note{
  These tests were derived under the assumption of continuous margins.
}
\seealso{
  \code{\link{cpTestU}()} for related tests based on
  U-statistics, \code{\link{cpTestFn}()} for a related test
  based on the multivariate
  empirical c.d.f., \code{\link{cpTestCn}()} for a related test based on
  the empirical copula.
}
\examples{
\dontrun{
require(copula)
n <- 100
k <- 50 ## the true change-point
u <- rCopula(k,gumbelCopula(1.5))
v <- rCopula(n-k,gumbelCopula(3))
x <- rbind(u,v)
cp <- cpTestRho(x)
cp
## estimated change-point
which(cp$rho == max(cp$rho))
}
}
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{ts}

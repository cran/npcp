\name{cpTestU}
\alias{cpTestU}
\title{Some tests for change-point detection based on U-statistics}
\description{
  Nonparametric tests for change-point detection particularly sensitive
  to changes in certain quantities that can be estimated using
  one-sample U-statistics of order two. Thus far, the quantities under
  consideration are the variance, Gini's mean difference and Kendall's
  tau (a generic mecanism for defining the U-statistic will be
  implemented in future releases).
  The observations can be serially independent or dependent (strongly mixing).
  Approximate p-values for the test statistic are obtained by means of a
  \emph{multiplier} approach or by estimating the asymptotic null
  distribution. Details can be found in first reference.
}
\usage{
cpTestU(x, statistic = c("kendall", "variance", "gini"),
        method = c("seq", "nonseq", "asym.var"),
        b = 1, weights = c("parzen", "bartlett"),
        N = 1000, init.seq = NULL)
}
\arguments{
  \item{x}{a data matrix whose rows are continuous observations.}
  \item{statistic}{a string specifying the statistic of interest;
    can be either \code{"kendall"} (Kendall's tau, in which case
    \code{ncol(x)} must be greater than one), \code{"variance"}
    or \code{"gini"} (the variance or Gini's mean difference, in which case
    \code{ncol(x)} must be equal to one).}
  \item{method}{a string specifying the method for computing the
    approximate p-value for the test statistic; can be either
    \code{"seq"} (the 'check' approach in the first reference),
    \code{"nonseq"} (the 'hat' approach in the first reference), or
    \code{"asym.var"} (the approach based on the estimation of the
    asymptotic null distriution of the test statistic described
    in the first reference). The 'seq' approach appears overall to lead to
    better behaved tests when \code{statistic == "kendall"}. More experiments
    are necessary for the other two statistics.}
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
  When \code{method} is either \code{"seq"} or \code{"nonseq"},
  the approximate p-value is computed as
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
  \item{u}{the values of the \code{nrow(x)-3} intermediate
    change-point statistics; the test statistic is defined as
    the maximum of those.}
  \item{b}{the value of parameter \code{b}.}
}
\references{
  A. B\enc{ü}{u}cher and I. Kojadinovic (2014), Dependent multiplier
  bootstraps for non-degenerate U-statistics under mixing conditions
  with applications, \url{http://arxiv.org/abs/1412.5875}.

  A. B\enc{ü}{u}cher and I. Kojadinovic (2014), A dependent multiplier
  bootstrap for the sequential empirical copula process under strong
  mixing, \emph{Bernoulli}, in press, \url{http://arxiv.org/abs/1306.3930}.
}
\note{
  A generic mecanism for defining the U-statistic will be
  implemented in future releases.
}
\seealso{
  \code{\link{cpTestFn}()} for a related test based on the multivariate
  empirical c.d.f., \code{\link{cpTestCn}()} for a related test based on
  the empirical copula, \code{\link{cpTestRho}()} for a related test based on
  Spearman's rho.
}
\examples{
\dontrun{
require(copula)
n <- 100
k <- 50 ## the true change-point
u <- rCopula(k,gumbelCopula(1.5))
v <- rCopula(n-k,gumbelCopula(3))
x <- rbind(u,v)
cp <- cpTestU(x)
cp
## estimated change-point
which(cp$u == max(cp$u))
}
}
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{ts}

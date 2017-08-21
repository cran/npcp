\name{cpDist}
\alias{cpDist}
\title{Test for Change-Point Detection in Possibly Multivariate Observations
  Sensitive to Changes in the Distribution Function}
\description{
  Nonparametric test for change-point detection based on the
  (multivariate) empirical distribution function. The observations can
  be continuous univariate or multivariate, and serially independent or
  dependent (strongly mixing).  Approximate p-values for the test
  statistics are obtained by means of a \emph{multiplier}
  approach. The first reference treats the serially
  independent case while details about the serially
  dependent case can be found in second and third references.
}
\usage{
cpDist(x, statistic = c("cvmmax", "cvmmean", "ksmax", "ksmean"),
       method = c("nonseq", "seq"), b = NULL,
       weights = c("parzen", "bartlett"),
       m = 5, L.method=c("max","median","mean","min"),
       N = 1000, init.seq = NULL, include.replicates = FALSE)
}
\arguments{
  \item{x}{a data matrix whose rows are continuous observations.}
  \item{statistic}{a string specifying the statistic whose value and
    p-value will be displayed; can be either \code{"cvmmax"} or
    \code{"cvmmean"} (the maximum or average of the \code{nrow(x)-1}
    intermediate Cram\enc{é}{e}r-von Mises statistics), or
    \code{"ksmax"} or \code{"ksmean"} (the maximum or average of the
    \code{nrow(x)-1} intermediate Kolmogorov-Smirnov statistics); see
    Section 3 in the first reference. The four statistics and the
    corresponding p-values are computed at each execution.}
  \item{method}{a string specifying the simulation method for
    generating multiplier replicates of the test statistic;
    can be either \code{"nonseq"} (the 'check' approach
    in the first reference) or \code{"seq"} (the 'hat' approach
    in the first reference). The 'check' approach appears to lead to
    better behaved tests and is recommended.}
  \item{b}{strictly positive integer specifying the value of the
    bandwidth parameter determining the serial dependence when
    generating dependent multiplier sequences using the 'moving average
    approach'; see Section 5 of the second reference. The
    value 1 will create i.i.d. multiplier
    sequences suitable for serially independent observations. If set to
    \code{NULL}, \code{b} will be estimated from \code{x} using the
    function \code{\link{bOptEmpProc}()}; see the procedure described in
    Section 5 of the second reference.}
 \item{weights}{a string specifying the kernel for creating the
   weights used in the generation of dependent multiplier sequences within the
   'moving average approach'; see Section 5 of the second reference.}
 \item{m}{a strictly positive integer specifying the number of points of the
   uniform grid on \eqn{(0,1)^d}{(0,1)^d} (where \eqn{d}{d} is
   \code{ncol(x)}) involved in the estimation of the bandwidth
   parameter; see Section 5 of the third reference. The number of
   points of the grid is given by \code{m^ncol(x)} so that \code{m} needs to be
   decreased as \eqn{d}{d} increases.}
 \item{L.method}{a string specifying how the parameter \eqn{L}{L} involved
   in the estimation of the bandwidth parameter is computed; see
   Section 5 of the second reference.}
 \item{N}{number of multiplier replications.}
 \item{init.seq}{a sequence of independent standard normal variates of
   length \code{N * (nrow(x) + 2 * (b - 1))} used to generate dependent
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
  \item{cvm}{the values of the \code{nrow(x)-1} intermediate
    Cram\enc{é}{e}r-von Mises change-point statistics.}
  \item{ks}{the values of the \code{nrow(x)-1} intermediate
    Kolmogorov-Smirnov change-point statistics.}
  \item{all.statistics}{the values of all four test statistics.}
  \item{all.p.values}{the corresponding p-values.}
  \item{b}{the value of parameter \code{b}.}
}
\references{
  M. Holmes, I. Kojadinovic and J-F. Quessy (2013), Nonparametric tests
  for change-point detection \enc{à}{a} la Gombay and Horv\enc{á}{a}th,
  \emph{Journal of Multivariate Analysis} \bold{115}, pages 16-32.

  A. B\enc{ü}{u}cher and I. Kojadinovic (2016), A dependent multiplier
  bootstrap for the sequential empirical copula process under strong
  mixing, \emph{Bernoulli} \bold{22:2}, pages 927-968,
  \url{http://arxiv.org/abs/1306.3930}.

  A. B\enc{ü}{u}cher, J.-D. Fermanian and I. Kojadinovic (2017),
  Combining cumulative sum change-point detection tests for assessing
  the stationarity of continuous univariate time series,
  \url{http://arxiv.org/}.
}
\note{
  Note that when the observations are continuous univariate and serially
  independent, independent realizations of the tests statistics under
  the null hypothesis of no change in the distribution can be
  obtained by simulation; see Section 4 in the first reference.
}
\seealso{
  \code{\link{cpCopula}()} for a related test based on the empirical
  copula, \code{\link{cpRho}()} for a related test based on
  Spearman's rho, \code{\link{cpTau}()} for a related test based on
  Kendall's tau, \code{\link{bOptEmpProc}()} for the function used to
  estimate \code{b} from \code{x} if \code{b = NULL}.
}
\examples{
## A univariate example
n <- 100
k <- 50 ## the true change-point
y <- rnorm(k)
z <- rexp(n-k)
x <- matrix(c(y,z))
cp <- cpDist(x, b = 1)
cp

## All statistics
cp$all.statistics
## Corresponding p.values
cp$all.p.values

## Estimated change-point
which(cp$cvm == max(cp$cvm))
which(cp$ks == max(cp$ks))

## A very artificial trivariate example
## with a break in the first margin
n <- 100
k <- 50 ## the true change-point
y <- rnorm(k)
z <- rnorm(n-k, mean = 2)
x <- cbind(c(y,z),matrix(rnorm(2*n), n, 2))
cp <- cpDist(x, b = 1)
cp

## All statistics
cp$all.statistics
## Corresponding p.values
cp$all.p.values

## Estimated change-point
which(cp$cvm == max(cp$cvm))
which(cp$ks == max(cp$ks))
}
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{ts}
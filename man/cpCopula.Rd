\name{cpCopula}
\alias{cpCopula}
\title{Test for Change-Point Detection in Multivariate Observations
  Sensitive to Changes in the Copula}
\description{
  Nonparametric test for change-point detection particularly sensitive
  to changes in the copula of multivariate continuous observations. The
  observations can be serially independent or dependent (strongly mixing).
  Approximate p-values for the test statistic are obtained by means of a
  \emph{multiplier} approach. Details can be found in the first reference.
}
\usage{
cpCopula(x, method = c("seq", "nonseq"), b = NULL,
         weights = c("parzen", "bartlett"), m = 5,
         L.method=c("max","median","mean","min"),
         N = 1000, init.seq = NULL, include.replicates = FALSE)
}
\arguments{
  \item{x}{a data matrix whose rows are multivariate continuous
    observations.}
  \item{method}{a string specifying the simulation method for
    generating multiplier replicates of the test statistic;
    can be either \code{"seq"} (the 'check' approach
    in the first reference) or \code{"nonseq"} (the 'hat' approach
    in the first reference). The 'check' approach appears to lead to
    better behaved tests in the case of samples of moderate size. The
    'hat' approach is substantially faster.}
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
    Cram\enc{é}{e}r-von Mises change-point statistics; the test
    statistic is defined as the maximum of those.}
  \item{b}{the value of parameter \code{b}.}
}
\references{
  A. B\enc{ü}{u}cher, I. Kojadinovic, T. Rohmer and J. Segers (2014),
  Detecting changes in cross-sectional dependence in multivariate
  time series, \emph{Journal of Multivariate Analysis} \bold{132}, pages 111-128,
  \url{http://arxiv.org/abs/1206.2557}.

  A. B\enc{ü}{u}cher and I. Kojadinovic (2016), A dependent multiplier
  bootstrap for the sequential empirical copula process under strong
  mixing, \emph{Bernoulli} \bold{22:2}, pages 927-968,
  \url{http://arxiv.org/abs/1306.3930}.
}
\note{
  These tests were derived under the assumption of continuous margins.
}
\seealso{
  \code{\link{cpRho}()} for a related test based on
  Spearman's rho, \code{\link{cpTau}()} for a related test based on
  Kendall's tau, \code{\link{cpDist}()} for a related test based
  on the multivariate empirical d.f., \code{\link{bOptEmpProc}()} for the
  function used to estimate \code{b} from \code{x} if \code{b = NULL}.
}
\examples{
\dontrun{
require(copula)
n <- 100
k <- 50 ## the true change-point
u <- rCopula(k, gumbelCopula(1.5))
v <- rCopula(n - k, gumbelCopula(3))
x <- rbind(u,v)
cp <- cpCopula(x, b = 1)
cp
## Estimated change-point
which(cp$cvm == max(cp$cvm))}
}
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{ts}

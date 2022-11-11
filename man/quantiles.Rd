\name{quantiles}
\alias{quantiles}
\title{Estimated Quantiles for the Open-end Nonparametric Sequential
  Change-Point Detection Tests}
\description{
  Estimated quantiles for the open-end nonparametric sequential
  change-point detection tests described in
  \code{\link{seqOpenEndCpMean}} and
  \code{\link{seqOpenEndCpDist}}. More details can be found in the
  references below.
}
\usage{
data("quantiles")
}
\format{
  \code{\link{list}} of 6 elements. The first 5 are arrays containing
  the estimated 90\%, 95\% and 99\% quantiles necessary for carrying out
  the sequential tests described in \code{\link{seqOpenEndCpMean}}. The
  last element is a list containing the estimated 90\%, 95\% and 99\%
  quantiles as well as other estimated parameters necessary for carrying
  out the sequential test described in \code{\link{seqOpenEndCpDist}}.
}
\references{
  J. G\enc{ö}{o}smann, T. Kley and H. Dette (2021), A new approach for
  open-end sequential change point monitoring, \emph{Journal of Time
  Series Analysis} \bold{42:1}, pages 63-84,
  \url{https://arxiv.org/abs/1906.03225}.

   M. Holmes and I. Kojadinovic (2021), Open-end nonparametric
  sequential change-point detection based on the retrospective CUSUM
  statistic, \emph{Electronic Journal of Statistics} \bold{15:1}, pages
  2288-2335, \doi{10.1214/21-EJS1840}.

  L. Horv\enc{á}{a}th, M. Hu\enc{š}{s}kov\enc{á}{a}, P. Kokoszka and
  J. Steinebach (2004). Monitoring changes in linear
  models. \emph{Journal of Statistical Planning and Inference}
  \bold{126}, pages 225-251.

  M. Holmes, I. Kojadinovic and A. Verhoijsen (2022), Multi-purpose
  open-end monitoring procedures for multivariate observations based on
  the empirical distribution function, 45 pages,
  \url{https://arxiv.org/abs/2201.10311}.
}
\examples{
data("quantiles")
str(quantiles)
}
\keyword{datasets}

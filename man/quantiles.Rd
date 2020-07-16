\name{quantiles}
\alias{quantiles}
\title{Estimated Quantiles for the Open-end Nonparametric Sequential
  Change-Point Detection Tests Sensitive to Changes in the Mean}
\description{
  Estimated quantiles for the open-end nonparametric sequential
  change-point detection tests in \code{\link{seqCpMean}} sensitive to
  changes in the mean. More details can be found in the references below.
}
\usage{
data("quantiles")
}
\format{
  \code{\link{list}} of 5 arrays containing the estimated 90\%, 95\%
  and 99\% quantiles necessary for carrying out the sequential tests
  in \code{\link{seqCpMean}} described in the references below.
}
%% \source{
%%   M. Holmes and I. Kojadinovic (2020), Open-end nonparametric sequential
%%   change-point detection based on the retrospective CUSUM statistic, 40
%%   pages, available soon on the arXiv.
%% }
\references{
  J. G\enc{ö}{o}smann, T. Kley and H. Dette (2020), A new approach for
  open-end sequential change point monitoring, 49 pages,
  \url{http://arxiv.org/abs/1906.03225}.

  M. Holmes and I. Kojadinovic (2020), Open-end nonparametric sequential
  change-point detection based on the retrospective CUSUM statistic, 41
  pages, available on the arXiv.

  L. Horv\enc{á}{a}th, M. Hu\enc{š}{s}kov\enc{á}{a}, P. Kokoszka and
  J. Steinebach (2004). Monitoring changes in linear
  models. \emph{Journal of Statistical Planning and Inference}
  \bold{126}, pages 225-251.
}
\examples{
data("quantiles")
str(quantiles)
}
\keyword{datasets}

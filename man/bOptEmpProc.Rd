\name{bOptEmpProc}
\alias{bOpt}
\alias{bOptEmpProc}
\title{Bandwidth Parameter Estimation}
\description{
  In the context of the standard CUSUM test based on the sample mean or
  in a particular empirical process setting, the following functions
  estimate the bandwidth parameter controlling the serial dependence
  when generating dependent multiplier sequences using the 'moving
  average approach'; see Section 5 of the third reference. The function
  function \code{bOpt()} is called in the functions
  \code{\link{cpMean}()}, \code{\link{cpVar}()}, \code{\link{cpGini}()},
  \code{\link{cpAutocov}()}, \code{\link{cpCov}()},
  \code{\link{cpTau}()} and \code{\link{detOpenEndCpMean}()} when \code{b} is
  set to \code{NULL}. The function function \code{bOptEmpProc()} is
  called in the functions \code{\link{cpDist}()},
  \code{\link{cpCopula}()}, \code{\link{cpAutocop}()},
  \code{\link{stDistAutocop}()} and \code{\link{simClosedEndCpDist}()} when
  \code{b} is set to \code{NULL}.
}
\usage{
bOpt(influ, weights = c("parzen", "bartlett"))

bOptEmpProc(x, m=5, weights = c("parzen", "bartlett"),
            L.method=c("max","median","mean","min"))
}
\arguments{
  \item{influ}{a numeric containing the relevant influence coefficients,
            which, in the case of the standard CUSUM test based on the
            sample mean, are simply the available observations;
            see also the last reference.}
  \item{x}{a data matrix whose rows are continuous observations.}
  \item{weights}{a string specifying the kernel for creating the
   weights used in the generation of dependent multiplier sequences within the
   'moving average approach'; see Section 5 of the third reference.}
  \item{m}{a strictly positive integer specifying the number of points of the
   uniform grid on \eqn{(0,1)^d}{(0,1)^d} (where \eqn{d}{d} is
   \code{ncol(x)}) involved in the estimation of the bandwidth
   parameter; see Section 5 of the third reference. The number of
   points of the grid is given by \code{m^ncol(x)} so that \code{m} needs to be
   decreased as \eqn{d}{d} increases.}
 \item{L.method}{a string specifying how the parameter \code{L} involved
   in the estimation of the bandwidth parameter is computed; see
   Section 5 of the third reference.}
}
\details{
  The implemented approach results from an adaptation of the procedure
  described in the first two references (see also the references
  therein). The use of theses functions in a context different from that
  considered in the third or fourth reference may not be meaningful.

  Acknowledgment: Part of the code of the function results from an
  adaptation of R code of C. Parmeter and J. Racine, itself an
  adaptation of Matlab code by A. Patton.
}
\value{A strictly positive integer.}

\references{
  D.N. Politis and H. White (2004), Automatic block-length selection for the
  dependent bootstrap, \emph{Econometric Reviews} \bold{23(1)}, pages 53-70.

  D.N. Politis, H. White and A.J. Patton (2004), Correction: Automatic
  block-length selection for the dependent bootstrap,
  \emph{Econometric Reviews} \bold{28(4)}, pages 372-375.

  A. B\enc{ü}{u}cher and I. Kojadinovic (2016), A dependent multiplier
  bootstrap for the sequential empirical copula process under strong
  mixing, \emph{Bernoulli} \bold{22:2}, pages 927-968,
  \url{https://arxiv.org/abs/1306.3930}.

  A. B\enc{ü}{u}cher and I. Kojadinovic (2016), Dependent multiplier
  bootstraps for non-degenerate U-statistics under mixing conditions
  with applications, \emph{Journal of Statistical Planning and
    Inference} \bold{170} pages 83-105, \url{https://arxiv.org/abs/1412.5875}.

}
%\note{}
\seealso{
  \code{\link{cpDist}()}, \code{\link{cpCopula}()},
  \code{\link{cpAutocop}()}, \code{\link{stDistAutocop}()},
  \code{\link{cpMean}()}, \code{\link{cpVar}()}, \code{\link{cpGini}()},
  \code{\link{cpAutocov}()}, \code{\link{cpCov}()},
  \code{\link{cpTau}()}, \code{\link{seqOpenEndCpMean}} and
  \code{\link{seqClosedEndCpDist}}.
}
%\examples{}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{ts}

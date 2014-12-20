\name{bOptEmpProc}
\alias{bOptEmpProc}
\title{Bandwidth parameter estimation}
\description{
  In a particular empirical process setting, estimates the
  bandwidth parameter controlling the serial dependence when
  generating dependent multiplier sequences using the 'moving average
  approach'; see Sections 5 and 6.1 of the third reference. This
  function is called in the functions \code{\link{cpTestFn}()} and
  \code{\link{cpTestCn}()} if \code{b} is set to \code{NULL}.
}
\usage{
bOptEmpProc(x, m=5, weights = c("parzen", "bartlett"),
            L.method=c("max","median","mean","min"))
}
\arguments{
  \item{x}{a data matrix whose rows are continuous observations.}
  \item{weights}{a string specifying the kernel for creating the
   weights used in the generation of dependent multiplier sequences within the
   'moving average approach'; see Section 6.1 of the third reference.}
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
  therein). The use of this function in a context different from that
  considered in the third reference may not be meaningful.

  Acknowledgment: Part of the code of the function results from an
  adaptation of R code of C. Parmeter and J. Racine, itself an
  adaptation of Matlab code by A. Patton.
}
\value{A strictly positive integer.}

\references{
  D.N. Politis and H. White (2004), Automatic block-length selection for the
  dependent bootstrap, \emph{Econometric Reviews} 23(1):53–70.

  D.N. Politis, H. White and A.J. Patton (2004), Correction: Automatic
  block-length selection for the dependent bootstrap,
  \emph{Econometric Reviews} 28(4):372-375.

  A. B\enc{ü}{u}cher and I. Kojadinovic (2014), A dependent multiplier
  bootstrap for the sequential empirical copula process under strong
  mixing, Bernoulli, in press, \url{http://arxiv.org/abs/1306.3930}.
}
%\note{}
\seealso{
  \code{\link{cpTestFn}()}, \code{\link{cpTestCn}()}.
}
%\examples{}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{ts}

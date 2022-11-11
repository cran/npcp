\name{selectPoints}

\alias{selectPoints}

\title{
  A point selection procedure for multivariate data
}
\description{
  Returns a matrix of `representative' points.
}
\usage{
selectPoints(x, r, kappa = 1.5, plot = FALSE)
}
\arguments{
  \item{x}{a numeric matrix with \code{d} columns whose rows represent
    multivariate observations.}
  \item{r}{integer specifying the size of an initial uniformly-spaced
  grid `on the probability scale'; an upper bound for the number of
  selected points is \code{r^d}.}
  \item{kappa}{numeric constant required to be strictly greater than one
  involved in the point selection procedure.}
  \item{plot}{logical used only if \code{d = 2} specifying whether a plot
should be produced.}
}
\details{
  The selection procedure is described in detail in Section 3.2 of the
  reference below. Set \code{plot = TRUE} for visual feedback and
  information on the minimum number of `surrounding' observations
  required for a grid point to be selected. The initial grid `on the
  probability scale' is in blue, while the points selected by the procedure
  are in red.
}
\value{
  a matrix with \code{d} columns whose rows are the selected points.
}
\references{
  M. Holmes, I. Kojadinovic, and A. Verhoijsen, Multi-purpose open-end
  monitoring procedures for multivariate observations based on the
  empirical distribution function, 45 pages,
  \url{https://arxiv.org/abs/2201.10311}.
}
\seealso{
  \code{selectPoints()} is used in \code{\link{detOpenEndCpDist}()}.
}
\examples{
## Generate data
set.seed(123)
x1 <- rnorm(1000, 0, 1)
x2 <- rnorm(1000, 0.7 * x1, sqrt((1 - 0.7^2)))
x <- cbind(x1, x2)

## Point selection
selectPoints(x, r = 3, kappa = 1.5, plot = TRUE)
selectPoints(x, r = 3, kappa = 4, plot = TRUE)

selectPoints(x, r = 5, kappa = 1.5, plot = TRUE)
selectPoints(x, r = 5, kappa = 4, plot = TRUE)
}
\keyword{nonparametric}
\keyword{multivariate}

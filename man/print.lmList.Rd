% $Id: print.lmList.Rd,v 1.2 2000/04/17 19:03:40 bates Exp $
\name{print.lmList}
\title{Print an lmList Object}
\usage{
print(x, pool, ...)
}
\alias{print.lmList}
\arguments{
 \item{x}{an object inheriting from class \code{lmList}, representing
   a list of fitted \code{lm} objects.
 }
 \item{pool}{an optional logical value indicating whether a pooled
   estimate of the residual standard error should be used. Default is
   \code{attr(object, "pool")}.
 }
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.
 }
}
\description{
  Information describing the individual \code{lm} fits corresponding to
  \code{object} is printed. This includes the estimated coefficients and
  the residual standard error.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}}

\examples{
data(Orthodont)
fm1 <- lmList(Orthodont)
print(fm1)
}
\keyword{models}

% $Id: print.summary.lmList.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{print.summary.lmList}
\title{Print a summary.lmList Object}
\usage{
\method{print}{summary.lmList}(x, ...)
}
\alias{print.summary.lmList}
\arguments{
 \item{x}{an object inheriting from class \code{summary.lmList},
   representing a summarized \code{lme} object.
 }
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.
 }
}
\description{
  Information summarizing the individual \code{lm} fitted objects
  corresponding to \code{x} is printed. This includes the estimated
  coefficients and their respective standard errors, t-values, and
  p-values. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{summary.lmList}}, \code{\link{lmList}}}

\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
print(summary(fm1))
}
\keyword{models}

% $Id: print.intervals.lme.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{print.intervals.lme}
\title{Print an intervals.lme Object}
\usage{
print(x, ...)
}
\alias{print.intervals.lme}
\arguments{
 \item{x}{an object inheriting from class \code{intervals.lme},
   representing a list of data frames with confidence intervals and
   estimates for the coefficients in the \code{lme} object that produced
   \code{x}.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  The individual components of \code{x} are printed.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{intervals.lme}}}

\examples{
library(nlme)
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
print(intervals(fm1))
}
\keyword{models}

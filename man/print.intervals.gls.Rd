% $Id: print.intervals.gls.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{print.intervals.gls}
\title{Print an intervals.gls Object}
\usage{
print(x, ...)
}
\alias{print.intervals.gls}
\arguments{
 \item{x}{an object inheriting from class \code{intervals.gls},
   representing a list of data frames with confidence intervals and
   estimates for the coefficients in the \code{gls} object that produced
   \code{x}.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  The individual components of \code{x} are printed.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{intervals.gls}}}

\examples{
library(nlme)
data(Ovary)
%\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
print(intervals(fm1))
%}
}
\keyword{models}

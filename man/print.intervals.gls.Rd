% $Id: print.intervals.gls.Rd,v 1.4 2000/07/03 18:22:48 bates Exp $
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
data(Ovary)
%\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
print(intervals(fm1))
%}
}
\keyword{models}

% $Id: formula.gls.Rd,v 1.7 2001/10/30 20:51:15 bates Exp $
\name{formula.gls}
\title{Extract gls Object Formula}
\usage{
\method{formula}{gls}(x, \dots)
}
\alias{formula.gls}
\arguments{
 \item{x}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function extracts the linear model formula
  associated with \code{x}. 
}
\value{
  a two-sided linear formula specifying the linear model used to
  obtain \code{x}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}}

\examples{
data(Ovary)
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
formula(fm1)
}
\keyword{models}

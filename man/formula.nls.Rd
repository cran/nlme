% $Id: formula.nls.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{formula.nls}
\title{Extract Model Formula from nls Object}
\usage{
formula(object)
}
\alias{formula.nls}
\arguments{
 \item{object}{an object inheriting from class \code{nls}, representing
   a non-linear least squares fit.}
}
\description{
  Returns the model used to fit \code{object}. 
}
\value{
  a formula representing the model used to obtain \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{nls}}, \code{\link{formula}}}
\examples{
data(Orange)
fm1 <- nls(circumference ~ A/(1+exp((B-age)/C)), Orange,
  start = list(A=160, B=700, C = 350))
formula(fm1)
}
\keyword{models}

% $Id: formula.corStruct.Rd,v 1.6 2001/10/30 20:51:15 bates Exp $
\name{formula.corStruct}
\title{Extract corStruct Object Formula}
\usage{
\method{formula}{corStruct}(x, \dots)
}
\alias{formula.corStruct}
\arguments{
 \item{x}{an object inheriting from class \code{corStruct}
   representing a correlation structure.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function extracts the formula associated with a
  \code{corStruct} object, in which the covariate and the grouping
  factor, if any is present, are defined. 
}
\value{
  an object of class \code{formula} specifying the covariate and the
  grouping factor, if any is present, associated with \code{x}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{formula}}}

\examples{
cs1 <- corCAR1(form = ~ Time | Mare)
formula(cs1)
}
\keyword{models}

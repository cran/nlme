% $Id: formula.corStruct.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{formula.corStruct}
\title{Extract corStruct Object Formula}
\usage{
formula(object)
}
\alias{formula.corStruct}
\arguments{
 \item{object}{an object inheriting from class \code{corStruct}
   representing a correlation structure.}
}
\description{
  This method function extracts the formula associated with a
  \code{corStruct} object, in which the covariate and the grouping
  factor, if any is present, are defined. 
}
\value{
  an object of class \code{formula} specifying the covariate and the
  grouping factor, if any is present, associated with \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{formula}}}

\examples{
library(nlme)
cs1 <- corCAR1(form = ~ Time | Mare)
formula(cs1)
}
\keyword{models}

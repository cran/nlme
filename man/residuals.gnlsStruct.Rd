% $Id: residuals.gnlsStruct.Rd,v 1.4 2000/07/03 18:22:48 bates Exp $
\name{residuals.gnlsStruct}
\title{Calculate gnlsStruct Residuals}
\usage{
residuals(object)
}
\alias{residuals.gnlsStruct}
\arguments{
 \item{object}{an object inheriting from class \code{gnlsStruct},
   representing a list of model components, such as
   \code{corStruct} and \code{varFunc} objects, and attributes
   specifying the underlying nonlinear model and the response variable.}
}
\description{
  The residuals for the nonlinear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the residuals for the nonlinear model represented by
  \code{object}. 
}
\author{Jose Pinheiro and Douglas Bates}
\note{
  This method function is primarily used inside \code{gnls} and 
  \code{residuals.gnls}.
}
\seealso{\code{\link{gnls}}, \code{\link{residuals.gnls}},
  \code{\link{fitted.gnlsStruct}} } 
\keyword{models}

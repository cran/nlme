% $Id: fitted.gnlsStruct.Rd,v 1.1 2000/03/17 22:21:22 saikat Exp $
\name{fitted.gnlsStruct}
\title{Calculate gnlsStruct Fitted Values}
\usage{
fitted(object)
}
\alias{fitted.gnlsStruct}
\arguments{
 \item{object}{an object inheriting from class \code{gnlsStruct},
   representing a list of model components, such as
   \code{corStruct} and \code{varFunc} objects, and attributes
   specifying the underlying nonlinear model and the response variable.}
}
\description{
  The fitted values for the nonlinear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the fitted values for the nonlinear model represented by
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}
\note{
  This method function is generally only used inside \code{gnls} and 
  \code{fitted.gnls}.
}

\seealso{\code{\link{gnls}}, \code{\link{fitted.gnls}},
  \code{\link{residuals.gnlsStruct}} } 

\keyword{models}

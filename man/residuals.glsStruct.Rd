% $Id: residuals.glsStruct.Rd,v 1.1 1999/10/13 00:50:12 saikat Exp $
\name{residuals.glsStruct}
\title{Calculate glsStruct Residuals}
\usage{
residuals(object, glsFit)
}
\alias{residuals.glsStruct}
\arguments{
 \item{object}{an object inheriting from class \code{glsStruct},
   representing a list of linear model components, such as
   \code{corStruct} and \code{varFunc} objects.}
 \item{glsFit}{an optional list with components \code{logLik}   
   (log-likelihood), \code{beta} (coefficients), \code{sigma} (standard
   deviation for error term), \code{varBeta} (coefficients' covariance
   matrix), \code{fitted} (fitted values), and \code{residuals}
   (residuals). Defaults to \code{attr(object, "glsFit")}.}
}
\description{
  The residuals for the linear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the residuals for the linear model represented by
  \code{object}. 
}
\author{Jose Pinheiro and Douglas Bates}
\note{
  This method function is primarily used inside \code{gls} and 
  \code{residuals.gls}.
}
\seealso{\code{\link{gls}}, \code{\link{residuals.gls}},
  \code{\link{fitted.glsStruct}} } 

\keyword{models}

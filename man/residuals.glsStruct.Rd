% $Id: residuals.glsStruct.Rd,v 1.5 2001/10/30 20:51:15 bates Exp $
\name{residuals.glsStruct}
\title{Calculate glsStruct Residuals}
\usage{
\method{residuals}{glsStruct}(object, glsFit, \dots)
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
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
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

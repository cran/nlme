% $Id: fitted.glsStruct.Rd,v 1.3 2000/07/03 18:22:46 bates Exp $
\name{fitted.glsStruct}
\title{Calculate glsStruct Fitted Values}
\usage{
fitted(object, glsFit)
}
\alias{fitted.glsStruct}
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
  The fitted values for the linear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the fitted values for the linear model represented by
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}
\note{
  This method function is generally only used inside \code{gls} and 
  \code{fitted.gls}.
}

\seealso{\code{\link{gls}}, \code{\link{fitted.gls}},
  \code{\link{residuals.glsStruct}} } 

\keyword{models}
